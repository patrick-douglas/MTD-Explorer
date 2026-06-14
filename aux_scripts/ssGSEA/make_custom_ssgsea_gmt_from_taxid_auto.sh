#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Build custom ssGSEA GMT from TaxID using NCBI Datasets + eggNOG
#
# Logic:
#   1. Download protein/GFF/GTF by taxid or assembly using NCBI Datasets
#   2. Check if downloaded IDs are compatible with host.gct
#   3. Rename representative proteins to host.gct gene IDs
#   4. Run eggNOG-mapper
#   5. Convert eggNOG GO annotations to GMT
#   6. Validate final GMT overlap with host.gct
# ------------------------------------------------------------

usage() {
cat << EOF
Usage:
  bash make_custom_ssgsea_gmt_from_taxid_auto.sh \\
    --taxid 6526 \\
    --host-gct /path/to/host.gct \\
    --outdir /path/to/custom_gmt \\
    --threads 20

Required:
  --taxid INT              NCBI taxid
  --host-gct FILE          MTD ssGSEA host.gct file
  --outdir DIR             Output directory for custom GMT

Optional:
  --assembly ACCESSION     Prefer a specific assembly, e.g. GCA_947242115.1
  --protein-fasta FILE     Use local protein FASTA instead of downloading from NCBI
  --annotation-gff FILE    Optional local GFF/GFF3/GTF for protein-to-gene mapping
  --threads INT            Default: nproc
  --eggnog-db DIR          Default: \$HOME/MTD/eggnog_db
  --workdir DIR            Default: OUTDIR/work_taxid_TAXID
  --min-overlap-genes INT  Minimum compatible genes before continuing. Default: 100
  --min-overlap-pct FLOAT  Minimum % of host.gct genes compatible. Default: 1
  --force                  Re-run eggNOG and overwrite intermediate files

Outputs:
  custom_taxid_<taxid>_eggNOG_GO.gmt
  custom_taxid_<taxid>_eggNOG_GO_gene_table.tsv
  custom_taxid_<taxid>_eggNOG_GO_summary.tsv
  custom_taxid_<taxid>_eggNOG_GO_overlap_report.txt
EOF
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

log() {
    echo "[INFO] $*"
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"
}

TAXID=""
ASSEMBLY=""
HOST_GCT=""
OUTDIR=""
PROTEIN_FASTA=""
ANNOTATION_GFF=""
THREADS="$(nproc)"
EGGNOG_DB="$HOME/MTD/eggnog_db"
WORKDIR=""
MIN_OVERLAP_GENES="100"
MIN_OVERLAP_PCT="1"
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --taxid)
            TAXID="$2"
            shift 2
            ;;
        --taxid=*)
            TAXID="${1#*=}"
            shift
            ;;
        --assembly)
            ASSEMBLY="$2"
            shift 2
            ;;
        --assembly=*)
            ASSEMBLY="${1#*=}"
            shift
            ;;
        --host-gct)
            HOST_GCT="$2"
            shift 2
            ;;
        --host-gct=*)
            HOST_GCT="${1#*=}"
            shift
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --outdir=*)
            OUTDIR="${1#*=}"
            shift
            ;;
        --protein-fasta)
            PROTEIN_FASTA="$2"
            shift 2
            ;;
        --protein-fasta=*)
            PROTEIN_FASTA="${1#*=}"
            shift
            ;;
        --annotation-gff)
            ANNOTATION_GFF="$2"
            shift 2
            ;;
        --annotation-gff=*)
            ANNOTATION_GFF="${1#*=}"
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --threads=*)
            THREADS="${1#*=}"
            shift
            ;;
        --eggnog-db)
            EGGNOG_DB="$2"
            shift 2
            ;;
        --eggnog-db=*)
            EGGNOG_DB="${1#*=}"
            shift
            ;;
        --workdir)
            WORKDIR="$2"
            shift 2
            ;;
        --workdir=*)
            WORKDIR="${1#*=}"
            shift
            ;;
        --min-overlap-genes)
            MIN_OVERLAP_GENES="$2"
            shift 2
            ;;
        --min-overlap-genes=*)
            MIN_OVERLAP_GENES="${1#*=}"
            shift
            ;;
        --min-overlap-pct)
            MIN_OVERLAP_PCT="$2"
            shift 2
            ;;
        --min-overlap-pct=*)
            MIN_OVERLAP_PCT="${1#*=}"
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[[ -n "$TAXID" ]] || die "Missing --taxid"
[[ -n "$HOST_GCT" ]] || die "Missing --host-gct"
[[ -n "$OUTDIR" ]] || die "Missing --outdir"

require_file "$HOST_GCT"

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    die "--threads must be a positive integer. Got: $THREADS"
fi

mkdir -p "$OUTDIR"

if [[ -z "$WORKDIR" ]]; then
    WORKDIR="$OUTDIR/work_taxid_${TAXID}"
fi

mkdir -p "$WORKDIR"

log "TaxID: $TAXID"
log "Assembly: ${ASSEMBLY:-auto}"
log "host.gct: $HOST_GCT"
log "Output directory: $OUTDIR"
log "Work directory: $WORKDIR"
log "eggNOG DB directory: $EGGNOG_DB"
log "Threads: $THREADS"

require_cmd wget
require_cmd gunzip
require_cmd tar
require_cmd unzip
require_cmd python3
require_cmd Rscript
require_cmd emapper.py

if [[ -z "$PROTEIN_FASTA" ]]; then
    require_cmd datasets
fi

# ------------------------------------------------------------
# eggNOG database download/check
# ------------------------------------------------------------

mkdir -p "$EGGNOG_DB"

download_gzip_if_missing() {
    local url="$1"
    local final_file="$2"
    local gz_file="$3"

    if [[ -s "$final_file" ]]; then
        log "Found existing eggNOG file, skipping: $final_file"
        return 0
    fi

    log "Downloading: $url"
    wget -c -O "$gz_file" "$url"

    [[ -s "$gz_file" ]] || die "Download failed: $gz_file"

    log "Decompressing: $gz_file"
    gunzip -f "$gz_file"

    [[ -s "$final_file" ]] || die "Expected file not created: $final_file"
}

download_gzip_if_missing \
    "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz" \
    "$EGGNOG_DB/eggnog.db" \
    "$EGGNOG_DB/eggnog.db.gz"

download_gzip_if_missing \
    "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz" \
    "$EGGNOG_DB/eggnog_proteins.dmnd" \
    "$EGGNOG_DB/eggnog_proteins.dmnd.gz"

if [[ -s "$EGGNOG_DB/eggnog.taxa.db" && -s "$EGGNOG_DB/eggnog.taxa.db.traverse.pkl" ]]; then
    log "Found existing eggNOG taxa files, skipping."
else
    log "Downloading eggNOG taxa database..."
    wget -c -O "$EGGNOG_DB/eggnog.taxa.tar.gz" \
        "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz"

    [[ -s "$EGGNOG_DB/eggnog.taxa.tar.gz" ]] || die "Download failed: $EGGNOG_DB/eggnog.taxa.tar.gz"

    log "Extracting eggNOG taxa database..."
    tar -zxf "$EGGNOG_DB/eggnog.taxa.tar.gz" -C "$EGGNOG_DB"
    rm -f "$EGGNOG_DB/eggnog.taxa.tar.gz"

    [[ -s "$EGGNOG_DB/eggnog.taxa.db" ]] || die "Missing after extraction: $EGGNOG_DB/eggnog.taxa.db"
fi

log "eggNOG DB is ready:"
ls -lh "$EGGNOG_DB"/eggnog.db "$EGGNOG_DB"/eggnog_proteins.dmnd "$EGGNOG_DB"/eggnog.taxa.db 2>/dev/null || true

# ------------------------------------------------------------
# Download NCBI protein/GFF/GTF if local protein FASTA not provided
# ------------------------------------------------------------

PROTEIN_LIST="$WORKDIR/protein_candidates.txt"
GFF_LIST="$WORKDIR/annotation_candidates.txt"

rm -f "$PROTEIN_LIST" "$GFF_LIST"
touch "$PROTEIN_LIST" "$GFF_LIST"

if [[ -n "$PROTEIN_FASTA" ]]; then
    require_file "$PROTEIN_FASTA"
    echo "$PROTEIN_FASTA" > "$PROTEIN_LIST"

    if [[ -n "$ANNOTATION_GFF" ]]; then
        require_file "$ANNOTATION_GFF"
        echo "$ANNOTATION_GFF" > "$GFF_LIST"
    fi

    log "Using local protein FASTA:"
    log "  $PROTEIN_FASTA"
else
    NCBI_ZIP="$WORKDIR/ncbi_taxid_${TAXID}.zip"
    NCBI_DIR="$WORKDIR/ncbi_taxid_${TAXID}"

    if [[ "$FORCE" == "1" ]]; then
        rm -rf "$NCBI_ZIP" "$NCBI_DIR"
    fi

    if [[ ! -s "$NCBI_ZIP" ]]; then
        if [[ -n "$ASSEMBLY" ]]; then
            log "Downloading NCBI dataset by assembly accession: $ASSEMBLY"
            datasets download genome accession "$ASSEMBLY" \
                --include protein,gff3,gtf,seq-report \
                --filename "$NCBI_ZIP"
        else
            log "Downloading NCBI dataset by taxid: $TAXID"
            log "This uses --tax-exact-match and --annotated to reduce accidental broad downloads."

            datasets download genome taxon "$TAXID" \
                --annotated \
                --tax-exact-match \
                --exclude-atypical \
                --include protein,gff3,gtf,seq-report \
                --filename "$NCBI_ZIP"
        fi
    else
        log "Found existing NCBI zip, skipping download:"
        log "  $NCBI_ZIP"
    fi

    require_file "$NCBI_ZIP"

    if [[ ! -d "$NCBI_DIR" ]]; then
        mkdir -p "$NCBI_DIR"
        unzip -q "$NCBI_ZIP" -d "$NCBI_DIR"
    else
        log "Found existing extracted NCBI directory, skipping unzip:"
        log "  $NCBI_DIR"
    fi

    find "$NCBI_DIR" -type f \( \
        -name "*.faa" -o \
        -name "*.faa.gz" -o \
        -name "*protein*.fa" -o \
        -name "*protein*.faa" -o \
        -name "*protein*.fasta" \
    \) | sort > "$PROTEIN_LIST"

    find "$NCBI_DIR" -type f \( \
        -name "*.gff" -o \
        -name "*.gff3" -o \
        -name "*.gff.gz" -o \
        -name "*.gff3.gz" -o \
        -name "*.gtf" -o \
        -name "*.gtf.gz" \
    \) | sort > "$GFF_LIST"
fi

if [[ ! -s "$PROTEIN_LIST" ]]; then
    die "No protein FASTA candidates found. Cannot run eggNOG."
fi

log "Protein FASTA candidates:"
cat "$PROTEIN_LIST"

if [[ -s "$GFF_LIST" ]]; then
    log "Annotation candidates:"
    cat "$GFF_LIST"
else
    log "No GFF/GTF annotation candidates found. Will try FASTA header mapping only."
fi

# ------------------------------------------------------------
# Build representative FASTA with headers renamed to host.gct gene IDs
# ------------------------------------------------------------

REP_FASTA="$WORKDIR/taxid_${TAXID}_representative_proteins_renamed_to_host_gene_ids.fa"
COMPAT_REPORT="$WORKDIR/taxid_${TAXID}_compatibility_report.tsv"

if [[ "$FORCE" == "1" ]]; then
    rm -f "$REP_FASTA" "$COMPAT_REPORT"
fi

if [[ -s "$REP_FASTA" && -s "$COMPAT_REPORT" ]]; then
    log "Representative compatible FASTA already exists, skipping compatibility mapping:"
    log "  $REP_FASTA"
else
    log "Testing compatibility and creating representative protein FASTA..."

    python3 - "$HOST_GCT" "$PROTEIN_LIST" "$GFF_LIST" "$REP_FASTA" "$COMPAT_REPORT" "$MIN_OVERLAP_GENES" "$MIN_OVERLAP_PCT" <<'PY'
import sys
import os
import re
import gzip
from collections import defaultdict

host_gct = sys.argv[1]
protein_list_file = sys.argv[2]
gff_list_file = sys.argv[3]
out_fasta = sys.argv[4]
compat_report = sys.argv[5]
min_overlap_genes = int(float(sys.argv[6]))
min_overlap_pct = float(sys.argv[7])

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def read_host_genes(gct):
    genes = []
    with open_maybe_gzip(gct) as f:
        for idx, line in enumerate(f, start=1):
            if idx <= 3:
                continue
            line = line.rstrip("\n\r")
            if not line:
                continue
            genes.append(line.split("\t")[0])
    return set(g for g in genes if g)

host_genes = read_host_genes(host_gct)

if not host_genes:
    raise SystemExit("[ERROR] No genes extracted from host.gct")

def tokenize(text):
    text = re.sub(r"[^A-Za-z0-9_.-]+", " ", str(text))
    toks = [x.strip() for x in text.split() if x.strip()]
    return toks

def host_candidates_from_text(text):
    out = []

    for tok in tokenize(text):
        candidates = [tok]

        for prefix in ("gene-", "rna-", "cds-", "protein-", "transcript-", "id-", "ID-", "Name-"):
            if tok.startswith(prefix):
                candidates.append(tok[len(prefix):])

        if "." in tok:
            candidates.append(tok.split(".")[0])
            candidates.append(".".join(tok.split(".")[:-1]))

        if "-" in tok:
            candidates.append(tok.split("-")[-1])

        for c in candidates:
            if c in host_genes and c not in out:
                out.append(c)

    return out

def parse_attrs(attr):
    d = defaultdict(list)

    for part in str(attr).split(";"):
        part = part.strip()
        if not part:
            continue

        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
            v = v.strip().strip('"')
        else:
            continue

        k = k.strip()
        vals = [x.strip() for x in re.split(r"[,]", v) if x.strip()]
        d[k].extend(vals)

    return d

def clean_node_id(x):
    x = str(x).strip()
    for prefix in ("gene-", "rna-", "cds-", "transcript-", "protein-"):
        if x.startswith(prefix):
            return x[len(prefix):]
    return x

def load_gff_protein_to_host(gff_paths):
    protein_to_host = {}
    node_parent = {}
    gene_node_to_host = {}

    raw_lines = []

    for path in gff_paths:
        if not path or not os.path.exists(path):
            continue

        try:
            with open_maybe_gzip(path) as f:
                for line in f:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n\r").split("\t")
                    if len(parts) < 9:
                        continue
                    raw_lines.append(parts)
        except Exception:
            continue

    # First pass: genes and parent graph
    for parts in raw_lines:
        feature = parts[2]
        attr = parts[8]
        attrs = parse_attrs(attr)

        ids = attrs.get("ID", [])
        parents = attrs.get("Parent", [])

        for i in ids:
            if parents:
                node_parent[i] = parents[0]
                node_parent[clean_node_id(i)] = clean_node_id(parents[0])

        if feature.lower() == "gene":
            cands = host_candidates_from_text(attr)
            if cands:
                for i in ids:
                    gene_node_to_host[i] = cands[0]
                    gene_node_to_host[clean_node_id(i)] = cands[0]

    def resolve_host_from_parent(node):
        seen = set()
        cur = node

        for _ in range(20):
            if not cur or cur in seen:
                return None
            seen.add(cur)

            if cur in gene_node_to_host:
                return gene_node_to_host[cur]

            cur_clean = clean_node_id(cur)
            if cur_clean in gene_node_to_host:
                return gene_node_to_host[cur_clean]

            cur = node_parent.get(cur) or node_parent.get(cur_clean)

        return None

    # Second pass: CDS/protein features
    for parts in raw_lines:
        feature = parts[2].lower()
        attr = parts[8]
        attrs = parse_attrs(attr)

        if feature not in ("cds", "protein", "polypeptide"):
            continue

        host = None
        cands = host_candidates_from_text(attr)

        if cands:
            host = cands[0]
        else:
            for p in attrs.get("Parent", []):
                host = resolve_host_from_parent(p)
                if host:
                    break

        if not host:
            continue

        protein_ids = []

        for key in ("protein_id", "Name", "ID"):
            for val in attrs.get(key, []):
                protein_ids.append(val)
                protein_ids.append(clean_node_id(val))

        for val in attrs.get("Dbxref", []):
            if ":" in val:
                db, acc = val.split(":", 1)
                if db.lower() in ("genbank", "refseq", "protein_id", "uniprotkb/swiss-prot", "uniprotkb/trembl"):
                    protein_ids.append(acc)

        for pid in protein_ids:
            pid = pid.strip()
            if pid:
                protein_to_host[pid] = host
                if "." in pid:
                    protein_to_host[pid.split(".")[0]] = host

    return protein_to_host

with open(protein_list_file) as f:
    protein_paths = [x.strip() for x in f if x.strip()]

with open(gff_list_file) as f:
    gff_paths = [x.strip() for x in f if x.strip()]

gff_map = load_gff_protein_to_host(gff_paths)

def parse_fasta(path):
    records = []
    header = None
    seq = []

    with open_maybe_gzip(path) as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())

    if header is not None:
        records.append((header, "".join(seq)))

    return records

def protein_id_from_header(header):
    return header.split()[0].strip()

def map_header_to_host_gene(header, protein_id):
    cands = host_candidates_from_text(header)

    if cands:
        return cands[0], "fasta_header"

    if protein_id in gff_map:
        return gff_map[protein_id], "gff_protein_id"

    if "." in protein_id and protein_id.split(".")[0] in gff_map:
        return gff_map[protein_id.split(".")[0]], "gff_protein_id_no_version"

    return None, "unmapped"

best = None
report_rows = []

for path in protein_paths:
    records = parse_fasta(path)

    gene_to_records = defaultdict(list)
    method_count = defaultdict(int)
    total_records = 0
    mapped_records = 0

    for header, seq in records:
        total_records += 1
        pid = protein_id_from_header(header)
        gene, method = map_header_to_host_gene(header, pid)

        method_count[method] += 1

        if gene and gene in host_genes:
            mapped_records += 1
            gene_to_records[gene].append((pid, header, seq, len(seq), method))

    mapped_genes = len(gene_to_records)
    pct = (mapped_genes / len(host_genes)) * 100 if host_genes else 0

    report_rows.append([
        path,
        str(total_records),
        str(mapped_records),
        str(mapped_genes),
        "%.4f" % pct,
        ";".join([f"{k}:{v}" for k, v in sorted(method_count.items())])
    ])

    if best is None or mapped_genes > best["mapped_genes"]:
        best = {
            "path": path,
            "gene_to_records": gene_to_records,
            "mapped_genes": mapped_genes,
            "pct": pct,
            "total_records": total_records,
            "mapped_records": mapped_records
        }

with open(compat_report, "w") as out:
    out.write("protein_fasta\ttotal_records\tmapped_records\tmapped_host_genes\tmapped_host_gene_pct\tmapping_methods\n")
    for row in report_rows:
        out.write("\t".join(row) + "\n")

if best is None:
    raise SystemExit("[ERROR] No protein FASTA could be parsed.")

print("[INFO] Best protein FASTA:", best["path"])
print("[INFO] Total records:", best["total_records"])
print("[INFO] Mapped records:", best["mapped_records"])
print("[INFO] Mapped host genes:", best["mapped_genes"])
print("[INFO] Mapped host gene pct: %.4f" % best["pct"])
print("[INFO] Compatibility report:", compat_report)

if best["mapped_genes"] < min_overlap_genes or best["pct"] < min_overlap_pct:
    print()
    print("[ERROR] Downloaded/provided files are not compatible enough with host.gct.")
    print("[ERROR] host.gct genes:", len(host_genes))
    print("[ERROR] mapped host genes:", best["mapped_genes"])
    print("[ERROR] mapped host gene pct: %.4f" % best["pct"])
    print("[ERROR] minimum required genes:", min_overlap_genes)
    print("[ERROR] minimum required pct:", min_overlap_pct)
    print()
    print("Possible causes:")
    print("  1. NCBI downloaded a different assembly than the one used for MTD host reference.")
    print("  2. host.gct uses Ensembl IDs but NCBI files use RefSeq/XP_/LOC IDs.")
    print("  3. The downloaded annotation lacks protein-to-gene IDs compatible with host.gct.")
    print()
    print("Try:")
    print("  --assembly GCA_or_GCF_ACCESSION")
    print("or:")
    print("  --protein-fasta /path/to/exact/proteome.fa.gz --annotation-gff /path/to/exact/annotation.gff3")
    raise SystemExit(1)

with open(out_fasta, "w") as out:
    for gene in sorted(best["gene_to_records"]):
        recs = best["gene_to_records"][gene]
        best_rec = sorted(recs, key=lambda x: x[3], reverse=True)[0]
        pid, header, seq, length, method = best_rec

        out.write(f">{gene}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

print("[OK] Representative FASTA written:", out_fasta)
print("[OK] Representative genes:", len(best["gene_to_records"]))
PY
fi

require_file "$REP_FASTA"

log "Compatible representative proteins:"
grep -c "^>" "$REP_FASTA" || true

log "Compatibility report:"
column -s $'\t' -t "$COMPAT_REPORT" || cat "$COMPAT_REPORT"

# ------------------------------------------------------------
# Run eggNOG-mapper
# ------------------------------------------------------------

EGGNOG_OUTDIR="$WORKDIR/eggnog_results"
EGGNOG_PREFIX="taxid_${TAXID}_custom_ssgsea"
EGGNOG_ANNOT="$EGGNOG_OUTDIR/${EGGNOG_PREFIX}.emapper.annotations"

mkdir -p "$EGGNOG_OUTDIR"

if [[ "$FORCE" == "1" ]]; then
    rm -f "$EGGNOG_OUTDIR/${EGGNOG_PREFIX}".emapper.*
fi

if [[ -s "$EGGNOG_ANNOT" ]]; then
    log "eggNOG annotation already exists, skipping emapper:"
    log "  $EGGNOG_ANNOT"
else
    log "Running eggNOG-mapper..."

    emapper.py \
        -i "$REP_FASTA" \
        --itype proteins \
        -m diamond \
        --data_dir "$EGGNOG_DB" \
        --cpu "$THREADS" \
        --go_evidence all \
        --override \
        -o "$EGGNOG_PREFIX" \
        --output_dir "$EGGNOG_OUTDIR"
fi

require_file "$EGGNOG_ANNOT"

# ------------------------------------------------------------
# Convert eggNOG GO annotations to GMT
# ------------------------------------------------------------

OUT_GMT="$OUTDIR/custom_taxid_${TAXID}_eggNOG_GO.gmt"
OUT_GENE_TABLE="$OUTDIR/custom_taxid_${TAXID}_eggNOG_GO_gene_table.tsv"
OUT_SUMMARY="$OUTDIR/custom_taxid_${TAXID}_eggNOG_GO_summary.tsv"
OUT_OVERLAP="$OUTDIR/custom_taxid_${TAXID}_eggNOG_GO_overlap_report.txt"

if [[ "$FORCE" == "1" ]]; then
    rm -f "$OUT_GMT" "$OUT_GENE_TABLE" "$OUT_SUMMARY" "$OUT_OVERLAP"
fi

log "Converting eggNOG GO annotations to GMT..."

Rscript - "$HOST_GCT" "$EGGNOG_ANNOT" "$OUT_GMT" "$OUT_GENE_TABLE" "$OUT_SUMMARY" "$OUT_OVERLAP" "$MIN_OVERLAP_GENES" "$MIN_OVERLAP_PCT" <<'RS'
args <- commandArgs(trailingOnly = TRUE)

host_gct <- args[1]
eggnog_file <- args[2]
out_gmt <- args[3]
out_gene_table <- args[4]
out_summary <- args[5]
out_overlap <- args[6]
min_overlap_genes <- as.numeric(args[7])
min_overlap_pct <- as.numeric(args[8])

min_genes <- 5
max_genes <- 500

cat("[INFO] host.gct:", host_gct, "\n")
cat("[INFO] eggNOG annotations:", eggnog_file, "\n")
cat("[INFO] Output GMT:", out_gmt, "\n")

host_df <- read.delim(
  host_gct,
  skip = 2,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

host_genes <- unique(host_df[[1]])
host_genes <- host_genes[!is.na(host_genes) & host_genes != ""]

cat("[INFO] Genes in host.gct:", length(host_genes), "\n")

lines <- readLines(eggnog_file, warn = FALSE)
header_idx <- grep("^#query", lines)

if (length(header_idx) == 0) {
  stop("[ERROR] Could not find eggNOG header line starting with #query.")
}

header <- strsplit(sub("^#", "", lines[header_idx[1]]), "\t", fixed = TRUE)[[1]]
data_lines <- lines[!grepl("^#", lines) & nchar(lines) > 0]

if (length(data_lines) == 0) {
  stop("[ERROR] No data lines found in eggNOG annotation file.")
}

eggnog <- read.delim(
  text = paste(data_lines, collapse = "\n"),
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(eggnog) != length(header)) {
  stop(
    "[ERROR] eggNOG data/header column mismatch. Header: ",
    length(header), " Data: ", ncol(eggnog)
  )
}

colnames(eggnog) <- header

query_col <- "query"

if (!query_col %in% colnames(eggnog)) {
  query_col <- colnames(eggnog)[1]
}

go_candidates <- c("GOs", "GO_terms", "go_terms", "GO")
go_col <- go_candidates[go_candidates %in% colnames(eggnog)]

if (length(go_col) == 0) {
  stop("[ERROR] Could not find GO column. Expected one of: ", paste(go_candidates, collapse = ", "))
}

go_col <- go_col[1]

cat("[INFO] eggNOG annotation rows:", nrow(eggnog), "\n")
cat("[INFO] Using query column:", query_col, "\n")
cat("[INFO] Using GO column:", go_col, "\n")

query_ids <- eggnog[[query_col]]
go_raw <- eggnog[[go_col]]

valid <- query_ids %in% host_genes &
  !is.na(go_raw) &
  go_raw != "" &
  go_raw != "-"

eggnog_go <- data.frame(
  gene_id = query_ids[valid],
  GO_raw = go_raw[valid],
  stringsAsFactors = FALSE
)

cat("[INFO] Rows with host gene ID and GO:", nrow(eggnog_go), "\n")
cat("[INFO] Unique host genes with GO before GMT filtering:", length(unique(eggnog_go$gene_id)), "\n")

go_gene_pairs <- list()
counter <- 1

for (i in seq_len(nrow(eggnog_go))) {
  gene <- eggnog_go$gene_id[i]
  terms <- unlist(strsplit(eggnog_go$GO_raw[i], ",", fixed = TRUE))
  terms <- trimws(terms)
  terms <- terms[grepl("^GO:[0-9]+$", terms)]

  if (length(terms) == 0) {
    next
  }

  for (go in terms) {
    go_gene_pairs[[counter]] <- c(go_id = go, gene_id = gene)
    counter <- counter + 1
  }
}

if (length(go_gene_pairs) == 0) {
  stop("[ERROR] No valid GO-gene pairs found.")
}

go_long <- as.data.frame(do.call(rbind, go_gene_pairs), stringsAsFactors = FALSE)
go_long <- unique(go_long)

cat("[INFO] GO-gene pairs:", nrow(go_long), "\n")
cat("[INFO] Unique host.gct genes with GO:", length(unique(go_long$gene_id)), "\n")
cat("[INFO] Unique GO terms before filtering:", length(unique(go_long$go_id)), "\n")

write.table(
  go_long,
  out_gene_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

go_terms <- sort(unique(go_long$go_id))

gmt_list <- list()
summary_rows <- list()
idx <- 1

for (go in go_terms) {
  genes <- sort(unique(go_long$gene_id[go_long$go_id == go]))
  n <- length(genes)

  if (n >= min_genes && n <= max_genes) {
    gmt_list[[idx]] <- list(
      go_id = go,
      description = "eggNOG_GO",
      genes = genes,
      n_genes = n
    )

    summary_rows[[idx]] <- data.frame(
      go_id = go,
      description = "eggNOG_GO",
      n_genes = n,
      stringsAsFactors = FALSE
    )

    idx <- idx + 1
  }
}

if (length(gmt_list) == 0) {
  stop("[ERROR] No GO terms passed the gene set size filters.")
}

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$n_genes, decreasing = TRUE), ]

write.table(
  summary_df,
  out_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

con <- file(out_gmt, open = "wt")

for (item in gmt_list) {
  line <- c(item$go_id, item$description, item$genes)
  writeLines(paste(line, collapse = "\t"), con)
}

close(con)

gmt_genes <- sort(unique(unlist(lapply(gmt_list, function(x) x$genes))))
overlap <- intersect(sort(unique(host_genes)), gmt_genes)
overlap_pct <- length(overlap) / length(unique(host_genes)) * 100

writeLines(
  c(
    paste0("host_gct_genes\t", length(unique(host_genes))),
    paste0("gmt_genes\t", length(gmt_genes)),
    paste0("overlap_genes\t", length(overlap)),
    paste0("overlap_pct\t", sprintf("%.4f", overlap_pct)),
    paste0("gmt_terms\t", length(gmt_list)),
    paste0("min_required_overlap_genes\t", min_overlap_genes),
    paste0("min_required_overlap_pct\t", min_overlap_pct)
  ),
  out_overlap
)

cat("[OK] GMT written:", out_gmt, "\n")
cat("[OK] GO-gene table:", out_gene_table, "\n")
cat("[OK] Summary:", out_summary, "\n")
cat("[OK] Overlap report:", out_overlap, "\n")
cat("[OK] Overlap genes:", length(overlap), "\n")
cat("[OK] Overlap pct:", sprintf("%.4f", overlap_pct), "\n")
cat("[OK] GMT terms:", length(gmt_list), "\n")

if (length(overlap) < min_overlap_genes || overlap_pct < min_overlap_pct) {
  stop("[ERROR] Final GMT overlap with host.gct is below threshold.")
}
RS

require_file "$OUT_GMT"
require_file "$OUT_OVERLAP"

log "Custom GMT generated:"
echo "  $OUT_GMT"

log "Overlap report:"
cat "$OUT_OVERLAP"

log "Done."
