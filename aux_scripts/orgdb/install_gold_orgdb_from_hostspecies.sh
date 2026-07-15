#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# install_gold_orgdb_from_hostspecies.sh
# v6: robust protein FASTA parsing + clean previous OrgDb artifacts before rebuild.
#
# Build a reference-matched OrgDb from the same GTF/GFF and
# protein FASTA used by MTD. This avoids NCBI/Entrez ID drift.
#
# Main modes:
#   1) --eggnog FILE
#   2) --protein-fasta FILE  -> representative FASTA -> eggNOG -> OrgDb
#
# The eggNOG DB is searched by default at $MTDIR/eggnog_db.
# If missing, the script tries to download it there.
# ============================================================

TAXID=""
HOSTSPECIES=""
GTF=""
GENOME=""
PROTEIN_FASTA=""
EGGNOG=""
EGGNOG_DB=""
LIB=""
BUILD_DIR=""
VERSION="0.3.0"
THREADS="20"
FORCE=0
SKIP_IF_COMPATIBLE=0
SKIP_ORGDB_BUILD=0
CONDA_R_ENV="${MTD_ORGDB_ENV:-MTD_orgdb}"
CONDA_MTD_ENV="MTD"
SYMBOL_MODE="gene_id"
GENE_PATTERN='(?:gene[:=]|gene_id[:=]|gene=)([A-Za-z0-9_.:-]+)'
ID_AS_GENE=0
MIN_PROTEIN_GENES=100
MIN_PROTEIN_GTF_PCT="5"
MIN_EGGNOG_GENES=100
MIN_EGGNOG_GTF_PCT="1"
MIN_EXISTING_ORGDB_PCT="50"
MIN_GENOME_SEQ_OVERLAP_PCT="1"

# Resolve MTDIR from this script location when possible.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Expected location: $MTDIR/aux_scripts/orgdb
MTDIR_DEFAULT="$(cd "$SCRIPT_DIR/../.." 2>/dev/null && pwd || true)"
if [[ ! -d "$MTDIR_DEFAULT" || "$(basename "$MTDIR_DEFAULT")" == "aux_scripts" ]]; then
  MTDIR_DEFAULT="/home/me/MTD"
fi
MTDIR="${MTDIR:-$MTDIR_DEFAULT}"

HOSTSPECIES_DEFAULT="$MTDIR/HostSpecies.csv"

OFFLINE_CACHE_PATH_FILE="$MTDIR/offlineCachePath"
OFFLINE_CACHE_ROOT=""

if [[ -s "$OFFLINE_CACHE_PATH_FILE" ]]; then
    IFS= read -r OFFLINE_CACHE_ROOT < "$OFFLINE_CACHE_PATH_FILE"
    OFFLINE_CACHE_ROOT="${OFFLINE_CACHE_ROOT%$'\r'}"
fi

if [[ -n "$OFFLINE_CACHE_ROOT" ]]; then
    EGGNOG_DB_DEFAULT="$OFFLINE_CACHE_ROOT/eggNOG/emapperdb-5.0.2"
else
    # Compatibility fallback for installations created before offlineCachePath.
    EGGNOG_DB_DEFAULT="$MTDIR/eggnog_db"
fi

LIB_DEFAULT="$MTDIR/custom_R_libs"
BUILD_DIR_DEFAULT="$MTDIR/build/orgdb_gold"

show_help() {
cat <<EOF
Usage:
  install_gold_orgdb_from_hostspecies.sh \\
    --taxid TAXID \\
    --hostspecies /path/to/HostSpecies.csv \\
    --gtf /path/to/annotation.gtf.gz \\
    --skip-orgdb-build \\
    [--genome /path/to/genome.fa] \\
    [--protein-fasta /path/to/proteins.fa.gz | --eggnog /path/to/file.emapper.annotations] \\
    [--lib $LIB_DEFAULT] \\
    [--build-dir $BUILD_DIR_DEFAULT] \\
    [--threads 20] \\
    [--version 0.3.0] \\
    [--skip-if-compatible] \\
    [--force]

Modes:
  1) --eggnog FILE
     Build OrgDb from existing eggNOG annotations. The eggNOG query IDs must be gene IDs.

  2) --protein-fasta FILE
     Create representative-protein-per-gene FASTA, run eggNOG in the MTD env,
     then build OrgDb in the R412 env.

Defaults:
  --hostspecies  $HOSTSPECIES_DEFAULT
  --eggnog-db    $EGGNOG_DB_DEFAULT
  --lib          $LIB_DEFAULT
  --build-dir    $BUILD_DIR_DEFAULT/TAXID

Safety checks:
  --genome checks FASTA sequence names against GTF seqnames.
  --protein-fasta is checked against GTF gene IDs before running eggNOG.
  eggNOG annotations are checked against GTF gene IDs before OrgDb creation.
  --skip-if-compatible tests an existing OrgDb from HostSpecies.csv and skips creation if compatible.

Notes:
  - This script avoids AnnotationForge::makeOrgPackageFromNCBI.
  - It builds the OrgDb using the IDs from the actual GTF/proteome reference.
  - SYMBOL defaults to gene_id so DEG_Anno_Plot.R with keyType='SYMBOL' works for Ensembl-like IDs.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxid) TAXID="$2"; shift 2 ;;
    --taxid=*) TAXID="${1#*=}"; shift ;;
    --hostspecies) HOSTSPECIES="$2"; shift 2 ;;
    --hostspecies=*) HOSTSPECIES="${1#*=}"; shift ;;
    --gtf|--gff) GTF="$2"; shift 2 ;;
    --gtf=*|--gff=*) GTF="${1#*=}"; shift ;;
    --genome) GENOME="$2"; shift 2 ;;
    --genome=*) GENOME="${1#*=}"; shift ;;
    --protein-fasta) PROTEIN_FASTA="$2"; shift 2 ;;
    --protein-fasta=*) PROTEIN_FASTA="${1#*=}"; shift ;;
    --eggnog) EGGNOG="$2"; shift 2 ;;
    --eggnog=*) EGGNOG="${1#*=}"; shift ;;
    --eggnog-db) EGGNOG_DB="$2"; shift 2 ;;
    --eggnog-db=*) EGGNOG_DB="${1#*=}"; shift ;;
    --lib) LIB="$2"; shift 2 ;;
    --lib=*) LIB="${1#*=}"; shift ;;
    --build-dir) BUILD_DIR="$2"; shift 2 ;;
    --build-dir=*) BUILD_DIR="${1#*=}"; shift ;;
    --version) VERSION="$2"; shift 2 ;;
    --version=*) VERSION="${1#*=}"; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --threads=*) THREADS="${1#*=}"; shift ;;
    --conda-r-env) CONDA_R_ENV="$2"; shift 2 ;;
    --conda-r-env=*) CONDA_R_ENV="${1#*=}"; shift ;;
    --conda-mtd-env) CONDA_MTD_ENV="$2"; shift 2 ;;
    --conda-mtd-env=*) CONDA_MTD_ENV="${1#*=}"; shift ;;
    --symbol-mode) SYMBOL_MODE="$2"; shift 2 ;;
    --symbol-mode=*) SYMBOL_MODE="${1#*=}"; shift ;;
    --gene-pattern) GENE_PATTERN="$2"; shift 2 ;;
    --gene-pattern=*) GENE_PATTERN="${1#*=}"; shift ;;
    --id-as-gene) ID_AS_GENE=1; shift ;;
    --skip-if-compatible) SKIP_IF_COMPATIBLE=1; shift ;;
    --min-existing-orgdb-pct) MIN_EXISTING_ORGDB_PCT="$2"; shift 2 ;;
    --min-existing-orgdb-pct=*) MIN_EXISTING_ORGDB_PCT="${1#*=}"; shift ;;
    --min-protein-genes) MIN_PROTEIN_GENES="$2"; shift 2 ;;
    --min-protein-genes=*) MIN_PROTEIN_GENES="${1#*=}"; shift ;;
    --min-protein-gtf-pct) MIN_PROTEIN_GTF_PCT="$2"; shift 2 ;;
    --min-protein-gtf-pct=*) MIN_PROTEIN_GTF_PCT="${1#*=}"; shift ;;
    --min-eggnog-genes) MIN_EGGNOG_GENES="$2"; shift 2 ;;
    --min-eggnog-genes=*) MIN_EGGNOG_GENES="${1#*=}"; shift ;;
    --min-eggnog-gtf-pct) MIN_EGGNOG_GTF_PCT="$2"; shift 2 ;;
    --min-eggnog-gtf-pct=*) MIN_EGGNOG_GTF_PCT="${1#*=}"; shift ;;
    --skip-orgdb-build) SKIP_ORGDB_BUILD=1 shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1"; show_help; exit 1 ;;
  esac
done

HOSTSPECIES="${HOSTSPECIES:-$HOSTSPECIES_DEFAULT}"
EGGNOG_DB="${EGGNOG_DB:-$EGGNOG_DB_DEFAULT}"
LIB="${LIB:-$LIB_DEFAULT}"
if [[ -z "$BUILD_DIR" ]]; then
  if [[ -n "$TAXID" ]]; then BUILD_DIR="$BUILD_DIR_DEFAULT/$TAXID"; else BUILD_DIR="$BUILD_DIR_DEFAULT"; fi
fi

if [[ -z "$TAXID" ]]; then echo "[ERROR] Missing --taxid"; show_help; exit 1; fi
if [[ ! -s "$HOSTSPECIES" ]]; then echo "[ERROR] HostSpecies.csv not found: $HOSTSPECIES"; exit 1; fi
if [[ ! -s "$GTF" ]]; then echo "[ERROR] GTF/GFF not found: $GTF"; exit 1; fi
if [[ -n "$GENOME" && ! -s "$GENOME" ]]; then echo "[ERROR] Genome FASTA not found: $GENOME"; exit 1; fi
if [[ -z "$EGGNOG" && -z "$PROTEIN_FASTA" ]]; then echo "[ERROR] Provide --eggnog or --protein-fasta"; show_help; exit 1; fi
if [[ -n "$PROTEIN_FASTA" && ! -s "$PROTEIN_FASTA" ]]; then echo "[ERROR] Protein FASTA not found: $PROTEIN_FASTA"; exit 1; fi
if [[ -n "$EGGNOG" && ! -s "$EGGNOG" ]]; then echo "[ERROR] eggNOG annotations not found: $EGGNOG"; exit 1; fi

PY_REP="$SCRIPT_DIR/make_gene_representative_fasta.py"
R_BUILD="$SCRIPT_DIR/build_gold_orgdb_from_gtf_eggnog.R"

if [[ ! -s "$PY_REP" ]]; then echo "[ERROR] Missing helper: $PY_REP"; exit 1; fi
if [[ ! -s "$R_BUILD" ]]; then echo "[ERROR] Missing helper: $R_BUILD"; exit 1; fi

# Ensure conda is usable even when launched from a non-interactive script.
if ! command -v conda >/dev/null 2>&1; then
  if [[ -s "$MTDIR/condaPath" ]]; then
    CONDA_BASE="$(head -n 1 "$MTDIR/condaPath")"
    # shellcheck source=/dev/null
    source "$CONDA_BASE/etc/profile.d/conda.sh"
  elif [[ -s "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  fi
fi
if ! command -v conda >/dev/null 2>&1; then echo "[ERROR] conda not found."; exit 1; fi

read -r GENUS SPECIES SCI_NAME < <(python3 - "$HOSTSPECIES" "$TAXID" <<'PY'
import csv, sys
path, taxid = sys.argv[1], sys.argv[2]
with open(path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if str(row.get('Taxon_ID','')).strip() == str(taxid):
            sci = row.get('Scientific_name','').strip()
            parts = sci.split()
            if len(parts) < 2:
                raise SystemExit(f'[ERROR] Scientific_name is not Genus species: {sci}')
            print(parts[0], parts[1], sci)
            raise SystemExit(0)
raise SystemExit(f'[ERROR] TaxID {taxid} not found in HostSpecies.csv')
PY
)

mkdir -p "$LIB" "$BUILD_DIR"
WORK_DIR="$BUILD_DIR/work_taxid_${TAXID}"
if [[ "$FORCE" == "1" ]]; then rm -rf "$WORK_DIR"; fi
mkdir -p "$WORK_DIR"

# OrgDb package name expected by AnnotationForge and/or listed in HostSpecies.csv.
ORGDB_FROM_CSV="$(python3 - "$HOSTSPECIES" "$TAXID" <<'PYORGDB'
import csv, sys
path, taxid = sys.argv[1], sys.argv[2]
with open(path, newline='') as f:
    for row in csv.DictReader(f):
        if str(row.get('Taxon_ID','')).strip() == str(taxid):
            print((row.get('OrgDb') or '').strip())
            raise SystemExit(0)
PYORGDB
)"
ORGDB_GUESS="org.${GENUS:0:1}${SPECIES}.eg.db"

remove_path_if_exists() {
  local p="$1"
  if [[ -e "$p" || -L "$p" ]]; then
    echo "[CLEAN] Removing previous OrgDb artifact: $p"
    rm -rf "$p"
  fi
}

clean_previous_orgdb_artifacts() {
  local pkg
  local seen=""
  for pkg in "$ORGDB_FROM_CSV" "$ORGDB_GUESS"; do
    [[ -z "$pkg" || "$pkg" == "NA" || "$pkg" == "None" ]] && continue
    case ":$seen:" in *":$pkg:"*) continue ;; esac
    seen="$seen:$pkg"

    remove_path_if_exists "$LIB/$pkg"
    remove_path_if_exists "$LIB/00LOCK-$pkg"
    remove_path_if_exists "$BUILD_DIR/$pkg"
    remove_path_if_exists "$BUILD_DIR/00LOCK-$pkg"
  done
}

GTF_GENES="$WORK_DIR/gtf_gene_ids.txt"
GTF_SEQNAMES="$WORK_DIR/gtf_seqnames.txt"
GENOME_SEQNAMES="$WORK_DIR/genome_seqnames.txt"
INPUT_REPORT="$WORK_DIR/input_compatibility_report.tsv"

printf "[INFO] TaxID: %s\n" "$TAXID"
printf "[INFO] Scientific name: %s\n" "$SCI_NAME"
printf "[INFO] GTF/GFF: %s\n" "$GTF"
printf "[INFO] Genome: %s\n" "${GENOME:-not provided}"
printf "[INFO] Protein FASTA: %s\n" "${PROTEIN_FASTA:-not provided}"
printf "[INFO] eggNOG DB: %s\n" "$EGGNOG_DB"
printf "[INFO] Build dir: %s\n" "$BUILD_DIR"
printf "[INFO] Install lib: %s\n" "$LIB"

# ------------------------------------------------------------
# Input inspection: GTF genes and seqnames
# ------------------------------------------------------------
python3 - "$GTF" "$GTF_GENES" "$GTF_SEQNAMES" <<'PY'
import gzip, re, sys
path, gene_out, seq_out = sys.argv[1:4]

def open_text(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def attr_val(attr, keys):
    for key in keys:
        m = re.search(r'(?:^|;)\s*' + re.escape(key) + r'=([^;]+)', attr)
        if m:
            return re.sub(r'^(gene:|ID=|ID=gene:)', '', m.group(1).strip().strip('"'))
        m = re.search(re.escape(key) + r' "([^"]+)"', attr)
        if m:
            return re.sub(r'^(gene:|ID=|ID=gene:)', '', m.group(1).strip())
    return None

genes, seqs = set(), set()
with open_text(path) as fh:
    for line in fh:
        if not line.strip() or line.startswith('#'): continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9: continue
        seqs.add(parts[0])
        gid = attr_val(parts[8], ['gene_id','ID','locus_tag','Name'])
        if gid: genes.add(gid)
with open(gene_out, 'w') as out:
    for g in sorted(genes): out.write(g+'\n')
with open(seq_out, 'w') as out:
    for s in sorted(seqs): out.write(s+'\n')
print(f'[INFO] GTF unique gene IDs: {len(genes)}')
print(f'[INFO] GTF unique seqnames: {len(seqs)}')
if len(genes) == 0:
    raise SystemExit('[ERROR] No gene IDs found in GTF/GFF. Check annotation file and attributes.')
PY

GTF_GENE_COUNT="$(wc -l < "$GTF_GENES" | tr -d ' ')"
GTF_SEQ_COUNT="$(wc -l < "$GTF_SEQNAMES" | tr -d ' ')"

# ------------------------------------------------------------
# Optional genome-vs-GTF seqname compatibility check
# ------------------------------------------------------------
GENOME_OVERLAP="NA"
GENOME_OVERLAP_PCT="NA"
if [[ -n "$GENOME" ]]; then
  python3 - "$GENOME" "$GENOME_SEQNAMES" <<'PY'
import gzip, re, sys
path, outpath = sys.argv[1:3]
def open_text(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
seqs=set()
with open_text(path) as fh:
    for line in fh:
        if line.startswith('>'):
            h=line[1:].strip()
            h=re.sub(r'^kraken:taxid\|[0-9]+\|', '', h)
            h=re.sub(r'\s+\[organism=[^]]+\]\s*$', '', h)
            seqs.add(h.split()[0])
with open(outpath,'w') as out:
    for s in sorted(seqs): out.write(s+'\n')
print(f'[INFO] Genome FASTA headers: {len(seqs)}')
PY
  GENOME_SEQ_COUNT="$(wc -l < "$GENOME_SEQNAMES" | tr -d ' ')"
  GENOME_OVERLAP="$(comm -12 <(sort -u "$GTF_SEQNAMES") <(sort -u "$GENOME_SEQNAMES") | wc -l | tr -d ' ')"
  GENOME_OVERLAP_PCT="$(awk -v ov="$GENOME_OVERLAP" -v n="$GTF_SEQ_COUNT" 'BEGIN{if(n>0) printf "%.4f", 100*ov/n; else print 0}')"
  echo "[INFO] Genome/GTF seqname overlap: $GENOME_OVERLAP / $GTF_SEQ_COUNT (${GENOME_OVERLAP_PCT}%)"
  if awk -v pct="$GENOME_OVERLAP_PCT" -v min="$MIN_GENOME_SEQ_OVERLAP_PCT" 'BEGIN{exit !(pct < min)}'; then
    echo "[ERROR] Genome FASTA and GTF/GFF seqnames barely overlap."
    echo "[ERROR] This usually means genome and annotation are from different assemblies/releases."
    echo "[ERROR] Report: $INPUT_REPORT"
    exit 1
  fi
fi

printf "metric\tvalue\n" > "$INPUT_REPORT"
printf "taxid\t%s\n" "$TAXID" >> "$INPUT_REPORT"
printf "scientific_name\t%s\n" "$SCI_NAME" >> "$INPUT_REPORT"
printf "gtf_gene_ids\t%s\n" "$GTF_GENE_COUNT" >> "$INPUT_REPORT"
printf "gtf_seqnames\t%s\n" "$GTF_SEQ_COUNT" >> "$INPUT_REPORT"
printf "genome_gtf_seqname_overlap\t%s\n" "$GENOME_OVERLAP" >> "$INPUT_REPORT"
printf "genome_gtf_seqname_overlap_pct\t%s\n" "$GENOME_OVERLAP_PCT" >> "$INPUT_REPORT"

# ------------------------------------------------------------
# Existing OrgDb compatibility limiter
# ------------------------------------------------------------
if [[ "$SKIP_IF_COMPATIBLE" == "1" && "$FORCE" != "1" ]]; then
  if [[ -n "$ORGDB_FROM_CSV" && "$ORGDB_FROM_CSV" != "NA" && "$ORGDB_FROM_CSV" != "None" ]]; then
    echo "[INFO] Existing OrgDb listed in HostSpecies.csv: $ORGDB_FROM_CSV"
    ORGDB_COMPAT_REPORT="$WORK_DIR/existing_orgdb_compatibility.tsv"
    ORGDB_COMPAT_R="$WORK_DIR/check_existing_orgdb_compatibility.R"

    # IMPORTANT: do not use `conda run Rscript - <<RS` here.
    # Some conda versions do not reliably forward stdin to Rscript, which can leave
    # ORGDB_COMPAT_REPORT missing and crash the downstream awk parser.
    cat > "$ORGDB_COMPAT_R" <<'RS'
suppressPackageStartupMessages(library(AnnotationDbi))

orgdb_pkg <- Sys.getenv('ORGDB_PKG')
gtf_genes_file <- Sys.getenv('GTF_GENES')
out_report <- Sys.getenv('OUT_REPORT')
lib <- Sys.getenv('LIB')

if (nzchar(lib)) .libPaths(unique(c(lib, .libPaths())))

gtf_genes <- readLines(gtf_genes_file, warn = FALSE)

ok <- requireNamespace(orgdb_pkg, quietly = TRUE)

if (!ok) {
  res <- data.frame(
    orgdb = orgdb_pkg,
    status = 'not_installed',
    keytype = NA_character_,
    gtf_genes = length(gtf_genes),
    orgdb_keys = NA_integer_,
    overlap = 0L,
    overlap_pct = 0,
    stringsAsFactors = FALSE
  )
} else {
  suppressPackageStartupMessages(library(orgdb_pkg, character.only = TRUE))
  db <- get(orgdb_pkg)
  res <- do.call(rbind, lapply(keytypes(db), function(k) {
    ks <- tryCatch(keys(db, keytype = k), error = function(e) character())
    ov <- length(intersect(gtf_genes, ks))
    data.frame(
      orgdb = orgdb_pkg,
      status = 'installed',
      keytype = k,
      gtf_genes = length(gtf_genes),
      orgdb_keys = length(ks),
      overlap = ov,
      overlap_pct = round(100 * ov / length(gtf_genes), 4),
      stringsAsFactors = FALSE
    )
  }))
  res <- res[order(res$overlap, decreasing = TRUE), ]
}

write.table(res, out_report, sep = '\t', quote = FALSE, row.names = FALSE)
RS

    if ! ORGDB_PKG="$ORGDB_FROM_CSV" GTF_GENES="$GTF_GENES" OUT_REPORT="$ORGDB_COMPAT_REPORT" LIB="$LIB" \
      conda run -n "$CONDA_R_ENV" Rscript "$ORGDB_COMPAT_R"; then
      echo "[WARNING] Existing OrgDb compatibility check failed. Continuing with gold OrgDb creation."
      rm -f "$ORGDB_COMPAT_REPORT"
    fi

    if [[ -s "$ORGDB_COMPAT_REPORT" ]]; then
      BEST_PCT="$(awk 'NR==2 {print $7}' "$ORGDB_COMPAT_REPORT")"
      BEST_KEYTYPE="$(awk 'NR==2 {print $3}' "$ORGDB_COMPAT_REPORT")"
      BEST_OVERLAP="$(awk 'NR==2 {print $6}' "$ORGDB_COMPAT_REPORT")"
      STATUS="$(awk 'NR==2 {print $2}' "$ORGDB_COMPAT_REPORT")"
      BEST_PCT="${BEST_PCT:-0}"
      BEST_KEYTYPE="${BEST_KEYTYPE:-NA}"
      BEST_OVERLAP="${BEST_OVERLAP:-0}"
      STATUS="${STATUS:-unknown}"
      echo "[INFO] Existing OrgDb status: $STATUS"
      echo "[INFO] Best keytype: $BEST_KEYTYPE"
      echo "[INFO] Best overlap: $BEST_OVERLAP genes (${BEST_PCT}%)"
      echo "[INFO] Existing OrgDb compatibility report: $ORGDB_COMPAT_REPORT"
      if [[ "$STATUS" == "installed" ]] && \
       awk -v pct="$BEST_PCT" -v min="$MIN_EXISTING_ORGDB_PCT" \
         'BEGIN { exit !(pct >= min) }'
    then
    echo "[OK] Existing OrgDb is compatible enough."
    echo "[INFO] The R OrgDb package will not be rebuilt."
    echo "[INFO] Continuing to ensure that representative proteins,"
    echo "[INFO] eggNOG annotations and ssGSEA resources are available."

    SKIP_ORGDB_BUILD=1
fi
      echo "[WARNING] Existing OrgDb is missing or not compatible enough. Building gold OrgDb."
    else
      echo "[WARNING] Existing OrgDb compatibility report was not created. Building gold OrgDb anyway."
    fi
  fi
fi

# ------------------------------------------------------------
# Clean previous OrgDb package/source artifacts for this species
# ------------------------------------------------------------
# At this point, either --skip-if-compatible was not used, --force was used,
# or the existing OrgDb was missing/incompatible. Since we are building a gold
# reference-matched OrgDb, remove previous installs/source dirs for this package
# from the custom library/build folder to avoid stale or duplicated artifacts.
if [[ "$SKIP_ORGDB_BUILD" == "1" ]]; then
    echo "[INFO] Keeping the existing OrgDb package and build artifacts."
else
    clean_previous_orgdb_artifacts
fi

# ------------------------------------------------------------
# eggNOG DB check/download
# ------------------------------------------------------------
# Notes:
# - download_eggnog_data.py in some eggNOG-mapper releases still points to
#   the old eggnogdb.embl.de host. Here we use the current eggnog5.embl.de
#   endpoint directly and download with resume + retries.
# - We download .gz/.tar.gz into temporary compressed files, validate them,
#   then decompress/extract. This avoids leaving truncated final DB files.
# ------------------------------------------------------------
EGGNOG_BASE_URLS=(
  "https://eggnog5.embl.de/download/emapperdb-5.0.2"
  "http://eggnog5.embl.de/download/emapperdb-5.0.2"
)

validate_gzip_file() {
  local f="$1"
  [[ -s "$f" ]] || return 1
  gzip -t "$f" >/dev/null 2>&1
}

validate_tar_gz_file() {
  local f="$1"
  [[ -s "$f" ]] || return 1
  tar -tzf "$f" >/dev/null 2>&1
}

robust_fetch() {
  local url="$1"
  local out="$2"
  local tries="${3:-50}"
  local wait="${4:-20}"

  mkdir -p "$(dirname "$out")"

  if command -v aria2c >/dev/null 2>&1; then
    echo "[DOWNLOAD] aria2c resume/retry: $url"
    aria2c \
      --continue=true \
      --max-connection-per-server=8 \
      --split=8 \
      --min-split-size=10M \
      --retry-wait="$wait" \
      --max-tries="$tries" \
      --timeout=60 \
      --connect-timeout=30 \
      --lowest-speed-limit=20K \
      --allow-overwrite=true \
      --auto-file-renaming=false \
      --dir "$(dirname "$out")" \
      --out "$(basename "$out")" \
      "$url"
  elif command -v wget >/dev/null 2>&1; then
    echo "[DOWNLOAD] wget resume/retry: $url"
    wget \
      --continue \
      --tries="$tries" \
      --waitretry="$wait" \
      --read-timeout=60 \
      --timeout=60 \
      --retry-connrefused \
      --user-agent="Mozilla/5.0" \
      -O "$out" \
      "$url"
  else
    echo "[ERROR] Neither aria2c nor wget was found. Install one of them first."
    return 1
  fi
}

fetch_from_any_base() {
  local rel="$1"
  local out="$2"
  local base

  for base in "${EGGNOG_BASE_URLS[@]}"; do
    echo "[INFO] Trying: ${base}/${rel}"
    if robust_fetch "${base}/${rel}" "$out"; then
      return 0
    fi
    echo "[WARNING] Download failed from: ${base}/${rel}"
  done

  return 1
}

ensure_eggnog_db() {
  mkdir -p "$EGGNOG_DB"

  if [[ -s "$EGGNOG_DB/eggnog.db" && -s "$EGGNOG_DB/eggnog_proteins.dmnd" && -s "$EGGNOG_DB/eggnog.taxa.db" ]]; then
    echo "[INFO] eggNOG DB found at: $EGGNOG_DB"
    return 0
  fi

  echo "[WARNING] eggNOG DB not complete at: $EGGNOG_DB"
  echo "[INFO] Robust manual downloader will use eggnog5.embl.de with resume + retries."
  echo "[INFO] Target directory: $EGGNOG_DB"

  (
    set -euo pipefail
    cd "$EGGNOG_DB"

    # eggnog.db
    if [[ ! -s eggnog.db ]]; then
      echo "[INFO] Missing: eggnog.db"
      if [[ -s eggnog.db.gz ]] && ! validate_gzip_file eggnog.db.gz; then
        echo "[WARNING] Existing eggnog.db.gz is corrupt/truncated. Removing it."
        rm -f eggnog.db.gz
      fi
      fetch_from_any_base "eggnog.db.gz" "$EGGNOG_DB/eggnog.db.gz"
      validate_gzip_file eggnog.db.gz || { echo "[ERROR] Downloaded eggnog.db.gz failed gzip validation."; exit 1; }
      echo "[PROCESS] Decompressing eggnog.db.gz"
      gunzip -f eggnog.db.gz
    else
      echo "[OK] Existing: eggnog.db"
    fi

    # eggnog_proteins.dmnd
    if [[ ! -s eggnog_proteins.dmnd ]]; then
      echo "[INFO] Missing: eggnog_proteins.dmnd"
      if [[ -s eggnog_proteins.dmnd.gz ]] && ! validate_gzip_file eggnog_proteins.dmnd.gz; then
        echo "[WARNING] Existing eggnog_proteins.dmnd.gz is corrupt/truncated. Removing it."
        rm -f eggnog_proteins.dmnd.gz
      fi
      fetch_from_any_base "eggnog_proteins.dmnd.gz" "$EGGNOG_DB/eggnog_proteins.dmnd.gz"
      validate_gzip_file eggnog_proteins.dmnd.gz || { echo "[ERROR] Downloaded eggnog_proteins.dmnd.gz failed gzip validation."; exit 1; }
      echo "[PROCESS] Decompressing eggnog_proteins.dmnd.gz"
      gunzip -f eggnog_proteins.dmnd.gz
    else
      echo "[OK] Existing: eggnog_proteins.dmnd"
    fi

    # eggnog.taxa.db is inside eggnog.taxa.tar.gz
    if [[ ! -s eggnog.taxa.db ]]; then
      echo "[INFO] Missing: eggnog.taxa.db"
      if [[ -s eggnog.taxa.tar.gz ]] && ! validate_tar_gz_file eggnog.taxa.tar.gz; then
        echo "[WARNING] Existing eggnog.taxa.tar.gz is corrupt/truncated. Removing it."
        rm -f eggnog.taxa.tar.gz
      fi
      fetch_from_any_base "eggnog.taxa.tar.gz" "$EGGNOG_DB/eggnog.taxa.tar.gz"
      validate_tar_gz_file eggnog.taxa.tar.gz || { echo "[ERROR] Downloaded eggnog.taxa.tar.gz failed tar.gz validation."; exit 1; }
      echo "[PROCESS] Extracting eggnog.taxa.tar.gz"
      tar -xzf eggnog.taxa.tar.gz
      rm -f eggnog.taxa.tar.gz
    else
      echo "[OK] Existing: eggnog.taxa.db"
    fi
  )

  if [[ ! -s "$EGGNOG_DB/eggnog.db" || ! -s "$EGGNOG_DB/eggnog_proteins.dmnd" || ! -s "$EGGNOG_DB/eggnog.taxa.db" ]]; then
    echo "[ERROR] eggNOG DB is still incomplete after robust download attempt."
    echo "Expected:"
    echo "  $EGGNOG_DB/eggnog.db"
    echo "  $EGGNOG_DB/eggnog_proteins.dmnd"
    echo "  $EGGNOG_DB/eggnog.taxa.db"
    exit 1
  fi

  echo "[OK] eggNOG DB is ready: $EGGNOG_DB"
}

# ------------------------------------------------------------
# Protein FASTA mode: representative proteins -> eggNOG
# ------------------------------------------------------------
if [[ -n "$PROTEIN_FASTA" ]]; then
  ensure_eggnog_db
  echo "[INFO] Protein FASTA mode enabled. Creating representative FASTA and running eggNOG."
  REP_FASTA="$WORK_DIR/taxid_${TAXID}_representative_proteins_gene_ids.fa"
  REP_REPORT="$WORK_DIR/taxid_${TAXID}_representative_proteins_report.tsv"
  EGGNOG_OUTDIR="$WORK_DIR/eggnog_results"
  EGGNOG_PREFIX="taxid_${TAXID}_gold_orgdb"
  EGGNOG="$EGGNOG_OUTDIR/${EGGNOG_PREFIX}.emapper.annotations"
  mkdir -p "$EGGNOG_OUTDIR"

  ID_ARG=()
  if [[ "$ID_AS_GENE" == "1" ]]; then ID_ARG=(--id-as-gene); fi

  UNMAPPED_REPORT="$WORK_DIR/taxid_${TAXID}_representative_unmapped_headers.tsv"

  conda run -n "$CONDA_MTD_ENV" python3 "$PY_REP" \
    --protein-fasta "$PROTEIN_FASTA" \
    --out-fasta "$REP_FASTA" \
    --report "$REP_REPORT" \
    --gene-pattern "$GENE_PATTERN" \
    --gene-list "$GTF_GENES" \
    --gtf "$GTF" \
    --debug-unmapped "$UNMAPPED_REPORT" \
    "${ID_ARG[@]}"

  REP_GENE_COUNT="$(grep -c '^>' "$REP_FASTA" || true)"
  REP_GTF_PCT="$(awk -v n="$REP_GENE_COUNT" -v total="$GTF_GENE_COUNT" 'BEGIN{if(total>0) printf "%.4f", 100*n/total; else print 0}')"
  echo "[INFO] Representative proteins matching GTF genes: $REP_GENE_COUNT / $GTF_GENE_COUNT (${REP_GTF_PCT}%)"
  printf "representative_protein_genes\t%s\n" "$REP_GENE_COUNT" >> "$INPUT_REPORT"
  printf "representative_protein_gtf_pct\t%s\n" "$REP_GTF_PCT" >> "$INPUT_REPORT"

  if awk -v n="$REP_GENE_COUNT" -v min="$MIN_PROTEIN_GENES" 'BEGIN{exit !(n < min)}'; then
    echo "[ERROR] Too few proteins could be mapped to GTF gene IDs."
    echo "[ERROR] Check if --protein-fasta and --gtf are from the same annotation/release."
    echo "[ERROR] Representative report: $REP_REPORT"
    exit 1
  fi
  if awk -v pct="$REP_GTF_PCT" -v min="$MIN_PROTEIN_GTF_PCT" 'BEGIN{exit !(pct < min)}'; then
    echo "[ERROR] Protein/GTF compatibility is too low: ${REP_GTF_PCT}%"
    echo "[ERROR] Representative report: $REP_REPORT"
    exit 1
  fi

  if [[ "$FORCE" == "1" ]]; then rm -f "$EGGNOG_OUTDIR/${EGGNOG_PREFIX}.emapper."*; fi

  conda run -n "$CONDA_MTD_ENV" emapper.py \
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

if [[ ! -s "$EGGNOG" ]]; then
  echo "[ERROR] eggNOG annotations file was not created/found: $EGGNOG"
  exit 1
fi

# ------------------------------------------------------------
# eggNOG-vs-GTF compatibility check before OrgDb creation
# ------------------------------------------------------------
EGGNOG_COMPAT_REPORT="$WORK_DIR/eggnog_gtf_compatibility.tsv"
python3 - "$EGGNOG" "$GTF_GENES" "$EGGNOG_COMPAT_REPORT" <<'PY'
import sys
ann, genes_file, report = sys.argv[1:4]
gtf_genes = set(x.strip() for x in open(genes_file) if x.strip())
header = None
rows = 0
with_go = set()
for line in open(ann, errors='replace'):
    line=line.rstrip('\n')
    if line.startswith('#query'):
        header=line[1:].split('\t')
        continue
    if not line or line.startswith('#'):
        continue
    rows += 1
    parts=line.split('\t')
    if header and len(parts)==len(header):
        d=dict(zip(header, parts))
        q=d.get('query', parts[0])
        gos=d.get('GOs') or d.get('GO_terms') or d.get('GO') or ''
    else:
        q=parts[0]
        gos=''
    if q in gtf_genes and gos and gos != '-':
        with_go.add(q)
pct = 100*len(with_go)/len(gtf_genes) if gtf_genes else 0
with open(report,'w') as out:
    out.write('metric\tvalue\n')
    out.write(f'eggnog_rows\t{rows}\n')
    out.write(f'gtf_genes\t{len(gtf_genes)}\n')
    out.write(f'eggnog_genes_with_go_matching_gtf\t{len(with_go)}\n')
    out.write(f'eggnog_gtf_pct\t{pct:.4f}\n')
print(f'[INFO] eggNOG rows: {rows}')
print(f'[INFO] eggNOG genes with GO matching GTF: {len(with_go)} / {len(gtf_genes)} ({pct:.4f}%)')
PY
EGGNOG_MATCH_GENES="$(awk '$1=="eggnog_genes_with_go_matching_gtf"{print $2}' "$EGGNOG_COMPAT_REPORT")"
EGGNOG_MATCH_PCT="$(awk '$1=="eggnog_gtf_pct"{print $2}' "$EGGNOG_COMPAT_REPORT")"
printf "eggnog_genes_with_go_matching_gtf\t%s\n" "$EGGNOG_MATCH_GENES" >> "$INPUT_REPORT"
printf "eggnog_gtf_pct\t%s\n" "$EGGNOG_MATCH_PCT" >> "$INPUT_REPORT"

if awk -v n="$EGGNOG_MATCH_GENES" -v min="$MIN_EGGNOG_GENES" 'BEGIN{exit !(n < min)}'; then
  echo "[ERROR] Too few eggNOG annotated genes match the GTF: $EGGNOG_MATCH_GENES"
  echo "[ERROR] Report: $EGGNOG_COMPAT_REPORT"
  exit 1
fi
if awk -v pct="$EGGNOG_MATCH_PCT" -v min="$MIN_EGGNOG_GTF_PCT" 'BEGIN{exit !(pct < min)}'; then
  echo "[ERROR] eggNOG/GTF compatibility is too low: ${EGGNOG_MATCH_PCT}%"
  echo "[ERROR] Report: $EGGNOG_COMPAT_REPORT"
  exit 1
fi

# ------------------------------------------------------------
# Build OrgDb in R412
# ------------------------------------------------------------
if [[ "$SKIP_ORGDB_BUILD" == "1" ]]; then
    echo "[INFO] Skipping R OrgDb package construction."
    echo "[INFO] Reference-matched eggNOG annotations remain available:"
    echo "  $EGGNOG"
else
    echo "[INFO] Building OrgDb from eggNOG file: $EGGNOG"

    FORCE_ARG=()

    if [[ "$FORCE" == "1" ]]; then
        FORCE_ARG=(--force)
    fi

    conda run -n "$CONDA_R_ENV" Rscript "$R_BUILD" \
        --gtf "$GTF" \
        --eggnog "$EGGNOG" \
        --lib "$LIB" \
        --build-dir "$BUILD_DIR" \
        --version "$VERSION" \
        --taxid "$TAXID" \
        --genus "$GENUS" \
        --species "$SPECIES" \
        --symbol-mode "$SYMBOL_MODE" \
        "${FORCE_ARG[@]}"
fi

cat <<EOF

[OK] Finished gold OrgDb build for $SCI_NAME
[INFO] Input compatibility report:
  $INPUT_REPORT
[INFO] eggNOG compatibility report:
  $EGGNOG_COMPAT_REPORT

Use it before running MTD_SE/MTD_SE3.sh:
  export R_LIBS_USER="$LIB:\${R_LIBS_USER:-}"
EOF
