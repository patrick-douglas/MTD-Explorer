#!/usr/bin/env python3
"""
make_gene_representative_fasta.py

Robustly create one representative protein per gene from a protein FASTA.

Why this exists:
- Some Ensembl/EnsemblGenomes protein FASTA headers contain gene:GENEID.
- Some releases/species encode only protein/transcript IDs in the FASTA header.
- Some IDs differ only by a final version suffix, e.g. ENSG000001.1 vs ENSG000001.

This script can parse gene IDs directly from FASTA headers and/or build a
protein/transcript -> gene map from the matching GTF/GFF.
The selected FASTA records are renamed to the GTF gene_id.
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Iterable, Optional, Set, Tuple


def open_text(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


def strip_version(x: Optional[str]) -> Optional[str]:
    if x is None:
        return None
    # Strip common final numeric version suffix from Ensembl-like IDs.
    # Examples: ENSAPLG00000012345.1 -> ENSAPLG00000012345
    return re.sub(r'\.\d+$', '', x)


def clean_id(x: Optional[str]) -> Optional[str]:
    if x is None:
        return None
    x = x.strip().strip(';,').strip('"').strip("'")
    for prefix in (
        'gene:', 'gene=', 'gene_id:', 'gene_id=', 'GeneID:',
        'transcript:', 'transcript=', 'transcript_id:', 'transcript_id=',
        'protein:', 'protein=', 'protein_id:', 'protein_id=',
        'ID=gene:', 'ID=transcript:', 'ID=protein:', 'ID=',
        'Parent=gene:', 'Parent=transcript:', 'Parent=',
    ):
        if x.startswith(prefix):
            x = x[len(prefix):]
    return x.strip().strip(';,').strip('"').strip("'")


def read_gene_list(path: Optional[str]) -> Tuple[Optional[Set[str]], Dict[str, str]]:
    if not path:
        return None, {}
    genes: Set[str] = set()
    aliases: Dict[str, str] = {}
    with open(path, 'rt') as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            gene = clean_id(line.split('\t')[0].split(',')[0])
            if not gene:
                continue
            genes.add(gene)
            aliases[gene] = gene
            aliases[strip_version(gene)] = gene
    return genes, aliases


def parse_attributes(attr: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    # GTF style: key "value";
    for key, val in re.findall(r'([A-Za-z0-9_.:-]+)\s+"([^"]+)"', attr):
        d[key] = val
    # GFF3 style: key=value;key2=value2
    for field in attr.split(';'):
        field = field.strip()
        if not field or '=' not in field:
            continue
        key, val = field.split('=', 1)
        key = key.strip()
        val = val.strip().split(',')[0]
        d.setdefault(key, val)
    return d


def attr_first(attrs: Dict[str, str], keys: Iterable[str]) -> Optional[str]:
    for key in keys:
        val = attrs.get(key)
        if val:
            return clean_id(val)
    return None


def load_gtf_maps(gtf_path: Optional[str], gene_aliases: Dict[str, str]) -> Dict[str, str]:
    """Return alias map: protein/transcript/gene IDs -> canonical GTF gene_id."""
    aliases: Dict[str, str] = {}
    if not gtf_path:
        return aliases

    n_lines = 0
    n_genes = 0
    n_alias = 0

    with open_text(gtf_path) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith('#'):
                continue
            parts = raw.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            n_lines += 1
            feature = parts[2]
            attrs = parse_attributes(parts[8])
            gid = attr_first(attrs, ['gene_id', 'gene', 'geneID', 'ID', 'locus_tag', 'Name'])
            if gid and gid.startswith('gene:'):
                gid = clean_id(gid)
            if not gid:
                continue

            canonical_gid = gene_aliases.get(gid) or gene_aliases.get(strip_version(gid)) or gid
            if canonical_gid:
                n_genes += 1
                for x in (gid, strip_version(gid), canonical_gid, strip_version(canonical_gid)):
                    if x:
                        aliases[x] = canonical_gid

            # Map transcript and protein identifiers to the canonical gene ID.
            candidate_keys = [
                'transcript_id', 'transcript', 'transcriptId',
                'protein_id', 'protein', 'proteinId',
                'ccds_id', 'Name', 'ID', 'Parent',
                'Derives_from', 'derives_from'
            ]
            for key in candidate_keys:
                val = attrs.get(key)
                if not val:
                    continue
                # GFF Parent can be transcript:ID or gene:ID. Multiple values possible.
                for item in str(val).split(','):
                    item = clean_id(item)
                    if not item:
                        continue
                    for x in (item, strip_version(item)):
                        if x:
                            aliases[x] = canonical_gid
                            n_alias += 1

    sys.stderr.write(f'[INFO] GTF/GFF mapping aliases loaded: {len(aliases)}\n')
    if len(aliases) == 0:
        sys.stderr.write('[WARNING] No aliases loaded from GTF/GFF. FASTA header gene parsing only will be used.\n')
    return aliases


def parse_fasta(path: str):
    header = None
    seq_chunks = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_chunks)


def first_token(header: str) -> str:
    return header.split()[0]


def extract_header_fields(header: str, pattern: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    token = clean_id(first_token(header))
    if token:
        fields['record_id'] = token

    # User-provided pattern first.
    if pattern:
        m = re.search(pattern, header)
        if m:
            fields['pattern'] = clean_id(m.group(1) if m.lastindex else m.group(0))

    # Generic Ensembl/GFF-ish header tokens: gene:ID, transcript:ID, protein_id=ID, etc.
    for label in ('gene', 'gene_id', 'transcript', 'transcript_id', 'protein', 'protein_id'):
        m = re.search(r'(?:^|[\s;])' + re.escape(label) + r'[:=]([^\s;]+)', header)
        if m:
            fields[label] = clean_id(m.group(1))

    # Some FASTA descriptions carry [gene=...] [protein_id=...]
    for label in ('gene', 'gene_id', 'protein_id', 'transcript_id'):
        m = re.search(r'\[' + re.escape(label) + r'=([^\]]+)\]', header)
        if m:
            fields[f'bracket_{label}'] = clean_id(m.group(1))

    return fields


def choose_gene_id(fields: Dict[str, str], gtf_aliases: Dict[str, str], gene_aliases: Dict[str, str], id_as_gene: bool) -> Tuple[Optional[str], str, Optional[str]]:
    # Candidate priority. gene-specific tokens first, then GTF alias lookup by transcript/protein/record ID.
    candidate_items = []
    if id_as_gene and fields.get('record_id'):
        candidate_items.append(('record_id_as_gene', fields['record_id']))
    for key in ('gene', 'gene_id', 'pattern', 'bracket_gene', 'bracket_gene_id'):
        if fields.get(key):
            candidate_items.append((key, fields[key]))
    for key in ('record_id', 'transcript', 'transcript_id', 'protein', 'protein_id', 'bracket_protein_id', 'bracket_transcript_id'):
        if fields.get(key):
            candidate_items.append((key, fields[key]))

    for source, val in candidate_items:
        if not val:
            continue
        val = clean_id(val)
        for x in (val, strip_version(val)):
            if not x:
                continue
            if x in gene_aliases:
                return gene_aliases[x], source, val
            if x in gtf_aliases:
                return gtf_aliases[x], source, val

    # No allowed gene list: accept gene-like fields directly.
    if not gene_aliases:
        for source, val in candidate_items:
            if source in ('gene', 'gene_id', 'pattern', 'record_id_as_gene') and val:
                return clean_id(val), source, val
    return None, 'unmapped', None


def wrap(seq: str, width: int = 60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--protein-fasta', required=True)
    ap.add_argument('--out-fasta', required=True)
    ap.add_argument('--report', required=True)
    ap.add_argument('--gene-pattern', default=r'(?:gene[:=]|gene_id[:=]|gene=)([A-Za-z0-9_.:-]+)',
                    help='Regex used to extract gene ID from FASTA header. First capture group is used.')
    ap.add_argument('--id-as-gene', action='store_true',
                    help='Use FASTA record ID as gene ID before trying other fields.')
    ap.add_argument('--gene-list', default=None,
                    help='Optional list of allowed canonical GTF gene IDs.')
    ap.add_argument('--gtf', default=None,
                    help='Optional matching GTF/GFF. Used to map protein/transcript IDs to gene IDs.')
    ap.add_argument('--debug-unmapped', default=None,
                    help='Optional output file with examples of unmapped FASTA headers.')
    ap.add_argument('--max-debug', type=int, default=50)
    args = ap.parse_args()

    allowed_genes, gene_aliases = read_gene_list(args.gene_list)
    gtf_aliases = load_gtf_maps(args.gtf, gene_aliases)

    best = {}
    isoform_count = defaultdict(int)
    source_count = defaultdict(int)
    total = 0
    header_extracted = 0
    mapped = 0
    unmapped_examples = []

    for header, seq in parse_fasta(args.protein_fasta):
        total += 1
        fields = extract_header_fields(header, args.gene_pattern)
        if len(fields) > 1 or ('record_id' in fields and len(fields) == 1):
            header_extracted += 1
        gene_id, source, raw_id = choose_gene_id(fields, gtf_aliases, gene_aliases, args.id_as_gene)
        if not gene_id:
            if len(unmapped_examples) < args.max_debug:
                unmapped_examples.append((header, ';'.join(f'{k}={v}' for k, v in fields.items())))
            continue
        if allowed_genes is not None and gene_id not in allowed_genes:
            # Should be rare because gene_aliases maps to allowed IDs, but keep safe.
            if len(unmapped_examples) < args.max_debug:
                unmapped_examples.append((header, f'mapped_to_not_allowed={gene_id};raw={raw_id};source={source}'))
            continue
        mapped += 1
        source_count[source] += 1
        isoform_count[gene_id] += 1
        seq_id = first_token(header)
        length = len(seq.replace('*', ''))
        if gene_id not in best or length > best[gene_id]['length']:
            best[gene_id] = {
                'gene_id': gene_id,
                'seq_id': seq_id,
                'length': length,
                'header': header,
                'seq': seq,
                'source': source,
                'raw_id': raw_id or '',
            }

    os.makedirs(os.path.dirname(os.path.abspath(args.out_fasta)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)

    with open(args.out_fasta, 'wt') as out:
        for gene_id in sorted(best):
            rec = best[gene_id]
            out.write(f'>{gene_id}\n')
            for line in wrap(rec['seq']):
                out.write(line + '\n')

    with open(args.report, 'wt') as rep:
        rep.write('gene_id\tselected_seq_id\tlength\tn_isoforms\tmapping_source\traw_mapped_id\toriginal_header\n')
        for gene_id in sorted(best):
            rec = best[gene_id]
            rep.write(f"{gene_id}\t{rec['seq_id']}\t{rec['length']}\t{isoform_count[gene_id]}\t{rec['source']}\t{rec['raw_id']}\t{rec['header']}\n")

    if args.debug_unmapped:
        with open(args.debug_unmapped, 'wt') as dbg:
            dbg.write('header\textracted_fields_or_reason\n')
            for h, fields in unmapped_examples:
                dbg.write(h.replace('\t', ' ') + '\t' + fields.replace('\t', ' ') + '\n')

    sys.stderr.write(f'[INFO] Total FASTA records: {total}\n')
    sys.stderr.write(f'[INFO] FASTA records with parseable header fields: {header_extracted}\n')
    sys.stderr.write(f'[INFO] Records mapped to GTF genes: {mapped}\n')
    sys.stderr.write(f'[INFO] Representative genes: {len(best)}\n')
    if source_count:
        sys.stderr.write('[INFO] Mapping sources:\n')
        for k, v in sorted(source_count.items(), key=lambda kv: (-kv[1], kv[0])):
            sys.stderr.write(f'  {k}: {v}\n')
    sys.stderr.write(f'[OK] Wrote representative FASTA: {args.out_fasta}\n')
    sys.stderr.write(f'[OK] Wrote report: {args.report}\n')
    if args.debug_unmapped:
        sys.stderr.write(f'[OK] Wrote unmapped header examples: {args.debug_unmapped}\n')

    if len(best) == 0:
        sys.stderr.write('[ERROR] No representative proteins were generated.\n')
        sys.stderr.write('[ERROR] This usually means protein FASTA and GTF are from different releases/assemblies, or headers need a custom parser.\n')
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
