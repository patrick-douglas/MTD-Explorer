#!/usr/bin/env python3
"""Convert locally and lawfully acquired KEGG KO list and GO OBO to MTD dictionaries.
No downloads. No external dependencies. Historical resources are not reproduced by
newer snapshots; use the recorded historical files for exact reproduction.
"""
import argparse
import hashlib
import re
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def insert(rows, key, value, pattern):
    if not re.fullmatch(pattern, key) or not value or any(c in value for c in '\t\r\n'):
        raise ValueError(f'Invalid ID/description for {key!r} (tabs/newlines are not supported)')
    if key in rows:
        raise ValueError(f'Duplicate ID: {key}')
    rows[key] = value


def parse_ko(path):
    rows = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        key, description = line.split('\t', 1)
        key = 'ko:' + (key[3:] if key.startswith('ko:') else key)
        insert(rows, key, description, r'ko:K[0-9]{5}')
    return rows


def parse_go(path):
    rows, term, version = {}, None, None
    def finish():
        if term is not None:
            # Use primary IDs only. Obsolete terms remain unmapped rather than
            # being assigned a different term automatically.
            if term.get('is_obsolete') != 'true':
                insert(rows, term.get('id', ''), term.get('name', ''), r'GO:[0-9]{7}')
    with path.open(encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('data-version: ') and version is None:
                version = line[len('data-version: '):]
            elif line.startswith('[') and line.endswith(']'):
                finish()
                term = {} if line == '[Term]' else None
            elif term is not None:
                for field in ('id', 'name', 'is_obsolete'):
                    if line.startswith(field + ': '):
                        term[field] = line[len(field)+2:]
        finish()
    if not version:
        raise ValueError('GO OBO data-version is missing; supply a versioned ontology.')
    return rows, version


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kegg-list', type=Path, required=True)
    p.add_argument('--kegg-release', required=True, help='Release reported by KEGG or dated retrieval identifier')
    p.add_argument('--go-obo', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    if any(c in a.kegg_release for c in '\t\r\n'):
        p.error('KEGG release must be a single TSV-safe value')
    ko = parse_ko(a.kegg_list)
    go, version = parse_go(a.go_obo)
    if not ko or not go:
        p.error('A dictionary is empty')
    a.output_dir.mkdir(parents=True, exist_ok=True)
    paths = [a.output_dir / n for n in ('koterms.txt', 'goterms.txt', 'dictionary_sources.tsv')]
    if any(x.exists() or x.is_symlink() for x in paths):
        p.error('Refusing to overwrite existing resources; use a new output directory')
    for path, rows in zip(paths, (ko, go)):
        with path.open('x', encoding='utf-8', newline='\n') as f:
            f.write('terms\n')
            for key in sorted(rows):
                f.write(f'{key}\t{rows[key]}\n')
    with paths[2].open('x', encoding='utf-8', newline='\n') as f:
        f.write('resource\tversion\tsource_sha256\tdictionary_sha256\tsource\ttransformation\n')
        f.write(f'KO\t{a.kegg_release}\t{digest(a.kegg_list)}\t{digest(paths[0])}\tKEGG list/ko; https://www.genome.jp/kegg/legal.html\tID prefix normalization; sorted two-column dictionary\n')
        f.write(f'GO\t{version}\t{digest(a.go_obo)}\t{digest(paths[1])}\tGene Ontology Consortium; https://geneontology.org/page/use-and-license\tPrimary non-obsolete GO IDs and names; sorted two-column dictionary; CC BY 4.0 source\n')
    print(f'Prepared {len(ko)} KO and {len(go)} GO entries in {a.output_dir}')


if __name__ == '__main__':
    main()
