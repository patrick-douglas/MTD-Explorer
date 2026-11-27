#!/usr/bin/env bash
# Fast offline tests using local data; R translation regression when Rscript is available.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
source "$ROOT/aux_scripts/resources/functional_resources.sh"
WORK=$(mktemp -d)
trap 'rm -rf -- "$WORK"' EXIT
unset MTD_HUMANN_DICTIONARY_DIR MTD_SSGSEA_DEFAULT_GMT
mtd_resolve_functional_resources "$ROOT"
mtd_validate_dictionaries
mtd_validate_legacy_gmt
[[ $(mtd_resource_hash "$HUMANN_DICTIONARY_DIR/koterms.txt") == "$MTD_LEGACY_KO_SHA256" ]]
[[ $(mtd_resource_hash "$HUMANN_DICTIONARY_DIR/goterms.txt") == "$MTD_LEGACY_GO_SHA256" ]]
ORIGINAL_DICT="$HUMANN_DICTIONARY_DIR"
mkdir -p "$WORK/local data"
cp -- "$ORIGINAL_DICT"/{koterms,goterms}.txt "$WORK/local data/"
cp -- "$SSGSEA_DEFAULT_GMT" "$WORK/local data/legacy.gmt"
cd "$WORK"
export MTD_HUMANN_DICTIONARY_DIR='local data'
export MTD_SSGSEA_DEFAULT_GMT='local data/legacy.gmt'
mtd_resolve_functional_resources "$ROOT"
mtd_validate_dictionaries
mtd_start_resource_manifest "$WORK/default.tsv" comparison default
mtd_resource_record "$WORK/default.tsv" ssGSEA_analysis_GMT "$SSGSEA_DEFAULT_GMT" default_MSigDB_c2_symbols
grep -Fq "$MTD_LEGACY_GMT_SHA256" "$WORK/default.tsv"
grep -Fq "$MTD_LEGACY_GO_SHA256" "$WORK/default.tsv"
printf 'custom_set\tdescription\tA\tB\n' > "$WORK/custom.gmt"
export MTD_SSGSEA_DEFAULT_GMT='missing.gmt'
mtd_resolve_functional_resources "$ROOT"
mtd_start_resource_manifest "$WORK/exploratory.tsv" exploratory default
grep -q 'not executed' "$WORK/exploratory.tsv"
mtd_start_resource_manifest "$WORK/auto.tsv" comparison auto
mtd_resource_record "$WORK/auto.tsv" ssGSEA_master_GMT "$WORK/custom.gmt" master
mtd_resource_record "$WORK/auto.tsv" ssGSEA_analysis_GMT "$WORK/custom.gmt" auto
mtd_start_resource_manifest "$WORK/custom.tsv" comparison "$WORK/custom.gmt"
mtd_resource_record "$WORK/custom.tsv" ssGSEA_analysis_GMT "$WORK/custom.gmt" custom_path
if mtd_validate_legacy_gmt; then echo 'FAIL: missing GMT accepted'; exit 1; fi
export MTD_SSGSEA_DEFAULT_GMT='custom.gmt'
mtd_resolve_functional_resources "$ROOT"
if mtd_validate_legacy_gmt; then echo 'FAIL: wrong GMT accepted'; exit 1; fi
printf 'terms\nko:K00001\tExample\nK00001\tDuplicate\n' > "$WORK/bad.txt"
if mtd_validate_dictionary "$WORK/bad.txt" KO; then echo 'FAIL: duplicate dictionary accepted'; exit 1; fi
if mtd_validate_dictionary "$WORK/missing.txt" GO; then echo 'FAIL: missing dictionary accepted'; exit 1; fi
printf 'K00001\tExample #1\nko:K00002\tSecond\n' > "$WORK/ko.tsv"
cat > "$WORK/go.obo" <<'OBO'
format-version: 1.2
data-version: releases/2026-01-01

[Term]
id: GO:0000001
name: mitochondrion inheritance

[Term]
id: GO:0000002
name: old term
is_obsolete: true

[Typedef]
id: unrelated
name: unrelated
OBO
python3 "$ROOT/aux_scripts/resources/prepare_functional_dictionaries.py" --kegg-list "$WORK/ko.tsv" --kegg-release test --go-obo "$WORK/go.obo" --output-dir "$WORK/new"
mtd_validate_dictionary "$WORK/new/koterms.txt" KO
mtd_validate_dictionary "$WORK/new/goterms.txt" GO
[[ $(wc -l < "$WORK/new/goterms.txt") == 2 ]]
export MTD_HUMANN_DICTIONARY_DIR="$WORK/new"
mtd_resolve_functional_resources "$ROOT"
mtd_start_resource_manifest "$WORK/new.tsv" exploratory auto
[[ -s "$WORK/new_dictionary_sources.tsv" ]]
echo 'PASS: historical hashes, local path with spaces, GMT mode preflight, errors, dictionary conversion and provenance'
# Run the unchanged translator twice with the same dictionaries through different paths.
if ! command -v Rscript >/dev/null 2>&1; then
    echo 'SKIP: Rscript unavailable (run in R412 before committing)'
    exit 0
fi
for d in baseline relocated; do
    mkdir -p "$WORK/$d"
    printf '# Gene Family\tS1_Abundance-RPKs\tS2_Abundance-RPKs\nK00001\t1.5\t0\nK00001|example\t2\t3\nK99999\t0\t1\nUNGROUPED|unclassified\t1\t2\n' > "$WORK/$d/ko.tsv"
    printf '# Gene Family\tS1_Abundance-RPKs\tS2_Abundance-RPKs\nGO:0000001\t1\t2\nGO:0000001|example\t2\t3\nGO:9999999\t0\t1\nUNGROUPED|unclassified\t1\t2\n' > "$WORK/$d/go.tsv"
    if [[ "$d" == baseline ]]; then dict="$ORIGINAL_DICT"; else dict="$WORK/local data"; fi
    Rscript "$ROOT/aux_scripts/analysis/humann/humann_ID_translation_adjusted.R" "$WORK/$d/ko.tsv" "$WORK/$d/go.tsv" "$dict"
done
for name in ko_translated.tsv ko_translated_id_map.tsv go_translated.tsv go_translated_id_map.tsv; do
    cmp "$WORK/baseline/$name" "$WORK/relocated/$name"
done
echo 'PASS: four R translation outputs byte-identical with relocated dictionaries'
