#!/usr/bin/env bash
# Local functional resource resolution. No network access; no data redistribution.
# Sourced by MTD_explorer.sh; errors are handled by the caller.

MTD_LEGACY_GMT_SHA256=06f3f970c7cba8c202b2792fb69c40f2c5f7afc1fd5c6743631057124210b236
MTD_LEGACY_GO_SHA256=3f0ed5ec2dbd1c7a4f38db10170349207c4953108cded7063733684404af3a74
MTD_LEGACY_KO_SHA256=88ad9fac54da282c34b1126f7fa3adb40b8d728659ab4cae406b8f4f25de6f9e

mtd_resource_error() { printf '[Functional resources] %s\n' "$*" >&2; return 1; }

mtd_resource_hash() {
    local result
    result=$(sha256sum -- "$1") || return 1
    # sha256sum prefixes escaped filenames with a backslash. Reject such paths
    # (and TSV-unsafe names) explicitly rather than recording ambiguous hashes.
    result=${result%% *}
    [[ "$result" =~ ^[0-9a-f]{64}$ ]] || return 1
    printf '%s\n' "$result"
}

mtd_resource_path() {
    local resolved
    resolved=$(realpath -m -- "$1") || return 1
    [[ "$resolved" != *$'\t'* && "$resolved" != *$'\n'* && "$resolved" != *$'\r'* ]] || return 1
    printf '%s\n' "$resolved"
}

mtd_resolve_functional_resources() {
    local root="$1"
    HUMANN_DICTIONARY_DIR=$(mtd_resource_path "${MTD_HUMANN_DICTIONARY_DIR:-$root/aux_scripts/analysis/humann}") || return 1
    SSGSEA_DEFAULT_GMT=$(mtd_resource_path "${MTD_SSGSEA_DEFAULT_GMT:-$root/Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt}") || return 1
}

mtd_validate_dictionary() {
    local file="$1" kind="$2"
    [[ -r "$file" && -s "$file" && -f "$file" ]] || { mtd_resource_error "Missing $kind dictionary: $file. See docs/getting-started/functional-resources.md or set MTD_HUMANN_DICTIONARY_DIR."; return 1; }
    LC_ALL=C awk -F '\t' -v kind="$kind" '
        { sub(/\r$/, "") }
        NR == 1 { if ($0 != "terms") bad=1; next }
        {
            key=$1
            if (kind == "KO") { sub(/^ko:/, "", key); valid=(key ~ /^K[0-9]+$/ && length(key)==6) }
            else { valid=(key ~ /^GO:[0-9]+$/ && length(key)==10) }
            if (NF != 2 || !valid || $2 == "" || seen[key]++) bad=1
        }
        END { exit (bad || NR < 2) ? 1 : 0 }
    ' "$file" || { mtd_resource_error "Invalid $kind dictionary (expected header terms, unique IDs and one description per row): $file"; return 1; }
}

mtd_validate_dictionaries() {
    mtd_validate_dictionary "$HUMANN_DICTIONARY_DIR/koterms.txt" KO &&
    mtd_validate_dictionary "$HUMANN_DICTIONARY_DIR/goterms.txt" GO
}

mtd_validate_legacy_gmt() {
    local hash
    [[ -r "$SSGSEA_DEFAULT_GMT" && -s "$SSGSEA_DEFAULT_GMT" && -f "$SSGSEA_DEFAULT_GMT" ]] || { mtd_resource_error "Legacy MSigDB GMT not found: $SSGSEA_DEFAULT_GMT. Obtain v7.5.1 under its terms or use --ssgsea-gmt auto|FILE. See docs/getting-started/functional-resources.md."; return 1; }
    hash=$(mtd_resource_hash "$SSGSEA_DEFAULT_GMT") || return 1
    [[ "$hash" == "$MTD_LEGACY_GMT_SHA256" ]] || { mtd_resource_error "Legacy GMT SHA-256 does not match the historical v7.5.1 file. Use --ssgsea-gmt FILE for an explicitly chosen different file."; return 1; }
}

mtd_resource_record() {
    local manifest="$1" role="$2" file="$3" identity="$4" hash
    file=$(mtd_resource_path "$file") || return 1
    hash=$(mtd_resource_hash "$file") || return 1
    case "$hash" in
        "$MTD_LEGACY_GO_SHA256") identity='historical MTD GO dictionary; release unknown' ;;
        "$MTD_LEGACY_KO_SHA256") identity='historical MTD KO dictionary; release unknown' ;;
        "$MTD_LEGACY_GMT_SHA256") identity='historical MSigDB C2 v7.5.1 symbols' ;;
    esac
    printf '%s\t%s\t%s\t%s\n' "$role" "$file" "$hash" "$identity" >> "$manifest"
}

mtd_start_resource_manifest() {
    local manifest="$1" mode="$2" requested="$3" sidecar="$HUMANN_DICTIONARY_DIR/dictionary_sources.tsv"
    printf 'resource\tpath\tsha256\tidentity\n' > "$manifest" || return 1
    mtd_resource_record "$manifest" HUMAnN_KO "$HUMANN_DICTIONARY_DIR/koterms.txt" 'user-provided KO dictionary; see source metadata if available' || return 1
    mtd_resource_record "$manifest" HUMAnN_GO "$HUMANN_DICTIONARY_DIR/goterms.txt" 'user-provided GO dictionary; see source metadata if available' || return 1
    if [[ -f "$sidecar" && -s "$sidecar" ]]; then
        mtd_resource_record "$manifest" dictionary_source_metadata "$sidecar" 'user-provided source metadata' || return 1
        cp -- "$sidecar" "${manifest%.tsv}_dictionary_sources.tsv" || return 1
    else
        rm -f -- "${manifest%.tsv}_dictionary_sources.tsv" || return 1
    fi
    if [[ "$mode" == exploratory ]]; then
        printf 'ssGSEA\tNA\tNA\tnot executed in exploratory mode\n' >> "$manifest" || return 1
    elif [[ "$requested" == default ]]; then
        mtd_validate_legacy_gmt || return 1
    fi
}
