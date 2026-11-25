# Runtime commands

The legacy test used the following command with the original MTD checkout at commit `d488faf929b12d810d5b46742975e535e6f48397`, the historical run directory `BGL`, and a fresh working directory `OUT`:

```bash
SSGSEA="$ORIGINAL/Tools/ssGSEA2.0"
cd "$OUT"
conda run --no-capture-output -n R412 Rscript "$SSGSEA/ssgsea-cli.R" \
    -i "$BGL/ssGSEA/host.gct" -o "$OUT/legacy_ssgsea" \
    -d "$SSGSEA/db/msigdb/c2.all.v7.5.1.symbols.gmt" \
    -y "$SSGSEA/config.yaml" -u 1 > legacy_ssgsea_execution.log 2>&1
echo "$?" > exit_status.txt
```

The expected outcome is exit status 1 with the explicit no-overlap diagnostic, not any arbitrary error. The original run was performed by the operator; its raw log is retained in `runtime/`. A fresh installation of the legacy tool may need different dependencies. The OrgDb query can be repeated using the supplied R script and a compatible installed custom package:

```bash
export R_LIBS_USER="$REPO/custom_R_libs${R_LIBS_USER:+:$R_LIBS_USER}"
mkdir -p "$OUT/orgdb"
conda run --no-capture-output -n R412 Rscript /path/to/query_orgdb.R \
    "$OUT/orgdb" "$ORIGINAL/HostSpecies.csv"
```

This query also checks the original catalog's zero-length BioMart dataset selection without calling the remote service.
