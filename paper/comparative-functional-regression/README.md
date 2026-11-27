# Comparative functional regression: Biomphalaria glabrata

This is an isolated software regression, not a biological reanalysis or a new end-to-end benchmark. The original comparative materials remain unchanged.

- Before: `af59f09ca6b83bd630cce7cf02f3a3c3726ba0eb`
- After: `5be26aa2e03938abb12fe00a37a08b0c0cef3777`
- Inputs: 15 libraries, three groups and three contrasts, using saved GO/KEGG abundance matrices.
- [Audit report](MTD_Bglabrata_functional_regression_audit.md)
- [Public regression archive](MTD_Bglabrata_functional_regression_public.tar.gz) (contains both versions, inputs excluding two external dictionaries, scripts, logs, historical comparison tables, selected native PDF figures and provenance).

The pre-patch merged DEG tables match historical files byte-for-byte. Corrected GO/KEGG profiles retain all raw rows and omit one non-functional row per matrix from statistical testing. Statistical results are not interchangeable; model assumptions were not independently validated. The archive contains HAllA-ready functional matrices, not a new HAllA execution. The standard taxon–host HAllA driver does not consume these functional matrices.

## Public archive and historical resources

The public archive omits only `inputs/goterms.txt` and `inputs/koterms.txt`. All 171 retained historical files are byte-identical to the original archive (SHA-256 `0617d1f247bc5f3d3abbd50f362de7ba0c0e5d5828b6679b9f5b2e8ffbb3e59d`). The audit report retains its historical archive name and checksum intentionally. The public archive includes a new file manifest and a manifest for the available original inputs/source files; the original manifests remain unchanged as historical evidence. See `README_PUBLIC_RESOURCES.md` inside the archive for local restoration and resource terms. Derived results retain their functional labels; redistribution conditions for derived materials still require separate review.

## Verify and inspect

From this directory:

```bash
sha256sum -c SHA256SUMS.txt
mkdir regression
tar -xzf MTD_Bglabrata_functional_regression_public.tar.gz -C regression
cd regression
sha256sum -c provenance/public_archive_files.sha256
sha256sum -c provenance/public_input_and_source.sha256
cat provenance/stage_status.tsv
```

The archived `provenance/run_script.sh` records the actual execution, with original server paths and environment names; adapt these to a local installation before rerunning. The new web PNGs are rendered from the corrected `after/hmn_genefamily_abundance_files/Nonhost_hmn_DEG/{GO,KEGG}/{PCA_color,heatmap_thumbnail}.pdf` files with only outer white margins cropped. No native PDF was altered.
