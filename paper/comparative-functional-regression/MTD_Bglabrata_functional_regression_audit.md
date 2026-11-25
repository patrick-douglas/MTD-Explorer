# MTD-Explorer — comparative functional regression audit

Date: 2026-11-25. Scope: software regression, not biological inference.

Source archive: `MTD_Bglabrata_functional_regression_sWfKySHL.tar.gz`.
SHA-256: `0617d1f247bc5f3d3abbd50f362de7ba0c0e5d5828b6679b9f5b2e8ffbb3e59d`.

## Execution and integrity

- Before: `af59f09ca6b83bd630cce7cf02f3a3c3726ba0eb`.
- After: `5be26aa2e03938abb12fe00a37a08b0c0cef3777`.
- All 10 logged subprocesses exited 0; no MaAsLin2 error markers; 14 input/source checksums verified independently.
- Same 15 libraries (5 per group), three contrasts, no extended metadata; DESeq2 1.34.0 and MaAsLin2 1.8.0.
- DESeq2 design: `~ group + transcriptome_size`. MaAsLin2: group and numeric transcriptome_size, TSS normalization and LOG transformation.
- Both pre-patch merged DEG CSVs match the historical outputs byte-for-byte.
- Both KEGG runs automatically used a local dispersion fit when the parametric fit was unsuitable. R package build-version warnings and an existing HAllA directory warning were also recorded.

## Translation and downstream matrices

| Metric | GO | KEGG |
|---|---:|---:|
| Original and translated rows | 3750 | 3384 |
| Eligible canonical features after filtering | 3749 | 3383 |
| Eligible features without dictionary description | 37 | 3 |
| Excluded non-functional rows | 1 | 1 |
| Max absolute raw/translated difference | 2.98e-8 | 2.24e-8 |
| HAllA input before | 3750 x 15 | 3384 x 15 |
| HAllA input after | 3749 x 15 | 3383 x 15 |
| Merged complete-case DEG CSV before | 2074 | 1088 |
| Merged complete-case DEG CSV after | 2073 | 1022 |

All raw IDs were preserved in the corrected mapping; no duplicate IDs. Numerical equivalence of translated abundance matrices assessed using relative tolerance 1e-12 and absolute tolerance 1e-10. The excluded row was `UNGROUPED|unclassified` in each matrix. HAllA matrices are finite with unique row IDs, 15 samples, and corrected inputs equal the corresponding transformed matrices numerically. HAllA itself was not run. The pre-patch original IDs were reconstructed by the original merge's canonical-key sorting and checked against all abundance vectors (no duplicate abundance vectors within canonical-key groups).

## DESeq2 comparison on shared canonical features

Significance criterion: adjusted p < 0.05 and absolute log2 fold change > 0.5, as used by the pipeline's volcano/heatmap threshold. Counts are computed from individual contrast tables, not the merged complete-case DEG tables. The excluded non-functional row is not counted on either side.

| Rank | Contrast | Before | After | Shared significant | Jaccard | Pearson LFC | Median absolute LFC difference | Maximum absolute LFC difference |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| GO | infected_vs_uninfected | 1744 | 1732 | 1732 | 0.993119 | 0.998994 | 0.011448 | 15.602473 |
| GO | infected_vs_infection_failed | 1355 | 1344 | 1344 | 0.991882 | 0.998433 | 0.009645 | 15.688368 |
| GO | uninfected_vs_infection_failed | 189 | 189 | 189 | 1.000000 | 0.999774 | 0.001172 | 3.666246 |
| KEGG | infected_vs_uninfected | 2085 | 2039 | 2038 | 0.976989 | 0.997840 | 0.182137 | 14.814713 |
| KEGG | infected_vs_infection_failed | 2026 | 1985 | 1984 | 0.978786 | 0.999901 | 0.155327 | 0.282394 |
| KEGG | uninfected_vs_infection_failed | 65 | 67 | 65 | 0.970149 | 0.995767 | 0.000000 | 14.858248 |

No fold-change direction reversals among features significant in both runs. Large individual changes and changes in adjusted-p NA patterns remain; outputs are not numerically interchangeable.

## MaAsLin2 comparison on matched feature / metadata / value records

All output feature names were mapped to original IDs using the old translation reconstruction and corrected sidecar, including R-style name sanitization. One-to-one joins were verified. The number of records includes both group and transcriptome-size associations.

| Rank | Reference | Common records | Pearson coefficient correlation | Median absolute coefficient difference | Maximum absolute coefficient difference | Before q <= 0.25 (all rows) | After q <= 0.25 (all rows) |
|---|---|---:|---:|---:|---:|---:|---:|
| GO | uninfected | 9774 | 0.999711 | 0.009776 | 0.243372 | 3696 | 3768 |
| GO | infection_failed | 9774 | 0.999644 | 0.009776 | 0.282833 | 3725 | 3789 |
| KEGG | uninfected | 9804 | 0.734442 | 0.007165 | 7.579393 | 3958 | 4082 |
| KEGG | infection_failed | 9804 | 0.739133 | 0.007165 | 6.449859 | 3969 | 3997 |

UNGROUPED accounts for 1.04–20.31% of summed raw GO abundance and 78.65–99.93% of summed raw KEGG abundance across libraries. Removing it changes the TSS denominator substantially, especially for KEGG. This is a plausible explanation for the observed changes; a label-only/filter-only factorial experiment was not performed.

## Disposition

Operational regression: PASS for the tested functional modules. Translation/ID preservation and filtering: PASS. Statistical output equivalence: NOT established and not expected given the changed feature universe. The corrected functional results must not be mixed with historical comparative outputs as if interchangeable. Retain historical benchmark and provenance. If the manuscript presents functional comparative results or associations derived from these matrices, regenerate affected downstream outputs and explicitly link figures/tables to the analytic commit. Full HAllA and Slurm execution are not validated by this archive. Statistical model validity is outside this regression: the implementation rounds HUMAnN RPKs to integers for DESeq2, which should be disclosed and independently assessed if inferential claims are made.
