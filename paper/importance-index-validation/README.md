# Taxon Prioritization Index analysis

This directory contains the inputs, executable analysis, resampled rankings and supporting outputs for Supplementary Figure S2 and Supplementary Table S9. It accompanies the post hoc assessment of the Taxon Prioritization Index (TPI) in 15 *Biomphalaria glabrata* RNA-seq libraries from [BioProject PRJNA1306560](https://www.ncbi.nlm.nih.gov/bioproject/1306560). The two input matrices contain the same 425 species-level features, represented as Bracken-estimated reads or as within-library percentages. For the complete original rankings and component values, see [Supplementary Table S7](../supplementary-tables/MTD_Explorer_Supplementary_Tables.xlsx).

## Read the outputs

- [Supplementary Figure S2](figures/Supplementary_Figure_S2_importance_robustness.pdf) and its [caption](figures/Supplementary_Figure_S2_caption.txt) summarize weight sensitivity and sampling stability. The figures directory also contains an editable SVG, a PNG and a CSV for each panel.
- [Supplementary Table S9](tables/Supplementary_Table_S9.xlsx) contains the scenario and sampling summaries. The [table legend](tables/Supplementary_Table_S9_legend.txt) describes the worksheets.
- [Protocol](protocol.json) specifies the weights, detection and taxon inclusion criteria, random seed and comparisons. The [protocol and script hashes recorded before execution](provenance/protocol_identity_before_execution.json) identify the analysis that generated these results.
- [Historical R source](provenance/historical_importance_index.R) contains the ranking implementation from [MTD Explorer commit 04d61d4](https://github.com/patrick-douglas/MTD-Explorer/blob/04d61d412955ac7b9142ca8ebcc4d2edeb6bac5d/aux_scripts/exploratory/MTD_exploratory_detected_microbiome_rank.R). The repository provides the [license for that revision](https://github.com/patrick-douglas/MTD-Explorer/blob/04d61d412955ac7b9142ca8ebcc4d2edeb6bac5d/LICENCE.md). The original protocol, source scripts, column labels and filenames retain the names used when the study was conducted.

## Reproduce the analysis

Work in a copy of the package if you want to keep the distributed files for checksum comparison. The scripts write new outputs under `results/`, `figures/` and `provenance/independent/`. Start from this directory. The original analysis used Python 3.12.14.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python src/run_validation.py
python src/plot_validation.py
python src/audit_historical.py
python src/audit_validation.py
```

The pinned dependencies are NumPy 2.3.5 and Matplotlib 3.10.7. The scripts read the supplied matrices and ranking tables from `data/`. The plotting script uses Arial if installed. The numerical results and panel tables are independent of the font. Execution times, metadata timestamps and compressed file headers can differ on regeneration. To check the scientific outputs already supplied, run the final two commands. The first checks the two original rankings. The second compares all 186 scenario results, all 60 leave-one-out rankings, 96 selected full bootstrap rankings and all 20,000 bootstrap top-20 retention and Jaccard results. Both scripts exit with a nonzero status when these checks fail.

## Interpret the comparisons

The TPI combines prevalence, mean abundance, maximum abundance and total abundance with weights of 0.45, 0.30, 0.15 and 0.10. The analysis uses the complete matrices from five libraries per group. It applies 79 distinct weight vectors, ten named models including the default, two changes to the detection threshold and two restrictions on the taxon universe for each representation.

The sample-based analysis uses 15 single-library omissions and 5,000 bootstrap draws with replacement, preserving five libraries from each group. Both abundance representations and the two models (default and removal of total abundance) use the same draws. Each omission or bootstrap ranking is compared with that model's full-data ranking. A separate field, `historical_default_top20_retained`, gives retention against the original default. An absent taxon contributes zero to its inclusion frequency. Its percentile ranks are calculated across draws in which it was eligible. These are empirical measures of stability in the deposited libraries.

## Locate individual draws and taxa

| Files | Contents |
| --- | --- |
| `data/*_matrix.tsv`, `data/*_historical_*.tsv`, `data/sample_metadata.csv` | Full input matrices, original rankings and abundance summaries, and library labels. |
| `results/baseline_reproduction.csv`, `results/weight_design.csv`, `results/scenario_metrics.csv` | Reproduction of the original rankings, all weight vectors and all scenario-level comparisons. |
| `results/bootstrap_draws.csv.gz`, `results/bootstrap_draw_indices.npy` | The sampled libraries in every draw. The NumPy array gives zero-based column indices in input sample order. |
| `results/bootstrap_metrics.csv.gz`, `results/loo_metrics.csv` | Metrics for every bootstrap draw and omission in each model and representation. |
| `results/*_taxon_summary.csv` | Per-feature eligibility, inclusion frequencies and rank summaries. |
| `results/*arrays*.npz` | Full taxon ranks and scores for fixed scenarios and bootstrap replicates, and full ranks for omissions. |
| `results/summary.json`, `provenance/independent/` | Software versions, input identities, numerical summaries and independent checks. |
| `provenance/workbook/prepare_spec.py`, `provenance/workbook/workbook_spec.json` | Tabular data and column definitions used to assemble Table S9 from the results. Run the script to regenerate the typed specification. |

The second axis of each rank array is the taxon order in the archive's `taxa` array. Rank 426 denotes a taxon ineligible for a given ranking. Its score is NaN in arrays that store scores. Bootstrap row 0 corresponds to replicate 1 in the CSV. For example:

```python
import numpy as np

with np.load("results/bootstrap_arrays_relative_percent_default.npz", allow_pickle=False) as archive:
    print(archive["ranks"].shape)
    print(archive["taxa"][:3], archive["ranks"][0, :3])
```
