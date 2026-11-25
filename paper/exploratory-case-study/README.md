# Exploratory case study: Danio rerio

This supplementary snapshot demonstrates MTD-Explorer software execution and repeatability in exploratory mode, not biological inference or validation of taxonomic identity. Data: [PRJNA1005695](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1005695), Rice et al. (2024), [doi:10.1371/journal.pbio.3002606](https://doi.org/10.1371/journal.pbio.3002606). Analytical revision: [5be26aa2e03938abb12fe00a37a08b0c0cef3777](https://github.com/patrick-douglas/MTD-Explorer/tree/5be26aa2e03938abb12fe00a37a08b0c0cef3777). Six paired-end libraries, one analytical group, no contrasts. Contaminant removal disabled intentionally (local deletion of conta_ls.txt); no negative controls included in the selected libraries.

## Contents and reproduction

- `outputs/`: 116 source tables from the definitive execution, identical to the scientific export. Includes Bracken matrices, TPI, diversity, QC, host counts and microbial functional tables with ID maps.
- `provenance/`: selected ENA metadata, FASTQ MD5s, input samplesheet, reference catalog, runtime versions, parameters, stage markers, GNU time record and recorded source checksums. Local home path prefixes in public provenance are normalized where applicable; these are not byte-identical copies of all original logs.
- `validation/`: independently compared pre-patch and definitive tables, host count differences, and the server-side local regression log. Historical hashes describe the earlier run; the complete historical inputs are not included in this snapshot.
- `validation/verify_outputs.py`: independently recompute TPI and check functional ID maps and metadata using the included tables (Python 3 standard library only).

From this directory, run `sha256sum -c SHA256SUMS.txt` (Linux) or `shasum -a 256 -c SHA256SUMS.txt` (macOS), then `python3 validation/verify_outputs.py`. See the [published example](https://patrick-douglas.github.io/MTD-Explorer/user-guide/danio-rerio/) for download, reference preparation, command and original figures. The archived wrapper is a record of the server configuration, not a turnkey portable installer. Adapt paths and use the analytical commit in a dedicated checkout. Raw FASTQs, full databases, full logs and large intermediate files are not included. Database descriptions and distribution checksums do not constitute a full database snapshot.

## Computational performance and limitations

GNU time reported 4 h 12 min 53 s, 113107988 KiB maximum RSS, 1089% CPU and zero swaps, exit status 0. Output inventory: 678 files, none empty, approximately 459 GiB allocated. Reference and Bracken distribution preparation were separate. Maximum RSS is not sampled aggregate process-tree memory; this workload is not a controlled speed comparison with the historical B. glabrata benchmark. Three libraries had zero reconstructed pathway abundance; this is not evidence of biological absence. TPI describes assigned taxonomic profiles, not confirmed organisms. The comparative GO/KEGG and Slurm end-to-end regressions remain separate tasks. No supplementary table or figure numbering is assigned here.
