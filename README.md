<p align="center">
  <a href="https://patrick-douglas.github.io/MTD-Explorer/">
    <img src="docs/assets/images/mtd-explorer-logo-1024.png" alt="MTD Explorer logo" width="150">
  </a>
</p>

<h1 align="center">MTD Explorer</h1>

<p align="center">
  <strong>Joint host transcriptome and microbiome analysis from bulk RNA-seq data.</strong>
</p>

<p align="center">
  <a href="https://patrick-douglas.github.io/MTD-Explorer/"><strong>Documentation</strong></a>
  |
  <a href="https://patrick-douglas.github.io/MTD-Explorer/getting-started/installation/">Installation</a>
  |
  <a href="https://patrick-douglas.github.io/MTD-Explorer/getting-started/quick-start/">Quick start</a>
  |
  <a href="https://patrick-douglas.github.io/MTD-Explorer/user-guide/">User guide</a>
  |
  <a href="https://patrick-douglas.github.io/MTD-Explorer/reference/">References</a>
</p>

<p align="center">
  <a href="https://patrick-douglas.github.io/MTD-Explorer/">
    <img alt="Documentation" src="https://img.shields.io/badge/docs-MTD%20Explorer-0b3d66">
  </a>
  <a href="https://github.com/patrick-douglas/MTD">
    <img alt="GitHub repository" src="https://img.shields.io/badge/repository-GitHub-black">
  </a>
</p>

---

## What is MTD Explorer?

**MTD Explorer** is an open-source workflow for joint analysis of host
transcriptome and microbiome/metatranscriptome signals from **bulk RNA-seq**
data.

It builds upon the original **MTD** pipeline developed by Fei Wu, Yao-Zhong Liu,
and Binhua Ling, and extends it with updated installation, host support,
taxonomic exploration, functional profiling, reproducibility summaries, and
documentation.

The original MTD pipeline should be cited when using MTD Explorer.

## Documentation

The full documentation is available here:

**https://patrick-douglas.github.io/MTD-Explorer/**

Start with:

- [Installation](https://patrick-douglas.github.io/MTD-Explorer/getting-started/installation/)
- [Verify installation](https://patrick-douglas.github.io/MTD-Explorer/getting-started/verify-installation/)
- [Quick start](https://patrick-douglas.github.io/MTD-Explorer/getting-started/quick-start/)
- [User guide](https://patrick-douglas.github.io/MTD-Explorer/user-guide/)
- [Academic references](https://patrick-douglas.github.io/MTD-Explorer/reference/)

## Main analysis layers

MTD Explorer can generate outputs for:

- host transcriptome analysis;
- host differential expression;
- microbiome taxonomic profiling;
- taxonomic exploratory summaries;
- Krona and GraPhlAn visualizations;
- microbiome group comparisons;
- HUMAnN-derived functional profiling;
- GO and KEGG functional summaries;
- ssGSEA host gene-set activity analysis;
- HAllA host-microbiome association analysis;
- reproducibility and methods reporting.

## Scope

MTD Explorer currently focuses on **bulk RNA-seq** host-microbiome analysis.

The single-cell workflow from the original MTD pipeline is not part of the
current MTD Explorer scope, tests, or claims.

## Quick repository setup

Clone the repository:

```bash
git clone https://github.com/patrick-douglas/MTD.git
cd MTD
```

Check the available command-line options:

```bash
bash MTD_explorer.sh --help
```

For installation and database setup, use the documentation:

https://patrick-douglas.github.io/MTD-Explorer/getting-started/installation/

## Output documentation

The documentation includes dedicated pages for major output folders:

- [Output files](https://patrick-douglas.github.io/MTD-Explorer/user-guide/output-files/)
- [Taxonomic exploratory outputs](https://patrick-douglas.github.io/MTD-Explorer/user-guide/taxonomic-exploratory-outputs/)
- [Taxonomic visualizations](https://patrick-douglas.github.io/MTD-Explorer/user-guide/taxonomic-visualizations/)
- [Host expression outputs](https://patrick-douglas.github.io/MTD-Explorer/user-guide/host-expression-outputs/)
- [Microbiome comparison outputs](https://patrick-douglas.github.io/MTD-Explorer/user-guide/microbiome-comparison-outputs/)
- [Functional profiling outputs](https://patrick-douglas.github.io/MTD-Explorer/user-guide/functional-profiling-outputs/)
- [ssGSEA outputs](https://patrick-douglas.github.io/MTD-Explorer/user-guide/ssgsea-outputs/)
- [HAllA integration outputs](https://patrick-douglas.github.io/MTD-Explorer/user-guide/halla-integration-outputs/)
- [Methods and reproducibility outputs](https://patrick-douglas.github.io/MTD-Explorer/user-guide/methods-reproducibility-outputs/)

## Citation

If you use MTD Explorer, cite the original MTD publication:

> Wu, F., Liu, Y.-Z., & Ling, B. (2022). MTD: a unique pipeline for host and
> meta-transcriptome joint and integrative analyses of RNA-seq data. *Briefings
> in Bioinformatics*, 23(3), bbac111. https://doi.org/10.1093/bib/bbac111

Also cite the individual tools and databases used in your analysis. See:

https://patrick-douglas.github.io/MTD-Explorer/reference/

In particular, if GNU Parallel is used in work leading to publication, check
the recommended citation with:

```bash
parallel --citation
```

## Acknowledgements

MTD Explorer is built upon the original MTD pipeline by Fei Wu, Yao-Zhong Liu,
and Binhua Ling.

The original repository is available at:

https://github.com/FEI38750/MTD

## License

Please see the repository license file for usage terms.
