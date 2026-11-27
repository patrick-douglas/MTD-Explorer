# Third-party software and data notices

Scope: tracked resources inspected in MTD Explorer commit `4cef8f66af21144d731df4365862311dc9b0a564` and an archived source comparison with MTD `d488faf929b12d810d5b46742975e535e6f48397`. This is a component inventory, not a blanket clearance of all redistribution rights. Existing notices remain in effect. Installed dependencies have their own licenses and are not relicensed by MTD Explorer.

## Bundled source code

| Component / paths | Terms and attribution | Source / license |
| --- | --- | --- |
| `Tools/KrakenTools/` | GPL v3; individual headers may allow later versions. Original copyright notices apply. | [Local license](Tools/KrakenTools/LICENSE); https://github.com/jenniferlu717/KrakenTools |
| `Tools/combine_bracken_outputs.py` | GPL-3.0-or-later, Jennifer Lu, copyright 2016–2017; header retained. | [GPL v3 text](licenses/GPL-3.0.txt); https://github.com/jenniferlu717/Bracken |
| `Installation/hisat2_extract_exons.py`, `Installation/hisat2_extract_splice_sites.py` | GPL-3.0-or-later, Daehwan Kim, copyright 2015; headers retained. | [GPL v3 text](licenses/GPL-3.0.txt); https://github.com/DaehwanKimLab/hisat2 |
| `Tools/graphlan/` upstream code | MIT; Nicola Segata and Curtis Huttenhower (2012). | [Local license](Tools/graphlan/license.txt); https://github.com/biobakery/graphlan |
| `Tools/export2graphlan/` | MIT; Francesco Asnicar, Nicola Segata and Curtis Huttenhower (2015). | [Local license](Tools/export2graphlan/license.txt); https://github.com/SegataLab/export2graphlan |
| `Tools/ssGSEA2.0/` software (not gene-set data) | BSD-style, Broad Institute (2018–2020); also includes cmap/cmapR BSD-3-Clause notices in `src/io.R` and `src/utils.R`. | [Local license](Tools/ssGSEA2.0/LICENSE.txt); https://github.com/broadinstitute/ssGSEA2.0 |
| Kraken2-derived scripts | MIT at the upstream revision of the license inspected; retain Derrick Wood notices. See detailed scope below. | [Upstream license](https://github.com/DerrickWood/kraken2/blob/master/LICENSE); [local license text](licenses/Kraken2-MIT.txt) |

The directly identifiable Kraken2-derived scripts include `Installation/download_genomic_library.sh`, `Installation/download_genomic_library_plasmid.sh`, `Installation/rsync_from_ncbi.pl`, `Installation/rsync_from_ncbi_bacteria.pl`, and `Installation/rsync_from_ncbi_offline.pl`, which carry Derrick Wood notices. The files `Installation/rsync_from_ncbi_2.pl`, `Installation/rsync_from_ncbi_archaea.pl` and `aux_scripts/Kraken2/kraken2-build-download-taxonomy` also warrant source-history review for retained/adapted upstream code. The exact originating Kraken2 revision and any missing notices must be confirmed before release. `Installation/mask_low_complexity.sh` is a documented replacement for the legacy masking script; its origin and any retained upstream code should likewise be checked. The presence of MTD Explorer additions inside third-party directories does not by itself establish their license or authorship; this inventory distinguishes upstream material from project-controlled material.

## Historical and local data resources — release review required

| Path | Provenance and terms | Status |
| --- | --- | --- |
| `Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt` | Byte-identical to the MTD snapshot; MSigDB C2 v7.5.1, copyright Broad Institute, MIT and Regents of the University of California as specified in [MSigDB terms](https://gsea-msigdb.org/gsea/msigdb_license_terms.jsp). The collection includes KEGG legacy and BioCarta sets with additional conditions. | Not included in new source snapshots; locally supplied under applicable terms. Historical redistribution clearance not established. |
| `aux_scripts/analysis/humann/goterms.txt` | Byte-identical to MTD's `goterms.txt`; GO term dictionary. GO Consortium data are generally [CC BY 4.0](https://geneontology.org/page/use-and-license). | Original release and precise derivation undocumented; no longer tracked in new snapshots as a precaution. GO Consortium attribution applies. |
| `aux_scripts/analysis/humann/koterms.txt` | Byte-identical to MTD's `koterms.txt`; KEGG KO descriptions. Copyright Kanehisa Laboratories; see [KEGG terms](https://www.genome.jp/kegg/legal.html). | Not included in new source snapshots; obtain locally through permitted channels. Historical redistribution clearance not established. |
| `Installation/M33262_SIVMM239.fa` | Byte-identical to the MTD snapshot; header identifies sequence accession M33262. See [NCBI record](https://www.ncbi.nlm.nih.gov/nuccore/M33262). | Confirm source/version and applicable data terms. |

MSigDB also requires acknowledgment of BioCarta (2000–2017) and Kanehisa Laboratories (1995–2017 for legacy KEGG sets), and compliance with the linked conditions. This acknowledgment is not a claim of redistribution authorization. The GMT is distinct from the ssGSEA software license.

Externally acquired databases and packages (including Kraken2 reference databases, Ensembl resources, eggNOG, UniRef, HUMAnN/MetaPhlAn resources, KEGG and GO) remain subject to source-specific conditions. This notice does not authorize redistribution of these resources or of third-party materials in derived documentation and supplementary archives. Check the release contents and all bundled archives separately before deposition.

The three historical functional data files above are no longer tracked in new snapshots. Existing local copies are not altered by the pipeline. This does not clear past Git history or other archives. See [local resource documentation](docs/getting-started/functional-resources.md) for acquisition, fingerprints and migration.
