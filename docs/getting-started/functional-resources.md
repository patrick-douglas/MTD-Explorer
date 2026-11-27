# Local functional resources

MTD Explorer does not redistribute the historical MSigDB GMT or HUMAnN GO/KO dictionaries in new source snapshots. Obtain resources under their respective terms, then prepare them once locally. The analysis performs no network requests to acquire these files. Both exploratory and comparative modes need the GO and KO dictionaries. Only comparative runs explicitly selecting `--ssgsea-gmt default` need the historical MSigDB collection; `auto` and a custom GMT do not use it.

## Existing installations: back up before updating

Git may remove formerly tracked data when updating to a version that no longer tracks it. **Before pulling the update**, save all three files outside the checkout (if present and you are entitled to retain them for local use):

```bash
cd ~/MTD-Explorer
mkdir -p ~/MTD-local-functional/humann ~/MTD-local-functional/msigdb
cp -n aux_scripts/analysis/humann/{koterms,goterms}.txt ~/MTD-local-functional/humann/
cp -n Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt ~/MTD-local-functional/msigdb/
```

After updating, use absolute local paths (or restore to the old, now ignored paths). Do not add the data back to Git. Keep a private backup and apply the applicable source terms.

```bash
export MTD_HUMANN_DICTIONARY_DIR="$HOME/MTD-local-functional/humann"
export MTD_SSGSEA_DEFAULT_GMT="$HOME/MTD-local-functional/msigdb/c2.all.v7.5.1.symbols.gmt"
```

Add these exports to your run script, including Slurm submission/orchestration environments if applicable. The paths are resolved before the driver changes directories. The legacy default paths remain `aux_scripts/analysis/humann/` and `Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt` relative to the repository if the environment variables are not set.

## New installations: GO and KEGG dictionaries

The translator expects two local tab-delimited files with a one-field header `terms`, followed by a unique ID and a description on every row. KO IDs may have a `ko:` prefix, GO IDs have a `GO:` prefix. An existing compatible, lawfully obtained dictionary can be used unchanged. Otherwise prepare a new pair from authorized source files. Review [KEGG terms](https://www.genome.jp/kegg/legal.html) and [KEGG API restrictions](https://www.genome.jp/kegg/rest/) before obtaining data; the KEGG REST API is restricted to eligible academic users for academic use. Obtain a KO list in `list/ko` format through a permitted method and record the release/retrieval date. No automatic download is performed by MTD Explorer. Obtain a dated `go-basic.obo` from the [GO ontology download/archive](https://geneontology.org/docs/download-ontology/) and record its release. GO source data are attributed to the Gene Ontology Consortium under [CC BY 4.0](https://geneontology.org/page/use-and-license).

```bash
python3 aux_scripts/resources/prepare_functional_dictionaries.py \
    --kegg-list /path/to/authorized/ko-list.tsv \
    --kegg-release 'recorded release or retrieval date' \
    --go-obo /path/to/versioned/go-basic.obo \
    --output-dir "$HOME/MTD-local-functional/humann-new"
export MTD_HUMANN_DICTIONARY_DIR="$HOME/MTD-local-functional/humann-new"
```

The converter retains KEGG descriptions, normalizes KO prefixes, extracts primary non-obsolete GO terms and names, and writes a source manifest. It does not overwrite existing dictionaries, remap obsolete/alternate GO IDs, or silently alter descriptions. It rejects incompatible input. The historical R translator is unchanged, including its existing handling of comment characters in descriptions. The newly generated dictionaries may differ in coverage and labels from the historical dictionaries; their equivalence is not assumed. For exact reproduction of historical output use the original resources with the hashes below.

## Legacy MSigDB (optional)

For `--ssgsea-gmt default`, obtain `c2.all.v7.5.1.symbols.gmt` through the [MSigDB release archive](https://www.gsea-msigdb.org/gsea/downloads_archive.jsp) under the [license terms](https://gsea-msigdb.org/gsea/msigdb_license_terms.jsp), including additional conditions for KEGG and BioCarta gene sets. The driver checks the exact historical SHA-256. If you choose a different GMT or release, explicitly pass `--ssgsea-gmt /absolute/path/to/file.gmt` instead. The ordinary default (`--ssgsea-gmt auto`) uses the host master GMT created by `Create_custom_host.sh` and is unchanged.

## Provenance and reproducibility

Historical resource fingerprints (from the MTD Explorer snapshot `9f0ce1953d166afefd817c1c3bc3306f0e7ffe00`):

| File | SHA-256 |
| --- | --- |
| `koterms.txt` | `88ad9fac54da282c34b1126f7fa3adb40b8d728659ab4cae406b8f4f25de6f9e` |
| `goterms.txt` | `3f0ed5ec2dbd1c7a4f38db10170349207c4953108cded7063733684404af3a74` |
| `c2.all.v7.5.1.symbols.gmt` | `06f3f970c7cba8c202b2792fb69c40f2c5f7afc1fd5c6743631057124210b236` |

The historical dictionary release dates are unknown. Each run writes `methods/functional_resources.tsv`, with absolute paths, SHA-256 and identity of the input dictionaries, and of the actual ssGSEA GMT (plus master GMT in auto mode). When available, `dictionary_sources.tsv` is copied as `methods/functional_resources_dictionary_sources.tsv`. Exploratory runs record ssGSEA as not executed. The manifest does not embed dictionary or gene-set data. Keep resource files immutable throughout the run. Resource changes can change labels and, for gene-set changes, enrichment results. Do not share third-party data or derived materials without checking their conditions.
