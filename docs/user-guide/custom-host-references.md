# Custom host references

[MTD Explorer][mtd-explorer] can build the host reference needed for a selected
species with:

```bash
bash Create_custom_host.sh
```

The current workflow has two modes:

1. **Automatic mode** for species already curated in `HostSpecies.csv`.
2. **Manual mode** for species that are not yet resolvable from
   `HostSpecies.csv`.

This means users do **not** always need to provide genome, GTF, and protein
files manually. If the requested NCBI Taxon ID is present in the curated
`HostSpecies.csv` with complete reference metadata, MTD Explorer can resolve the
genome, GTF, and protein FASTA automatically.

## Quick use

### Automatic mode: species already in HostSpecies.csv

For a species with complete curated metadata in `HostSpecies.csv`, run only:

```bash
cd ~/MTD

bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id>
```

Example:

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id 59463
```

In automatic mode, `Create_custom_host.sh` reads `HostSpecies.csv`, resolves the
genome FASTA, GTF, and protein FASTA URLs, downloads or reuses them through the
persistent MTD cache, and builds the host reference.

### Manual mode: species not resolvable from HostSpecies.csv

Use manual mode when the Taxon ID is absent from `HostSpecies.csv`, or when the
row does not contain complete automatic reference metadata.

```bash
cd ~/MTD

bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --scientific-name "Genus species" \
  --genome /path/to/genome.fa.gz \
  --gtf-file /path/to/annotation.gtf.gz \
  --protein-fasta /path/to/proteins.fa.gz
```

Manual inputs can be local paths or URLs:

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id 6526 \
  --scientific-name "Biomphalaria glabrata" \
  --genome "https://server/species.dna.toplevel.fa.gz" \
  --gtf-file "https://server/species.annotation.gtf.gz" \
  --protein-fasta "https://server/species.pep.all.fa.gz"
```

!!! important

    Do not mix partial automatic and partial manual inputs. Use either
    automatic mode with only `--ncbi-taxon-id`, or manual mode with
    `--genome`, `--gtf-file`, and `--protein-fasta`.

## What HostSpecies.csv does

`HostSpecies.csv` is the host reference registry used by MTD Explorer.

It is not just a species list. It also stores the metadata required to resolve
automatic references.

For automatic mode, the most important fields are:

| Column | Purpose |
| --- | --- |
| `Taxon_ID` | NCBI Taxon ID requested with `--ncbi-taxon-id` |
| `Scientific_name` | Scientific name of the requested host |
| `Reference_Taxon_ID` | Taxon ID of the genome reference actually used |
| `Reference_Scientific_name` | Scientific name of the reference species |
| `Ensembl_name` | Ensembl-style species identifier |
| `Ensembl_division` | Ensembl division used to locate files |
| `Ensembl_release` | Ensembl release used for the reference |
| `Assembly` | Genome assembly name |
| `Genome_URL` | Genome FASTA URL |
| `GTF_URL` | Gene annotation GTF URL |
| `Pep_URL` | Protein FASTA URL |
| `Reference_status` | Whether the row is complete enough for automatic resolution |
| `OrgDb` | Custom or generated R OrgDb package associated with the Taxon ID, when available |
| `kegg` or KEGG-related field | Optional KEGG species code/name used only when species-specific KEGG outputs are possible |

A species can appear in `HostSpecies.csv` but still fail automatic resolution if
its reference metadata is incomplete or its `Reference_status` is not marked as
complete.

## Check whether a species is curated

Search by scientific name:

```bash
cd ~/MTD

grep -i "Carollia" HostSpecies.csv
grep -i "Sturnus" HostSpecies.csv
```

Search by NCBI Taxon ID:

```bash
grep ",59463," HostSpecies.csv
```

For a more careful check, inspect the full row:

```bash
awk -F',' '$1 == "59463" {print}' HostSpecies.csv
```

The automatic mode requires complete reference metadata, especially
`Genome_URL`, `GTF_URL`, and `Pep_URL`.

## Automatic mode behavior

When automatic mode is used, MTD Explorer:

1. reads the requested `Taxon_ID` from `HostSpecies.csv`;
2. validates that the row has complete curated reference metadata;
3. resolves `Genome_URL`, `GTF_URL`, and `Pep_URL`;
4. stores or reuses the files in the persistent installation cache;
5. builds the host [Kraken 2][kraken2] database;
6. builds [HISAT2][hisat2] indexes;
7. builds a [Magic-BLAST][magic-blast] database;
8. creates reference-matched functional host resources;
9. optionally builds a custom [OrgDb][orgdb] package.

The persistent cache is normally read from:

```bash
~/MTD/offlineCachePath
```

Manual override is still possible:

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --offline-folder /path/to/cache
```

## Where to download public reference files

Manual custom-host mode requires three coordinated files:

| File needed by MTD Explorer | Common public file type |
| --- | --- |
| `--genome` | genome FASTA, usually `.fa.gz`, `.fna.gz`, or `.fasta.gz` |
| `--gtf-file` | gene annotation GTF, usually `.gtf.gz` |
| `--protein-fasta` | protein FASTA, usually `.pep.all.fa.gz`, `.faa.gz`, or protein `.fa.gz` |

The safest rule is: **download genome FASTA, GTF, and protein FASTA from the
same source, assembly, and release whenever possible.**

Do not mix an [NCBI][ncbi] genome with an [Ensembl][ensembl] GTF or a protein
file from another release unless you have checked chromosome names, gene IDs,
and annotation compatibility.

### Recommended source order

| Situation | Recommended source |
| --- | --- |
| Vertebrate host species curated by Ensembl | [Ensembl][ensembl-downloads] |
| Non-vertebrate eukaryotes such as metazoa, fungi, plants, and protists | [Ensembl Genomes][ensembl-genomes-downloads] |
| Species not available in Ensembl or when a RefSeq/GenBank assembly is preferred | [NCBI Datasets][ncbi-datasets] |
| Human or mouse with GENCODE-style annotation needs | [GENCODE][gencode] |
| UCSC-only assemblies or browser-based reference needs | [UCSC Genome Browser downloads][ucsc-downloads] |

### Ensembl

[Ensembl][ensembl-downloads] is usually the best source for vertebrate species
when available, because the genome FASTA, GTF, and protein FASTA are organized
by release and species.

Typical URL pattern:

```text
https://ftp.ensembl.org/pub/release-<release>/fasta/<ensembl_species>/dna/
https://ftp.ensembl.org/pub/release-<release>/gtf/<ensembl_species>/
https://ftp.ensembl.org/pub/release-<release>/fasta/<ensembl_species>/pep/
```

Common file names:

```text
<Species>.<Assembly>.dna.toplevel.fa.gz
<Species>.<Assembly>.<release>.gtf.gz
<Species>.<Assembly>.pep.all.fa.gz
```

Example download structure:

```bash
mkdir -p ~/MTD_host_refs/example_species
cd ~/MTD_host_refs/example_species

wget -c "<genome_fasta_url>"  -O genome.fa.gz
wget -c "<gtf_url>"           -O annotation.gtf.gz
wget -c "<protein_fasta_url>" -O proteins.fa.gz
```

Then use manual mode:

```bash
cd ~/MTD

bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --scientific-name "Genus species" \
  --genome ~/MTD_host_refs/example_species/genome.fa.gz \
  --gtf-file ~/MTD_host_refs/example_species/annotation.gtf.gz \
  --protein-fasta ~/MTD_host_refs/example_species/proteins.fa.gz
```

### Ensembl Genomes

[Ensembl Genomes][ensembl-genomes-downloads] covers many non-vertebrate
organisms, including metazoa, plants, fungi, protists, and bacteria.

Typical URL pattern:

```text
https://ftp.ensemblgenomes.ebi.ac.uk/pub/<division>/release-<release>/fasta/<ensembl_species>/dna/
https://ftp.ensemblgenomes.ebi.ac.uk/pub/<division>/release-<release>/gtf/<ensembl_species>/
https://ftp.ensemblgenomes.ebi.ac.uk/pub/<division>/release-<release>/fasta/<ensembl_species>/pep/
```

The `<division>` part is usually one of:

```text
metazoa
plants
fungi
protists
bacteria
```

Use Ensembl Genomes in the same way as Ensembl: pick one release, one assembly,
and download the matching genome FASTA, GTF, and protein FASTA.

### NCBI Datasets

[NCBI Datasets][ncbi-datasets] is often the easiest option when a species is not
available in Ensembl or when you want a specific [RefSeq][ncbi-refseq] or
[GenBank][ncbi-genbank] assembly.

The [NCBI Datasets command-line tool][ncbi-datasets-cli] can download genome,
GTF, and protein files in one package.

Example by taxon name:

```bash
mkdir -p ~/MTD_host_refs/ncbi_example
cd ~/MTD_host_refs/ncbi_example

datasets download genome taxon "Genus species" \
  --reference \
  --include genome,gtf,protein \
  --filename ncbi_dataset.zip

unzip ncbi_dataset.zip
```

Example by assembly accession:

```bash
datasets download genome accession GCF_000000000.0 \
  --include genome,gtf,protein \
  --filename ncbi_dataset.zip

unzip ncbi_dataset.zip
```

Find the downloaded files:

```bash
find ncbi_dataset -type f \
  \( -name "*genomic.fna" -o -name "*genomic.fna.gz" \
     -o -name "*.gtf" -o -name "*.gtf.gz" \
     -o -name "*protein.faa" -o -name "*protein.faa.gz" \)
```

Then pass the matching files to `Create_custom_host.sh`:

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --scientific-name "Genus species" \
  --genome /path/to/*genomic.fna.gz \
  --gtf-file /path/to/*.gtf \
  --protein-fasta /path/to/*protein.faa
```

!!! warning

    Some NCBI assemblies may provide GFF3 but not GTF. MTD Explorer expects a
    GTF file for the host build. Prefer assemblies with GTF available, or
    convert GFF3 to GTF only if you know how to validate the result.

### NCBI Assembly FTP

For advanced users, [NCBI Assembly][ncbi-assembly] pages also link to FTP
folders for individual assemblies.

Common files in an NCBI assembly FTP folder include:

```text
*_genomic.fna.gz
*_genomic.gtf.gz
*_protein.faa.gz
```

When available, these can be used directly with manual mode.

### GENCODE

[GENCODE][gencode] is mainly useful for human and mouse analyses when a
GENCODE-specific reference is desired.

GENCODE provides GTF annotation and FASTA files. For MTD Explorer, make sure
you can obtain a compatible genome FASTA and protein FASTA for the same genome
build and annotation release.

### UCSC Genome Browser downloads

[UCSC Genome Browser downloads][ucsc-downloads] provides sequence and
annotation downloads for genome assemblies shown in the UCSC Genome Browser.

Use UCSC carefully for custom host creation, because MTD Explorer needs a
coordinated genome FASTA, GTF, and protein FASTA. For many non-model organisms,
Ensembl, Ensembl Genomes, or NCBI Datasets are usually easier and safer.

### File compatibility checklist

Before running manual mode, check:

```text
[ ] genome FASTA, GTF, and protein FASTA are from the same source
[ ] genome FASTA, GTF, and protein FASTA are from the same release
[ ] assembly name matches across files
[ ] chromosome/scaffold names in the GTF match the genome FASTA
[ ] protein FASTA gene/transcript IDs can be connected to the GTF
[ ] files are compressed or uncompressed formats accepted by Create_custom_host.sh
```

A quick chromosome-name check:

```bash
zcat genome.fa.gz 2>/dev/null | grep '^>' | head
zcat annotation.gtf.gz 2>/dev/null | awk '$0 !~ /^#/ {print $1; exit}'
```

If the genome uses names such as `1`, `2`, `3`, but the GTF uses `chr1`,
`chr2`, `chr3`, the files may not be compatible without preprocessing.

## Manual mode behavior

Manual mode is used when the Taxon ID is not resolvable from `HostSpecies.csv`.

In this mode, MTD Explorer requires:

- `--ncbi-taxon-id`;
- `--scientific-name`;
- `--genome`;
- `--gtf-file`;
- `--protein-fasta`.

When the Taxon ID is absent, MTD Explorer can append a minimal
`MANUAL_CUSTOM` entry to `HostSpecies.csv`. A backup of the CSV is created
before modification.

The manual row allows the custom host to be tracked by Taxon ID and scientific
name, but the reference files still come from the command-line inputs.

## Why protein FASTA is required

The protein FASTA is not optional for a complete host build.

MTD Explorer uses it to generate:

- representative proteins per gene;
- [eggNOG][eggnog] annotations;
- a reference-matched GO master GMT;
- optional custom [OrgDb][orgdb] resources.

This is why the genome FASTA, GTF, and protein FASTA should come from the same
source and release whenever possible.

## KEGG limitation and workaround

Some host species do not have a species-specific [KEGG][kegg] code or KEGG
organism name.

This is a limitation of KEGG species-level mapping, not a failure of the host
reference build.

In MTD Explorer, the host build does **not** depend only on a KEGG species code.
The custom host workflow builds reference-matched functional resources from the
protein FASTA and GTF using [eggNOG][eggnog]. It then creates a GO master GMT
that can be used for host gene-set analysis, including [ssGSEA][ssgsea].

Practical interpretation:

| Situation | What happens |
| --- | --- |
| KEGG code/name exists | KEGG-specific host outputs may be available when mappings are available |
| KEGG code/name is missing | KEGG-specific enrichment may be limited or absent |
| GO/eggNOG annotation is available | GO-based and custom GMT-based host gene-set analyses can still be generated |
| Protein FASTA and GTF are poor or mismatched | functional annotation quality may be incomplete |

So, lack of a KEGG name does not necessarily mean the host is unsupported. It
means KEGG-specific interpretation should be treated carefully, while
GO/eggNOG-derived resources provide the main fallback route.

## Useful options

### Skip only the R OrgDb build

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --skip-orgdb
```

This skips only construction of the custom R OrgDb package. The workflow still
generates representative proteins, eggNOG annotations, and the master GMT used
by ssGSEA.

### Force rebuild

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --force-orgdb
```

Use this when a compatible OrgDb appears to exist but you want to rebuild it.

### Clean previous host files only

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --clean-only
```

This removes previous files for that Taxon ID and exits before rebuilding.

### Keep previous files

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --no-clean
```

Use this only when you intentionally do not want to remove earlier files for the
same Taxon ID.

## Verify the host reference

After the build finishes, run:

```bash
bash MTD_check_installation.sh --full
```

Then run MTD Explorer using the same Taxon ID:

```bash
bash ~/MTD/MTD_explorer.sh \
  -i /path/to/samplesheet.csv \
  -o /path/to/output_directory \
  -h <host_taxon_id>
```

## What the host reference enables

A successful host reference setup enables or improves:

- host read mapping;
- host read removal before microbiome classification;
- host gene count matrix generation;
- host DEG analysis;
- host PCA and heatmaps;
- host GO/KEGG functional interpretation when mappings are available;
- ssGSEA host gene-set analysis;
- host-microbiome integration through [HAllA][halla].

## Common problems

### TaxID is not present in HostSpecies.csv

Use manual mode:

```bash
bash Create_custom_host.sh \
  --ncbi-taxon-id <host_taxon_id> \
  --scientific-name "Genus species" \
  --genome /path/to/genome.fa.gz \
  --gtf-file /path/to/annotation.gtf.gz \
  --protein-fasta /path/to/proteins.fa.gz
```

### TaxID is present, but automatic mode fails

Check whether the row has complete values for:

- `Genome_URL`;
- `GTF_URL`;
- `Pep_URL`;
- `Reference_status`.

If these fields are incomplete, use manual mode or update the curated
`HostSpecies.csv` row.

### Genome, GTF, and protein files are from different releases

Avoid this. Mismatched releases can reduce mapping quality, gene counting
accuracy, and functional annotation compatibility.

### KEGG outputs are missing

Check whether the host has a KEGG code/name and whether KEGG mappings were
available.

If not, use GO/eggNOG-derived results and the custom master GMT as the main
functional interpretation layer.

## Recommended reporting

In manuscripts or methods sections, report:

- host species scientific name;
- NCBI Taxon ID;
- whether automatic or manual custom host mode was used;
- genome assembly name;
- annotation source;
- annotation release;
- protein FASTA source;
- whether a custom OrgDb was built or skipped;
- whether KEGG species-specific mapping was available;
- MTD Explorer version or commit;
- date of host reference creation.

Also cite the original [MTD][mtd-paper] paper and the genome/annotation source.

[mtd-explorer]: https://github.com/patrick-douglas/MTD-Explorer
[mtd-paper]: https://doi.org/10.1093/bib/bbac111
[kraken2]: https://ccb.jhu.edu/software/kraken2/index.shtml
[hisat2]: https://daehwankimlab.github.io/hisat2/
[magic-blast]: https://ncbi.github.io/magicblast/
[eggnog]: http://eggnog-mapper.embl.de/
[orgdb]: https://bioconductor.org/packages/AnnotationForge/
[kegg]: https://www.genome.jp/kegg/
[ssgsea]: https://gsea-msigdb.github.io/ssGSEA-gpmodule/v10/
[halla]: https://github.com/biobakery/halla

[ensembl-downloads]: https://www.ensembl.org/info/data/ftp/index.html
[ensembl-genomes-downloads]: https://ftp.ensemblgenomes.ebi.ac.uk/pub/
[ncbi-datasets]: https://www.ncbi.nlm.nih.gov/datasets/
[ncbi-datasets-cli]: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
[ncbi-refseq]: https://www.ncbi.nlm.nih.gov/refseq/
[ncbi-genbank]: https://www.ncbi.nlm.nih.gov/genbank/
[ncbi-assembly]: https://www.ncbi.nlm.nih.gov/assembly/
[gencode]: https://www.gencodegenes.org/
[ucsc-downloads]: https://hgdownload.cse.ucsc.edu/downloads.html
[ensembl]: https://www.ensembl.org/
[ensembl-genomes]: https://ensemblgenomes.org/
[ncbi]: https://www.ncbi.nlm.nih.gov/
