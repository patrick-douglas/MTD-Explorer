# Custom microbiome references

[MTD Explorer][mtd-explorer] uses a [Kraken 2][kraken2] database for microbial
taxonomic classification and [Bracken][bracken] files for abundance
re-estimation.

The microbiome reference is independent from the host reference.

| Reference type | Main purpose | Main MTD Explorer option |
| --- | --- | --- |
| Host reference | host read mapping, host gene counts, host functional analyses | `-h`, `--kraken-host-db` |
| Microbiome reference | microbial classification and abundance estimation | `--kraken-micro-db` |

Use this page when you want to create or use a custom microbial database instead
of the default microbiome database.

The generic helper script is:

```bash
Create_custom_micro.sh
```

## When to create a custom microbiome database

Create a custom microbiome database when:

- the default microbiome database is not specific enough for the project;
- the study focuses on a targeted group such as viruses, bacteria, fungi,
  archaea, or protozoa;
- local genomes or curated FASTA files must be added;
- public microbial libraries need to be combined with project-specific
  sequences;
- a cleaner database is needed after removing contaminant, host, or unresolved
  TaxID records.

A custom microbiome database can be broad, such as bacteria plus viruses plus
fungi, or targeted, such as a virus-only database.

!!! important

    MTD Explorer classifies reads only against the database provided with
    `--kraken-micro-db`. If the database contains only viral sequences, the
    microbiome output should be interpreted as a viral profiling result, not as
    a complete bacteria-fungi-virus microbiome profile.

## Quick examples

### Broad microbial database

```bash
cd ~/MTD

bash Create_custom_micro.sh \
  --db-name kraken2DB_micro_custom \
  --libraries bacteria,viral,archaea,fungi,protozoa \
  --threads 20 \
  --bracken-read-len 75
```

### Virus-only database

```bash
cd ~/MTD

bash Create_custom_micro.sh \
  --db-name kraken2DB_viral_custom \
  --libraries viral \
  --threads 20 \
  --bracken-read-len 75
```

### Add a custom FASTA with TaxIDs already in the headers

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_project_micro \
  --add-fasta /path/to/custom_sequences.fa \
  --threads 20 \
  --bracken-read-len 75
```

### Add a custom FASTA where all records belong to one TaxID

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_project_viral \
  --add-fasta-with-taxid /path/to/project_viral_sequences.fa:10239 \
  --threads 20 \
  --bracken-read-len 75
```

In this example, all records are assigned to NCBI TaxID `10239`, the viral root
TaxID. Prefer species-level or genus-level TaxIDs when the sequence identity is
known more precisely.

### Mixed public and custom database

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_micro_plus_project \
  --libraries bacteria,viral \
  --add-fasta-with-taxid /path/to/project_virus.fa:10239 \
  --threads 20 \
  --bracken-read-len 75
```

## Use the custom database in MTD Explorer

After the database is built, point MTD Explorer to it:

```bash
bash ~/MTD/MTD_explorer.sh \
  -i /path/to/samplesheet.csv \
  -o /path/to/output_directory \
  -h <host_taxon_id> \
  --kraken-micro-db "$HOME/MTD/kraken2DB_micro_custom" \
  --bracken-read-len 75
```

The `--bracken-read-len` value should match the read length used when building
the Bracken files with `Create_custom_micro.sh`.

## What Create_custom_micro.sh does

The script can:

- create an empty [Kraken 2][kraken2] database directory;
- download or reuse [NCBI Taxonomy][ncbi-taxonomy];
- add public Kraken 2 libraries;
- add custom FASTA files;
- rewrite FASTA headers with `kraken:taxid|TAXID` when requested;
- build the Kraken 2 database;
- build [Bracken][bracken] files for the selected read length;
- validate that the expected database files exist;
- write an `MTD_custom_micro_manifest.txt` file inside the database directory.

A completed Kraken 2 database should contain:

```text
hash.k2d
opts.k2d
taxo.k2d
```

A completed Bracken database should also contain:

```text
database75mers.kmer_distrib
```

The number in `database75mers.kmer_distrib` changes with the Bracken read
length.

## Public Kraken 2 libraries

`Create_custom_micro.sh` can call `kraken2-build --download-library` for public
Kraken 2 libraries.

Common options include:

| Library | Typical use |
| --- | --- |
| `bacteria` | bacterial profiling |
| `viral` | viral profiling |
| `archaea` | archaeal profiling |
| `fungi` | fungal profiling |
| `protozoa` | protozoan profiling |
| `plasmid` | plasmid sequences |
| `UniVec_Core` | vector/contaminant screening |

Example:

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_micro_custom \
  --libraries bacteria,viral,archaea,fungi,protozoa,UniVec_Core \
  --threads 20 \
  --bracken-read-len 75
```

## Custom FASTA files and TaxIDs

[Kraken 2][kraken2] needs taxonomic information to classify reads.

For custom FASTA files, there are two common cases.

### Case 1: headers are already TaxID-aware

If the FASTA headers already include `kraken:taxid|TAXID`, use:

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_custom_taxid_headers \
  --add-fasta /path/to/sequences_with_taxids.fa
```

Example FASTA header:

```text
>sequence_001|kraken:taxid|10239
```

### Case 2: all records should receive the same TaxID

If all records in a FASTA belong to the same taxon, use:

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_custom_single_taxid \
  --add-fasta-with-taxid /path/to/sequences.fa:10239
```

The script will create a prepared FASTA inside the database directory with
`kraken:taxid|10239` added to each header.

!!! warning

    Avoid assigning a broad parent TaxID when a more precise species-level or
    genus-level TaxID is known. Broad TaxIDs can make downstream profiles less
    informative.

## FASTA list file

For many FASTA files, create a list:

```text
/path/to/bacteria.fa
/path/to/project_virus.fa	10239
/path/to/fungus.fa,4751
```

Then run:

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_micro_from_list \
  --fasta-list fasta_inputs.txt \
  --threads 20 \
  --bracken-read-len 75
```

Lines with only a FASTA path are added directly. Lines with a FASTA path plus
TaxID are rewritten with `kraken:taxid|TAXID`.

## Example: expanded viral database

During MTD Explorer development, a custom viral database was built as a targeted
example.

The example database combined:

- [RefSeq Viral][refseq-viral];
- [Virus-Host DB][virus-host-db];
- [C-RVDB v32][rvdb];
- [NCBI Taxonomy][ncbi-taxonomy].

The cleaned final FASTA used for the viral database contained:

```text
1,333,150 viral sequences
22,384 non-viral records excluded
0 unresolved TaxIDs
0 Homo sapiens TaxID 9606 records
```

The final database name used in testing was:

```text
kraken2DB_viral_extended_rvdb32_clean
```

Using the generic script, a cleaned viral FASTA can be turned into a database
with:

```bash
cd ~/MTD

bash Create_custom_micro.sh \
  --db-name kraken2DB_viral_extended_rvdb32_clean \
  --add-fasta "$HOME/MTD/viral_extended_rvdb32_clean.fa" \
  --threads 20 \
  --bracken-read-len 75
```

If the cleaned FASTA does not already contain TaxID-aware headers, use:

```bash
bash Create_custom_micro.sh \
  --db-name kraken2DB_viral_extended_rvdb32_clean \
  --add-fasta-with-taxid "$HOME/MTD/viral_extended_rvdb32_clean.fa:10239" \
  --threads 20 \
  --bracken-read-len 75
```

The earlier viral helper scripts used during development were:

```text
build_refseq_viral_taxid_map.py
build_virushost_taxid_map.py
build_nonredundant_viral_fasta.py
```

These were specific to the expanded viral database. `Create_custom_micro.sh` is
the generic database builder.

## Rebuild Bracken only

Sometimes the Kraken 2 database is already built, but the Bracken file is
missing or was built for the wrong read length.

Use:

```bash
bash Create_custom_micro.sh \
  --output-db "$HOME/MTD/kraken2DB_micro_custom" \
  --rebuild-bracken-only \
  --threads 20 \
  --bracken-read-len 75
```

This keeps the existing Kraken 2 database and rebuilds only:

```text
database75mers.kmer_distrib
```

Use the same read length in MTD Explorer:

```bash
--bracken-read-len 75
```

## Validate an existing custom database

```bash
bash Create_custom_micro.sh \
  --output-db "$HOME/MTD/kraken2DB_micro_custom" \
  --validate-only \
  --bracken-read-len 75
```

Manual checks:

```bash
DB="$HOME/MTD/kraken2DB_micro_custom"

ls -lh "$DB"/hash.k2d "$DB"/opts.k2d "$DB"/taxo.k2d
ls -lh "$DB"/database75mers.kmer_distrib

kraken2-inspect --db "$DB" | head
```

## Common problems

### Bracken file is missing

Rebuild Bracken with the same read length used by MTD Explorer:

```bash
bash Create_custom_micro.sh \
  --output-db "$DB" \
  --rebuild-bracken-only \
  --threads 20 \
  --bracken-read-len 75
```

Then run MTD Explorer with:

```bash
--bracken-read-len 75
```

### The database is too narrow

A virus-only database cannot detect bacterial or fungal taxa. Use a broader
database if the biological question requires broader microbiome profiling.

### Many reads are unclassified

Possible causes include:

- the database is too narrow for the sample type;
- important microbial groups are missing;
- custom FASTA files lack correct TaxIDs;
- references are too distant from organisms in the sample;
- read quality or trimming is too strict.

### Unexpected taxa appear

Check whether the custom FASTA headers have correct TaxIDs. Also check whether
sequences were assigned to broad parent TaxIDs instead of precise taxa.

### Host or human sequences appear in the microbial database

Remove host or contaminant sequences before building the database. In the
expanded viral database example, records assigned to `Homo sapiens` TaxID 9606
were excluded.

## Recommended reporting

In manuscripts or methods sections, report:

- database name;
- database purpose, such as broad microbiome or virus-only profiling;
- public Kraken 2 libraries used;
- custom FASTA sources;
- release or download date;
- filtering rules;
- number of final sequences when available;
- Kraken 2 version;
- Bracken version;
- Kraken 2 k-mer length;
- Bracken read length;
- MTD Explorer version or commit.

[mtd-explorer]: https://github.com/patrick-douglas/MTD-Explorer
[kraken2]: https://ccb.jhu.edu/software/kraken2/
[bracken]: https://ccb.jhu.edu/software/bracken/
[ncbi-taxonomy]: https://www.ncbi.nlm.nih.gov/taxonomy
[refseq-viral]: https://www.ncbi.nlm.nih.gov/genome/viruses/
[virus-host-db]: https://www.genome.jp/virushostdb/
[rvdb]: https://rvdb.dbi.udel.edu/
