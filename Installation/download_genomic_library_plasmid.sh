#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Download specific genomic libraries for use with Kraken 2.
# Supported libraries were chosen based on support from NCBI's FTP site
#   in easily obtaining a good collection of genomic data.  Others may
#   be added upon popular demand.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

# Kraken2 supplies the database path through KRAKEN2_DB_NAME.
# It may already be an absolute path, so it must not be prefixed
# with the MTD installation directory.
if [[ -z "${KRAKEN2_DB_NAME:-}" ]]; then
    echo "Error: KRAKEN2_DB_NAME is not defined." >&2
    exit 1
fi

if [[ "$KRAKEN2_DB_NAME" = /* ]]; then
    KRAKEN2_DB_DIR="$KRAKEN2_DB_NAME"
else
    KRAKEN2_DB_DIR="$(readlink -m -- "$KRAKEN2_DB_NAME")"
fi

LIBRARY_DIR="$KRAKEN2_DB_DIR/library"

NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
FTP_SERVER="https://$NCBI_SERVER"
RSYNC_SERVER="rsync://$NCBI_SERVER"
THIS_DIR=$PWD

library_name="$1"
ftp_subdir=$library_name
library_file="library.fna"
if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
  library_file="library.faa"
fi

function download_file() {
  file="$1"
  if [ -n "$KRAKEN2_USE_FTP" ]
  then
    wget ${FTP_SERVER}${file}
  else
    rsync --no-motd ${RSYNC_SERVER}${file} .
  fi
}

case $library_name in
  "archaea" | "bacteria" | "viral" | "fungi" | "plant" | "human" | "protozoa")
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    rm -f assembly_summary.txt
    remote_dir_name=$library_name
    if [ "$library_name" = "human" ]; then
      remote_dir_name="vertebrate_mammalian/Homo_sapiens"
    fi
    if ! download_file "/genomes/refseq/$remote_dir_name/assembly_summary.txt"; then
      1>&2 echo "Error downloading assembly summary file for $library_name, exiting."
      exit 1
    fi
    if [ "$library_name" = "human" ]; then
      grep "Genome Reference Consortium" assembly_summary.txt > x
      mv -f x assembly_summary.txt
    fi
    rm -rf all/ library.f* manifest.txt rsync.err
    rsync_from_ncbi.pl assembly_summary.txt
    scan_fasta_file.pl $library_file >> prelim_map.txt
    ;;
"plasmid")
    mkdir -p "$LIBRARY_DIR/plasmid"
    cd "$LIBRARY_DIR/plasmid"

    rm -f \
        library.fna \
        library.faa \
        manifest.txt \
        prelim_map.txt \
        plasmid.* \
        .listing

    local_download_dir="${MTD_KRAKEN2_PLASMID_CACHE:-}"

    if [[ -z "$local_download_dir" ]]; then
        echo "Error: MTD_KRAKEN2_PLASMID_CACHE is not defined." >&2
        echo "The MTD installer must export the plasmid cache path." >&2
        exit 1
    fi

    if [[ ! -d "$local_download_dir" ]]; then
        echo "Error: plasmid cache directory not found:" >&2
        echo "  $local_download_dir" >&2
        exit 1
    fi

    echo "Using cached plasmid files from:"
    echo "  $local_download_dir"

    mapfile -d '' plasmid_files < <(
        find "$local_download_dir" \
            -maxdepth 1 \
            -type f \
            -name '*.genomic.fna.gz' \
            -size +0c \
            -print0 |
        sort -z
    )

    if (( ${#plasmid_files[@]} == 0 )); then
        echo "Error: no non-empty *.genomic.fna.gz files found in:" >&2
        echo "  $local_download_dir" >&2
        exit 1
    fi

    echo "Plasmid files found: ${#plasmid_files[@]}"

    for plasmid_file in "${plasmid_files[@]}"; do
        if ! gzip -t "$plasmid_file"; then
            echo "Error: invalid gzip file:" >&2
            echo "  $plasmid_file" >&2
            exit 1
        fi
    done

    printf '%s\n' "${plasmid_files[@]##*/}" > manifest.txt

    if [[ ! -s manifest.txt ]]; then
        echo "Error: manifest.txt was not created." >&2
        exit 1
    fi

    : > "$library_file"

    current_file=0
    total_files="${#plasmid_files[@]}"

    for plasmid_file in "${plasmid_files[@]}"; do
        current_file=$((current_file + 1))

        printf '[%d/%d] Processing %s\n' \
            "$current_file" \
            "$total_files" \
            "$(basename "$plasmid_file")"

        if ! gzip -cd -- "$plasmid_file" >> "$library_file"; then
            echo "Error while extracting:" >&2
            echo "  $plasmid_file" >&2
            rm -f "$library_file"
            exit 1
        fi
    done

    if [[ ! -s "$library_file" ]]; then
        echo "Error: $library_file is empty." >&2
        exit 1
    fi

    if ! grep -m 1 -q '^>' "$library_file"; then
        echo "Error: no FASTA header found in $library_file." >&2
        exit 1
    fi

    if ! scan_fasta_file.pl "$library_file" > prelim_map.txt; then
        echo "Error: scan_fasta_file.pl failed." >&2
        exit 1
    fi

    if [[ ! -s prelim_map.txt ]]; then
        echo "Error: prelim_map.txt is empty." >&2
        exit 1
    fi

    echo "Plasmid library prepared successfully."
    echo "Library FASTA: $LIBRARY_DIR/plasmid/$library_file"
    echo "Manifest:      $LIBRARY_DIR/plasmid/manifest.txt"
    echo "Preliminary map: $LIBRARY_DIR/plasmid/prelim_map.txt"
    ;;
  "nr" | "nt")
    protein_lib=0
    if [ "$library_name" = "nr" ]; then
      protein_lib=1
    fi
    if (( protein_lib == 1 )) && [ -z "$KRAKEN2_PROTEIN_DB" ]; then
      1>&2 echo "$library_name is a protein database, and the Kraken DB specified is nucleotide"
      exit 1
    fi
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    rm -f $library_name.gz
    echo -n "Downloading $library_name database from server... "
    download_file "/blast/db/FASTA/$library_name.gz"
    echo "done."
    echo -n "Uncompressing $library_name database..."
    gunzip $library_name.gz
    mv $library_name $library_file
    echo "done."
    echo -n "Parsing $library_name FASTA file..."
    # The nr/nt files tend to have non-standard sequence IDs, so
    # --lenient is used here.
    scan_fasta_file.pl --lenient $library_file >> prelim_map.txt
    echo "done."
    ;;
  "UniVec" | "UniVec_Core")
    if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
      1>&2 echo "$library_name is for nucleotide databases only"
      exit 1
    fi
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    echo -n "Downloading $library_name data from server... "
    download_file "/pub/UniVec/$library_name"
    echo "done."
    # 28384: "other sequences"
    special_taxid=28384
    echo -n "Adding taxonomy ID of $special_taxid to all sequences... "
    sed -e "s/^>/>kraken:taxid|$special_taxid|/" $library_name > library.fna
    scan_fasta_file.pl library.fna > prelim_map.txt
    echo "done."
    ;;
  *)
    1>&2 echo "Unsupported library.  Valid options are: "
    1>&2 echo "  archaea bacteria viral fungi plant protozoa human plasmid"
    1>&2 echo "  nr nt UniVec UniVec_Core"
    exit 1
    ;;
esac

if [ -n "$KRAKEN2_MASK_LC" ]; then
  echo -n "Masking low-complexity regions of downloaded library..."
  mask_low_complexity.sh .
  echo " done."
fi

