#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Reads an assembly_summary.txt file, which indicates taxids and FTP paths for
# genome/protein data.  Performs the download of the complete genomes from
# that file, decompresses, and explicitly assigns taxonomy as needed.

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use Getopt::Std;
use List::Util qw/max/;

my $PROG = basename $0;

# Specify the original directory where your downloaded .gz files are stored
my $local_download_dir = "/media/me/4TB_BACKUP_LBN/Compressed/MTD/Kraken2DB_micro/library/bacteria/all/";

# Specify the new directory where you want to copy the files
my $new_download_dir = getcwd() . "/all";

# Create the new directory if it doesn't exist
make_path($new_download_dir) unless -d $new_download_dir;

# Copy all files from the original directory to the new directory
opendir(my $dh, $local_download_dir) || die "Can't opendir $local_download_dir: $!";
while (my $file = readdir($dh)) {
    next if $file =~ /^\./;  # skip hidden files/directories
    my $source = "$local_download_dir/$file";
    my $destination = "$new_download_dir/$file";
    next unless -f $source;

# Skip files that have already been copied completely.
if (-f $destination && -s $destination == -s $source) {
    next;
}

copy($source, $destination)
    or die "$PROG: Copy failed from $source to $destination: $!\n";
}
closedir($dh);

# Update $local_download_dir to the new directory
$local_download_dir = $new_download_dir;

my $ignored_non_genome = 0;
my $missing_genomes = 0;

# Manifest hash maps filenames (keys) to taxids (values)
my %manifest;
while (<>) {
    next if /^#/;
    chomp;
    my @fields = split /\t/;
    my ($taxid, $asm_level, $ftp_path) = @fields[5, 11, 19];
    # Possible TODO - make the list here configurable by user-supplied flags
    next unless grep {$asm_level eq $_} ("Complete Genome", "Chromosome");
    next if $ftp_path eq "na";  # Skip if no provided path

    my $filename = basename($ftp_path);

# Ignore malformed or non-assembly FTP entries such as "identical".
# Valid NCBI assembly directories normally begin with GCF_ or GCA_.
if ($filename !~ /^GC[AF]_\d+\.\d+/) {
    $ignored_non_genome++;
    next;
}

# Different NCBI assembly names require different local filename forms.
# In particular, some assembly directory names already end in "_genomic".
my @candidate_paths = (
    "$local_download_dir/$filename.gz",
    "$local_download_dir/$filename.fna.gz",
    "$local_download_dir/${filename}_genomic.fna.gz",
);

my ($local_path) = grep { -e $_ } @candidate_paths;

if (defined $local_path) {
    $manifest{$local_path} = $taxid;
}
else {
    $missing_genomes++;

    print STDERR "$PROG: Valid bacterial genome file not found locally:\n";
    print STDERR "  $_\n" for @candidate_paths;
    print STDERR "  Skipping this genome.\n";
    }
}

open MANIFEST, ">", "manifest.txt"
    or die "$PROG: can't write manifest: $!\n";
print MANIFEST "$_\n" for keys %manifest;
close MANIFEST;

print STDERR "Manifest verification complete:\n";
print STDERR "  Found valid genomes: " . scalar(keys %manifest) . "\n";
print STDERR "  Ignored non-genome entries: $ignored_non_genome\n";
print STDERR "  Missing valid genome files: $missing_genomes\n";
#print STDERR "  Not found: " . (58075 - scalar(keys %manifest)) . "\n";
#print STDERR "  Replaced: 0\n";

print STDERR "Step 1/2: Processing locally downloaded files\n";
my $output_file = "library.fna";
open OUT, ">", $output_file
    or die "$PROG: can't write $output_file: $!\n";
my $projects_added = 0;
my $sequences_added = 0;
my $ch_added = 0;
my $ch = "bp";
my $max_out_chars = 0;
for my $in_filename (keys %manifest) {
    my $taxid = $manifest{$in_filename};
    open IN, "gunzip -c $in_filename |" or die "$PROG: can't read $in_filename: $!\n";
    while (<IN>) {
        if (/^>/) {
            s/^>/>kraken:taxid|$taxid|/;
            $sequences_added++;
        } else {
            $ch_added += length($_) - 1;
        }
        print OUT;
    }
    close IN;
    # Remove the unlink command to preserve the original files
    # unlink $in_filename;
    $projects_added++;
    my $out_line = progress_line($projects_added, scalar keys %manifest, $sequences_added, $ch_added) . "...";
    $max_out_chars = max(length($out_line), $max_out_chars);
    my $space_line = " " x $max_out_chars;
    print STDERR "\r$space_line\r$out_line" if -t STDERR;
}
close OUT;
print STDERR " done.\n" if -t STDERR;

print STDERR "All files processed, cleaning up...\n";
print STDERR " done, library complete.\n";

sub progress_line {
    my ($projs, $total_projs, $seqs, $chs) = @_;
    my $line = "Processed ";
    $line .= ($projs == $total_projs) ? "$projs" : "$projs/$total_projs";
    $line .= " project" . ($total_projs > 1 ? 's' : '') . " ";
    $line .= "($seqs sequence" . ($seqs > 1 ? 's' : '') . ", ";
    my $prefix;
    my @prefixes = qw/k M G T P E/;
    while (@prefixes && $chs >= 1000) {
        $prefix = shift @prefixes;
        $chs /= 1000;
    }
    if (defined $prefix) {
        $line .= sprintf '%.2f %s%s)', $chs, $prefix, $ch;
    } else {
        $line .= "$chs $ch)";
    }
    return substr($line, 0, 79);
}

