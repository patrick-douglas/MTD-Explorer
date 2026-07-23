#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use Getopt::Std;
use List::Util qw/max/;

my $PROG = basename $0;

# Local folder where archaeal .gz files are stored
my $local_download_dir = "/media/me/4TB_BACKUP_LBN/Compressed/MTD/Kraken2DB_micro/library/archaea/all/";

# Folder where Kraken2 will read/copy the files
my $new_download_dir = getcwd() . "/all";

make_path($new_download_dir) unless -d $new_download_dir;

opendir(my $dh, $local_download_dir) || die "Can't opendir $local_download_dir: $!";
while (my $file = readdir($dh)) {
    next if $file =~ /^\./;
    next unless $file =~ /\.gz$/;

    my $source = "$local_download_dir/$file";
    my $destination = "$new_download_dir/$file";

    copy($source, $destination) or die "Copy failed: $source -> $destination: $!";
}
closedir($dh);

$local_download_dir = $new_download_dir;

my %manifest;

while (<>) {
    next if /^#/;
    chomp;

    my @fields = split /\t/;
    my ($taxid, $asm_level, $ftp_path) = @fields[5, 11, 19];

    next unless grep { $asm_level eq $_ } ("Complete Genome", "Chromosome");
    next if $ftp_path eq "na";

    my $filename = basename($ftp_path);

    my $local_path = "$local_download_dir/$filename.gz";
    my $alt_local_path = "$local_download_dir/${filename}_genomic.fna.gz";

    if (-e $local_path) {
        $manifest{$local_path} = $taxid;
    } elsif (-e $alt_local_path) {
        $manifest{$alt_local_path} = $taxid;
    } else {
        print STDERR "$PROG: Local file $local_path or $alt_local_path not found. Skipping.\n";
    }
}

open MANIFEST, ">", "manifest.txt"
    or die "$PROG: can't write manifest: $!\n";

print MANIFEST "$_\n" for keys %manifest;
close MANIFEST;

print STDERR "Manifest verification complete:\n";
print STDERR "  Found: " . scalar(keys %manifest) . "\n";

print STDERR "Step 1/2: Processing locally downloaded archaeal files\n";

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

    open IN, "gunzip -c $in_filename |"
        or die "$PROG: can't read $in_filename: $!\n";

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

    $projects_added++;

    my $out_line = progress_line(
        $projects_added,
        scalar keys %manifest,
        $sequences_added,
        $ch_added
    ) . "...";

    $max_out_chars = max(length($out_line), $max_out_chars);
    my $space_line = " " x $max_out_chars;

    print STDERR "\r$space_line\r$out_line" if -t STDERR;
}

close OUT;

print STDERR " done.\n" if -t STDERR;
print STDERR "All archaeal files processed, cleaning up...\n";
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
