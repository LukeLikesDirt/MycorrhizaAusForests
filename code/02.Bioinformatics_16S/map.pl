#!/usr/bin/perl                                                             

use warnings;
use strict;

die "Three arguments needed: fasta1, uc, and fasta2
 fasta1:  original fasta1 file       
 uc:      mapping
 fasta2:  fasta2 file with accepted sequences
\n" unless scalar @ARGV > 2;

my ($fasta1, $uc, $fasta2) = @ARGV;

# read fasta2 file with accepted sequences                                      

my %accepted = ();

open(F2, $fasta2) || die "FATAL ERROR:\n Unable to open fasta2 ($fasta2)\n";
while (<F2>) {

    if (/^>([^ ;]+)/) {
     	$accepted{$1} = 1;
    }
}
close F2;

# read uc file with mapping                                                     

open(UC, $uc) || die "FATAL ERROR:\n Unable to open uc ($uc)\n";

while (<UC>) {
    chomp;
    my @uc_file_columns = split /\t/;

    my $query_name;
    if ($uc_file_columns[8] =~ /^([^ ;*]+)/)  {
     	$query_name = $1;
    }

    my $target_name;
    if ($uc_file_columns[9] =~ /^([^ ;*]+)/)  {
     	$target_name = $1;
    }

    if ((defined $target_name) && ($accepted{$target_name}) && (defined $query_name)) {
        $accepted{$query_name} = 1;
    }
}
close UC;

# read original fasta1 file                                                     

my $ok = 0;

open(F1, $fasta1)  || die "FATAL ERROR:\n Unable to open fasta1 ($fasta1)\n";

while (<F1>) {

    if (/^>([^ ;]+)/) {
     	$ok = $accepted{$1};
    }
    print if $ok;
}

close F1;