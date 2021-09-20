#!/usr/bin/env perl

use Bio::DB::Fasta;

my $fastaFile = shift;
my $queryFile = shift;
my $gffFile = shift;
my $outFile = shift;

my @contigs = ();

open (IN, $queryFile);
while (<IN>){
    next if (/^[^>]/);
    chomp;
    my @f = split(/>/);
#    print $f[1], "\n";
    push @contigs, $f[1];
}
close(IN);

my %bin = map { $_ => 1 } @contigs;
 
my $db = Bio::DB::Fasta->new( $fastaFile );
open (OUT, '>', $outFile);
open (IN, $gffFile);
while (<IN>){
    next if (/^\#/);
    chomp;
    my @f = split(/\t/);
    if ($f[2] eq "CDS" && exists $bin{$f[0]}) { #0: contig 1:source 2:feature 3:start 4:end 5:. 6:sense 7:score 8:attribute
        my ($seq) = ($f[8] =~ /ID=([^;]+);/);
        my $sequence = $db->seq($seq);
        if  (!defined( $sequence )) {
            print STDERR "Sequence $seq not found. \n";
            next;
        }   
        print OUT ">$seq\n", "$sequence\n";        
   }
}
close(IN);
close(OUT);
