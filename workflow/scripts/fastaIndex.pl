#!/usr/bin/env perl

use Bio::DB::Fasta;

my $fastaFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );
