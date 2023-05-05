#! /usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;


unless (scalar @ARGV == 2) {
	die "\nIncorrect command line arguments.\n".
		"Usage : process_repeat_library_for_10x_v3.pl <Repeat library name> <Repeat library fasta> \n\n";
}

my $repeatLibName = shift @ARGV;
my $sequence_file = shift @ARGV;

###################################
# 1. Prepare output files

#get base name to create new file name
$sequence_file =~ /(.+)\.fa/;

# use captured pattern to create output file name
my $gtf_out = $1.".gtf";
my $fa_out = $1.".clean.fa";

open (OUTGTF,'>',$gtf_out) or die "Couldn't write to $gtf_out file.: $!\n";
open (OUTFA,'>',$fa_out) or die "Couldn't write to $fa_out file.: $!\n";
###################################


###################################
# 2. Go through file to see components

my $inseq = Bio::SeqIO->new('-file' => "<$sequence_file",'-format' => 'fasta' ) ;

my $gid = 0;

# go through the fasta names
while (my $seq_obj = $inseq->next_seq ) {

	my $scaffold_id  = $seq_obj->id;
	my $scaffold_seq = $seq_obj->seq;
	my $desc         = $seq_obj->desc();
	my $scaffold_length = length($scaffold_seq);
	
	print $desc."\n";
	
	## seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
	## source - name of the program that generated this feature, or the data source (database or project name)
	## feature - feature type name, e.g. Gene, Variation, Similarity
	## start - Start position of the feature, with sequence numbering starting at 1.
	## end - End position of the feature, with sequence numbering starting at 1.
	## score - A floating point value.
	## strand - defined as + (forward) or - (reverse).
	## frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
	## attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
	
	my $seqname         =   $scaffold_id."_".$desc   ;
	
	# clean up repeat name if needed (if comment or slash present in name)
	$seqname            =~ s/\#/_/g;
	$seqname            =~ s/\//_/g;
	$seqname            =~ s/-//g;
	$seqname            =~ s/\s/_/g;

	my $source          =   $repeatLibName  ;
	my @features        =   ("gene","transcript","exon")   ;
	my $start           =   1  ;
	my $end             =   $scaffold_length  ;
	my $score           =   "."   ;
	my $strand          =   "+"   ;
	my $frame           =   "."   ;
	my $attribute_gene  =   qq/gene_id "Repeat_$gid"; gene_version "1"; gene_name "$seqname"; gene_source "$source"; gene_biotype "Transposon";/   ;
	my $attribute_trans =   qq/gene_id "Repeat_$gid"; gene_version "1"; transcript_id "Repeat_transcript_$gid"; transcript_version "1"; gene_name "$seqname"; gene_source "$source"; gene_biotype "Transposon"; transcript_name "$seqname"; transcript_source "$source"; transcript_biotype "protein_coding"; tag "basic"; transcript_support_level "1";/   ;
    my $attribute_exon  =   qq/gene_id "Repeat_$gid"; gene_version "1"; transcript_id "Repeat_transcript_$gid"; transcript_version "1"; exon_number "1"; gene_name "$seqname"; gene_source "$source"; gene_biotype "Transposon"; transcript_name "$seqname"; transcript_source "$source"; transcript_biotype "protein_coding"; exon_id "Repeat_exon_1"; exon_version "1"; tag "basic"; transcript_support_level "1";/   ;
	
	print OUTGTF join("\t",($scaffold_id,$source,"gene",       $start,$end,$score,$strand,$frame,$attribute_gene))."\n";
	print OUTGTF join("\t",($scaffold_id,$source,"transcript", $start,$end,$score,$strand,$frame,$attribute_trans))."\n";
	print OUTGTF join("\t",($scaffold_id,$source,"exon",       $start,$end,$score,$strand,$frame,$attribute_exon))."\n";
	
	$scaffold_seq = uc($scaffold_seq); #demask

	print OUTFA ">".$scaffold_id."\n";
	print OUTFA $scaffold_seq."\n";

	
	# increment before looping
	++$gid;

}# end while

close OUTGTF;
close OUTFA;

exit;