#!/usr/bin/perl

use warnings;
use strict;

#### inspired by https://www.biostars.org/p/224372/

my $input  = "GCF_001465895.1_Nfu_20140520_genomic_CLEAN_MT_exon.filtered.gtf";

open(GTF,$input);  # Input file

my %gene;  # hash to store gene->transcript->exon information
my %genepos; # hash to store gene's start and end position
my %transpos; # hash to store transcript's start and end position
my %geneinfo; # hash to store information details of gene like chr, strand etc

# Parse input file and create all the data structures in the form of hash
while (my $line = <GTF>){
	
	my @linedata = get_line_data ($line);

	# only derive gene from transcript
	unless ($linedata[2] eq "transcript") {
		next;
	}
	
	# Fetch gene id and name
	# NC_029649.1	Gnomon	transcript	4898	59338	.	+	.	transcript_id "rna-XM_015955279.1"; gene_id "gene-LOC107382895"; gene_name "LOC107382895"; Dbxref "GeneID:107382895,Genbank:XM_015955279.1"; Name "XM_015955279.1"; Note "The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred CDS: added 121 bases not found in genome assembly"; exception "annotated by transcript or proteomic data"; gbkey "mRNA"; gene "LOC107382895"; inference "similar to RNA sequence (same species):INSD:GAIB01200069.1"; model_evidence "Supporting evidence includes similarity to: 4 ESTs%2C 12 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 2 samples with support for all annotated introns"; partial "true"; product "matrix-remodeling-associated protein 5-like%2C transcript variant X1"; start_range ".,4898"; transcriptID "XM_015955279.1"; gene_biotype "protein_coding"; CDS_Dbxref "GeneID:107382895,Genbank:XP_015810765.1"; CDS_Name "XP_015810765.1"; CDS_gbkey "CDS"; CDS_product "matrix-remodeling-associated protein 5-like isoform X1"; protein_id "XP_015810765.1"; locus "RLOC_00000001";
	$line =~ m/gene_id \"(.+)\"; gene_name \"(.+)\"; Dbxref/;
	my $geneid=$1;
	#print $geneid."\n";
	
	my $genename=$2;
	#print $genename."\n";
	
	# target formating
	#chr1	HAVANA	gene	3073253	3074322	.	+	.	gene_id "ENSMUSG00000102693.1"; gene_type "TEC"; gene_name "4933401J01Rik"; level 2; havana_gene "OTTMUSG00000049935.1";

	my $gene_gtf =  join("\t",(@linedata[0..1],"gene",@linedata[3..7],qq(gene_id "$geneid"; gene_type "protein_coding"; gene_name "$genename";)))  ;
	print $gene_gtf."\n";
	#push @{$geneinfo{$geneid}},@temp[0,1,6];
}

close GTF;

exit;



###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}
