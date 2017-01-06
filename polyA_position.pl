#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
die "perl $0 <collapsed.gff> <collapsed.rep.fq>  <flnc alignment result>\n" unless(@ARGV==3);
#gene	strand	aligned reads	num sites	locations
#GBSCAFFOLD17453.23	+	1	0
my $input=shift;
#orginal_gene_id "PB.1";
#C30067912       PacBio  transcript      148     1097    .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
#C30067912       PacBio  exon    148     1097    .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
my $input2=shift;
#PB.1.1|C30067912:148-1097(+)|c1/f1p0/951
my $input3=shift;
#c1/f1p0/951     0       C30067912       148     3       105M1I845M

open IN,$input;
my %hash;
while(<IN>){
	chomp;
	my @temp=split /\t/;
	if($temp[2]=~/transcript/){
		my $transcript_id=$1 if($temp[8]=~/transcript_id "(\S+)"/);
		my $orginal_gene_id=$1 if($temp[8]=~/orginal_gene_id "(\S+)"/);
		$hash{$transcript_id}->[0]=$temp[4];
		$hash{$transcript_id}->[1]=$temp[6];
		$hash{$transcript_id}->[2]=$orginal_gene_id;
	}
}
close IN;
my %relation;
open IN2,$input2;
while(<IN2>){
	chomp;
	my @temp=split /\|/;
	my $transcript_id=$temp[0];
	$transcript_id=~s/^@//;
	my $reads=$temp[-1];
#	print "$transcript_id\t$reads\n";
	$relation{$reads}=$transcript_id;
	<IN2>;
	<IN2>;
	<IN2>;
	
}
close IN2;
open IN3,$input3;
while(<IN3>){
	chomp;
	next if(/^@/);
	my @temp=split;
	my $cigar=$temp[5];
	my $S_num;
	if($cigar=~/(\d+)[S|H]/){
		$S_num=$1;
	}else{
		$S_num=0;
	}
	next if(!exists $relation{$temp[0]});
	if($S_num<10 && exists $hash{$relation{$temp[0]}}){
		print $relation{$temp[0]},"\t",$hash{$relation{$temp[0]}}->[1],"\t",$temp[2],"\t",$hash{$relation{$temp[0]}}->[0],"\t",$hash{$relation{$temp[0]}}->[2],"\n";
	}
}
close IN3;
