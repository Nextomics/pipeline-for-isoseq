#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
my $input=shift; #all.besthit_longest.collapsed.rep.fq
#@PB.1.1|chrA01:2669-5604(-)|c1/f1p0/2210
my $input2=shift;# all.besthit_longest.result.sam
#c1/f1p0/2210    0       chrA01  2669    40      509M80N137M82N189M94N96M99N195M64N67M81N243M75N209M82N634M
my $input3=shift;#all.besthit_longest.result.bed
#chrA01  2668    3177    c1/f1p0/2210    40      +
my %hash;
open IN,$input;
while(<IN>){
	chomp;
	$_=~s/^@//;
	my @temp=split /\|/;
	$hash{$temp[2]}->[0]=$_;
	<IN>;<IN>;<IN>;
}
close IN;
open IN2,$input2;
while(<IN2>){
	chomp;
	next if(/^@/);
	my @temp=split;
	if(exists $hash{$temp[0]}){
#		$hash{$temp[0]}->[1]=$temp[5];
		my @cigar=split /\d+N/,$temp[5];
#		print Dumper \@cigar;
		
		for(my $i=0;$i<=$#cigar;$i++){

			$hash{$temp[0]}->[$i+1]=$cigar[$i];
		}
	}
}
#print Dumper \%hash;
open IN3,$input3;
my %reads;
while(<IN3>){
	chomp;
	my @temp=split;
	#chrA01  2668    3177    c1/f1p0/2210    40      +
	$reads{$temp[3]}->[0]++;
	my $length=(split /\//,$temp[3])[-1];
	my ($score,$start,$end,$gap_string);
#	$reads{$temp[3]}->[1]=0 if(!defined $reads{$temp[3]}->[1]);
#	$start=$reads{$temp[3]}->[1]+1;
	if(exists $hash{$temp[3]}){
		my $cigar=$hash{$temp[3]}->[$reads{$temp[3]}->[0]];
		my @M=($cigar=~/(\d+)M/g);
		my @I=($cigar=~/(\d+)I/g);
		my @D=($cigar=~/(\d+)D/g);
		$cigar=~s/(\d+)([SMIDH])/$2$1/g;
		my @all_tag=($cigar=~/([MID]\d+)/g);
		my $M_len=&sum(\@M);
		my $I_len=&sum(\@I);
		my $D_len=&sum(\@D);
#		my $end=$start+$M_len+$I_len-1;
		my $score=100-$I_len-$D_len;
		if($temp[5] eq '+'){
			if(!defined $reads{$temp[3]}->[1]){
				$reads{$temp[3]}->[1]=0;
#				if($all_tag[0]=~/[SH](\d+)/){
				if($cigar=~/^[SH](\d+)/){
					$reads{$temp[3]}->[1]=$1;
					
				}
			}
#			$all_tag[0]=~s/[SH](\d+)//g;
#			$all_tag[-1]=~s/[SH](\d+)//g;
			$start=$reads{$temp[3]}->[1]+1;
			$end=$start+$M_len+$I_len-1;
			$gap_string=join " ",@all_tag;
			print "$temp[0]\tcotton\tcDNA_match\t$temp[1]\t$temp[2]\t$score\t$temp[5]\t.\tID=$hash{$temp[3]}->[0].path1;Name=$hash{$temp[3]}->[0];Target=$hash{$temp[3]}->[0] $start $end;Gap=$gap_string\n";
		}else{
			if(!defined $reads{$temp[3]}->[1]){
				$reads{$temp[3]}->[1]=$length+1;
#				if($all_tag[0]=~/[SH](\d+)/){
				if($cigar=~/^[SH](\d+)/){
					$reads{$temp[3]}->[1]=$length-$1+1;
				}
			}
#			$all_tag[0]=~s/[SH](\d+)//g;
#			$all_tag[-1]=~s/[SH](\d+)//g;
			$start=$reads{$temp[3]}->[1]-1;
			$end=$start-$M_len-$I_len+1;
			@all_tag=reverse @all_tag;
			$gap_string=join " ",@all_tag;
			print "$temp[0]\tcotton\tcDNA_match\t$temp[1]\t$temp[2]\t$score\t$temp[5]\t.\tID=$hash{$temp[3]}->[0].path1;Name=$hash{$temp[3]}->[0];Target=$hash{$temp[3]}->[0] $end $start;Gap=$gap_string\n";
		}
		$reads{$temp[3]}->[1]=$end;
	}
}
#print Dumper \%reads;
close IN3;
sub sum{
	my $array_p=shift;
	my $sum=0;
	foreach my $num(@$array_p){
		$sum+=$num;
	}
	return $sum;
}
