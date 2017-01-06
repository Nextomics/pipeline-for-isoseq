#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Soft;
use FindBin qw($Bin $Script);
use lib $Bin;
my $config="$Bin/program.config.txt";
our $samtools=parse_config($config,"samtools");
my $bamfile0=shift;
my $bamfile1=shift;
my $phasefile=shift;
my $vcf_file=shift;
my $prefix=shift;
my %block;
open IN,$phasefile;
my $block_name;
my @merge;
while(<IN>){
	chomp;
	my @temp=split;
	if(/^PS/){
		$block_name="$temp[1]\t$temp[2]\t$temp[3]";
		push @merge,[$temp[1],$temp[2],$temp[3]];
	}
	if(/^M1/){
		$block{$temp[3]}=$block_name;
	}
}
my %new_merge;
my ($new_chr,$new_start,$new_end);

if($#merge==0){
	$new_chr=$merge[0]->[0];
	$new_start=1;
	$new_end=$merge[0]->[2]+10000;
	$new_merge{"$merge[0]->[0]\t$merge[0]->[1]\t$merge[0]->[2]"}="$new_chr\t$new_start\t$new_end";

}else{
	for(my $i=0;$i<=$#merge;$i++){
		if($i==0){
			$new_chr=$merge[$i]->[0];
			$new_start=1;
			$new_end=int (($merge[$i]->[2]+$merge[$i+1]->[1])/2);
		}elsif($i>0 && $i<$#merge){
			$new_chr=$merge[$i]->[0];
			$new_start=int (($merge[$i-1]->[2]+$merge[$i]->[1])/2)+1;
			$new_end=int (($merge[$i]->[2]+$merge[$i+1]->[1])/2);
		}elsif($i==$#merge){
			$new_chr=$merge[$i]->[0];
			$new_start=int (($merge[$i-1]->[2]+$merge[$i]->[1])/2)+1;
			$new_end=$merge[$i]->[2]+10000;
	
		}
		$new_merge{"$merge[$i]->[0]\t$merge[$i]->[1]\t$merge[$i]->[2]"}="$new_chr\t$new_start\t$new_end";
	}
}
close IN;
#$block{$temp[3]}=$block_name;
foreach my $pos (keys %block){
	$block{$pos}=$new_merge{$block{$pos}};
}
#print Dumper \%block;
open VCF,$vcf_file;
my %result;
while(<VCF>){
	chomp;
	next if(/^#/);
	my @temp=split;
	if(exists $block{$temp[1]}){
		if($temp[-2]!~/0\/0/ && $temp[-2]!~/\.\/\./  ){
			$result{$block{$temp[1]}}->[0]++;
		}
		if($temp[-1]!~/0\/0/ && $temp[-1]!~/\.\/\./ ){		
			$result{$block{$temp[1]}}->[1]++;
		}
	}	
}
close VCF;
#print Dumper \%result;
system ("$samtools view -H $bamfile0 > $prefix/aln0.sam");
system ("$samtools view -H $bamfile1 > $prefix/aln1.sam");
foreach my $block_region (keys %result){
	open TEMP,">temp.bed.txt";
	print TEMP $block_region;
	$result{$block_region}->[0]=0 if(!defined $result{$block_region}->[0]);
	$result{$block_region}->[1]=0 if(!defined $result{$block_region}->[1]);
	if($result{$block_region}->[0]<$result{$block_region}->[1]){
		
		system ("$samtools view -L temp.bed.txt $bamfile0 >> $prefix/aln0.sam");
		system ("$samtools view -L temp.bed.txt $bamfile1 >> $prefix/aln1.sam");
	}else{
		system ("$samtools view -L temp.bed.txt $bamfile0 >> $prefix/aln1.sam");
		system ("$samtools view -L temp.bed.txt $bamfile1 >> $prefix/aln0.sam");
	}	
	close TEMP;

#	print $block_region,"\t",$result{$block_region},"\n";
}
$bamfile0=~s/bam/adjust.bam/;
$bamfile1=~s/bam/adjust.bam/;
system ("$samtools view -bS $prefix/aln0.sam > $bamfile0 ");
system ("$samtools view -bS $prefix/aln1.sam > $bamfile1 ");
system ("rm $prefix/aln0.sam");
system ("rm $prefix/aln1.sam");
