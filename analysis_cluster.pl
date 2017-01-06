#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
my $input=shift;#group
my $input2=shift;#bam
my $input3=shift; #fq
#@c1/f1p1/811 isoform=c1;full_length_coverage=1;non_full_length_coverage=1;isoform_length=811
#PB.1.1  m160618_205236_42221_c100964032550000001823218007011634_s1_p0/161883/1010_59_CCS
#m160620_192404_42221_c100963832550000001823218007011606_s1_p0/15372/1317_58_CCS 0       C30057488       1       40      506S370M1D383M
my %hash;
open IN,$input;
my $cluster_num;
my %cluster;
while(<IN>){
	chomp;
	$cluster_num++;
	my @temp=split;
	my @reads = split /,/,$temp[1];
	my $reads_num=@reads;
	my $max_len=0;
	my $max_reads;
	foreach my $reads(@reads){
		my $temp_length=&cal_pb_reads_len($reads);
		if($max_len < $temp_length){
			$max_len = $temp_length;
			$max_reads=$reads;
		}
	}
	$hash{$max_reads}="c$cluster_num/f$reads_num"."p0/$max_len";
	$cluster{$cluster_num}=$temp[0];
}
close IN;
open OUT2,">id_relation.txt";
foreach my $key (sort {$a<=>$b} keys %cluster){
	print OUT2 "$key\t$cluster{$key}\n";
}
close OUT2;
open IN2,$input2;
while(<IN2>){
	chomp;
	if(/^@/){
		print $_,"\n";
	}else{
		my @temp=split;
		#m160620_192404_42221_c100963832550000001823218007011606_s1_p0/15372/1317_58_CCS 0       C30057488       1       40      506S370M1D383M
		if(exists $hash{$temp[0]}){
			$temp[0]=$hash{$temp[0]};
			my $line=join "\t",@temp;
			print $line,"\n";
		}
	}
}
close IN2;
open IN3,$input3;
open OUT,">chose.fq";
while(my $head=<IN3>){
	my $seq=<IN3>;
	my $line=<IN3>;
	my $line2=<IN3>;
	#@c1/f1p1/811 isoform=c1;full_length_coverage=1;non_full_length_coverage=1;isoform_length=811
	#@m160606_111855_42221_c100964102550000001823218007011637_s1_p0/11/3114_58_CCS
	my $reads_name=$1 if($head=~/^@(\S+)/);
	if(exists $hash{$reads_name}){
		my ($cluster,$fl,$nfl,$len)=($1,$2,$3,$4) if($hash{$reads_name}=~/(c\d+)\/f(\d+)p(\d+)\/(\d+)/);
		print OUT "\@$hash{$reads_name} isoform=$cluster;full_length_coverage=$fl;non_full_length_coverage=$nfl;isoform_length=$len\n";
		print OUT "$seq$line$line2";
	}	
}	
sub cal_pb_reads_len{
	my $reads=shift;
	my ($num1,$num2)=($1,$2) if($reads=~/(\d+)_(\d+)_CCS/);
	my $length = abs ($num1 -$num2);
}

