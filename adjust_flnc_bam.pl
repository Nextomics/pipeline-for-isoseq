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
#my $csv_list="/share1/RNA/liangfan/mianhua_HZND_isoseq/pacbio_data/consensus/all_flnc.csv.list";
my $samlist=shift; #sam.list
my $result_dir=shift; #"/share1/RNA/liangfan/mianhua_HZND_isoseq/phase/split_bam/chr";
#my $fa_chr_len="/share1/RNA/liangfan/mianhua_HZND_isoseq/database/cotton.draft.genome.fa.len";
my $sam_file=shift; #"/share1/RNA/liangfan/mianhua_HZND_isoseq/phase/split_bam/all_merge_flnc.sam";
my @sam=`cat $samlist`;
chomp @sam;
#my @csv_list=`cat $csv_list`;
#chomp @csv_list;
our %hash;
foreach my $sam_file(@sam){
	my $prefix=(split /\//,$sam_file)[0];
	my $vcf_file="$result_dir/split/$prefix/$prefix.raw.vcf";
	my $bam_0_file="$result_dir/split/$prefix/$prefix.0.adjust.bam";
	my $bam_1_file="$result_dir/split/$prefix/$prefix.1.adjust.bam";
	my $bam_chimera="$result_dir/split/$prefix/$prefix.chimera.bam";
	my $phase_file="$result_dir/split/$prefix/$prefix.phased";
	my $snp_num_0=0;
	my $snp_num_1=1;
	&read_bam($bam_0_file,$snp_num_0);
	&read_bam($bam_1_file,$snp_num_1);
	&read_bam($bam_chimera,3);
	
	close IN;
#0 for eq 1 for other 3 for chimera 2 for not	
}
#print Dumper \%hash;

open SAM,$sam_file;
while(<SAM>){
	chomp;
	if(/^@/){
		print $_,"\n";
	}else{
		my @temp=split;
		if(exists $hash{$temp[0]}){
			if($hash{$temp[0]}->[0] eq $temp[2] && $hash{$temp[0]}->[1] eq $temp[3]){
				print $_,"\n";
			}
		}else{
			print $_,"\n";
		}
			
	}

}
close SAM;

sub overhang_len{
	my $string=shift;
    my @temp=split /\s+/,$string;
    my $cigar=$temp[5];
    my ($head,$tail)=(0,0);
    $head=$1 if($temp[5]=~/^(\d+)[S|H]/);
    $tail=$1 if($temp[5]=~/(\d+)[S|H]$/);
	my $len=$head+$tail;
	return $len;
}
sub read_bam{
    my $bam_file=shift;
    my $num=shift;

    open IN,"$samtools view -h $bam_file|";
    while(<IN>){
        chomp;
        next if(/^@/);
        my @temp=split;
        my $overhang=&overhang_len($_);
        if(!exists $hash{$temp[0]}){
            @{$hash{$temp[0]}}=($temp[2],$temp[3],$num,$overhang);
        }else{
            if($num<$hash{$temp[0]}->[2]){
                @{$hash{$temp[0]}}=($temp[2],$temp[3],$num,$overhang);
            }elsif($num==$hash{$temp[0]}->[2] && $overhang < $hash{$temp[0]}->[3]){
                @{$hash{$temp[0]}}=($temp[2],$temp[3],$num,$overhang);
            }

        }
    }
    close IN;
}
