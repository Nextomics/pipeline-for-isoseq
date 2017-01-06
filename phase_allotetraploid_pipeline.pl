#! /usr/bin/perl -w
=head1 NAME:
	phase_allotetraploid_pipeline.pl -- A pipeline for phasing of alloterraploid plants sequenced by SMRT platform.

=head1 USAGE:
perl phase_allotetraploid_pipeline.pl
	--flnc Flnc reads generated from smrtanalysis
	--gmap_genome_directory  gmap genome directory
	--gmap_genome_database gmap genome database
	--reference_fasta referece fasta
	--outdir oudir
=cut
use strict;
use FindBin qw($Bin $Script);
use lib $Bin;
use Soft;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
#liangfan@nextomics.org
our ($flnc,$g_dir,$g_database,$ref_fasta);
our $program_start_time=`date`;
our $outdir;
$program_start_time=~s/\s+/_/g;
my $cmd="perl $0 ".join " ",@ARGV;
GetOptions(
	"flnc:s"=>\$flnc,
	"gmap_genome_directory:s"=>\$g_dir,
	"gmap_genome_database:s"=>\$g_database,
	"reference_fasta:s"=>\$ref_fasta,
	"outdir:s"=>\$outdir,
);

die `pod2text $0` unless (defined $flnc && defined $g_dir && defined $g_database && defined $outdir && defined $ref_fasta);


my $config="$Bin/program.config.txt";
#get program
our $gmap=parse_config($config,"gmap");
our $samtools=parse_config($config,"samtools");
our $bcftools=parse_config($config,"bcftools");
#step 0 change dir
#&check_dir($outdir);
$outdir= abs_path($outdir);
$g_dir=abs_path($g_dir);
$flnc=abs_path($flnc);
$ref_fasta=abs_path($ref_fasta);
&creat_log($cmd,"CMD");
chdir $outdir;

#$cmd="cd $outdir\n";
#step 1 mapping
$cmd="$gmap -D $g_dir -d $g_database -f samse -n 2 -t 24 $flnc > $outdir/flnc.sam 2> $outdir/gmap.log\n";
$cmd.="$samtools sort $outdir/flnc.sam $outdir/flnc.sort\n";
$cmd.="$samtools view -h $outdir/flnc.sort.bam > $outdir/flnc.sort.sam\n";
$cmd.="rm $outdir/flnc.sort.bam\n";

#&monitor_jobs ($cmd,"STEP1 MAPPING");

#step 2 prepare bam
&check_dir("split");
chdir "split";
$cmd="perl $Bin/prepare_bam_for_phase.pl $outdir/flnc.sort.sam";

&monitor_jobs ($cmd,"STEP2 PREPARING");

#step 3 phasing
my $sam_list="$outdir/split/sam.list";
#index reference
#&check_dir("$outdir/ref");
#$cmd = "ln -s $ref_fasta $outdir/ref/$ref_fasta_file_name\n";
#$cmd.="samtools faidx $outdir/ref/$ref_fasta_file_name\n";
$cmd=&phase_cmd($sam_list);
&monitor_jobs ($cmd,"STEP3 PHASING");

#step 4 adjust alignments
$cmd="perl $Bin/adjust_flnc_bam.pl $sam_list $outdir $outdir/flnc.sort.sam > $outdir/flnc_adjust.sam";
&monitor_jobs ($cmd,"STEP4 ADJUSTALIGNMENTS");

#step 5 collapsed and consensus
#/home/xiaoyh/software/python/bin/collapse_isoforms_by_sam.py -c 0.95 --input fastq/all.fastq --fq -s all.sort.besthit.sam -o all
#method 1 (longest)
#perl chose_longest.pl all.collapsed.group.txt  all.sort.besthit.sam fastq/all.fastq > all.besthit_longest.result.sam
#method 2 (dagcon)
#/home/xiaoyh/software/python/bin/collapse_isoforms_by_sam.py -c 0.95 --input chose.fq --fq -s all.besthit_longest.result.sam -o all.besthit_longest
#
#step 6 quiver

#step 8 alternative splicing
sub phase_cmd{
	my $list=shift;
	my @list=`cat $list`;
	chomp @list;
	my $cmd;
	foreach my $sam_file(@list){
		my $bam_file=$sam_file;
		$bam_file=~s/sam/bam/;
		my $prefix=(split /\//,$sam_file)[0];
#		if(!defined $cmd){
#			$cmd="$samtools phase -Q 5 -b $outdir/split/$prefix/$prefix $outdir/split/$sam_file > $outdir/split/$prefix/$prefix.phased\n";
#		}else{
#			$cmd.="$samtools phase -Q 5 -b  $outdir/split/$prefix/$prefix $outdir/split/$sam_file > $outdir/split/$prefix/$prefix.phased\n";
#		}
		##call snp
		if(!defined $cmd){
			$cmd="$samtools view -bS $outdir/split/$sam_file > $bam_file\n";
		}else{
			$cmd.="$samtools view -bS $outdir/split/$sam_file > $bam_file\n";
		}
		$cmd.="$samtools phase -Q 5 -b  $outdir/split/$prefix/$prefix $outdir/split/$bam_file > $outdir/split/$prefix/$prefix.phased\n";
		$cmd.="$samtools mpileup -t DP,SP -m 2 -ugf $ref_fasta $outdir/split/$prefix/$prefix.0.bam $outdir/split/$prefix/$prefix.1.bam | $bcftools call -vmO v -o $outdir/split/$prefix/$prefix.raw.vcf\n";
		##ajust phased bam
		$cmd.="perl $Bin/adjust_phased_bam.pl $outdir/split/$prefix/$prefix.0.bam $outdir/split/$prefix/$prefix.1.bam $outdir/split/$prefix/$prefix.phased $outdir/split/$prefix/$prefix.raw.vcf $outdir/split/$prefix\n";
	}
	return $cmd;

}


sub monitor_jobs{
	my $cmd=shift;
	my $step=shift;
	&creat_log($step,"START");
	my @unit_cmd=split /\n/,$cmd;
	foreach my $unit_cmd (@unit_cmd){
		my $check=system ($unit_cmd);
		&creat_log($unit_cmd,"CMD");
		if($check!=0){
			&creat_log($unit_cmd,"ERROR");
			&creat_log($step,"ERROR");
			die "$step error\n";
		}
	}
	&creat_log($step,"DONE");
	
}
sub check_dir{
	my $dir=shift;
	if(!-d $dir){
		system ("mkdir $dir");
	}else{
		die "$dir is exists, we will not override it\n";
	}
}
sub creat_log{
	my $string=shift;
	my $step=shift;
	my $status=shift;
	open OUT,">>$outdir/workflow_$program_start_time.log";
	my $date=`date`;
	$date=~s/\n//;
	my @string=split /\n/,$string;
	foreach my $cmd (@string){
		print OUT "[INFO] $date [$step] $cmd\n";	
	}
}
