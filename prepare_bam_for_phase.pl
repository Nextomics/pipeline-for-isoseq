#!/usr/bin/perl -w
use strict;
use FileHandle;
my $input=shift;
if($input=~/gz$/){
	open IN,"zcat $input|";
}else{
	open IN,$input;
}
our %reads;
my @sam;
while(<IN>){
	chomp;
	my $head;
	if(/^#/ || /^@/){
		my $fh = FileHandle->new(">>head.txt");
		print $fh $_,"\n";
		$fh->close;	
	}else{
		my @temp=split;
		next if($temp[2]=~/\*/);
		&check_dir($temp[2]);
		
		if(!-f "$temp[2].sam"){
			`cp head.txt $temp[2]/$temp[2].sam`;
		}
		
		my $fh = FileHandle->new(">>$temp[2]/$temp[2].sam");
		####
		my $overhang=&overhang_len($_);
		if(exists $reads{$temp[2]}{$temp[0]}){
			if($reads{$temp[2]}{$temp[0]}->[0]<$overhang){
				$reads{$temp[2]}{$temp[0]}->[0]=$overhang;
				$reads{$temp[2]}{$temp[0]}->[1]=$temp[3];
			}
		}else{
			$reads{$temp[2]}{$temp[0]}->[0]=$overhang;
			$reads{$temp[2]}{$temp[0]}->[1]=$temp[3];
		}
		####
		my $line=&recal_MQ($_);
	
		print $fh $line,"\n";
		$fh->close;
	}

}
close IN;
open LIST,">sam.list";
foreach my $chr (keys %reads){
	my $bam_file="$chr/$chr.sam";
	print LIST "$chr/$chr.besthit.sam\n";
	&read_sam($bam_file);
}
close LIST;
sub recal_MQ{
	my $string=shift;
	my @temp=split /\t/,$string;
	my $mq=$temp[4];
	if($mq==0){
		$temp[4]=40;
	}
	my $newline=join "\t",@temp;
	return $newline;
}

sub read_sam{
    my $sam_file=shift;

    open IN,$sam_file;
	my $outsam=$sam_file;
	$outsam=~s/sam/besthit.sam/;
	open SAM,">$outsam";
    while(<IN>){
        chomp;

        if(/^@/){
			print SAM $_,"\n";
		}else{
	        my @temp=split;
			if(exists $reads{$temp[2]}{$temp[0]}){
				if($reads{$temp[2]}{$temp[0]}->[1]==$temp[3]){
					print SAM $_,"\n";
				}
			}
		}
    }
    close IN;
	close SAM;

}
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
sub check_dir{
    my $dir=shift;
    if(!-d $dir){
        system ("mkdir $dir");
    }
}
