package  Soft;
use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(parse_config);
##parse the software.config file, and check the existence of each software
################################################

sub parse_config{
        my ($config,$soft)=@_;
        open IN,$config || die;
        my %ha;
        while(<IN>){
                chomp;
                next if /^#|^$/;
                s/\s+//g;
                my @p=split/=/,$_;
                $ha{$p[0]}=$p[1];
        }
        close IN;
        if(exists $ha{$soft}){
                if(-e $ha{$soft}){
                        return $ha{$soft};
                }else{
                         die "\nConfig Error: $soft wrong path in $config\n";
                }
        }else{
                die "\nConfig Error: $soft not set in $config\n";
        }
}
1;              
__END__         
