### program list
##smrtanalysis (http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)
## gmap (http://research-pub.gene.com/gmap/)
## phase_allotetraploid_pipeline.pl  (attached)
## collapse_isoforms_by_sam.py from pbtranscript-tofu (https://github.com/PacificBiosciences/cDNA_primer)
## analysis_cluster.pl in-house Perl script (attached)
## bed2cDNA_match.pl in-house Perl script (attached)
## samtools (http://www.htslib.org/)
## bedtools (https://github.com/arq5x/bedtools2/)
## alternative_splice.py developed by ourselves (attached)
## polyA_position.pl (attached)
## fusion_finder.py (https://github.com/PacificBiosciences/cDNA_primer/blob/master/pbtranscript-tofu/pbtranscript/pbtools/pbtranscript/fusion_finder.py)
## PLEKModelling.py from plek (https://sourceforge.net/projects/plek/files/)
## PLEK.py from plek (https://sourceforge.net/projects/plek/files/)

 
1. Quality control
2. Classification
#step 1 QC and Classification, protocol version="2.3.0" id="RS_IsoSeq.1", default  parameters.
smrtpipe.py --distribute  --params=settings.xml --output=outputdir xml:input.xml 2> smrtpipe.stderr 1> smrtpipe.stdout
3. Clustering

#  step 1 mapping and phasing
perl phase_allotetraploid_pipeline.pl –flnc flnc.fastq --gmap_genome_directory database/ --gmap_genome_database databasename –outdir ./result --reference_fasta ref.fasta

#step 2 doing isoform-level-cluster according to alignments. 
python collapse_isoforms_by_sam.py -c 0.90 -i 0.90 --input flnc.fastq --fq -s flnc.sort.sam -o all

#step 3 consensus, each cluster generate one consensus sequence.
perl analysis_cluster.pl all.collapsed.group.txt flnc.sort.sam flnc.fastq  > flnc.best.sort.sam

#step 4 doing isoform-level-cluster again.
python collapse_isoforms_by_sam.py -c 0.90 -i 0.90 --input chose.fq --fq -s flnc.best.sort.sam -o all.consensus

#step 5 convert bam format to gff format.
samtools view -bS all.consensus.collapsed.rep.fq.sam > all.consensus.bam

bedtools bamtobed  -split -i all.consensus.bam > all.consensus.bed

perl bed2cDNA_match.pl all.consensus.collapsed.rep.fq all.consensus.collapsed.rep.fq.sam > all.consensus.cDNA_match.gff


4. Transcriptome analysis

#step 1 alternative splicing analysis.
python alternative_splice.py -i all.consensus.cDNA_match.gff -g ref.gtf -f ref.fasta -o ./ -os -as -ats T -op

#step 2 alternative polyadenylation analysis.
perl polyA_position.pl all.consensus.collapsed.gff all.consensus.collapsed.rep.fq flnc.sort.sam > transcript_polyA.result

#step 3 finding fusion gene.
python fusion_finder.py --input flnc.fastq --fq -s flnc.sort.sam -o ./fusion 

#step 4 finding non-coding RNA.
python PLEKModelling.py -lncRNA high_quality_lncRNA.fa -prefix species -mRNA mRNA.fasta

python PLEK.py  -fasta flnc.fasta -out lncRNA.predicted -thread 10 -range species.range -model species.model -k 4
