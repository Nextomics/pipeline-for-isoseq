<ol>
<li>Quality control</li>
<li>
<p>Classification
step 1 QC and Classification, protocol version=&quot;2.3.0&quot; id=&quot;RS_IsoSeq.1&quot;, default  parameters.
smrtpipe.py --distribute  --params=settings.xml --output=outputdir xml:input.xml 2&gt; smrtpipe.stderr 1&gt; smrtpipe.stdout</p>
</li>
<li>
<p>Clustering
step 1 mapping and phasing
perl phase<em>allotetraploid</em>pipeline.pl –flnc flnc.fastq --gmap<em>genome<em>directory database/ --gmap</em>genome</em>database databasename –outdir ./result --reference<em>fasta ref.fasta
step 2 doing isoform-level-cluster according to alignments. 
python collapse</em>isoforms<em>by<em>sam.py -c 0.90 -i 0.90 --input flnc.fastq --fq -s flnc.sort.sam -o all
step 3 consensus, each cluster generate one consensus sequence.
perl analysis</em>cluster.pl all.collapsed.group.txt flnc.sort.sam flnc.fastq  &gt; flnc.best.sort.sam
step 4 doing isoform-level-cluster again.
python collapse<em>isoforms<em>by</em>sam.py -c 0.90 -i 0.90 --input chose.fq --fq -s flnc.best.sort.sam -o all.consensus
step 5 convert bam format to gff format.
samtools view -bS all.consensus.collapsed.rep.fq.sam &gt; all.consensus.bam
bedtools bamtobed  -split -i all.consensus.bam &gt; all.consensus.bed
perl bed2cDNA</em>match.pl all.consensus.collapsed.rep.fq all.consensus.collapsed.rep.fq.sam &gt; all.consensus.cDNA</em>match.gff</p>
</li>
<li>
<p>Transcriptome analysis
step 1 alternative splicing analysis.
python alternative<em>splice.py -i all.consensus.cDNA</em>match.gff -g ref.gtf -f ref.fasta -o ./ -os -as -ats T -op
step 2 alternative polyadenylation analysis.
perl polyA<em>position.pl all.consensus.collapsed.gff all.consensus.collapsed.rep.fq flnc.sort.sam &gt; transcript<em>polyA.result
step 3 finding fusion gene.
python fusion</em>finder.py --input flnc.fastq --fq -s flnc.sort.sam -o ./fusion 
step 4 finding non-coding RNA.
python PLEKModelling.py -lncRNA high</em>quality_lncRNA.fa -prefix species -mRNA mRNA.fasta
python PLEK.py  -fasta flnc.fasta -out lncRNA.predicted -thread 10 -range species.range -model species.model -k 4</p>
</li>
</ol>
<p>5.program list
smrtanalysis (http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)
gmap (http://research-pub.gene.com/gmap/)
phase<em>allotetraploid</em>pipeline.pl  (attached)
collapse<em>isoforms<em>by</em>sam.py from pbtranscript-tofu (https://github.com/PacificBiosciences/cDNA</em>primer)
analysis<em>cluster.pl in-house Perl script (attached)
bed2cDNA</em>match.pl in-house Perl script (attached)
samtools (http://www.htslib.org/)
bedtools (https://github.com/arq5x/bedtools2/)
alternative_splice.py developed by ourselves (attached)
polyA<em>position.pl (attached)
fusion</em>finder.py (https://github.com/PacificBiosciences/cDNA<em>primer/blob/master/pbtranscript-tofu/pbtranscript/pbtools/pbtranscript/fusion</em>finder.py)
PLEKModelling.py from plek (https://sourceforge.net/projects/plek/files/)
PLEK.py from plek (https://sourceforge.net/projects/plek/files/)</p>

