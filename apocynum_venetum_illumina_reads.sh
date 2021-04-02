!/bin/bash

echo "DOWNLOAD RAW READS"
screen -L fastq-dump --defline-qual "+" --split-files --gzip --clip SRR7360366
fastq-stats SRR7360366_1.fastq.gz > venetumRAWR1.fastq.gz_stats.txt
fastq-stats SRR7360366_2.fastq.gz > venetumRAWR2.fastq.gz_stats.txt

echo "PROCESS READS"
fastq-mcf -o venetumR1.fastq.gz -o venetumR2.fastq.gz /data/99_genomicclass/00_shared_data/Genomics/IlluminaAdapters_V2.fasta <(gunzip -c SRR7360366_1.fastq.gz;) <(gunzip -c SRR7360366_2.fastq.gz;)
fastq-stats venetumR1.fastq.gz > venetumR1.fastq.gz_stats.txt
fastq-stats venetumR2.fastq.gz > venetumR2.fastq.gz_stats.txt

echo "MAP READS (BOWTIE 2)"
bowtie2-build ref.fa refv
time bowtie2 -p 8 -x refv -1 venetumR1.fastq.gz -2 venetumR2.fastq.gz | samtools view -Sb -F 4 -o venetumBT2.bam -
samtools sort -o venetumBT2.bam venetumBT2.bam
samtools index venetumBT2.bam
samtools faidx ref.fa
bedtools genomecov -d -ibam venetumBT2.bam -g ref.fa.fai > venetumBT2.covbed.txt
samtools fastq venetumBT2.bam -1 venetumBT2R1.fq -2 venetumBT2R2.fq
fastq-stats venetumBT2R1.fq > venetumBT2R1.fq_stats.txt
fastq-stats venetumBT2R2.fq > venetumBT2R2.fq_stats.txt
freebayes -v venetumBT2.vcf -b venetumBT2.bam -f ref.fa
grep -v "#" venetumBT2.vcf | wc -l > venetumBT2_total_variants.txt

echo "MAP READS (BWA)"
bwa index ref.fa
time bwa mem ref.fa venetumR1.fastq.gz venetumR2.fastq.gz -t 8 | samtools view -Sb -F 4 -o venetumBWA.bam -
samtools sort -o venetumBWA.bam venetumBWA.bam
samtools index venetumBWA.bam
bedtools genomecov -d -ibam venetumBWA.bam -g ref.fa.fai > venetumBWA.covbed.txt
samtools fastq venetumBWA.bam -1 venetumBWAR1.fq -2 venetumBWAR2.fq
fastq-stats venetumBWAR1.fq > venetumBWAR1.fq_stats.txt
fastq-stats venetumBWAR2.fq > venetumBWAR2.fq_stats.txt
freebayes -v venetumBWA.vcf -b venetumBWA.bam -f ref.fa
grep -v "#" venetumBWA.vcf | wc -l > venetumBWAtotal_variants.txt

echo "ASSEMBLE READS"
time abyss-pe name="abyssBT2venetum" k=63 in="venetumBT2R1.fq venetumBT2R2.fq"
FastaSeqStats -i abyssBT2venetum-scaffolds.fa > abyssBT2venetum-scaffolds.fa_stats.txt
union -filter abyssBT2venetum-scaffolds.fa > abyssBT2venetum.fasta

time abyss-pe name="abyssBWAvenetum" k=63 in="venetumBWAR1.fq venetumBWAR2.fq"
FastaSeqStats -i abyssBWAvenetum-scaffolds.fa > abyssBWAvenetum-scaffolds.fa_stats.txt
union -filter abyssBWAvenetum-scaffolds.fa > abyssBWAvenetum.fasta
