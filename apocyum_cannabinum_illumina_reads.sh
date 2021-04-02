!/bin/bash

echo "DOWNLOAD RAW READS"
screen -L fastq-dump --defline-qual "+" --split-files --gzip --clip SRR6425639
fastq-stats SRR6425639_1.fastq.gz > cannabinumRAWR1.fastq.gz_stats.txt
fastq-stats SRR6425639_2.fastq.gz > cannabinumRAWR2.fastq.gz_stats.txt

echo "PROCESS READS"
fastq-mcf -o cannabinumR1.fastq.gz -o cannabinumR2.fastq.gz /data/99_genomicclass/00_shared_data/Genomics/IlluminaAdapters_V2.fasta <(gunzip -c SRR6425639_1.fastq.gz;) <(gunzip -c SRR6425639_2.fastq.gz;)
fastq-stats cannabinumR1.fastq.gz > cannabinumR1.fastq.gz_stats.txt
fastq-stats cannabinumR2.fastq.gz > cannabinumR2.fastq.gz_stats.txt

echo "MAP READS (BOWTIE 2)"
bowtie2-build ref.fa refc
time bowtie2 -p 8 -x refc -1 cannabinumR1.fastq.gz -2 cannabinumR2.fastq.gz | samtools view -Sb -F 4 -o cannabinumBT2.bam -
samtools sort -o cannabinumBT2.bam cannabinumBT2.bam
samtools index cannabinumBT2.bam
samtools faidx ref.fa
bedtools genomecov -d -ibam cannabinumBT2.bam -g ref.fa.fai > cannabinumBT2.covbed.txt
samtools fastq cannabinumBT2.bam -1 cannabinumBT2R1.fq -2 cannabinumBT2R2.fq
fastq-stats cannabinumBT2R1.fq > cannabinumBT2R1.fq_stats.txt
fastq-stats cannabinumBT2R2.fq > cannabinumBT2R2.fq_stats.txt
freebayes -v cannabinumBT2.vcf -b cannabinumBT2.bam -f ref.fa
grep -v "#" cannabinumBT2.vcf | wc -l > cannabinumBT2_total_variants.txt

echo "MAP READS (BWA)"
bwa index ref.fa
time bwa mem ref.fa cannabinumR1.fastq.gz cannabinumR2.fastq.gz -t 8 | samtools view -Sb -F 4 -o cannabinumBWA.bam -
samtools sort -o cannabinumBWA.bam cannabinumBWA.bam
samtools index cannabinumBWA.bam
bedtools genomecov -d -ibam cannabinumBWA.bam -g ref.fa.fai > cannabinumBWA.covbed.txt
samtools fastq cannabinumBWA.bam -1 cannabinumBWAR1.fq -2 cannabinumBWAR2.fq
fastq-stats cannabinumBWAR1.fq > f1_stats.txt
fastq-stats cannabinumBWAR2.fq > f2_stats.txt
freebayes -v cannabinumBWA.vcf -b cannabinumBWA.bam -f ref.fa
grep -v "#" avn.vcf | wc -l > cannabinumBWA_total_variants.txt

echo "ASSEMBLE READS"
time abyss-pe name="abyssBT2cannabinum" k=63 in="cannabinumBT2R1.fq cannabinumBT2R2.fq"
FastaSeqStats -i abyssBT2cannabinum-scaffolds.fa > abyssBT2cannabinum-scaffolds_stats.txt
union -filter abyssBT2cannabinum-scaffolds.fa > abyssBT2cannabinum.fasta

time abyss-pe name="abyssBWAcannabinum" k=63 in="cannabinumBWAR1.fq cannabinumBWAR2.fq"
FastaSeqStats -i abyssBWAcannabinum-scaffolds.fa > abyssBWAcannabinum-scaffolds_stats.txt
union -filter abyssBWAcannabinum-scaffolds.fa > abyssBWAcannabinum.fasta
