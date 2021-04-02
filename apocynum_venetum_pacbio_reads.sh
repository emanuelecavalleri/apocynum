!/bin/bash

echo "DOWNLOAD RAW READS"
screen -L fastq-dump --defline-qual "+" --split-files --gzip --clip SRR7360364
fastq-stats SRR7360364_1.fastq.gz > venetumRAW.fastq.gz_stats.txt
# Select 200X of PacBio reads. The average read size is 7.9kb. Genome size ~ 159kb x 200 = 31,800kb / 7.9kb = 4,025
seqtk sample -s1000 SRR7360364_1.fastq.gz 4025 > venetumRAW.fq

echo "MAP READS (NGMLR)"
time ngmlr -t 8 -r ref.fa -q venetumRAW.fq | samtools view -Sb -F 4 -o venetumNGMLR.bam -
samtools sort -o venetumNGMLR.bam venetumNGMLR.bam
bedtools genomecov -d -ibam venetumNGMLR.bam -g ref.fa.fai > venetumNGMLR.covbed.txt
samtools fastq venetumNGMLR.bam > venetumNGMLR.fq
fastq-stats venetumNGMLR.fq > venetumNGMLR.fq_stats.txt

echo "MAP READS (MINIMAP2)"
minimap2 -d ref.mmi ref.fa
time minimap2 -t 8 -a ref.mmi venetumRAW.fq | samtools view -Sb -F 4 -o venetumMM2.bam -
samtools sort -o venetumMM2.bam venetumMM2.bam
bedtools genomecov -d -ibam venetumMM2.bam -g ref.fa.fai > venetumMM2.covbed.txt
samtools fastq venetumMM2.bam > venetumMM2.fq
fastq-stats venetumMM2.fq > venetumMM2.fq_stats.txt

echo "ASSEMBLE READS"
time canu -d venetumCHL_canuNGMLR -p venetumCHL genomeSize=159k -pacbio-raw venetumNGMLR.fq
FastaSeqStats -i venetumCHL_canuNGMLR/venetumCHL.contigs.fasta > venetumCHL.contigs.fasta_stats.txt
union -filter venetumCHL_canuNGMLR/venetumCHL.contigs.fasta > venetumCHL.fasta

time canu -d venetumCHL_canuMM2 -p venetumCHL genomeSize=159k -pacbio-raw venetumMM2.fq
FastaSeqStats -i venetumCHL_canuMM2/venetumCHL.contigs.fasta > venetumCHL.contigs.fasta_stats.txt
union -filter venetumCHL_canuMM2/venetumCHL.contigs.fasta > venetumCHL.fasta
