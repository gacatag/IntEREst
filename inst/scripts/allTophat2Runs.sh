#!/bin/bash
#
# Run this in your mapping folder. Read scripts/readme.txt in IntEREst package.
#indexPATH=~/software/sequencing/bowtie2/indexes/hg19NoColor/hg19
#MAPPINGPATH=~/analysis/zrsr2/mapping
#RAWDATAPATH=~/projects/data/zrsr2
cd $MAPPINGPATH
mkdir $MAPPINGPATH/SRR1691633/
tophat2 -o $MAPPINGPATH/SRR1691633/ -p 6 $indexPATH $RAWDATAPATH/SRR1691633_1.fastq $RAWDATAPATH/SRR1691633_2.fastq
samtools sort $MAPPINGPATH/SRR1691633/accepted_hits.bam $MAPPINGPATH/SRR1691633/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691633/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691634/
tophat2 -o $MAPPINGPATH/SRR1691634/ -p 6 $indexPATH $RAWDATAPATH/SRR1691634_1.fastq $RAWDATAPATH/SRR1691634_2.fastq
samtools sort $MAPPINGPATH/SRR1691634/accepted_hits.bam $MAPPINGPATH/SRR1691634/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691634/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691635/
tophat2 -o $MAPPINGPATH/SRR1691635/ -p 6 $indexPATH $RAWDATAPATH/SRR1691635_1.fastq $RAWDATAPATH/SRR1691635_2.fastq
samtools sort $MAPPINGPATH/SRR1691635/accepted_hits.bam $MAPPINGPATH/SRR1691635/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691635/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691636/
tophat2 -o $MAPPINGPATH/SRR1691636/ -p 6 $indexPATH $RAWDATAPATH/SRR1691636_1.fastq $RAWDATAPATH/SRR1691636_2.fastq
samtools sort $MAPPINGPATH/SRR1691636/accepted_hits.bam $MAPPINGPATH/SRR1691636/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691636/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691637/
tophat2 -o $MAPPINGPATH/SRR1691637/ -p 6 $indexPATH $RAWDATAPATH/SRR1691637_1.fastq $RAWDATAPATH/SRR1691637_2.fastq
samtools sort $MAPPINGPATH/SRR1691637/accepted_hits.bam $MAPPINGPATH/SRR1691637/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691637/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691638/
tophat2 -o $MAPPINGPATH/SRR1691638/ -p 6 $indexPATH $RAWDATAPATH/SRR1691638_1.fastq $RAWDATAPATH/SRR1691638_2.fastq
samtools sort $MAPPINGPATH/SRR1691638/accepted_hits.bam $MAPPINGPATH/SRR1691638/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691638/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691639/
tophat2 -o $MAPPINGPATH/SRR1691639/ -p 6 $indexPATH $RAWDATAPATH/SRR1691639_1.fastq $RAWDATAPATH/SRR1691639_2.fastq
samtools sort $MAPPINGPATH/SRR1691639/accepted_hits.bam $MAPPINGPATH/SRR1691639/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691639/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691640/
tophat2 -o $MAPPINGPATH/SRR1691640/ -p 6 $indexPATH $RAWDATAPATH/SRR1691640_1.fastq $RAWDATAPATH/SRR1691640_2.fastq
samtools sort $MAPPINGPATH/SRR1691640/accepted_hits.bam $MAPPINGPATH/SRR1691640/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691640/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691641/
tophat2 -o $MAPPINGPATH/SRR1691641/ -p 6 $indexPATH $RAWDATAPATH/SRR1691641_1.fastq $RAWDATAPATH/SRR1691641_2.fastq
samtools sort $MAPPINGPATH/SRR1691641/accepted_hits.bam $MAPPINGPATH/SRR1691641/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691641/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691642/
tophat2 -o $MAPPINGPATH/SRR1691642/ -p 6 $indexPATH $RAWDATAPATH/SRR1691642_1.fastq $RAWDATAPATH/SRR1691642_2.fastq
samtools sort $MAPPINGPATH/SRR1691642/accepted_hits.bam $MAPPINGPATH/SRR1691642/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691642/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691643/
tophat2 -o $MAPPINGPATH/SRR1691643/ -p 6 $indexPATH $RAWDATAPATH/SRR1691643_1.fastq $RAWDATAPATH/SRR1691643_2.fastq
samtools sort $MAPPINGPATH/SRR1691643/accepted_hits.bam $MAPPINGPATH/SRR1691643/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691643/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691644/
tophat2 -o $MAPPINGPATH/SRR1691644/ -p 6 $indexPATH $RAWDATAPATH/SRR1691644_1.fastq $RAWDATAPATH/SRR1691644_2.fastq
samtools sort $MAPPINGPATH/SRR1691644/accepted_hits.bam $MAPPINGPATH/SRR1691644/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691644/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691645/
tophat2 -o $MAPPINGPATH/SRR1691645/ -p 6 $indexPATH $RAWDATAPATH/SRR1691645_1.fastq $RAWDATAPATH/SRR1691645_2.fastq
samtools sort $MAPPINGPATH/SRR1691645/accepted_hits.bam $MAPPINGPATH/SRR1691645/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691645/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691646/
tophat2 -o $MAPPINGPATH/SRR1691646/ -p 6 $indexPATH $RAWDATAPATH/SRR1691646_1.fastq $RAWDATAPATH/SRR1691646_2.fastq
samtools sort $MAPPINGPATH/SRR1691646/accepted_hits.bam $MAPPINGPATH/SRR1691646/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691646/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691647/
tophat2 -o $MAPPINGPATH/SRR1691647/ -p 6 $indexPATH $RAWDATAPATH/SRR1691647_1.fastq $RAWDATAPATH/SRR1691647_2.fastq
samtools sort $MAPPINGPATH/SRR1691647/accepted_hits.bam $MAPPINGPATH/SRR1691647/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691647/accepted_hits.sorted.bam

mkdir $MAPPINGPATH/SRR1691648/
tophat2 -o $MAPPINGPATH/SRR1691648/ -p 6 $indexPATH $RAWDATAPATH/SRR1691648_1.fastq $RAWDATAPATH/SRR1691648_2.fastq
samtools sort $MAPPINGPATH/SRR1691648/accepted_hits.bam $MAPPINGPATH/SRR1691648/accepted_hits.sorted
samtools index $MAPPINGPATH/SRR1691648/accepted_hits.sorted.bam
