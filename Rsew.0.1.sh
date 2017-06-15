#!usr/bin/sh 
# convert fastq files into sorted bam files with tophat2 and samtools packages.
# Renyi Wu, 6-15-2017
set -x
pwd

#user input
FQ_EXT="fq.gz" #fastq file extension, or ".fq.gz", ".fastq.gz"

GENE_ANNO="/home/Administrator/Documents/mm9.2_for_bismark/Mus_musculus_UCSC_mm9.2/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf" #GENE_annotation file
GENE_REF="/home/Administrator/Documents/mm9.2_for_bismark/Mus_musculus_UCSC_mm9.2/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome" # Gene reference files


#create necessary folders
mkdir fastq  tophat_out bam bam_sorted cuff 

#Move all .gz files to fastq\ folder
mv *.gz ./fastq/


#run tophat2 for all fastq files, then move "accepted_hits.bam" files to "bam" folder with correlated names
for fqfile in ./fastq/*
do
basefn=$(basename $fqfile $FQ_EXT)
mkdir -pv tophat_out/$basefn
tophat2 -G $GENE_ANNO -p 8 -o "./tophat_out/$basefd" $GEN_REF $fqfile #$HOME_DIR"BAM_Files" 
mv ./tophat_out/$basefn/accepted_hits.bam ./bam/$basefn.bam
done

for bamfile in ./bam/*.bam
do
basefn=$(basename $bamfile .bam)
samtools sort -@ 8 $bamfile ./bam_sorted/$basefn.sorted
done
#Run samtools index, optional
mkdir -pv bam_sorted/index
for bamefilesort in ./bam_sorted/*.sorted.bam
do
samtools index $bamefilesort
cp ./bam_sorted/*.bai ./bam_sorted/index/
done
ls -lh ./bam_sorted/index/




