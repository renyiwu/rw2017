#!/bin/bash 
# convert fastq files into sorted bam files with tophat2 and samtools packages.
# Renyi Wu,    v0.3-added fastqc on 6-15-20.;
cd "$(dirname "$0")" #or cd "${0%/*}" or cd "${0%/*}"
#Run FastQC,optional
mkdir FastQC_reports_auto #make directory for reports, optional
fastqc --extract --outdir=./FastQC_reports_auto/ *.gz
#echo "done"
#echo "counting reads"
date >fastqc_reads_for_all.txt
for i in ./FastQC_reports_auto/*/fastqc_data.txt
do
grep "Filename" "$i" >>fastqc_reads_for_all.txt
grep "Total Sequences" "$i" >>fastqc_reads_for_all.txt
grep ">>Sequence Duplication Levels" "$i" >>fastqc_reads_for_all.txt
grep "#Total Deduplicated Percentage" "$i" >>fastqc_reads_for_all.txt
grep "#Duplication Level" "$i" >>fastqc_reads_for_all.txt
#grep "^1\t" "$i" >>fastqc_reads_for_all.txt
awk '/^1\t/{n++}n==5{print; exit}' "$i" >>fastqc_reads_for_all.txt
echo "\n" >>fastqc_reads_for_all.txt
done
#read -n1 -r -p "Press any key to quit" key
