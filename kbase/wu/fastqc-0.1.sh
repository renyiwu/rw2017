#!/bin/bash 
# convert fastq files into sorted bam files with tophat2 and samtools packages.
# Renyi Wu,    v0.3-added fastqc on 6-15-20.;
cd "$(dirname "$0")" #or cd "${0%/*}" or cd "${0%/*}"


#Run FastQC,optional
mkdir FastQC_reports_auto #make directory for reports, optional
fastqc --extract --outdir=./FastQC_reports_auto/ *.gz
echo "done"
read -n1 -r -p "Press any key to quit" key
