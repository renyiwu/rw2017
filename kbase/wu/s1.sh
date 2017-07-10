#!/bin/bash 
# convert fastq files into sorted bam files with tophat2 and samtools packages.
# Renyi Wu,    v0.3-added fastqc on 6-15-20.;
cd "$(dirname "$0")" #or cd "${0%/*}" or cd "${0%/*}"
# redirect stdout/stderr to a file
#exec &> log.txt
# OR else to redirect only stdout use: 
exec 3>&1 4>&2
exec >> logfile.txt 2>&1
printf "new\n\n"
date
echo "current in: ${0%/*} , with these filefolder(s)"
ls -hp
date
sh s2.sh
echo "job done."
date
exec >&3 2>&4
date
printf "Job done. \nPress any key to quit."
read -n1 a





