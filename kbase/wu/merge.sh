#Merge .bam files with samtools merge function 

cd /bam/file/dir/
samtools merge -@ 8 RW04.sorted.bam RW041.sorted.bam RW042.sorted.bam #set 8 thread, out.bam in1.bam in2. bam [...]

