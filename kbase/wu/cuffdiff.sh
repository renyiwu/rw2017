#find differentially expressed genes with cuffdiff
#
#
cd ./bam_sorted/
GENE_ANNO="/home/administrator/Documents/mm9.2_for_bismark/Mus_musculus_UCSC_mm9.2/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
#
cuffdiff -p 12 -o ../cuff/ -L control,TPA,TPA1,CA,FX,CDDO,mITC $GENE_ANNO RW03.sorted.bam RW04.sorted.bam RW01.sorted.bam RW05.sorted.bam RW06.sorted.bam RW07.sorted.bam RW02.sorted.bam
