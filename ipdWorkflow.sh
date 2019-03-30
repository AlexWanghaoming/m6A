#!/bin/bash

## make bax fofn
fileTag=(`ls *.bax.h5 | cut -d "-" -f 2 | uniq`)
for i in ${fileTag[@]};do
{
	touch $i.fofn 
	for j in `ls *$i*.bax.h5`;do
		absPath=`pwd`/$j
		echo $absPath >> $i.fofn
	done
}
done
## multi-threads run baxFofn 
for baxFofn in `ls *.fofn`;do
{
	NAME=${baxFofn}
	bamPrefix=${NAME%.*}
	bax2bam --fofn $baxFofn --subread -o $bamPrefix
} &
done
wait 

## merge multi bam files
touch bamList.fofn
for bam in ls *.bam;do
	bamabsPath=`pwd`/$bam
	echo $bamabsPath >> bamList.fofn
done
samtools merge -b bamList.fofn mergeOut.bam -@ 16

## run blasr and ipdSummary
FASTA=`ls *.fasta`
blasr mergeOut.bam $FASTA --bam --out blasr_mapping.bam --nproc 16
samtools sort ./blasr_mapping.bam -o blasr_mapping_sorted.bam -@ 16
pbindex blasr_mapping_sorted.bam
samtools index blasr_mapping_sorted.bam
ipdSummary blasr_mapping_sorted.bam --reference $FASTA --identify m6A,m4C --methylFraction --gff basemods.gff --csv basemods.csv -j 16 --minCoverage 15


