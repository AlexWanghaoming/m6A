#!/bin/bash
# used for calculate methylation fraction
# Author: Wang Haoming
GFF_fileName=$1
echo "ipdSummary result file: $GFF_fileName"
m6ASiteNum=`awk '$3~/m6A/' $1 | wc -l | cut -f 1`
echo "Total 6mA site: $m6ASiteNum"
m4CSiteNum=`awk '$3~/m4C/' $1 | wc -l | cut -f 1`
echo "Total m4C site: $m4CSiteNum"

FASTA_FileName=$2
AT_baseCount=`grep -v "^>" $2 | grep -o "[AaTt]" | wc -l | cut -f 1`
CG_baseCount=`grep -v "^>" $2 | grep -o "[CcGg]" | wc -l | cut -f 1`

m6A_ratio=`awk 'BEGIN{printf"%.2f%\n",('$m6ASiteNum'/'$AT_baseCount')*100}'`
echo "6mA methylation ratio: $m6A_ratio"

m4C_ratio=`awk 'BEGIN{printf"%.2f%\n",('$m4CSiteNum'/'$CG_baseCount')*100}'`
echo "m4C methylation ratio: $m4C_ratio"
