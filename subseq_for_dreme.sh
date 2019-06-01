#!/bin/bash
# extract +-5nt around 6mA sites
awk -F ";|=" 'BEGIN{i=1}{print ">"i"\n"substr($4,16,11);i++}' m6A.gff > m6A_context.fasta
# run dreme
dreme-py3 -o m6A_motif -p m6A_context.fasta
