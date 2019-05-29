#!/bin/bash
OUTPUT=$(pwd)/read_assembly_ITS
mkdir $OUTPUT

#raconloop
for j in 0 1 2 3 4 5 6 7 8 9 10
  do
    echo ''
    echo ${j}
    echo ''
    minimap2 ${OUTPUT}/consensus_${j}.fasta $1 > ${OUTPUT}/consensus_${j}.paf
    racon $1 ${OUTPUT}/consensus_${j}.paf ${OUTPUT}/consensus_${j}.fasta > ${OUTPUT}/consensus_$((j+1)).fasta
  done



OUTPUT=$(pwd)/read_assembly_TEF
mkdir $OUTPUT


#raconloop
for j in 0 1 2 3 4 5 6 7 8 9 10
  do
    echo ''
    echo ${j}
    echo ''
    minimap2 ${OUTPUT}/consensus_${j}.fasta $2 > ${OUTPUT}/consensus_${j}.paf
    racon -e 0.5 $2 ${OUTPUT}/consensus_${j}.paf ${OUTPUT}/consensus_${j}.fasta > ${OUTPUT}/consensus_$((j+1)).fasta
  done