#!/bin/bash

names=(
"sample"
"reads_in_star"
"av_read_length"
"av_mapped_length"
"align_rate"
)

header=`printf ",%s" "${names[@]}" | sed "s/^.//" && printf "\n"`
echo $header > qc.csv

samples=($(ls -d ../../data/star/*/))

for i in {1..36}; do 
    
    sample=$(basename ${samples[${i}-1]})
    echo $sample

    values=()
    values+=($sample)

    ifn_st="../../data/star/${sample}/*.final.out"

    values+=(`grep "Number of input reads" ${ifn_st} | 
        sed 's/[^0-9.]//g'`)
    values+=(`grep "Average input read length" ${ifn_st} | 
        sed 's/[^0-9.]//g'`)
    values+=(`grep "Average mapped length" ${ifn_st} | 
        sed 's/[^0-9.]//g'`)
    values+=(`grep "Uniquely mapped reads %" ${ifn_st} | 
        sed 's/[^0-9.]//g'`)

    line=`printf ",%s" "${values[@]}" | sed "s/^.//" && printf "\n"`
    echo $line >> qc.csv

done
