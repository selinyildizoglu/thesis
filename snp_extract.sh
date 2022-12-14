#!/bin/bash

while getopts S:r:d:o:a:s:c:t:vi flag
do
    case "${flag}" in
        S) sample=${OPTARG};;
        r) ref=${OPTARG};;
        d) dfd=${OPTARG};;
        o) ofd=${OPTARG};;
        a) afd=${OPTARG};;
        s) sfd=${OPTARG};;
        t) maf_th=${OPTARG};;
        c) min_cov=${OPTARG};;
        v) verb="-v";;
        i) inf="-i";;
    esac
done

### Data folder

if [[ -z "$dfd" ]]; then
    echo "Warning: data folder not specified" >&2
    dfd=$(realpath "$(pwd)/../../data/")
    echo -e "\tAssuming ${dfd}" >&2
fi

if [[ ! -d "$dfd" ]]; then 
    echo "Error: data folder does not exist" >&2
    echo -e "\tExiting" >&2; exit 1
fi

### Alignments

if [[ -z "$afd" ]]; then
    afd="${dfd}/bam_rdna/"
else
    afd="${dfd}/${afd}/"
fi

if [[ ! -d "$afd" ]]; then 
    echo "Error: Alignment folder does not exist" >&2
    echo -e "\tExiting" >&2; exit 1
fi


### SNPs

if [[ -z "$sfd" ]]; then
    echo "Error: SNPs file or folder name missing" >&2
    echo -e "\tExiting" >&2; exit 1
fi

if [[ ! -f "$sfd" ]]; then
    if [[ ! -d "$sfd" ]]; then
        echo "Error: SNPs file or folder does not exist" >&2
        echo -e "\tExiting" >&2; exit 1
    fi
fi

### Reference

if [[ -z "$ref" ]]; then
    ref="BK000964.3_looped_3008"
fi

### Output

if [[ -z "$ofd" ]]; then
    echo "Warning: Output folder name missing" >&2
    ofd="${dfd}snp_extract/"
    echo -e "\tAssuming ${ofd}" >&2
else
    ofd="${dfd}${ofd}/"
fi

### Sample

if [[ -z "$sample" ]]; then
    echo "Error: Root name missing" >&2
    echo -e "\tExiting" >&2; exit 1
fi

### Parameters

if [[ -z "$maf_th" ]]; then
    maf_th="0.01"
fi

if [[ -z "$min_cov" ]]; then
    min_cov="10"
fi


module load -s python/3.6.3
source ~/envs/bisvar/bin/activate

python /data/Blizard-Rakyan/Rakyan_Lab_Files/Analysis_Pipelines/var/snp_extract/snp_extract.py \
    -o $ofd \
    -a $afd \
    -s $sfd \
    -c $min_cov \
    -t $maf_th \
    -r $ref \
    ${verb} \
    $sample 

deactivate

module unload python
