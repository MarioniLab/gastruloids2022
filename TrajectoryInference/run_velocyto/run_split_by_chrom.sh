#!/usr/local/bin/bash

while getopts i:r:n flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        r) reffile=${OPTARG};;
        n) chromnum=true;;
    esac
done

mkdir -p "$indir"lsf_output/
mkdir -p "$indir""chrom_split/"
readarray -t chromosomes < <(awk -F$'\t' '/^[^#]/ { print $1 }' "$reffile" | sort | uniq)
for chrom in "${chromosomes[@]}"
do
    if [ "$chromnum" = true ] && [ "${chrom:0:3}" = "chr" ] && [ $chrom != "chrM" ] ; then
        chrom="${chrom:3}"
    elif [ "$chromnum" = true ] && [ $chrom = "chrM" ] ; then
        chrom="MT"
    fi
    mkdir -p "$indir"chrom_split/"$chrom"/
    bsub -M 25000 -o "$indir"lsf_output/splitchrom_output_"$chrom".txt -e "$indir"lsf_output/splitchrom_error_"$chrom".txt /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/scvelo/split_by_chrom.sh -i "$indir" -c "$chrom"
done