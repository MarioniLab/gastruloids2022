#!/usr/local/bin/bash

while getopts i:r:m:n flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        r) reffile=${OPTARG};;
        m) maskfile=${OPTARG};;
        n) chromnum=true;;
    esac
done

mkdir -p "$indir"lsf_output/
mkdir -p "$indir"velocyto/
readarray -t chromosomes < <(awk -F$'\t' '/^[^#]/ { print $1 }' "$reffile" | sort | uniq)
for chrom in "${chromosomes[@]}"
do
    if [ "$chromnum" = true ] && [ "${chrom:0:3}" = "chr" ] && [ $chrom != "chrM" ] ; then
        chrom="${chrom:3}"
    elif [ "$chromnum" = true ] && [ $chrom = "chrM" ] ; then
        chrom="MT"
    fi
    mkdir -p "$indir"velocyto/"$chrom"/
    bsub -M 50000 -o "$indir"lsf_output/velocyto_output_"$chrom".txt -e "$indir"lsf_output/velocyto_error_"$chrom".txt velocyto run -b "$indir"outs/filtered_feature_bc_matrix/barcodes.tsv.gz -m "$maskfile" -o "$indir"velocyto/"$chrom"/ "$indir"chrom_split/"$chrom"/cellsorted_possorted_genome_bam_*.bam "$reffile"
done
