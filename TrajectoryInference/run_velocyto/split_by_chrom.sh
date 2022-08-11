#!/usr/local/bin/bash

while getopts i:c: flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        c) chrom=${OPTARG};;
    esac
done

samtools view "$indir""outs/possorted_genome_bam.bam" "$chrom" -b > "$indir"chrom_split/"$chrom"/possorted_genome_bam_"$chrom".bam
samtools index "$indir"chrom_split/"$chrom"/possorted_genome_bam_"$chrom".bam "$indir"chrom_split/"$chrom"/possorted_genome_bam_"$chrom".bam.bai
samtools sort -m 20000M -O BAM -o "$indir"chrom_split/"$chrom"/cellsorted_possorted_genome_bam_"$chrom".bam "$indir"chrom_split/"$chrom"/possorted_genome_bam_"$chrom".bam
samtools index "$indir"chrom_split/"$chrom"/cellsorted_possorted_genome_bam_"$chrom".bam "$indir"chrom_split/"$chrom"/cellsorted_possorted_genome_bam_"$chrom".bam.bai
