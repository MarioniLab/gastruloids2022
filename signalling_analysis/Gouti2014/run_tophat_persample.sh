#!/usr/local/bin/bash

# conda activate /nfs/research/marioni/Leah/software/conda_envs/alignment_env
#bsub -M 1000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/other_data/gouti2014/run_tophat_persample.sh -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/other_data/gouti2014 -r /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/other_data/gouti2014/E-MTAB-2268.sdrf.txt -o /hps/nobackup/marioni/Leah

while getopts i:r:o: flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        r) reffile=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done

samples=($(tail --lines=+2 $reffile | awk -F '\t' '/^Header/ {next} ; {print $25}'))

bgadd -L 40 /tophat

mkdir -p ${outdir}
mkdir -p ${outdir}/tophat_out/
mkdir -p ${outdir}/tophat_out/outfiles/

project=gastruloid_scRNAseq

for s in "${samples[@]}"; do
  mkdir -p ${outdir}/tophat_out/${s}/
  cd $indir/
  bsub \
  -g /tophat \
  -o ${outdir}/tophat_out/outfiles/${s}.txt \
  -M 10000 \
  tophat -o $outdir/tophat_out/${s}/ . $indir/${s}.fastq.gz
done