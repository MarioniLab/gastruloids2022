#!/usr/local/bin/bash

# conda activate /nfs/research/marioni/Leah/software/conda_envs/htseq_env
#bsub -M 1000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/other_data/gouti2014/run_htseq_persample.sh -i /hps/nobackup/marioni/Leah/tophat_out/ -r /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/other_data/gouti2014/E-MTAB-2268.sdrf.txt -g /nfs/research/marioni/common/references/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/other_data/gouti2014/

while getopts i:r:o:g: flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        r) reffile=${OPTARG};;
        o) outdir=${OPTARG};;
        g) gtffile=${OPTARG};;
    esac
done

samples=($(tail --lines=+2 $reffile | awk -F '\t' '/^Header/ {next} ; {print $25}'))

bgadd -L 40 /htseq

mkdir -p ${outdir}
mkdir -p ${outdir}htseq_out/
mkdir -p ${outdir}htseq_out/outfiles/

project=gastruloid_scRNAseq

for s in "${samples[@]}"; do
  cd ${indir}
  bsub \
  -g /htseq \
  -o ${outdir}htseq_out/outfiles/${s}.txt \
  -M 5000 \
  htseq-count -c ${outdir}htseq_out/${s}.mtx ${indir}${s}/accepted_hits.bam ${gtffile}
done