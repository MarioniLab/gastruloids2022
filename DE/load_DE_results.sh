#!/usr/local/bin/bash

# conda activate basic_renv
# cd /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline_final/DE/
# bsub -M 1000 ./load_DE_results.sh

bgadd -L 500 /DEfull

outbase="/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/DE_out/"

ct_seq=( 'All' 'Anterior Primitive Streak' 'Cardiopharyngeal Mesoderm' 'Caudal Epiblast' 'Caudal Neurectoderm' 'Early Anterior PSM' 'Early Nascent Mesoderm' 'Early Neurectoderm' 'Early Posterior PSM' 'Early Spinal Cord' 'Endothelium' 'Epiblast' 'Head Mesoderm' 'Late Anterior PSM' 'Late Nascent Mesoderm' 'Late Neurectoderm' 'Late Posterior PSM' 'Late Spinal Cord' 'Mature Endoderm' 'Mature Somites' 'NMPs' 'Neurons' 'Notochord' 'PGCs' 'Primitive Streak' 'Somites' )

for ct1 in "${ct_seq[@]}"
do
    for ct2 in "${ct_seq[@]}"
    do
        bsub -g /DEfull -o "$outbase"outfiles/loadfiles_"${ct1}"_"${ct2}".txt -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline_final/DE/load_DE_results.R -o "${outbase}" -c "${ct1}" -C "${ct2}"
    done
done