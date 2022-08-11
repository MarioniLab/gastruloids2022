#!/usr/local/bin/bash

# conda activate /nfs/research/marioni/Leah/software/conda_envs/deconvolution
# cd /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/
# bsub -M 1000 ./load_DA_simulation_files.sh

bgadd -L 500 /hetero_sims_load

outbase="/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/20220628/"

#DA_meths=( 'poisson' 'NB' 'Deconv' )
DA_meths=( 'Deconv' )
ct_seq=( '' 'Anterior Primitive Streak' 'Cardiopharyngeal Mesoderm' 'Caudal Epiblast' 'Caudal Neurectoderm' 'Early Anterior PSM' 'Early Nascent Mesoderm' 'Early Neurectoderm' 'Early Posterior PSM' 'Early Spinal Cord' 'Endothelium' 'Epiblast' 'Head Mesoderm' 'Late Anterior PSM' 'Late Nascent Mesoderm' 'Late Neurectoderm' 'Late Posterior PSM' 'Late Spinal Cord' 'Mature Endoderm' 'Mature Somites' 'NMPs' 'Neurons' 'Notochord' 'PGCs' 'Primitive Streak' 'Somites' )

for da_meth in "${DA_meths[@]}"
do
    for ct in "${ct_seq[@]}"
    do
        bsub -g /hetero_sims_load -o "$outbase"outfiles/loadfiles_"${da_meth}"_"${ct}".txt -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/load_ct_DA_files.R -o "${outbase}" -c "${ct}" -d "${da_meth}"
    done
done