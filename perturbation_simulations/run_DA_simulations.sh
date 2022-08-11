#!/usr/local/bin/bash

# conda activate /nfs/research/marioni/Leah/software/conda_envs/deconvolution
# cd /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/
# bsub -M 1000 ./run_DA_simulations.sh

bgadd -L 500 /hetero_sims

now=$(date '+%Y%m%d')
#outbase="/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/""$now""/"
outbase="/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/20220628/"
#mkdir -p $outbase
#mkdir -p "$outbase""/outfiles"

#simtypes=('poisson' 'gampois')
simtypes=()
#DA_meths=('poisson' 'NB')
DA_meths=()
for simtype in "${simtypes[@]}"
do
    for da_meth in "${DA_meths[@]}"
    do
        bsub -g /hetero_sims -o "$outbase"outfiles/general_sim_"${simtype}"_"${da_meth}".txt -M 5000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_general.R -f /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_functions.R -s "${simtype}" -d "${da_meth}" -p uniform -m unimodal -o "${outbase}"
    done
done


#tps=( 'd4' 'd4.5' )
tps=( 'd4.5' )
simtypes=( 'full' )
DA_meths=( 'poisson' 'NB' 'Deconv' )
IoD_seq=(  )
IoD_contr_seq=( 500 )
norg_seq=( 1 5 10 15 20 25 30 )
nsamps_seq=( 2 3 4 5 6 7 8 9 10 )
ct_seq=( '' 'Anterior Primitive Streak' 'Cardiopharyngeal Mesoderm' 'Caudal Epiblast' 'Caudal Neurectoderm' 'Early Anterior PSM' 'Early Nascent Mesoderm' 'Early Neurectoderm' 'Early Posterior PSM' 'Early Spinal Cord' 'Endothelium' 'Epiblast' 'Head Mesoderm' 'Late Anterior PSM' 'Late Nascent Mesoderm' 'Late Neurectoderm' 'Late Posterior PSM' 'Late Spinal Cord' 'Mature Endoderm' 'Mature Somites' 'NMPs' 'Neurons' 'Notochord' 'PGCs' 'Primitive Streak' 'Somites' )
fc_seq=( 0 0.25 0.5 0.75 )
ntrials=100

for simtype in "${simtypes[@]}"
do
    for da_meth in "${DA_meths[@]}"
    do
        for tp in "${tps[@]}"
        do
            PCA_rds=/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/"${tp}"_pca.RDS
            
            for IoD_contr in "${IoD_contr_seq[@]}"
            do
                for norg in "${norg_seq[@]}"
                do
                    for nsamps in "${nsamps_seq[@]}"
                    do
                        for (( t=1; t<=$ntrials; t++ ))
                        do
                            for fc in "${fc_seq[@]}"
                            do
                                for ct in "${ct_seq[@]}"
                                do
                                    if [ "${simtype}" == 'gampois' ] ; then
                                        for IoD in "${IoD_seq[@]}"
                                        do
                                            bsub -g /hetero_sims -q short -o "$outbase"outfiles/trial"${t}"_gastr_sim_"${simtype}"_"${tp}"_"${da_meth}"_"${ct}"_"${fc}"_"${IoD}"_"${nsamps}"_"${norg}"_"${IoD_contr}".txt -M 5000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_gastruloid.R -g /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/sample_metadata_final2.txt.gz -f /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_functions.R -s "${simtype}" -t "${tp}" -d "${da_meth}" -o "${outbase}"trial"${t}"_ -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/plotting_settings.R -i "${IoD}" -c "${IoD_contr}" -O "${norg}" -n "${nsamps}" -T 1 -P "${PCA_rds}" -F "${fc}" -C "${ct}"
                                        done
                                    bsub -g /hetero_sims -o "$outbase"outfiles/trial"${t}"_gastr_sim_"${simtype}"_"${tp}"_"${da_meth}"_"${ct}"_"${fc}"_"${nsamps}"_"${norg}"_"${IoD_contr}".txt -M 5000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_gastruloid.R -g /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/sample_metadata_final2.txt.gz -f /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_functions.R -s "${simtype}" -t "${tp}" -d "${da_meth}" -o "${outbase}"trial"${t}"_ -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/plotting_settings.R -c "${IoD_contr}" -O "${norg}" -n "${nsamps}" -T 1 -P "${PCA_rds}" -F "${fc}" -C "${ct}"
                                    fi
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done