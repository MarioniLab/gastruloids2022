#!/usr/local/bin/bash

# conda activate basic_renv
# cd /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline_final/DE/
# bsub -o submitting_outfile1.txt -M 1000 ./run_all_DEG_gastrtype_tp.sh -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline_final/DE/run_DEG_gastrtype_tp.R -f '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline_final/DE/DEG_functions.R' -d '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/sce_final.rds' -o '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/DE_out/' -g '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/feature_metadata_final.txt.gz'

while getopts s:d:f:o:g: flag
do
    case "${flag}" in
        s) script=${OPTARG};;
        d) sce=${OPTARG};;
        f) functions=${OPTARG};;
        o) outbase=${OPTARG};;
        g) genemeta=${OPTARG};;
    esac
done

mkdir -p "$outbase"outfiles/

bgadd -L 4000 /DEfull

#celltypes=("All" 'Mature Somites')
celltypes=("All" "PGCs" "Epiblast" "Primitive Streak" "Anterior Primitive Streak" "Mature Endoderm" "Notochord" "Early Nascent Mesoderm" "Late Nascent Mesoderm" "Head Mesoderm" "Cardiopharyngeal Mesoderm" "Endothelium" "Early Posterior PSM" "Late Posterior PSM" "Early Anterior PSM" "Late Anterior PSM" "Somites" "Mature Somites" "Caudal Epiblast" "NMPs" "Caudal Neurectoderm" "Early Neurectoderm" "Late Neurectoderm" "Early Spinal Cord" "Late Spinal Cord" "Neurons")
#gastrtypes=("All" 'mesodermal')
gastrtypes=('All' 'mesodermal' 'neural' 'intermediate')
#tps=("All" 'd3')
tps=('All' 'd3' 'd3.5' 'd4' 'd4.5' 'd5')
for iC in "${!celltypes[@]}"
do
    for ic in "${!celltypes[@]}"
    do
        for iG in "${!gastrtypes[@]}"
        do
            for ig in "${!gastrtypes[@]}"
            do
                for iT in "${!tps[@]}"
                do
                    for it in "${!tps[@]}"
                    do
                        if ( [ "$iC" -le "$ic" ] && [ "$iG" -le "$ig" ] && [ "$iT" -le "$it" ] ) && ( [ $(("$it" - "$iT")) -le 1 ] || [ "${tps[$iT]}" == 'All' ] ) && ! ( [ "$iC" -eq "$ic" ] && [ "$iG" -eq "$ig" ] && [ "$iT" -eq "$it" ] ) ; then # Make sure we don't duplicate, make sure subsequent timepoints, and make sure not all equal
                            #bsub -g /DEfull -o "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt -M 20000 "$script" -f "$functions" -s "$sce" -o "$outbase" -M "$genemeta" -T "${tps[$iT]}" -t "${tps[$it]}" -G "${gastrtypes[$iG]}" -g "${gastrtypes[$ig]}" -C "${celltypes[$iC]}" -c "${celltypes[$ic]}"
                            if grep -q "TERM_MEMLIMIT"  "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt; then
                                rm -f "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt
                                bsub -g /DEfull -o "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt -M 100000 "$script" -f "$functions" -s "$sce" -o "$outbase" -M "$genemeta" -T "${tps[$iT]}" -t "${tps[$it]}" -G "${gastrtypes[$iG]}" -g "${gastrtypes[$ig]}" -C "${celltypes[$iC]}" -c "${celltypes[$ic]}"
                            #elif grep -q "Error" "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt; then # checks for both an error, and if it just wasn't run yet
                                #rm -f "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt
                                #bsub -g /DEfull -o "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt -M 20000 "$script" -f "$functions" -s "$sce" -o "$outbase" -M "$genemeta" -T "${tps[$iT]}" -t "${tps[$it]}" -G "${gastrtypes[$iG]}" -g "${gastrtypes[$ig]}" -C "${celltypes[$iC]}" -c "${celltypes[$ic]}"
                            #elif ! grep -q "Successfully completed" "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt; then # checks for both an error, and if it just wasn't run yet
                                #rm -f "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt
                                #bsub -g /DEfull -o "$outbase"outfiles/"${tps[$iT]}"_"${gastrtypes[$iG]}"_"${celltypes[$iC]}"_vs_"${tps[$it]}"_"${gastrtypes[$ig]}"_"${celltypes[$ic]}".txt -M 20000 "$script" -f "$functions" -s "$sce" -o "$outbase" -M "$genemeta" -T "${tps[$iT]}" -t "${tps[$it]}" -G "${gastrtypes[$iG]}" -g "${gastrtypes[$ig]}" -C "${celltypes[$iC]}" -c "${celltypes[$ic]}"
                            #else
                            fi
                        fi
                    done
                done
            done
        done
    done
done
