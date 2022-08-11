#!/usr/local/bin/bash

# run in scvelo conda environment
# call using: /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/scvelo/run_merge_chroms.sh -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/scvelo/merge_chroms.py

while getopts s: flag
do
    case "${flag}" in
        s) script=${OPTARG};;
    esac
done

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/EB_d3/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/EB_d3/sample_metadata.txt.gz

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/EB_d4/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/EB_d4/sample_metadata.txt.gz

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/EB_d5/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/EB_d5/sample_metadata.txt.gz

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/ESC/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/ESC/sample_metadata.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/gastr_d3/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastr_d3/sample_metadata.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/gastr_d4/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastr_d4/sample_metadata.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/human_gastr/exp1_gastr_72h/SLX20571_SITTG10_human_gastr_Yang_72h_MULTIseq_A/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/human_gastr/scRNAseq/sample_metadata_nominmito.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/human_gastr/exp1_gastr_72h/SLX20571_SITTH10_human_gastr_Yang_72h_MULTIseq_B/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/human_gastr/scRNAseq/sample_metadata_nominmito.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp1_d5/CellRanger6_AGAAACTC_10x_MULTI_gastr_d5_sample1/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp1_d5/sample1_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp1_d5/CellRanger6_ACTTAGGA_10x_MULTI_gastr_d5_sample2/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp1_d5/sample2_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp1_d5/CellRanger6_AAAGAAGA_10x_MULTI_gastr_d5_sample3/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp1_d5/sample3_sample_metadata.txt.gz

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp2A_d3_d3.5/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp2A_d3_d3.5/sample_metadata.txt.gz

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp2C_d3_d3.5/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp2C_d3_d3.5/sample_metadata.txt.gz

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp2D_d3_d3.5/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp2D_d3_d3.5/sample_metadata.txt.gz

#bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3A_d4_d4.5/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp3A_d4_d4.5/sample_metadata.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3B_d4_d4.5/d4_B2_sample13/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp3B_d4_d4.5/sample2_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3B_d4_d4.5/d4_B3_sample14/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp3B_d4_d4.5/sample3_sample_metadata.txt.gz


bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3C_d4_d4.5/d4_C1_sample15/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp3C_d4_d4.5/sample1_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3C_d4_d4.5/d4_C2_sample16/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp3C_d4_d4.5/sample2_sample_metadata.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp4_d3/SITTA12_gastr_d3_MULTIseq_3prime_RNA_A/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp4_d3/sampleA_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp4_d3/SITTB12_gastr_d3_MULTIseq_3prime_RNA_B/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp4_d3/sampleB_sample_metadata.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTC12_gastr_d3_5_d4_MULTIseq_3prime_RNA_A/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp5_d3.5_d4/sampleA_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTD12_gastr_d3_5_d4_MULTIseq_3prime_RNA_B/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp5_d3.5_d4/sampleB_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTE12_gastr_d3_5_d4_MULTIseq_3prime_RNA_C/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp5_d3.5_d4/sampleC_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTF12_gastr_d3_5_d4_MULTIseq_3prime_RNA_D/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp5_d3.5_d4/sampleD_sample_metadata.txt.gz

bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp6_d4.5_d5/SITTE10_gastr_d4_5_d5_MULTIseq_3prime_RNA_A/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp6_d4.5_d5/sampleA_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp6_d4.5_d5/SITTF10_gastr_d4_5_d5_MULTIseq_3prime_RNA_B/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp6_d4.5_d5/sampleB_sample_metadata.txt.gz
bsub -M 15000 "$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp6_d4.5_d5/SITTG10_gastr_d4_5_d5_MULTIseq_3prime_RNA_C/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp6_d4.5_d5/sampleC_sample_metadata.txt.gz

