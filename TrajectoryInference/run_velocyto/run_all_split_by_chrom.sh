#!/usr/local/bin/bash

# run in scvelo conda environment
# call using: /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/scvelo/run_all_split_by_chrom.sh -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/scvelo/run_split_by_chrom.sh -r /nfs/research/marioni/common/references/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf

while getopts s:r: flag
do
    case "${flag}" in
        s) script=${OPTARG};;
        r) reffile=${OPTARG};;
    esac
done

"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/EB_d3/ -r "$reffile" -n

"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/EB_d4/ -r "$reffile" -n

"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/EB_d5/ -r "$reffile" -n

"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/ESC/ -r "$reffile" -n

"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/gastr_d3/ -r "$reffile" -n

"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/gastr_d4/ -r "$reffile" -n

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/human_gastr/exp1_gastr_72h/SLX20571_SITTG10_human_gastr_Yang_72h_MULTIseq_A/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/human_gastr/exp1_gastr_72h/SLX20571_SITTH10_human_gastr_Yang_72h_MULTIseq_B/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp1_d5/CellRanger6_AGAAACTC_10x_MULTI_gastr_d5_sample1/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp1_d5/CellRanger6_ACTTAGGA_10x_MULTI_gastr_d5_sample2/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp1_d5/CellRanger6_AAAGAAGA_10x_MULTI_gastr_d5_sample3/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp2A_d3_d3.5/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp2C_d3_d3.5/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp2D_d3_d3.5/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3A_d4_d4.5/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3B_d4_d4.5/d4_B2_sample13/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3B_d4_d4.5/d4_B3_sample14/ -r "$reffile"


#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3C_d4_d4.5/d4_C1_sample15/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp3C_d4_d4.5/d4_C2_sample16/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp4_d3/SITTA12_gastr_d3_MULTIseq_3prime_RNA_A/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp4_d3/SITTB12_gastr_d3_MULTIseq_3prime_RNA_B/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTC12_gastr_d3_5_d4_MULTIseq_3prime_RNA_A/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTD12_gastr_d3_5_d4_MULTIseq_3prime_RNA_B/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTE12_gastr_d3_5_d4_MULTIseq_3prime_RNA_C/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp5_d3.5_d4/SITTF12_gastr_d3_5_d4_MULTIseq_3prime_RNA_D/ -r "$reffile"

#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp6_d4.5_d5/SITTE10_gastr_d4_5_d5_MULTIseq_3prime_RNA_A/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp6_d4.5_d5/SITTF10_gastr_d4_5_d5_MULTIseq_3prime_RNA_B/ -r "$reffile"
#"$script" -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp6_d4.5_d5/SITTG10_gastr_d4_5_d5_MULTIseq_3prime_RNA_C/ -r "$reffile"


