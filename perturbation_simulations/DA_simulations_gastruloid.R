#!/usr/bin/env Rscript

# /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_gastruloid.R -g /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/plot_df_joined.txt.gz -P /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/d4.5_pca.RDS -f /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/DA_simulations_functions.R -s full -t d4.5 -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/ -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/plotting_settings.R -F 0.5 -C "Early Spinal Cord" -c 500 

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(DirichletReg))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(abind))


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-g", "--g_meta"), type="character", default=NULL, 
                help="gastruloid metadata", metavar="character"),
    make_option(c("-P", "--pca"), type="character", default=NULL, 
                help="pca rds", metavar="character"),
    make_option(c("-f", "--functions"), type="character", default=NULL, 
                help="DA simulation functions", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="plotting settings for gastr_type annotations", metavar="character"),
    make_option(c("-s", "--sim_type"), type="character", default="poisson", 
                help="simulation type (full, poisson, gampois, multinomial, or a gastruloid class)", metavar="character"),
    make_option(c("-t", "--timepoint"), type="character", default="d5", 
                help="timepoint", metavar="character"),
    make_option(c("-o", "--outbase"), type="character", default=NULL, 
                help="basic outpath", metavar="character"),
    #make_option(c("-p", "--param_tovary"), type="character", default=NULL, 
    #            help="which paramter to vary (either ng or ns)", metavar="character"),
    make_option(c("-d", "--DA_meth"), type="character", default=NULL, 
                help="which method to use for DA (poisson, NB, or deconv)", metavar="character"),
    make_option(c("-F", "--fc_perturb"), type="double", default=0.5, 
                help="fold change of the perturbed celltype", metavar="double"),
    make_option(c("-C", "--ct_perturb"), type="character", default=NULL, 
                help="cell type to be perturbed", metavar="character"),
    make_option(c("-i", "--IoD"), type="double", default=NA, 
                help="IoD of gamma prior", metavar="double"),
    make_option(c("-c", "--IoD_contr"), type="double", default=NULL, 
                help="IoD for how much each gastruloid contributes to sample", metavar="double"),
    make_option(c("-O", "--norg"), type="integer", default=NULL, 
                help="how many organoids per sample", metavar="integer"),
    make_option(c("-n", "--nsamps"), type="integer", default=NULL, 
                help="how many samples per condition", metavar="integer"),
    make_option(c("-T", "--ntrials"), type="integer", default=100, 
                help="how many trials", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

source(opts$functions)
source(opts$settings)
g_meta <- fread(opts$g_meta)

if (!is.null(opts$pca)) {
    pca <- readRDS(opts$pca)
} else {
    pca <- do_pca(df_ref=g_meta, tp=opts$timepoint, io=io)
}

####################
## Run Simulation ##
####################

fc_ct = rep(1, length(pca$ct_ordering))
if (!is.null(opts$ct_perturb)) fc_ct[pca$ct_ordering == opts$ct_perturb] <- opts$fc_perturb

if ((opts$DA_meth == 'NB') & (opts$nsamps == 1)) {
    stop('Need more than 1 sample for NBinom')
} else {
    out <- sim_DA_analysis_gastruloid(g_meta,
                                      fc_ct_list = list(fc_ct),
                                      tp=opts$timepoint,
                                      gastr_type_subset=NULL,
                                      cells_per_samp=10000,
                                      #n_samps_percond=c(3),
                                      n_samps_percond=c(opts$nsamps),
                                      #ng_seq=1:50,
                                      ng_seq=c(opts$norg),
                                      #IoD_seq = 2^seq(1:10),
                                      IoD_seq = c(opts$IoD),
                                      #IoD_contr_seq = c(0, 2^seq(0, 5, 2)),
                                      IoD_contr_seq = c(opts$IoD_contr),
                                      DA_meth = opts$DA_meth,
                                      trials=opts$ntrials,
                                      sim_type=opts$sim_type,
                                      norm.method='TMM',
                                      pca=pca,
                                      io=io
                                     )
}

fwrite(out$out, paste0(opts$outbase, 'gastr_sim_', opts$sim_type, '_', opts$timepoint, '_', opts$DA_meth, '_', opts$ct_perturb, '_', opts$fc_perturb, '_', opts$IoD, '_', opts$nsamps, '_', opts$norg, '_', opts$IoD_contr, '_out.txt.gz'))
fwrite(out$IoDs, paste0(opts$outbase, 'gastr_sim_', opts$sim_type, '_', opts$timepoint, '_', opts$DA_meth, '_', opts$ct_perturb, '_', opts$fc_perturb, '_', opts$IoD, '_', opts$nsamps, '_', opts$norg, '_', opts$IoD_contr, '_IoDout.txt.gz'))
