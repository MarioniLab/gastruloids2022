#!/usr/bin/env Rscript

# /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline/evaluating/run_DEG.R -f /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline/processing/DEG_functions.R -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/sce.RDS -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/DE_out/new/ -g /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gene_metadata.txt.gz -T d5 -t d4.5 -G meso -g neural -C "Spinal Cord" -c "Spinal Cord"

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-f", "--deg_functions"), type="character", default=NULL, 
                help="path to DEG functions", metavar="character"),
    make_option(c("-s", "--sce"), type="character", default=NULL, 
                help="path to the sce file", metavar="character"),
    make_option(c("-M", "--gene_metadata"), type="character", default=NULL, 
                help="path to the gene_metadata file", metavar="character"),
    make_option(c("-T", "--tp1"), type="character", default=NULL, 
                help="timepoint 1", metavar="character"),
    make_option(c("-t", "--tp2"), type="character", default=NULL, 
                help="timepoint 2 (can be same as tp1)", metavar="character"),
    make_option(c("-G", "--gastr_type1"), type="character", default=NULL, 
                help="gastruloid type to choose at timepoint 1", metavar="character"),
    make_option(c("-g", "--gastr_type2"), type="character", default=NULL, 
                help="gastruloid type to choose at timepoint 2", metavar="character"),
    make_option(c("-C", "--ct1"), type="character", default=NULL, 
                help="cell type to choose at timepoint 1", metavar="character"),
    make_option(c("-c", "--ct2"), type="character", default=NULL, 
                help="cell type to choose at timepoint 2", metavar="character"),
    make_option(c("-o", "--outbase"), type="character", default=NULL, 
                help="basic outpath", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

###############
## Load Data ##
###############
 
source(opts$deg_functions)

sce <- readRDS(opts$sce)
colData(sce)$gastr_type[is.na(colData(sce)$gastr_type)] <- 'unknown'

gene_metadata <- fread(opts$gene_metadata)
gene_metadata <- gene_metadata[ens_id%in%rownames(sce)]
gene_metadata <- setnames(gene_metadata, c("symbol"), c("gene"))

if (opts$tp1 == 'All') opts$tp1 <- unique(sce$timepoint)
if (opts$tp2 == 'All') opts$tp2 <- unique(sce$timepoint)
if (opts$gastr_type1 == 'All') opts$gastr_type1 <- unique(sce$gastr_type)
if (opts$gastr_type2 == 'All') opts$gastr_type2 <- unique(sce$gastr_type)
if (opts$ct1 == 'All') opts$ct1 <- unique(sce$cluster)
if (opts$ct2 == 'All') opts$ct2 <- unique(sce$cluster)

sce$group <- 'None'
sce$group[(colData(sce)$timepoint %in% opts$tp1) & (colData(sce)$cluster %in% opts$ct1) & (colData(sce)$gastr_type %in% opts$gastr_type1)] <- 'group1'
sce$group[(colData(sce)$timepoint %in% opts$tp2) & (colData(sce)$cluster %in% opts$ct2) & (colData(sce)$gastr_type %in% opts$gastr_type2)] <- 'group2'
if ((any((colData(sce)$timepoint %in% opts$tp1) & (colData(sce)$cluster %in% opts$ct1) & (colData(sce)$gastr_type %in% opts$gastr_type1) & (colData(sce)$timepoint %in% opts$tp2) & (colData(sce)$cluster %in% opts$ct2) & (colData(sce)$gastr_type %in% opts$gastr_type2)))) {
    sce$group[(colData(sce)$timepoint %in% opts$tp1) &
              (colData(sce)$cluster %in% opts$ct1) &
              (colData(sce)$gastr_type %in% opts$gastr_type1) &
              (colData(sce)$timepoint %in% opts$tp2) &
              (colData(sce)$cluster %in% opts$ct2) &
              (colData(sce)$gastr_type %in% opts$gastr_type2)
             ] <- 'None'
    if ((length(opts$tp1) == 1) & (length(opts$tp2) > 1) &
        (((length(opts$gastr_type1) > 1) & (length(opts$gastr_type2) > 1)) | ((length(opts$gastr_type1) == 1) & (length(opts$gastr_type2) == 1))) &
        (((length(opts$ct1) > 1) & (length(opts$ct2) > 1)) | ((length(opts$ct1) == 1) & (length(opts$ct2) == 1)))
       ) {
        sce$group[(colData(sce)$timepoint %in% opts$tp1) & (colData(sce)$cluster %in% opts$ct1) & (colData(sce)$gastr_type %in% opts$gastr_type1)] <- 'group1'
    }
    if ((length(opts$tp2) == 1) & (length(opts$tp1) > 1) &
        (((length(opts$gastr_type1) > 1) & (length(opts$gastr_type2) > 1)) | ((length(opts$gastr_type1) == 1) & (length(opts$gastr_type2) == 1))) &
        (((length(opts$ct1) > 1) & (length(opts$ct2) > 1)) | ((length(opts$ct1) == 1) & (length(opts$ct2) == 1)))
       ) {
        sce$group[(colData(sce)$timepoint %in% opts$tp2) & (colData(sce)$cluster %in% opts$ct2) & (colData(sce)$gastr_type %in% opts$gastr_type2)] <- 'group2'
    }
    if ((length(opts$gastr_type1) == 1) & (length(opts$gastr_type2) > 1) &
        (((length(opts$tp1) > 1) & (length(opts$tp2) > 1)) | ((length(opts$tp1) == 1) & (length(opts$tp2) == 1))) &
        (((length(opts$ct1) > 1) & (length(opts$ct2) > 1)) | ((length(opts$ct1) == 1) & (length(opts$ct2) == 1)))
       ) {
        sce$group[(colData(sce)$timepoint %in% opts$tp1) & (colData(sce)$cluster %in% opts$ct1) & (colData(sce)$gastr_type %in% opts$gastr_type1)] <- 'group1'
    }
    if ((length(opts$gastr_type2) == 1) & (length(opts$gastr_type1) > 1) &
        (((length(opts$tp1) > 1) & (length(opts$tp2) > 1)) | ((length(opts$tp1) == 1) & (length(opts$tp2) == 1))) &
        (((length(opts$ct1) > 1) & (length(opts$ct2) > 1)) | ((length(opts$ct1) == 1) & (length(opts$ct2) == 1)))
       ) {
        sce$group[(colData(sce)$timepoint %in% opts$tp2) & (colData(sce)$cluster %in% opts$ct2) & (colData(sce)$gastr_type %in% opts$gastr_type2)] <- 'group2'
    }
    if ((length(opts$ct1) == 1) & (length(opts$ct2) > 1) &
        (((length(opts$gastr_type1) > 1) & (length(opts$gastr_type2) > 1)) | ((length(opts$gastr_type1) == 1) & (length(opts$gastr_type2) == 1))) &
        (((length(opts$tp1) > 1) & (length(opts$tp2) > 1)) | ((length(opts$tp1) == 1) & (length(opts$tp2) == 1)))
       ) {
        sce$group[(colData(sce)$timepoint %in% opts$tp1) & (colData(sce)$cluster %in% opts$ct1) & (colData(sce)$gastr_type %in% opts$gastr_type1)] <- 'group1'
    }
    if ((length(opts$ct2) == 1) & (length(opts$ct1) > 1) &
        (((length(opts$gastr_type1) > 1) & (length(opts$gastr_type2) > 1)) | ((length(opts$gastr_type1) == 1) & (length(opts$gastr_type2) == 1))) &
        (((length(opts$tp1) > 1) & (length(opts$tp2) > 1)) | ((length(opts$tp1) == 1) & (length(opts$tp2) == 1)))
       ) {
        sce$group[(colData(sce)$timepoint %in% opts$tp2) & (colData(sce)$cluster %in% opts$ct2) & (colData(sce)$gastr_type %in% opts$gastr_type2)] <- 'group2'
    }
}

if (((sum(sce$group == 'group1')>0) & (sum(sce$group == 'group2')>0)) &
    (any(colData(sce)$experiment[colData(sce)$group=='group1'] %in% colData(sce)$experiment[colData(sce)$group=='group2'])) &
    (sum(colData(sce)$group=='group1')>2) &
    (sum(colData(sce)$group=='group2')>2)
   ) {
    sce <- sce[,sce$group %in% c('group1', 'group2')]

    sce$batch <- colData(sce)$experiment
    groups <- unique(sce$group)
    sce$group <- factor(sce$group, levels = groups)
    out <- do_DE(sce, gene_metadata, min_detection_rate_per_group=0.1, max_detection_rate_per_group=1, min.logFC=2, threshold_fdr=0.1)

    if (length(opts$tp1) > 1) opts$tp1 <- 'All'
    if (length(opts$tp2) > 1) opts$tp2 <- 'All'
    if (length(opts$gastr_type1) > 1) opts$gastr_type1 <- 'All'
    if (length(opts$gastr_type2) > 1) opts$gastr_type2 <- 'All'
    if (length(opts$ct1) > 1) opts$ct1 <- 'All'
    if (length(opts$ct2) > 1) opts$ct2 <- 'All'

    fwrite(out, paste0(opts$outbase, opts$tp1, '_', opts$gastr_type1, '_', opts$ct1, '_vs_', opts$tp2, '_', opts$gastr_type2, '_', opts$ct2, '_out.csv.gz'))
    
    out_to.plot <- out
    out_to.plot[p.value==1]$p.value <- max(out_to.plot[p.value!=1]$p.value, na.rm=TRUE)
    out_to.plot[p.value==0]$p.value <- min(out_to.plot[p.value!=0]$p.value, na.rm=TRUE)
    setorder(out_to.plot, -sig, padj_fdr, na.last=T)
    p <- gg_volcano_plot(out_to.plot, top_genes=10)
    p <- p + ggtitle(paste0(opts$tp1, ' ', opts$gastr_type1, ' ', opts$ct1, ' <--> ', opts$tp2, ' ', opts$gastr_type2, ' ', opts$ct2))

    ggsave(paste0(opts$outbase, opts$tp1, '_', opts$gastr_type1, '_', opts$ct1, '_vs_', opts$tp2, '_', opts$gastr_type2, '_', opts$ct2, "_volcano.pdf"),
           plot = p,
           width = 10,
           height = 5,
           units = "in",
           device='pdf'
          )
}