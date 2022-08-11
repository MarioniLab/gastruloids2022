#!/usr/bin/env Rscript

# conda activate basic_renv
# bsub -M 200000 -q bigmem /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/pseudobulk.R -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_QCed.rds -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(Matrix.utils))
suppressPackageStartupMessages(library(optparse))


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-s", "--srat"), type="character", default=NULL, 
                help="path to the srat file", metavar="character"),
    make_option(c("-i", "--iters"), type="integer", default=10000, 
                help="the number of times to simulate the gastruloids", metavar="number"),
    make_option(c("-o", "--outbase"), type="character", default="All", 
                help="basic outpath", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);


###############
## Load Data ##
###############
 
#mapping <- readRDS(opts$mapping)
#day <- "d5"

srat <- readRDS(opts$srat)
srat <- subset(srat, subset = (MULTI_class != 'Negative') & (MULTI_class != 'Doublet') & (MULTI_class != 'unidentifiable') & (tech=="MULTI"))
exps <- unique(srat@meta.data$experiment)
days <- unique(srat@meta.data$timepoint)
gastruloids <- unique(srat@meta.data$MULTI_class)
for (exp in exps) {
    
    for (day in days) {
        
        if (any((srat@meta.data$experiment == exp) & (srat@meta.data$timepoint == day))) {
        
            srat_sub <- subset(srat, subset = (timepoint==day) & (experiment==exp))

            # Create single cell experiment object
            sce <- SingleCellExperiment(assays = list(counts = srat_sub@assays$RNA@counts[rownames(srat_sub@assays$RNA@meta.features)[!(srat@assays$RNA@meta.features$cc_gene)],]), 
                                       colData = srat_sub@meta.data)

            # Aggregate the counts per sample_id and cluster_id

            # Subset metadata to only include the cluster and sample IDs to aggregate across
            groups <- colData(sce)[, c("MULTI_class")]

            # Aggregate across cluster-sample groups
            pb <- aggregate.Matrix(t(counts(sce)), 
                                   groupings = groups, fun = "sum") 

            write.csv(as.matrix(pb), file=paste0(opts$outbase, day, "_", exp, "_pb.csv"))

            pb_random <- list()
            for (i in 1:opts$iters) {

                # Aggregate across cluster-sample groups randomly
                pb_random[[i]] <- aggregate.Matrix(t(counts(sce)), 
                                       groupings = sample(groups), fun = "sum")

            }
            saveRDS(pb_random, paste0(opts$outbase, day, "_", exp, "_pb_simulation.csv"))
            
        }
        
    }

}






































