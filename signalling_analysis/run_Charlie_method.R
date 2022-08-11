#!/usr/bin/env Rscript

# conda activate basic_renv
# /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/d3_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -d {day} -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/{day}_{exp}_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/{day_exp}/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d3_exp4_d3_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp4_d3 -d d3 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d3_exp4_d3_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d3_exp4_d3/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d3.5_exp5_d3.5_d4_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp5_d3.5_d4 -d d3.5 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d3.5_exp5_d3.5_d4_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d3.5_exp5_d3.5_d4/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4_exp5_d3.5_d4_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp5_d3.5_d4 -d d4 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4_exp5_d3.5_d4_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d4_exp5_d3.5_d4/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4_exp3B_d4_d4.5_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp3B_d4_d4.5 -d d4 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4_exp3B_d4_d4.5_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d4_exp3B_d4_d4.5/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4.5_exp3B_d4_d4.5_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp3B_d4_d4.5 -d d4.5 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4.5_exp3B_d4_d4.5_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d4.5_exp3B_d4_d4.5/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4.5_exp6_d4.5_d5_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp6_d4.5_d5 -d d4.5 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d4.5_exp6_d4.5_d5_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d4.5_exp6_d4.5_d5/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d5_exp6_d4.5_d5_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp6_d4.5_d5 -d d5 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d5_exp6_d4.5_d5_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d5_exp6_d4.5_d5/
#bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/run_Charlie_method.R -c /nfs/research/marioni/Leah/data/cell_cycle_genes/all_mouse_cc.txt.gz -p /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d5_exp1_d5_pb.csv -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/seurat_forCarine.rds -e exp1_d5 -d d5 -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/settings.R -S /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gastruloid_pseudobulk/d5_exp1_d5_pb_simulation.csv -w /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/signalling_eval/CharlieMethod/get_wgcna.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/CharlieMethodOut/d5_exp1_d5/


# Adpated from: https://gitlab.ebi.ac.uk/petsalakilab/phenotype_networks/-/blob/master/NAFLD/lm_module.R
# Barker et al., bioRxiv 2021 (https://www.biorxiv.org/content/10.1101/2021.02.11.430597v2)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(flashClust))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sctransform))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scWGCNA))
suppressPackageStartupMessages(library(gprofiler2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(enrichR))


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-c", "--cc_genes"), type="character", default=NULL, 
                help="path to cell cycle gene table", metavar="character"),
    make_option(c("-p", "--pb_data"), type="character", default=NULL, 
                help="path to the pseudobulked table", metavar="character"),
    make_option(c("-m", "--srat"), type="character", default=NULL, 
                help="path to mapping output RDS", metavar="character"),
    make_option(c("-e", "--experiment"), type="character", default=NULL, 
                help="experiment to use", metavar="character"),
    make_option(c("-d", "--day"), type="character", default=NULL, 
                help="differentiation day to use", metavar="character"),
    make_option(c("-s", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-o", "--outbase"), type="character", default=NULL, 
                help="out directory", metavar="character"),
    make_option(c("-S", "--sims"), type="character", default=NULL, 
                help="RDS containing the output of the random pseudobulking", metavar="character"),
    make_option(c("-w", "--wgcna_functions"), type="character", default=NULL, 
                help="path to the get_wgcna.R script", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);


###############
## Load Data ##
###############

if (!dir.exists(opts$outbase)) {
    dir.create(opts$outbase, recursive = TRUE)
}

source(opts$settings)
source(opts$wgcna_functions)
ccgenes <- fread(opts$cc_genes)

#### Pseudobulk Data ####

GeneXData <- read.csv(opts$pb_data,
                      row.names=1
                     ) %>%
    .[!(rownames(.) %in% c('Negative', 'Doublet')),] %>%
    as.matrix %>%
    t %>%
    .[!(rownames(.) %in% ccgenes$ensembl_gene_id),]
GeneXData_cpm <- apply(GeneXData,2, function(x) {(x/sum(x))*1000000})
gene.IDs <- rownames(GeneXData_cpm)
# format properly for WGCNA
GeneXDiff <- GeneXData_cpm
GeneXData <- as.data.frame(t(GeneXData_cpm))
# LOG2 Transformation
logXdat <- log2(GeneXData + 1)
                       
# Run this to check if there are gene outliers
gsg = goodSamplesGenes(logXdat, verbose = 3)

if (!gsg$allOK) {
    if (sum(!gsg$goodGenes)>0) {
        message(paste("Removing genes:", paste(names(logXdat)[!gsg$goodGenes], collapse= ", ")));
    }
    if (sum(!gsg$goodSamples)>0) {
        message(paste("Removing samples:", paste(rownames(logXdat)[!gsg$goodSamples], collapse=", ")))
    }
  logXdat= logXdat[gsg$goodSamples, gsg$goodGenes]
}

#### Single Cell Data ####

srat <- readRDS(opts$srat)[colnames(logXdat),]
exps <- unique(srat@meta.data$experiment)
exps <- exps[!(exps %in% c("exp2A_d3_d3.5", "exp2C_d3_d3.5"))]

#atlas.meta <- fread("/nfs/research/marioni/Leah/gastrulation10x/data/Leah/atlas_metadata.txt.gz")[(stripped==FALSE) & (doublet==FALSE)]

rownames(srat@meta.data) <- srat@meta.data$cell
srat@meta.data <- srat@meta.data[Cells(srat),]

srat <- subset(srat, experiment==opts$experiment)
if (!(is.null(opts$day))) {
    srat <- subset(srat, timepoint==opts$day)
}


##############################
## Preprocess and Plot srat ##
##############################

srat <- NormalizeData(srat, verbose=FALSE)
srat@assays$RNA@scale.data <- as.matrix(srat@assays$RNA@data, verbose=FALSE)
srat <- FindVariableFeatures(srat, verbose=FALSE)
srat <- RunPCA(srat, npcs = 50, verbose=FALSE)    
srat <- FindNeighbors(srat, dims = 1:50, reduction = "pca")
srat <- FindClusters(srat, resolution = 1)
srat <- RunUMAP(srat, dims = 1:50, verbose=FALSE)

p <- DimPlot(srat, reduction = "umap")
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_cluster_umap.pdf"),
       plot = p,
       width = 7,
       height = 5,
       units = "in",
       device='pdf'
      )

#p <- DimPlot(srat, reduction = "umap", group.by="celltype.mapped", cols = celltype_colours)
#ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_celltype.mapped_umap.pdf"),
#       plot = p,
#       width = 15,
#       height = 5,
#       units = "in",
#       device='pdf'
#      )

p <- DimPlot(srat, reduction = "umap", group.by="batch")
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_sample_umap.pdf"),
       plot = p,
       width = 7,
       height = 5,
       units = "in",
       device='pdf'
      )

p <- DimPlot(srat, reduction = "umap", group.by="Phase")
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_Phase_umap.pdf"),
       plot = p,
       width = 7,
       height = 5,
       units = "in",
       device='pdf'
      )

p <- DimPlot(srat, reduction = "umap", group.by = 'cluster', cols=cluster_colours)
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_newcelltype_umap.pdf"),
       plot = p,
       width = 15,
       height = 5,
       units = "in",
       device='pdf'
      )


#################
## Run scWGCNA ##
#################

tmp <- srat
tmp = Seurat::FindNeighbors(srat,
                     reduction = "pca",
                     dims = 1:20,
                     k.param = 10)
nn.matrix = tmp@graphs[[paste0(Seurat::DefaultAssay(tmp),"_nn")]]

seeds <- 0.2

# To keep count
my.seeds = list()
seeds.count = c()

# Here, to check which set of randomly selected seed we'll use
if (is.numeric(0.2)) {

message("Choosing seeds")

for (i in 1:50) {
  seed.set = c()

  #go trhough each cluster or ID
  for (cluster in levels(Seurat::Idents(tmp)) ) {

    # How many seeds? According to proportion
    n.seeds = floor(table(Seurat::Idents(tmp))[[cluster]]/(1/seeds))
    # Choose the seeds
    seed.set=c(seed.set, sample(rownames(tmp@meta.data)[Seurat::Idents(tmp) == cluster], n.seeds) )
  }

  # Keep count of the seeds
  my.seeds[[i]] = seed.set
  rm(seed.set)

  #How many cells would be aggragated using this particular set of seeds?
  seeds.count = c(seeds.count, length(which(Matrix::colSums(nn.matrix[my.seeds[[i]],]) > 0)) )

}

# Choose the one with the highest count
seeds = my.seeds[[which.max(seeds.count)]]

}

#Subset the nn matrix, to only keep the seeds
nn.matrix = nn.matrix[seeds, ]
#Only keep the cells that would be aggregated
nn.matrix = nn.matrix[,Matrix::colSums(nn.matrix) > 0 ]

message(ncol(nn.matrix)," out of ", ncol(tmp)," Cells will be agreggated into ",nrow(nn.matrix)," Pseudocells")

#create a data frame, to keep record of cells and pseudocells. Non-asigned cells are 00
my.pseudocells = data.frame(pseudocell = rep("00",nrow(tmp@meta.data)), row.names = rownames(tmp@meta.data))

my.pseudocells$pseudocell = as.character(my.pseudocells$pseudocell)

# First, each seed is assigned to iself
my.pseudocells[seeds,1] = seeds
# We remove the seeds from the universe to choose from
nn.matrix = nn.matrix[,-which(colnames(nn.matrix) %in% seeds)]
# Which cells are up for grabs?
remaining.cells = colnames(nn.matrix)

message("Assining pseudocells")

while (length(remaining.cells) > 0) {
    #we go seed by seed, starting by the poorest seeds
    for (s in rownames(nn.matrix)[order(Matrix::rowSums(nn.matrix[,remaining.cells, drop = F]))]) {
        # If the seed still has options to choose from
        if ( sum(nn.matrix[s, remaining.cells]) > 0 & length(remaining.cells) > 0) {
            # take one of the seed nn
            c = sample(which(nn.matrix[s, remaining.cells] == 1),1)
            # Assign that nn to the seed
            my.pseudocells[remaining.cells[c],1] = s
            # Remove the nn from the universe to choose from
            remaining.cells = remaining.cells[-c]
        }
        # If all cells are assigned, stop
        if (length(remaining.cells) == 0) {break}
    }
}

# Make a new seurat object, using only the cells we assigned
ps.seurat = subset(tmp, cells = rownames(my.pseudocells)[my.pseudocells$pseudocell != "00"])
# Add a new metadata column, with the pseudocell (seeds) info and set the identity
ps.seurat = Seurat::AddMetaData(object = ps.seurat, metadata = my.pseudocells, col.name = "pseudo.ident")
ps.seurat = Seurat::SetIdent(ps.seurat, value = "pseudo.ident")

message("Aggregating pseudocell expression")

# Aggregate the cells using seurat. get the average expression
ps.seurat = Seurat::AverageExpression(object = ps.seurat,
                                      return.seurat = T,
                                      verbose = F,
                                      group.by="pseudo.ident")

# Rescue the original cluster of each pseudocell!
ps.seurat@meta.data$orig.cluster = Seurat::Idents(tmp)[match(rownames(ps.seurat@meta.data), rownames(tmp@meta.data))]

# Add our records of cells to the misc slot
ps.seurat@misc[["peudocells"]] = my.pseudocells

gnames <- gconvert(query = VariableFeatures(srat), 
                         organism = "mmusculus", 
                         target="ENSG", 
                         mthreshold = Inf, 
                         filter_na = FALSE)[,c("input","name")]

params = list(
  #Parse the parameters
  data = ps.seurat,
  sc.data = srat,
  gene.names = gnames,
  project.name = "gastr_d3_signalling",
  sp = "Mm",
  cells = FALSE,
  features = VariableFeatures(srat),
  reduction = "umap",
  dir = opts$outbase,
  is.pseudocell = TRUE,
  GO=FALSE,
  min.cell = 10
)

# placeholder for the variables and options we will get from an upper level pipeline
project_name = params$project.name
p.Wdata = params$sc.data
gnames = params$gene.names
gnames[,2]= as.character(gnames[,2])
rownames(gnames) = gnames[,1]
Wdata = params$data
# The subset in form of cell identities, F if using the whole sample
my.subset = params$cells
# Variable genes, or genes to use, F if the genes will be calculated in the script
my.vargenes = params$features
my.dr = params$reduction
my.sp = params$sp
is.pseudocell = params$is.pseudocell
my.min.cell = params$min.cell

nonex = which(apply(p.Wdata@assays$RNA@counts, 1, function(x) {length(which(x >0))}) < my.min.cell)

# First get rid of non-expressed genes
p.Wdata = subset(p.Wdata, features = rownames(p.Wdata@assays$RNA)[-nonex])
# Find the variable genes
p.Wdata=Seurat::SCTransform(p.Wdata, variable.features.n=6000, verbose=FALSE)
#Seurat::VariableFeaturePlot(p.Wdata, assay = "RNA")
Expr = Seurat::VariableFeatures(object = p.Wdata, assay = "SCT")
#Expr = my.vargenes
#if(length(nonex)>0) {Expr <- rownames(p.Wdata@assays$RNA)[-nonex]} else {Expr <- rownames(p.Wdata@assays$RNA)}

#### Let's check which genes are expressed in only one cell (less than 3 actually). To avoid one-celled modules. This also applies to genes that have in general very low expression, and one (three) outlier(s) with high expression
my.oc=which(apply(p.Wdata@assays$RNA@counts[Expr,], 1, function(x) { length(which(x > (max(x)/3))) } ) < 3)
if (length(my.oc) > 0) {
  print(paste0("The following genes were only highly expressed in less than 3 cells :  ",
               Expr[my.oc] ))
  Expr=Expr[-my.oc]
}

# We had a problem, where some genes end up not expressed in the pseudocells. I fixed this in this way: 
datExpr=Wdata@assays$RNA@counts[Expr,]
datExpr = datExpr[which(Matrix::rowSums(datExpr)>0),]
if (length(which(Matrix::rowSums(datExpr)==0))<1) {
    print(paste0("The following variable genes were not found expressed in the pseudocell object:  ", names(which(Matrix::rowSums(datExpr)==0))))
}
Expr = rownames(datExpr)

datExpr=Wdata@assays$RNA@counts[Expr,]

# This is my addition because afaik the don't log normalise???? Oddly, the average counts are stable w.r.t sequencing depth???
datExpr <- t(as.matrix(log2(datExpr + 1)))

enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = "signed", corFnc = "bicor", dataIsExpr = TRUE)

par(mfrow = c(1,2));
cex1 = 0.9;
#  Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# These are the scale-free topology indexes
indexes = (-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
# If we don't have aany index above 0.75, we stop the script
if ( !any( indexes > 0.75 ) ) {
  print("The scale-free topology index didn't reach 0.75 with any of the chosen powers, please consider changing the set of genes or cells provided")
  # quit(save = "no", 1, F)
}
if ( !any( indexes > 0.4 ) ) {
  print("The scale-free topology index didn't reach 0.4 with any of the chosen powers, the script will not continue, since it might fail")
  
  knitr::knit_exit()
}
# Take the smalles power that gives us an index over 0.9, or the highest index if we don't reach 0.9
if ( any( indexes > 0.9 ) ) {
    my.power = sft$fitIndices$Power[min(which(indexes > 0.9))]
} else { my.power = sft$fitIndices$Power[which.max(indexes)] }

print(paste0("Or power is ", my.power))

# Running WGCNA iteratively
my.Clnumber = 20
change = 0
genesets=list()
nonsig = 1
while(nonsig != 0) {
  
  my.adjacency =WGCNA::adjacency(datExpr,type = "signed", power = my.power, corFnc = "bicor")
  
  # Turn adjacency into topological overlap (high overlap if they share the same "neighborhood")
  TOM=WGCNA::TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = my.power, corType = "bicor")
  
  #Put the names in the tree
  colnames(TOM) <- gnames[colnames(datExpr),2]
  
  rownames(TOM) <- gnames[colnames(datExpr),2]
  
  #Make it a distance
  dissTOM = 1-TOM
  
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average")
  
  
  # Here I calculate the cutting height. Using the same formula and approach that the WGCNA package uses for the automatic function
  nMerge = length(geneTree$height) # The whole height of the tree
  refQuantile = 0.05 # What's the quantile that we want to exclude
  refMerge = round(nMerge * refQuantile) 
  refHeight = geneTree$height[refMerge]
  cutheight = signif(0.99 * (max(geneTree$height) - refHeight) + refHeight,4)
  
  # We construct THE TABLE that will help us make decisions
  # Min cluster sizes, from 7 to 30
  x=seq(7,30,1)
  # The height, up and down from the calculated height. We expect some "No module detected"
  y=seq(cutheight-0.0005,cutheight + 0.0005,0.0001)
  
  # The actual dataframe
  w=data.frame()
  # Populate, with i=min cluster size, j=cutting height, z=total number of clusters, z.1.=what's the first cluster? 0 is grey 1 is something else,
  # z.1.'=what's the size of the first cluster?
  for (i in x) {
    for (j in y) {
      sink("aux")
      z=table(dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = i, deepSplit = T, cutHeight = j, verbose = 0))
      sink(NULL)
      v=data.frame(i,j,dim(z),names(z[1]),unname(z[1]))
      w=rbind(w,v)
    }
  }
  
  # The height is then the one where we have the least number of genes in the first cluster
  my.height = w$j[which(w$unname.z.1..==min(w$unname.z.1..))]
  
  # Since different heights can give us the minimum grey size, we chose the computed height, if present, or the highest one.
  if (cutheight %in% my.height) {
    my.Clsize = w[which(w$j == cutheight),]
  } else { my.Clsize = w[which(w$j == max(my.height)),] }
  
  
  # This is to know, if we're looking for a minimum of cluster numbers
  # If we still have a lot of genes, we don't want to limit the number of clusters
  # if ( ((dim(datExpr)[2]) / length(Expr)) > 0.6 ) {change = 0}
  
  # If this is the first iteration after 0.6 of the genes are gone, we assign the number of clusters (and an extra for the grey in the case)
  if (change == 2) {
    my.Clnumber = length(table(dynamicColors)) + 1
  }
  # Count another iteration
  change = change + 1
  
  # If we don't have a gray cluster anymore, then we subset for those rows, and set a new height. ONLY if we get the same amount of clusters!
  if (any(w$names.z.1.. == 1)) { #any combination gives us no grey
    
    my.Clsize = w[which(w$names.z.1.. == 1),] # Take all combinations that gives us no grey
    
    if (any(my.Clsize$dim.z. >= (my.Clnumber -1) )) { # If there is any giving us the determined amount or more
      my.Clsize = my.Clsize[which(my.Clsize$dim.z. >= (my.Clnumber - 1) ),,drop=F] # Subset for those
    } else { my.Clsize = my.Clsize[which(my.Clsize$dim.z. == max(my.Clsize$dim.z.)),,drop=F] } # Or for the highest
    
    # Take the ones with the smallest number of clusters
    my.Clsize = my.Clsize[which(my.Clsize$dim.z. == min(my.Clsize$dim.z.)),,drop=F]
    # Take the one with the highest min cluster size
    my.Clsize = my.Clsize[which(my.Clsize$i == max(my.Clsize$i)),,drop=F]
    
    if (cutheight %in% my.Clsize$j) { # if original computed height is in,
      my.height = cutheight # take it
    } else { my.height = max(my.Clsize$j) } # Otherwise, the highest height
    
    my.Clsize = max(my.Clsize$i)
    
  }
  
    # Subset the table again, for those sizes that will gives the same number of clusters or more. IF NONE, use the highest number
  if (!any(w$names.z.1.. == 1)){
    if (any(my.Clsize$dim.z. >= my.Clnumber)) {
      my.Clsize = my.Clsize[which(my.Clsize$dim.z. >= my.Clnumber),,drop=F]
    } else {
      my.Clsize = my.Clsize[which(my.Clsize$dim.z. == max(my.Clsize$dim.z.)),,drop=F]}
    
   # Take the ones with the smallest number of clusters
    my.Clsize = my.Clsize[which(my.Clsize$dim.z. == min(my.Clsize$dim.z.)),,drop=F]
    # Take the one with the highest min cluster size
    my.Clsize = my.Clsize[which(my.Clsize$i == max(my.Clsize$i)),,drop=F]
    
    if (cutheight %in% my.Clsize$j) { # if original computed height is in,
      my.height = cutheight # take it
    } else { my.height = max(my.Clsize$j) } # Otherwise, the highest height
    
    my.Clsize = max(my.Clsize$i)
  }
  
  # If we still have more than 60% of the genes, we just use the min size of 15 regardles
  if ( change < 3 ) {
    my.Clsize = 15
    my.height = cutheight
  }
  
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = my.Clsize, deepSplit = T, cutHeight = my.height)
  
  #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, minClusterSize = minModuleSize)
  
  # Convert numeric lables into colors
  dynamicColors = WGCNA::labels2colors(dynamicMods)
  
  # Calculate eigengenes
  MEList = WGCNA::moduleEigengenes(as.matrix(datExpr), colors = dynamicColors)
  
  MEs = MEList$eigengenes
  
  # Calculate the module membership
  geneModuleMembership = as.data.frame(WGCNA::signedKME(datExpr, MEs))
  MMPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
  
  # We're gonna make a list, where we keep the genes that are significantly associated with each module
  x=c()
  xy=list()
  
  # We also need a vector with all the dynamic colors
  dcols = 1:length(levels(as.factor(dynamicColors)))
  
  # Getting rid of the grey module
  grey.genes = length(which(dynamicColors == "grey"))
  if (any(levels(as.factor(dynamicColors)) == "grey")) {
    dcols = dcols[-which(levels(as.factor(dynamicColors)) == "grey")]
  }
  
  # Run the loop to get the genes
  for (i in dcols) {
    modGenes = rownames(MMPvalue)[which(dynamicColors==levels(as.factor(dynamicColors))[i] & MMPvalue[,i]<0.01)]
    x=c(x,modGenes)
    xy[[i]]=modGenes
    #print(paste0(levels(as.factor(dynamicColors))[i]," ",length(modGenes),
    #" of ", length(which(dynamicColors==levels(as.factor(dynamicColors))[i]))))
    #print(gnames[modGenes,2])
  }
  
  # Make a new list, where we keep ALL the gens thar are left from the iteration, that will be used to make the new object. To keep track
  genesets[[length(genesets)+1]] = colnames(datExpr)
  
  # Give me a message saying how many genes are gone this time
  cat( paste0( grey.genes, " genes not assigned to any module.", '\n',
                 length(which(!(colnames(datExpr)%in%x))) - grey.genes, " genes excluded due to significance."))
  # Save this also, cause if it's 0 then we stop the whole thing
  nonsig = length(which(!(colnames(datExpr)%in%x)))
  
  # If it ain't 0, subset the dynamic colors and the expression data
  if (length(which(!(colnames(datExpr)%in%x))) != 0) {
    dynamicColors=dynamicColors[-which(!(colnames(datExpr)%in%x))]
    datExpr=datExpr[,-(which(!(colnames(datExpr)%in%x)))]
  } 
}

# We calculate here the expression of genes using the single-cell data
p.MEList = MEList
raw.datExpr = p.Wdata@assays$RNA@data[colnames(datExpr),]
raw.datExpr = t(as.matrix(raw.datExpr))
raw.MEList = WGCNA::moduleEigengenes(raw.datExpr, colors = dynamicColors)
p.MEList = raw.MEList

## Modules of co-expression

# We can see what are the expression levels of our co-expression modules in the single cells. Using the dimensionality reduction provided. Default is a tSNE  

options(repr.plot.width = 20, repr.plot.height = 40)
xx=list()
yy=c(levels(as.factor(dynamicColors)))
for (i in 1:length(yy)) {
  toplot = data.frame(Seurat::Embeddings(p.Wdata[[my.dr]]))
  xx[[i]] = ggplot2::ggplot(toplot[order(p.MEList$averageExpr[,i]),],
                   ggplot2::aes_string(x=colnames(toplot)[1], y=colnames(toplot)[2])) +
    ggplot2::geom_point(ggplot2::aes_string(color=p.MEList$averageExpr[order(p.MEList$averageExpr[,i]),i]), size=0.5) +
    ggplot2::scale_size(range = c(1, 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position="none") +
    ggplot2::ggtitle(colnames(p.MEList$eigengenes)[i]) +
    ggplot2::scale_colour_gradientn(colours = c("gray90", "gray90", yy[i], yy[i])) +
    ggplot2::labs(colour=levels(as.factor(dynamicColors))[i])
}
p <- gridExtra::grid.arrange(grobs=xx, ncol=4)
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_module_umaps.pdf"),
       plot = p,
       width = 20,
       height = 40,
       units = "in",
       device='pdf'
      )

WGCNA_data = list()
WGCNA_data[["datExpr"]] = datExpr
WGCNA_data[["dynamicMods"]] = dynamicMods
WGCNA_data[["dynamicColors"]] = dynamicColors
WGCNA_data[["MEList"]] = MEList
WGCNA_data[["MEs"]] = MEs
WGCNA_data[["modGenes"]] = xy
WGCNA_data[["MEList.sc"]] = p.MEList
WGCNA_data[["genesets"]] = genesets
WGCNA_data[["TOM"]]= TOM
WGCNA_data[["adjacency"]]= my.adjacency
saveRDS(WGCNA_data, file = paste0(opts$outbase, opts$day, "_", opts$exp, "_WGCNA_data_", str_replace_all(as.character(Sys.Date()), '-', ''), ".rds"))

#Recalculate MEs for gastruloids
logXdat <- logXdat[,colnames(datExpr)]
MEList_gastr = moduleEigengenes(logXdat, dynamicColors)
MEs_gastr = orderMEs(MEList_gastr$eigengenes)
gastr_cor <- cor(t(MEs_gastr))
gastr_dist <- as.dist(1 - gastr_cor)
gastr_tree <- hclust(gastr_dist, method="complete")
my_order <- colnames(gastr_cor)[gastr_tree$order]
gastr_cor <- cor(t(MEs_gastr))[my_order,my_order]

ME_cor <- cor(MEs)
ME_dist <- as.dist(1 - ME_cor)
ME_tree <- hclust(ME_dist, method="complete")
ME_order <- colnames(ME_cor)[ME_tree$order]
ME_cor <- cor(MEs)[ME_order,ME_order]

pb_sim_list <- readRDS(opts$sims)

var_emperical <- matrix(0L, nrow=length(pb_sim_list), ncol = ncol(MEList_gastr$varExplained))
colnames(var_emperical) <- colnames(MEList_gastr$eigengenes)
for (i in 1:length(pb_sim_list)) {
    GeneXData_sim <- pb_sim_list[[i]] %>% as.matrix %>% t
    GeneXData_sim <- apply(GeneXData_sim,2, function(x) { (x/sum(x))*1000000 })
    # format properly for WGCNA
    GeneXData_sim <- as.data.frame(t(GeneXData_sim))
    # LOG2 Transformation
    logXdat_sim <- log2(GeneXData_sim + 1)
    logXdat_sim <- logXdat_sim[,colnames(datExpr)]
    MEList_sim <- moduleEigengenes(logXdat_sim, dynamicColors)
    var_sim <- MEList_sim$varExplained
    colnames(var_sim) <- colnames(MEList_sim$eigengenes)
    var_emperical[i,] <- as.matrix(var_sim)[1,colnames(var_emperical)]
}

percentiles <- list()
xx <- list()
for (i in 1:ncol(var_emperical)) {
    percentile <- ecdf(var_emperical[,i])
    percentiles[[colnames(var_emperical)[i]]] <- percentile(MEList_gastr$varExplained[1,i])
    xx[[i]] <- ggplot(as.data.table(var_emperical), aes_string(x=colnames(var_emperical)[i])) +
      geom_density() +
      ggtitle(paste0(colnames(var_emperical)[i], ', percentile=', percentiles[[colnames(var_emperical)[i]]])) +
      geom_vline(xintercept=MEList_gastr$varExplained[1,i]) +
      theme_bw()
}
p <- gridExtra::grid.arrange(grobs=xx[order(unlist(percentiles))], ncol=4)
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_variance_percentiles.pdf"),
       plot = p,
       width = 13,
       height = 25,
       units = "in",
       device='pdf'
      )

#strong.correlated.mods.names <- colnames(MEList_gastr$eigengenes)[as.numeric(MEList_gastr$varExplained[1,])/as.numeric(MEList$varExplained[1,]) > 1.5]
strong.correlated.mods.names <- names(percentiles[(1-unlist(percentiles)) < (0.05/ncol(var_emperical))])
strong.correlated.mods <- MEs_gastr[my_order,colnames(MEs_gastr) %in% strong.correlated.mods.names]
textMatrix= paste(signif(strong.correlated.mods, 2))

ME_cor <- cor(MEs[,strong.correlated.mods.names])
ME_dist <- as.dist(1 - ME_cor)
ME_tree <- hclust(ME_dist, method="complete")
pdf(file = paste0(opts$outbase, opts$day, "_", opts$exp, "_sig_ME_dendrogram.pdf"),
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches
plot(ME_tree, hang = -1)
dev.off()

ME_order <- colnames(ME_cor)[ME_tree$order]
ME_cor <- cor(MEs)[ME_order,ME_order]
pdf(file = paste0(opts$outbase, opts$day, "_", opts$exp, "_ME_corrplot.pdf"),
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
corrplot(ME_cor, col= blueWhiteRed(50), col.lim= c(-1,1))
dev.off()

gastr_cor <- cor(t(MEs_gastr[,strong.correlated.mods.names]))
gastr_dist <- as.dist(1 - gastr_cor)
gastr_tree <- hclust(gastr_dist, method="complete")
pdf(file = paste0(opts$outbase, opts$day, "_", opts$exp, "_sig_gastruloid_dendrogram.pdf"),
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches
plot(gastr_tree, hang = -1)
dev.off()

my_order2 <- colnames(cor(t(MEs_gastr[,strong.correlated.mods.names])))[gastr_tree$order]
gastr_cor <- cor(t(MEs_gastr[,strong.correlated.mods.names]))[my_order2,my_order2]
pdf(file = paste0(opts$outbase, opts$day, "_", opts$exp, "_gastruloid_corrplot.pdf"),
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
corrplot(gastr_cor, col= blueWhiteRed(50), col.lim= c(-1,1))
dev.off()

p <- ggplot(data = melt(as.matrix(strong.correlated.mods[my_order2,ME_order])), aes(x=Var2, y = Var1, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low=rgb(0.05, 0.55, 1.00), mid=rgb(1,1,1), high=rgb(1.0, 0.2, 0), midpoint=0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_sigME_gastruloid_hm.pdf"),
       plot = p,
       width = 10,
       height = 10,
       units = "in",
       device='pdf'
      )


########################################################
### GET EVERYTHING NECESSARY FOR GENE SET ENRICHMENT ###
########################################################

genes <- colnames(logXdat)
gene_mod <- cbind(genes, module=dynamicColors)

moduleColors <- dynamicColors

for (k in (1:length(moduleColors))) {
  moduleColors[k] <- paste("ME",moduleColors[k],sep="")
}
gene_mod <- cbind(genes, module=moduleColors)
gene_mod <- as.data.frame(gene_mod)

modules <- unique(moduleColors)
modulek <- modules[1]
genesprmod <- list(gene_mod[gene_mod$module == modulek,"genes"])
for (k in (2:length(modules))) {
  modulek <- modules[k]
  genesprmod <- append(genesprmod,list(gene_mod[gene_mod$module == modulek,"genes"]))
}
names(genesprmod)<-modules
write.table(genesprmod[modulek],file = paste0(opts$outbase, opts$day, "_", opts$exp, "_ALLgenesprmodule.tab"), row.names = F, append = F,sep="\t",quote = F)
for (i in names(genesprmod)[!names(genesprmod) %in% modulek]) {
  write.table(genesprmod[i],file = paste0(opts$outbase, opts$day, "_", opts$exp, "_ALLgenesprmodule.tab"), row.names = F, append = T,sep="\t",quote = F)
}

to.plot <- data.frame(lengths(genesprmod))
colnames(to.plot) <- c('ngenes')
to.plot$module <- rownames(to.plot)
to.plot$f <- substr(to.plot$module, 3, nchar(to.plot$module))
to.plot <- as.data.table(to.plot)

p <- ggplot(to.plot, aes(x=module, y=ngenes, fill=f)) +
    geom_bar(stat="identity", colour='black') +
    scale_fill_identity() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_module_ngenes.pdf"),
       plot = p,
       width = 6,
       height = 4,
       units = "in",
       device='pdf'
      )

p <- ggplot(to.plot, aes(x=module, y=ngenes, fill=(module %in% strong.correlated.mods.names))) +
    geom_bar(stat="identity", colour='black') +
    #scale_fill_identity() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_module_ngenes_sig.pdf"),
       plot = p,
       width = 6,
       height = 4,
       units = "in",
       device='pdf'
      )

to.plot2 <- to.plot[module %in% strong.correlated.mods.names]
p <- ggplot(to.plot2, aes(x=module, y=ngenes, fill=f)) +
    geom_bar(stat="identity", colour='black') +
    #scale_fill_manual(values = to.plot2$f) +
    scale_fill_identity() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_sig_module_ngenes.pdf"),
       plot = p,
       width = 6,
       height = 4,
       units = "in",
       device='pdf'
      )


###########################
### GENE SET ENRICHMENT ###
###########################

set.seed(Sys.time())
seed<-sample(1:20000, 1, replace=F)
print(seed) #1446 #19607
set.seed(seed)

wgcnas <- get_wgcna(path_wgcna = paste0(opts$outbase, opts$day, "_", opts$exp, "_ALLgenesprmodule.tab"), 
                    organism = "mmusculus"
                   )
correlations_table <- strong.correlated.mods
wgcna.split <- wgcnas#[correlations_table$modules]
wgcna_df <- stack(wgcna.split)
colnames(wgcna_df) <- c("Genes", "Module Name")

#randomised modules for significance testing 
mods <- names(wgcnas)
gene_list <- as.character(unlist(wgcnas, recursive=FALSE))
len_vec <- as.numeric(lapply(wgcnas, function(x) { length(x) } ))
prob_vec <- len_vec/sum(len_vec)
ss <- sample(1:length(prob_vec),
             size=length(gene_list),
             replace=TRUE,
             prob=prob_vec
            )
shuffled_modules <- setNames(split(gene_list,ss), mods)

wgcna.split <- wgcnas#[uniq_mod]

#Enrichment of TFs using the following databases

#TRRUST_Transcription_Factors_2019

# Remove TF_Perturbations_Followed_by_Expression because of how unspecific it is (same TF is detected in many modules)
# Remove ChEA_2016 because of how unspecific it is (many TFs are enriched in a given module)
#tf_db<- c( "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "ENCODE_TF_ChIP-seq_2015")
tf_db<- c(#"ChEA_2016", # grep 'Mouse'; ChEA contains results from transcription factor ChIP-seq studies, extracted from supporting materials of publications. Each entry has the transcription factor, PMID, cell type and organism. A peak at the promoter of gene was detected in each study.
          "TRANSFAC_and_JASPAR_PWMs", # grep '(mouse)'; Using the PWM from TRANSFAC and JASPAR the binding motifs were detected at the promoter of gene.
          "ENCODE_TF_ChIP-seq_2015", # grep 'mm9'; Processed ChIP-seq data from the ENCODE project to associate detected peaks near genes.
          "TRRUST_Transcription_Factors_2019", # grep 'mouse'; Current version of TRRUST contains 8,444 and 6,552 TF-target regulatory relationships of 800 human TFs and 828 mouse TFs, respectively.
          #"ARCHS4_TFs_Coexp", # THIS IS ONLY HUMAN; Top 300 genes from ARCHS4 that are co-expressed with transcription factors.
          #"Enrichr_Submissions_TF-Gene_Coocurrence", # DOESN'T DIFFERENTIATE MOUSE AND HUMAN; All gene co-occurrence adjacency matrix was created from ~300,000 gene sets submitted to Enrichr by non-frequent users. The genes in each set are those most commonly co-occur with the transcription factors, which are the labels of each set
          
          # NOTE: FOR NOW I DON'T CARE ABOUT POS/NEG REGULATION
          #"TF_Perturbations_Followed_by_Expression", #check 'MOUSE' and neg_perturb --> DOWN or pos_perturb --> UP; based on GEO: [TF] [Perturbation Type, see neg and pos above] [Organism] [GEO accession] CREEDSID GENE [???] [UP/DOWN]
          "TF-LOF_Expression_from_GEO" #check 'mouse' and 'lof'-->'down', 'gof'-->'up'
         )

df.names<-data.frame()
for (module in names(wgcna.split)) {
  print(module)
  enriched_prizes<- enrichr(wgcna.split[module][[1]], tf_db)
    
  enriched_prizes[["ChEA_2016"]] <- enriched_prizes[["ChEA_2016"]][grepl('Mouse', enriched_prizes[["ChEA_2016"]][["Term"]]),]
  enriched_prizes[["TRANSFAC_and_JASPAR_PWMs"]] <- enriched_prizes[["TRANSFAC_and_JASPAR_PWMs"]][grepl('(mouse)', enriched_prizes[["TRANSFAC_and_JASPAR_PWMs"]][["Term"]]),]
  enriched_prizes[["ENCODE_TF_ChIP-seq_2015"]] <- enriched_prizes[["ENCODE_TF_ChIP-seq_2015"]][grepl('mm9', enriched_prizes[["ENCODE_TF_ChIP-seq_2015"]][["Term"]]),]
  enriched_prizes[["TRRUST_Transcription_Factors_2019"]] <- enriched_prizes[["TRRUST_Transcription_Factors_2019"]][grepl('mouse', enriched_prizes[["TRRUST_Transcription_Factors_2019"]][["Term"]]),]
  #enriched_prizes[["TF_Perturbations_Followed_by_Expression"]] <- enriched_prizes[["TF_Perturbations_Followed_by_Expression"]][grepl('MOUSE', enriched_prizes[["TF_Perturbations_Followed_by_Expression"]][["Term"]]),]
  #perturb_expr_terms <- enriched_prizes[["TRRUST_Transcription_Factors_2019"]][["Term"]]
  #idx <- grepl(paste(neg_perturb, collapse='|'), perturb_expr_terms) & grepl('DOWN', perturb_expr_terms)
  enriched_prizes[["TF-LOF_Expression_from_GEO"]] <- enriched_prizes[["TF-LOF_Expression_from_GEO"]][grepl('mouse', enriched_prizes[["TF-LOF_Expression_from_GEO"]][["Term"]]),]
  #enriched_prizes[["TF-LOF_Expression_from_GEO"]] <- enriched_prizes[["TF-LOF_Expression_from_GEO"]][grepl("lof.*down|gof.*up",enriched_prizes[["TF-LOF_Expression_from_GEO"]][["Term"]]),]
    
  enrichdt<-rbindlist(lapply(1:length(enriched_prizes), function(x){ setDT(enriched_prizes[[x]])[, id:=names(enriched_prizes)[x]]}), use.names=TRUE, fill=TRUE)
  enrichdt<-enrichdt[enrichdt$Adjusted.P.value < 0.05,]
  enrichdt<-enrichdt[order(enrichdt$Adjusted.P.value),]
  if (dim(enrichdt)[1] != 0) {
    row.to.add <- data.frame(module, 
                             enrichdt$Term,
                             enrichdt$id, 
                             enrichdt$Overlap, 
                             enrichdt$P.value,
                             enrichdt$Adjusted.P.value, 
                             enrichdt$Odds.Ratio, 
                             enrichdt$Genes)
    df.names<-rbind(df.names, row.to.add)
  }
}


df.names$enrichdt.Term <- unlist(lapply(strsplit(df.names$enrichdt.Term, " "), head, n = 1L))
sup.table <- df.names
sup.table$minuslog10pvalue <- -log10(sup.table$enrichdt.Adjusted.P.value)
sup.table <- as.data.table(sup.table)

fwrite(sup.table, file=paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_TFs.txt"))

dots_tf<-ggplot(data = sup.table[enrichdt.Adjusted.P.value<0.05,], aes(x=module, y = reorder(enrichdt.Term, -minuslog10pvalue), color = enrichdt.Odds.Ratio, size = minuslog10pvalue)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_tf_allmods.pdf"),
       plot = dots_tf,
       width = 15,
       height = 15,
       units = "in",
       device='pdf'
      )

sup.table.tmp <- sup.table[(enrichdt.Adjusted.P.value<0.05) & (module %in% colnames(strong.correlated.mods)),]
sup.table.tmp$module <- factor(sup.table.tmp$module, levels=ME_order)
dots_tf<-ggplot(data = sup.table.tmp, aes(x=module, y = reorder(enrichdt.Term, -minuslog10pvalue), color = enrichdt.Odds.Ratio, size = minuslog10pvalue)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_tf_sigmods.pdf"),
       plot = dots_tf,
       width = 15,
       height = 15,
       units = "in",
       device='pdf'
      )

dots_tf<-ggplot(data = sup.table[(enrichdt.Adjusted.P.value<0.05) & !(module %in% colnames(strong.correlated.mods)),], aes(x=module, y = reorder(enrichdt.Term, -minuslog10pvalue), color = enrichdt.Odds.Ratio, size = minuslog10pvalue)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_tf_nonsigmods.pdf"),
       plot = dots_tf,
       width = 15,
       height = 15,
       units = "in",
       device='pdf'
      )

colnames(sup.table) <- c("Module Name", "TF", "Source Database", "Overlap", "P.value", "Adjusted P.value", "Odds Ratio", "Gene List", "-log10(P)")
#write.csv(x = sup.table, file = "~/cell_shapes/manuscript/figures/cell shape figures/supplementary/STable2.csv")

split.names <- df.names %>%
  group_by(module)
split.names <- group_split(split.names)
names(split.names) <- sort(unique(df.names$module))

reactome.df<-data.frame()
for (module in names(split.names)) {
  if (  length(unique(split.names[module][[1]]$enrichdt.Term)) < 0.1) {
    next
  }
  #enriched_prizes <- enrichr(split.names[module][[1]]$enrichdt.Term, "KEGG_2019_Mouse")
  #enriched_prizes <- enrichr(split.names[module][[1]]$enrichdt.Term, "Reactome_2016")
  enriched_prizes <- enrichr(split.names[module][[1]]$enrichdt.Term, "WikiPathways_2019_Mouse")
  #enriched_prizes <- enrichr(names(table(split.names[module][[1]]$enrichdt.Term))[which(sort(table(split.names[module][[1]]$enrichdt.Term)) > 0)]
  #  , "Reactome_2016")
  enrichdt <- rbindlist(lapply(1:length(enriched_prizes), function(x) { setDT(enriched_prizes[[x]])[, id:=names(enriched_prizes)[x]]}), use.names=TRUE, fill=TRUE)
  enrichdt <- enrichdt[enrichdt$Adjusted.P.value < 0.05,]
  enrichdt <- enrichdt[order(enrichdt$Adjusted.P.value),]
  if (dim(enrichdt)[1] != 0) {
    row.to.add <- data.frame(module, 
                             enrichdt$Term,
                             enrichdt$id, 
                             enrichdt$Overlap, 
                             enrichdt$P.value,
                             enrichdt$Adjusted.P.value, 
                             enrichdt$Odds.Ratio, 
                             enrichdt$Genes)
    reactome.df <- rbind(reactome.df, row.to.add)
  }
}
reactome.df$minuslog10pvalue <- -log2(reactome.df$enrichdt.Adjusted.P.value)

reactome.df$REACTOME_ID <- unlist(lapply(strsplit(reactome.df$enrichdt.Term, " "), tail, n = 1L))
dot.plot <- reactome.df
if (any(is.infinite(dot.plot$enrichdt.Odds.Ratio))) {
    dot.plot[is.infinite(dot.plot$enrichdt.Odds.Ratio),]$enrichdt.Odds.Ratio <- 10
}
dot.plot <- as.data.table(dot.plot)

dots <- ggplot(data = dot.plot[enrichdt.Adjusted.P.value<0.05,], aes(x=module, y = reorder(enrichdt.Term, -minuslog10pvalue), color = log2(enrichdt.Odds.Ratio), size = minuslog10pvalue)) + 
  geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_pathway_allmods.pdf"),
       plot = dots,
       width = 15,
       height = 15,
       units = "in",
       device='pdf'
      )

dot.plot.tmp <- dot.plot[(enrichdt.Adjusted.P.value<0.05) & (module %in% colnames(strong.correlated.mods)),]
dot.plot.tmp$module <- factor(dot.plot.tmp$module, levels=ME_order)
dots <- ggplot(data = dot.plot.tmp, aes(x=module, y = reorder(enrichdt.Term, -minuslog10pvalue), color = log2(enrichdt.Odds.Ratio), size = minuslog10pvalue)) + 
  geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_pathway_sigmods.pdf"),
       plot = dots,
       width = 15,
       height = 15,
       units = "in",
       device='pdf'
      )

dots <- ggplot(data = dot.plot[(enrichdt.Adjusted.P.value<0.05) & !(module %in% colnames(strong.correlated.mods)),], aes(x=module, y = reorder(enrichdt.Term, -minuslog10pvalue), color = log2(enrichdt.Odds.Ratio), size = minuslog10pvalue)) + 
  geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_pathway_nonsigmods.pdf"),
       plot = dots,
       width = 15,
       height = 15,
       units = "in",
       device='pdf'
      )

fwrite(dot.plot, file=paste0(opts$outbase, opts$day, "_", opts$exp, "_dots_pathway.txt"))




























































