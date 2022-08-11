#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(SeuratDisk))

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to the settings file", metavar="character"),
    make_option(c("-f", "--functions"), type="character", default=NULL, 
                help="path to the functions file", metavar="character"),
    make_option(c("-I", "--input.seurat"), type="character", default=NULL, 
                help="path to input seurat file", metavar="character"),
    make_option(c("-i", "--input.metadata"), type="character", default=NULL, 
                help="path to input sample metadata file", metavar="character"),
    make_option(c("-G", "--input.gene.metadata"), type="character", default=NULL, 
                help="path to input gene metadata file", metavar="character"),
    make_option(c("-p", "--plot.outbase"), type="character", default=NULL, 
                help="path to directory in which to save plots", metavar="character"),
    make_option(c("-O", "--output.seurat"), type="character", default=NULL, 
                help="path to output seurat file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-g", "--output.gene.metadata"), type="character", default=NULL, 
                help="path to output gene metadata file", metavar="character")#,
    #make_option(c("-h", "--output.h5ad"), type="character", default=NULL, 
    #            help="path to output h5ad file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

source(opts$settings)
source(opts$functions)

message("Loading data")

# Load metadata
meta_query <- fread(opts$input.metadata)
meta_query <- meta_query[(pass_QC==TRUE) & !(hybrid_call2) & ((experiment != "exp1_d5") | (MULTI_class != "Doublet")) & (!(experiment %in% c("exp2A_d3_d3.5", "exp2C_d3_d3.5")))]

# Load SingleCellExperiment
sce_query <- Seurat::as.SingleCellExperiment( readRDS(opts$input.seurat)[,meta_query$cell] )

stopifnot(all(meta_query$cell == colnames(sce_query)))

message("Filtering genes...")

mart = useMart('ensembl', dataset='mmusculus_gene_ensembl', host='www.ensembl.org')
go_cellcycle = getBM(
    attributes = c('ensembl_gene_id','external_gene_name'), 
    filters    = 'go', 
    values     = c('GO:0007049', # cell cycle
                   'GO:0007059', # chromosome segregation
                   'GO:0006266', # DNA ligation
                   'GO:0000819', # sister chromatid segregation
                   'GO:0032508', # DNA duplex unwinding
                   'GO:1905462', # regulation of DNA duplex unwinding
                   'GO:0030261', # chromosome condensation
                   'GO:0010389', # regulation of G2/M transition of mitotic cell cycle
                   'GO:0072686', # mitotic spindle'
                   'GO:0006335', # DNA replication-dependent nucleosome assembly
                   'GO:0007052', # mitotic spindle organization
                   'GO:0000281', # mitotic cytokinesis
                   'GO:0000132', # establishment of mitotic spindle orientation
                   'GO:0000786', # nucleosome
                   'GO:0000278', # mitotic cell cycle
                   'GO:0007094', # mitotic spindle assembly checkpoint signaling
                   'GO:0045930', # negative regulation of mitotic cell cycle
                   'GO:0000086', # G2/M transition of mitotic cell cycle
                   'GO:0045740', # positive regulation of DNA replication
                   'GO:0051782', # negative regulation of cell division
                   'GO:0006270', # DNA replication initiation
                   'GO:0006261', # DNA-dependent DNA replication
                   'GO:0007051', # spindle organization
                   'GO:0007098', # centrosome cycle
                   'GO:0006260', # DNA replication
                   'GO:0071897', # DNA biosynthetic process
                   'GO:0006275', # regulation of DNA replication
                   'GO:0031297', # replication fork processing
                   'GO:0022402'#, # cell cycle process
                  ),
    mart       = mart)
go_cellcycle <- do.call("rbind",
                        list(go_cellcycle,
                             c('ENSMUSG00000020330', 'Hmmr'), 
                             c('ENSMUSG00000035293', 'G2e3'), # literally called "G2/M phase-specific E3 ubiquitin-protein ligase"
                             c('ENSMUSG00000044783', 'Hjurp'), # the protein is associated with all the above go terms (nearly), this ENS ID has no GO terms.
                             c('ENSMUSG00000022033', 'Pbk'), # no GO terms, but "Seems to be active only in mitosis"
                             c('ENSMUSG00000027203', 'Dut'), # produces dUMP, the immediate precursor of thymidine nucleotides
                             c('ENSMUSG00000020737', 'Jpt1') # Plays a role in the regulation of cell cycle, no GO terms
                            )
                       )
sce_query  <- sce_query[!(rownames(sce_query) %in% go_cellcycle$ensembl_gene_id),]


message("MNN mapping...")

meta_query$sample <- meta_query$batch
meta_query$batch <- meta_query$experiment
meta_query$timepoint.demultiplexed <- meta_query$timepoint
meta_query$timepoint <- meta_query$timepoint.base
correct_sce <- mapNoAtlas(sce_query, meta_query, c("d5", 'd4.5_d5', "d4.5", 'd4_d4.5', "d4", "d3.5_d4", "d3.5", "d3_d3.5", "d3"), k = io$k, npcs = io$npcs, return.list = FALSE)
meta_query$batch <- meta_query$sample

message("Done\n")

message("QCing...")

UMAP <- umap(correct_sce$corrected, random_state=2402)
tmp <- data.frame(UMAP$layout)
tmp$cell <- rownames(tmp)
plot_df_joined <- merge(tmp, meta_query, by=c("cell")) %>% as.data.table

p <- ggplot(data=plot_df_joined, mapping = aes(x=-X1, y=-X2, colour=timepoint)) +
  geom_point(size=0.1, alpha=0.5) +
  # ggrastr::geom_point_rast(size=io$dot_size, alpha=io$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_timepoint_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=Phase)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_phase_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=percent.mt)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_viridis(trans='log2') +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_percentmito_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=percent.ribo)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_viridis() +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_precentribo_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=nFeature_RNA)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_viridis() +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_nFeatureRNA_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=nCount_RNA)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_viridis(trans='log2') +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_nCountRNA_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=hybrid_call)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_hybridcall_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=hybrid_map)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_viridis() +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_hybridmap_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=batch)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_batch_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

p <- ggplot(data=plot_df_joined, mapping = aes(x=X1, y=X2, colour=experiment)) +
  geom_point(size=0.5, alpha=1) +
  # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic()
ggsave(paste0(opts$plot.outbase, "first_experiment_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

# Load data
srat <- readRDS(opts$input.seurat)
srat <- subset(srat, cells = rownames(correct_sce$corrected))
umap_embeddings <- as.matrix(plot_df_joined[,c("X1", "X2")])
rownames(umap_embeddings) <- plot_df_joined$cell
srat@meta.data <- as.data.frame(plot_df_joined[cell %in% Cells(srat)])
rownames(srat@meta.data) <- srat@meta.data$cell
srat@meta.data <- srat@meta.data[Cells(srat),]

srat@reductions[["pca"]] <- CreateDimReducObject(
  embeddings = correct_sce$corrected[Cells(srat),],
  loadings = metadata(correct_sce$pca)$rotation,
  key = "PC",
  assay = "RNA"
)
srat@reductions[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embeddings[Cells(srat),],
  key = "UMAP",
  assay = "RNA"
)

DefaultAssay(srat) <- "RNA"

VariableFeatures(srat) <- rownames(srat@reductions$pca@feature.loadings)

genemeta <- fread(opts$input.gene.metadata)
genemeta <- genemeta[ens_id %in% rownames(srat)]
genemeta[,cc_gene := (ens_id %in% go_cellcycle$ensembl_gene_id)]
genemeta[,gastr_only_hvg := (ens_id %in% rownames(srat@reductions$pca@feature.loadings))]
rownames(genemeta) <- genemeta$ens_id
srat[["RNA"]] <- AddMetaData(
    object = srat[["RNA"]],
    metadata = genemeta
)

sce <- as.SingleCellExperiment(srat)
sce <- multiBatchNorm(sce, batch=colData(sce)$batch)
srat@assays$RNA@data <- logcounts(sce)

srat <- FindNeighbors(srat, dims = 1:50, k.param = io$k)
srat <- FindClusters(srat, resolution = 1)
p <- DimPlot(srat, reduction = "umap", label = TRUE)
ggsave(paste0(opts$plot.outbase, "first_cluster_UMAP.pdf"), plot = p, device="pdf", width = 6, height = 5, units = "in")

message("Done\n")

message("Integrating final dataset...")
meta_query$sample <- meta_query$batch
meta_query$batch <- meta_query$experiment
#meta_query <- meta_query[cell %in% srat@meta.data[(!(srat@meta.data$`RNA_snn_res.1` %in% c("20"))) & (srat@meta.data$hybrid_map < 1/3),]$cell]
meta_query <- meta_query[cell %in% srat@meta.data[(!(srat@meta.data$`RNA_snn_res.1` %in% c("18"))),]$cell]
#sce_query <- sce_query[,srat@meta.data[(!(srat@meta.data$`RNA_snn_res.1` %in% c("20"))) & (srat@meta.data$hybrid_map < 1/3),]$cell]
sce_query <- sce_query[,srat@meta.data[(!(srat@meta.data$`RNA_snn_res.1` %in% c("18"))),]$cell]
correct_sce <- mapNoAtlas(sce_query, meta_query, c("d5", 'd4.5_d5', "d4.5", 'd4_d4.5', "d4", "d3.5_d4", "d3.5", "d3_d3.5", "d3"), k = io$k, npcs = io$npcs, return.list = FALSE)
meta_query$timepoint <- meta_query$timepoint.demultiplexed
meta_query$batch <- meta_query$sample

UMAP <- umap(correct_sce$corrected, random_state=2402)
tmp <- data.frame(UMAP$layout)
tmp$cell <- rownames(tmp)
plot_df_joined <- merge(tmp, meta_query, by=c("cell")) %>% as.data.table

# Load data
srat <- readRDS(opts$input.seurat)
srat <- subset(srat, cells = rownames(correct_sce$corrected))
umap_embeddings <- as.matrix(plot_df_joined[,c("X1", "X2")])
rownames(umap_embeddings) <- plot_df_joined$cell
srat@meta.data <- as.data.frame(plot_df_joined[cell %in% Cells(srat)])
rownames(srat@meta.data) <- srat@meta.data$cell
srat@meta.data <- srat@meta.data[Cells(srat),]

srat@reductions[["pca"]] <- CreateDimReducObject(
  embeddings = correct_sce$corrected[Cells(srat),],
  loadings = metadata(correct_sce$pca)$rotation,
  key = "PC",
  assay = "RNA"
)
srat@reductions[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embeddings[Cells(srat),],
  key = "UMAP",
  assay = "RNA"
)

DefaultAssay(srat) <- "RNA"

VariableFeatures(srat) <- rownames(srat@reductions$pca@feature.loadings)

genemeta <- fread(opts$input.gene.metadata)
genemeta <- genemeta[ens_id %in% rownames(srat)]
genemeta[,cc_gene := (ens_id %in% go_cellcycle$ensembl_gene_id)]
genemeta[,gastr_only_hvg := (ens_id %in% rownames(srat@reductions$pca@feature.loadings))]
rownames(genemeta) <- genemeta$ens_id
srat[["RNA"]] <- AddMetaData(
    object = srat[["RNA"]],
    metadata = genemeta
)

sce <- as.SingleCellExperiment(srat)
sce <- multiBatchNorm(sce, batch=colData(sce)$batch)
srat@assays$RNA@data <- logcounts(sce)

srat@meta.data <- srat@meta.data[Cells(srat),]

saveRDS(srat, opts$output.seurat)
fwrite(srat@meta.data, opts$output.metadata, quote=F, na="NA", sep="\t")
fwrite(genemeta, opts$output.gene.metadata, quote=F, na="NA", sep="\t")

#SaveH5Seurat(g_srat,
#             filename = str_replace_all(opts$output.h5ad, 'h5ad', 'h5Seurat')
#            )
#Convert(str_replace_all(opts$output.h5ad, 'h5ad', 'h5Seurat'), dest = "h5ad")



