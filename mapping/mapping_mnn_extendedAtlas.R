#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MouseGastrulationData))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Seurat))

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-A", "--atlas.RDS"), type="character", default=NULL, 
                help="path to the mouse gastrulation 10x atlas SCE RDS", metavar="character"),
    make_option(c("-a", "--atlas.metadata"), type="character", default=NULL, 
                help="path to the mouse gastrulation 10x atlas metdata txt", metavar="character"),
    make_option(c("-Q", "--query.seurat"), type="character", default=NULL, 
                help="path to the query to be mapped (a Seurat file)", metavar="character"),
    make_option(c("-q", "--query.metadata"), type="character", default=NULL, 
                help="path to the sample metadata file for the query", metavar="character"),
    make_option(c("-G", "--input.gene.metadata"), type="character", default=NULL, 
                help="path to input gene metadata file", metavar="character"),
    make_option(c("-b", "--query_batches"), type="character", default=NULL, 
                help="name of batch in query", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-f", "--functions"), type="character", default=NULL, 
                help="path to file with functions used", metavar="character"),
    make_option(c("-O", "--output.RDS"), type="character", default=NULL, 
                help="path to output RDS file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-g", "--output.gene.metadata"), type="character", default=NULL, 
                help="path to output gene metadata file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);


source(opts$settings)
source(opts$functions)

if (is.null(opts$query_batches)) opts$query_batches = "all"
opts$query_batches <- str_remove(opts$query_batches, 'MULTI/')


################
## Load atlas ##
################

message("Loading atlas...")

meta_atlas <- fread(opts$atlas.metadata, sep="\t", na="NA", quote=F) %>% as.data.table
#meta_atlas <- meta_atlas[meta_atlas$doublet==FALSE,]
#meta_atlas <- meta_atlas[meta_atlas$stripped==FALSE,]

# Load SingleCellExperiment
sce_atlas  <- readRDS(opts$atlas.RDS)[,meta_atlas$cell]
sce_atlas <- multiBatchNorm(sce_atlas, batch=meta_atlas$sample, preserve.single=TRUE)


################
## Load query ##
################

message("Loading query...")

# Load metadata
meta_query <- fread(opts$query.metadata) %>% as.data.table
meta_query <- meta_query[(pass_QC==TRUE) & !(hybrid_call) & ((experiment != "exp1_d5") | (MULTI_class != "Doublet")) & (!(experiment %in% c("exp2A_d3_d3.5", "exp2C_d3_d3.5")))]
if (opts$query_batches != 'all') {
    meta_query <- meta_query[(experiment %in% opts$query_batches)]
}

# Load SingleCellExperiment
sce_query <- Seurat::as.SingleCellExperiment( readRDS(opts$query.seurat)[,meta_query$cell] )


#####################
## Intersect genes ## 
#####################

message("Filtering genes...")

genemeta <- fread(opts$input.gene.metadata)

tmp1 <- rownames(sce_atlas)[rownames(sce_atlas) %in% genemeta[!(cc_gene)]$ens_id]
tmp2 <- rownames(sce_query)[rownames(sce_query) %in% genemeta[!(cc_gene)]$ens_id]
genes <- intersect(tmp1, tmp2)
sce_query  <- sce_query[genes,]
sce_atlas <- sce_atlas[genes,]

#########
## Map ##
#########

meta_query$batch <- meta_query$experiment
meta_query$timepoint <- meta_query$timepoint.base
tp_order <- c("d5", 'd4.5_d5', "d4.5", 'd4_d4.5', "d4", "d3.5_d4", "d3.5", "d3_d3.5", "d3")
tp_order <- tp_order[tp_order %in% unique(meta_query$timepoint)]
mapping <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_query, map_meta = meta_query,
  tp_order = tp_order,
  k = io$k, npcs = io$npcs
)
genemeta[,new_hvg := (ens_id %in% rownames(metadata(mapping$pca)$rotation))]
genemeta <- setnames(genemeta, c("new_hvg"), c(paste0("ExtendedAtlas_", opts$query_batches, "_hvg")))


mapping$mapping <- get_meta.extended(correct_atlas = mapping$correct_atlas,
                                     atlas_meta = meta_atlas,
                                     correct_map = mapping$correct_map,
                                     map_meta = meta_query,
                                     k_map = io$k
                                    )

message("Computing mapping scores...") 
for (i in seq(from = 1, to = io$k)) {
mapping$closest.cells[[i]]     <- sapply(mapping$mapping, function(x) x$cells.mapped[i])
mapping$celltypes.mapped[[i]]  <- sapply(mapping$mapping, function(x) x$celltypes.mapped[i])
mapping$cellstages.mapped[[i]] <- sapply(mapping$mapping, function(x) x$stages.mapped[i])
}  
multinomial.prob <- getMappingScore(mapping)
message("Done\n")

message("Writing output...") 
ct <- sapply(mapping$mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
st <- sapply(mapping$mapping, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
cm <- sapply(mapping$mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
mapping$mapping <- data.frame(
  cell                      = names(mapping$mapping), 
  celltype.mapped.extended  = unlist(ct),
  stage.mapped.extended     = unlist(st),
  closest.cell.extended     = unlist(cm),
  stringsAsFactors=FALSE)

mapping$mapping <- cbind(mapping$mapping,multinomial.prob)
message("Done\n")

meta <- fread(paste0(opts$query.metadata)) %>%
  merge(mapping$mapping %>% as.data.table, by="cell", all.x=TRUE)


#######################
## Assign Descendant ##
#######################

message("Assigning likely descendants...")
descendantmap <- get_meta.descendant(mapping$correct_atlas[meta_atlas$cell,],
                                     meta_atlas,
                                     mapping$correct_map[meta_query$cell,],
                                     meta_query
                                    )

descendantmap.out <- list()
for (i in seq(from = 1, to = io$k)) {
    descendantmap.out$closest.cells.descendant[[i]]     <- sapply(descendantmap, function(x) x$cells.mapped.descendant[i])
    descendantmap.out$celltypes.mapped.descendant[[i]]  <- sapply(descendantmap, function(x) x$celltypes.mapped.descendant[i])
    descendantmap.out$cellstages.mapped.descendant[[i]] <- sapply(descendantmap, function(x) x$stages.mapped.descendant[i])
}  
tmp.multinomial.prob <- getMappingScore.descendant(descendantmap.out)

ct <- sapply(descendantmap, function(x) x$celltype.mapped.descendant); is.na(ct) <- lengths(ct) == 0
st <- sapply(descendantmap, function(x) x$stage.mapped.descendant); is.na(st) <- lengths(st) == 0
cm <- sapply(descendantmap, function(x) x$cells.mapped.descendant[1]); is.na(cm) <- lengths(cm) == 0
descendantmap.out$mapping <- data.frame(
  cell                                = names(descendantmap), 
  celltype.mapped.descendant.extended = unlist(ct),
  stage.mapped.descendant.extended    = unlist(st),
  closest.cell.descendant.extended    = unlist(cm),
  stringsAsFactors           = FALSE)

descendantmap.out$mapping <- cbind(descendantmap.out$mapping, tmp.multinomial.prob)
descendantmap.mapping <- descendantmap.out$mapping

meta <- merge(meta, descendantmap.mapping, by="cell", all=TRUE)

mapping$descendantmap.out <- descendantmap.out
message("Done\n")


###############################
## Reassign Mesodermal Cells ##
###############################

#meta[,celltype.mapped.base := celltype.mapped]
#atlas.mesocells <- meta_atlas[(celltype_old %in% c("Paraxial mesoderm", "Somitic mesoderm", "Intermediate mesoderm")) & (stage == "E8.5")]$cell

#gastr.mesocells <- intersect(meta[(celltype.mapped %in% c("Paraxial mesoderm", "Somitic mesoderm", "Intermediate mesoderm"))]$cell,
#                                       rownames(mapping$correct_map)
#                                      )
#mesomap <- get_meta.meso(mapping$correct_atlas[atlas.mesocells,],
#                         meta_atlas[cell %in% atlas.mesocells],
#                         mapping$correct_map[gastr.mesocells,],
#                         meta[cell %in% gastr.mesocells]
#                        )
#
#mesomap.out <- list()
#for (i in seq(from = 1, to = io$k)) {
#    mesomap.out$closest.cells.meso[[i]]     <- sapply(mesomap, function(x) x$cells.mapped[i])
#    mesomap.out$celltypes.mapped.meso[[i]]  <- sapply(mesomap, function(x) x$celltypes.mapped[i])
#    mesomap.out$cellstages.mapped.meso[[i]] <- sapply(mesomap, function(x) x$stages.mapped[i])
#}  
#tmp.multinomial.prob <- getMappingScore.meso(mesomap.out)

#ct <- sapply(mesomap, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
#st <- sapply(mesomap, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
#cm <- sapply(mesomap, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
#mesomap.out$mapping <- data.frame(
#  cell                 = names(mesomap), 
#  celltype.mapped.meso = unlist(ct),
#  stage.mapped.meso    = unlist(st),
#  closest.cell.meso    = unlist(cm),
#  stringsAsFactors=FALSE)

#mesomap.out$mapping <- cbind(mesomap.out$mapping,tmp.multinomial.prob)
#mesomap.mapping <- mesomap.out$mapping

#meta <- merge(meta, mesomap.mapping, by="cell", all=TRUE)
#meta[!is.na(celltype.mapped.meso), celltype.mapped := celltype.mapped.meso]

#mapping$mesomap.out <- mesomap.out

##########
## Save ##
##########

message("Saving Output...") 
saveRDS(mapping, opts$output.RDS)
fwrite(meta, opts$output.metadata, sep="\t", na="NA", quote=F)
fwrite(genemeta, opts$output.gene.metadata, quote=F, na="NA", sep="\t")
message("Done\n")