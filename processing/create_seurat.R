#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scds))
suppressPackageStartupMessages(library(dplyr))

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-i", "--inputdir"), type="character", default=NULL, 
                help="directory of the input files", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-O", "--output.seurat"), type="character", default=NULL, 
                help="path to output seurat file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-g", "--output.gene.metadata"), type="character", default=NULL, 
                help="path to output gene metadata file", metavar="character"),
    make_option(c("-b", "--output.barcodes"), type="character", default=NULL, 
                help="path to output gene metadata file", metavar="character"),
    make_option(c("-p", "--subset.proteincoding"), type="character", default=NULL, 
                help="path to proteincoding genes, if subsetting", metavar="character"),
    make_option(c("-e", "--experiment"), type="character", default=NULL, 
                help="specify experiment", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default=NULL, action="store",
                help="specify sample", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (!is.null(opts$sample)) { opts$sample <-substr(opts$sample,1,nchar(opts$sample)-1) }

if (is.null(opts$inputdir)){
    print_help(opt_parser)
    stop("An input directory must be supplied.n", call.=FALSE)
} else if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.seurat)) {
    print_help(opt_parser)
    stop("A directory to a Seurat output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.metadata)) {
    print_help(opt_parser)
    stop("A directory to a sample metadata output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.gene.metadata)) {
    print_help(opt_parser)
    stop("A directory to a cleaned gene metadata file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.barcodes)) {
    print_help(opt_parser)
    stop("A directory to a cleaned barcode output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$experiment)) {
    print_help(opt_parser)
    stop("The experiment must be specified.n", call.=FALSE)
}

##########################
## Source settings file ##
##########################

source(opts$settings)

if ((!is.null(io$qc)) && (io$qc != TRUE)){
    stop("In the settings file, qc must be NULL or TRUE.n", call.=FALSE)
} else if ((!is.null(io$subset.proteincoding)) && (!is.string(io$subset.proteincoding))){
    stop("In the settings file, subset.proteincoding must be NULL or a path.n", call.=FALSE)
}else if ((!is.null(io$scaledata)) && (io$scaledata != TRUE)){
    stop("In the settings file, scaledatav must be NULL or TRUE.n", call.=FALSE)
}


##############################
## Load and merge data sets ##
##############################

message("Loading datasets")

mtx <- list()
cell.info <- list()
gene.info <- list()
  
file.list <- list.files(paste0(opts$inputdir, '/', opts$sample), full.names = TRUE)
barcode.loc <- file.list[grepl("barcodes.tsv.gz", file.list)]
barcode.loc <- barcode.loc[!grepl("unfiltered", barcode.loc)]
gene.loc <- file.list[grepl("features.tsv.gz", file.list)]
gene.loc <- gene.loc[!grepl("unfiltered", gene.loc)]
matrix.loc <- file.list[grepl("matrix.mtx.gz", file.list)]
matrix.loc <- matrix.loc[!grepl("unfiltered", matrix.loc)]

if ((length(barcode.loc) != 1) | (length(gene.loc) != 1) | (length(matrix.loc) != 1)) {
  next
}

# Load cell metadata
cell.info <- read.delim(barcode.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
colnames(cell.info) <- c("barcode")

if (grepl('exp2|exp3', opts$experiment)) {
tmp <- strsplit(basename(barcode.loc),"_")[[1]]
abrv_exp <- strsplit(opts$experiment,"/")[[1]][2]
cell.info$experiment <- abrv_exp
if (tmp[1] == 'd3') {
    cell.info$timepoint <- 'd3_d3.5'
} else if (tmp[1] == 'd4') {
    cell.info$timepoint <- 'd4_d4.5'
} else {
    cell.info$timepoint <- 'unknown'
}
cell.info$organoid <- "gastr"
cell.info$tech <- "MULTI"
cell.info$class <- paste0("gastr_",cell.info$timepoint)
cell.info$lane <- tmp[4]
cell.info$batch <- paste0(abrv_exp,"_",opts$sample)
} else if (grepl('exp4|exp5|exp6', opts$experiment)) {
tmp <- strsplit(basename(barcode.loc),"_")[[1]]
abrv_exp <- strsplit(opts$experiment,"/")[[1]][2]
cell.info$experiment <- abrv_exp
if ((tmp[3] == 'd3') & (tmp[4] != '5')) {
    cell.info$timepoint <- 'd3'
} else if ((tmp[3] == 'd3') & (tmp[4] == '5')) {
    cell.info$timepoint <- 'd3.5_d4'
} else if (tmp[3] == 'd4') {
    cell.info$timepoint <- 'd4.5_d5'
} else {
    cell.info$timepoint <- 'unknown'
}
cell.info$organoid <- "gastr"
cell.info$tech <- "MULTI"
cell.info$class <- paste0("gastr_",cell.info$timepoint)
cell.info$lane <- tmp[10]
cell.info$batch <- paste0(abrv_exp,"_",opts$sample)
} else if (grepl('exp1', opts$experiment, fixed = TRUE)) {
tmp <- strsplit(basename(barcode.loc),"_")[[1]]
abrv_exp <- strsplit(opts$experiment,"/")[[1]][2]
cell.info$experiment <- abrv_exp
cell.info$timepoint <- tmp[4]
cell.info$organoid <- "gastr"
cell.info$tech <- "MULTI"
cell.info$class <- paste0("gastr_",tmp[4])
cell.info$lane <- tmp[6]
cell.info$batch <- paste0(abrv_exp,"_",opts$sample)
} else {
cell.info$experiment <- opts$experiment
cell.info$timepoint <- strsplit(opts$experiment,"_")[[1]][2]
cell.info$organoid <- strsplit(opts$experiment,"_")[[1]][1]
cell.info$tech <- "10X"
cell.info$class <- opts$experiment
cell.info$lane <- NA
cell.info$batch <- opts$experiment
}

# Load gene metadata (note we could just load this once)
gene.info <- read.delim(gene.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")[,c(1,2)]
colnames(gene.info) <- c("ens_id","symbol")

# Load matrix
mtx <- Matrix::readMM(matrix.loc)
rownames(mtx) <- gene.info$symbol
colnames(mtx) <- cell.info$barcode

# Remove cells with few UMIs 
foo <- Matrix::colSums(mtx)>io$minUMIs
cells.to.keep <- colnames(mtx)[foo]
mtx <- mtx[,cells.to.keep]

# Subset cell metadata
cell.info <- as.data.table(cell.info)[barcode %in% colnames(mtx),]

# bind gene names and remove human alignments
rownames(gene.info) <- NULL
gene.info <- unique(gene.info)

# Concatenate cell  metadata
cell.info$cell <- paste("cell",1:nrow(cell.info),cell.info$batch,sep="_")
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
colnames(mtx) <- cell.info$cell

################
## Processing ##
################

message("Processing datasets...")

# Optionally subset protein-coding genes
if (!is.null(opts$subset.proteincoding)){
    genes <- fread(opts$subset.proteincoding)[,ens_id]
    genes <- genes[genes %in% mouse.genes]
    mouse.genes <- mouse.genes[mouse.genes %in% genes]
    mtx <- mtx[mouse.genes,]
}

########
## QC ##
########

if (io$qc == TRUE) {

    message("QCing...")

    srat <- CreateSeuratObject(mtx, meta.data = cell.info)
    srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
    ribo.genes <- c(grep(pattern = "^Rpl", x = rownames(srat), value = TRUE),grep(pattern = "^Rps", x = rownames(srat), value = TRUE))
    srat[["percent.ribo"]] <- PercentageFeatureSet(srat, features = ribo.genes)

    md <- as.data.table(srat@meta.data)

    foo <- ifelse(is.null(opts$sample), io$min_nFeature_RNA[[opts$experiment]], io$min_nFeature_RNA[[opts$experiment]][[opts$sample]])
    bar <- ifelse(is.null(opts$sample), io$min_nCount_RNA[[opts$experiment]], io$min_nCount_RNA[[opts$experiment]][[opts$sample]])
    baz <- ifelse(is.null(opts$sample), io$max_percent.mt[[opts$experiment]], io$max_percent.mt[[opts$experiment]][[opts$sample]])
    bam <- ifelse(is.null(opts$sample), io$min_percent.mt[[opts$experiment]], io$min_percent.mt[[opts$experiment]][[opts$sample]])
    srat_sub <- subset(srat, subset = nFeature_RNA > foo & nCount_RNA > bar & percent.mt > bam & percent.mt < baz)

    sce_tmp <- as.SingleCellExperiment(srat_sub)
    sce_tmp <- cxds_bcds_hybrid(sce_tmp, estNdbl=TRUE)
    md <- merge(as.data.table(md), as.data.table(colData(sce_tmp)), all.x=TRUE)
    md[,pass_QC := FALSE]
    md[cell %in% Cells(srat_sub), pass_QC := TRUE]
    rownames(md) <- md$cell
    srat_sub@meta.data <- md[pass_QC==TRUE]
    rownames(srat_sub@meta.data) <- srat_sub@meta.data$cell
    
    srat_sub <- NormalizeData(srat_sub)
    srat_sub <- FindVariableFeatures(srat_sub)
    srat_sub@assays$RNA@scale.data <- as.matrix(srat_sub@assays$RNA@data)
    srat_sub <- RunPCA(srat_sub, verbose=FALSE)
    srat_sub <- FindNeighbors(srat_sub, k.param=io$k, reduction = "pca")
    rownames(srat_sub@meta.data) <- srat_sub@meta.data$cell
    srat_sub <- CellCycleScoring(srat_sub,
                                 s.features = as.data.table(gene.info)[ens_id %in% fread(io$s.genes.ensid, header = FALSE)$V1]$symbol,
                                 g2m.features = as.data.table(gene.info)[ens_id %in% fread(io$g2m.genes.ensid, header = FALSE)$V1]$symbol,
                                 set.ident = FALSE
                                )

    tmp <- data.table(hybrid_map = rowSums(srat_sub@graphs$RNA_nn[,md[hybrid_call==TRUE]$cell])/io$k,
                      cell = rownames(srat_sub@graphs$RNA_snn))
    md <- merge(md, tmp, by = "cell", all.x=TRUE)
    md[hybrid_call==TRUE, hybrid_map := NA]
    md[,hybrid_call2 := TRUE]
    #md[hybrid_map<=((sum(md[pass_QC==TRUE]$hybrid_call==TRUE))/(dim(md[pass_QC==TRUE])[1])), hybrid_call2 := FALSE]
    md[hybrid_map <= 1/3, hybrid_call2 := FALSE]
    md[pass_QC==FALSE, hybrid_call2 := NA]
    rownames(md) <- md$cell

    cell.info <- merge(md, srat_sub@meta.data[,c("cell", "S.Score", "G2M.Score", "Phase")], by="cell", all.x=TRUE)
    rownames(cell.info) <- cell.info$cell
    rm(sce_tmp, srat_sub, srat, md)

}

############
## Seurat ##
############

# Create seurat object

message("Creating seurat...")

rownames(mtx) <- gene.info$ens_id
srat <- CreateSeuratObject(mtx, meta.data = cell.info)


##########################
## SingleCellExperiment ##
##########################

# sce <- as.SingleCellExperiment(srat)

##########
## Save ##
##########

saveRDS(srat, opts$output.seurat)
fwrite(cell.info, opts$output.metadata, quote=F, na="NA", sep="\t")
fwrite(gene.info, opts$output.gene.metadata, quote=F, na="NA", sep="\t")

if (io$qc == TRUE) {

    message("Saving updated barcodes...")

    fwrite(list(cell.info[(pass_QC==TRUE) & (!(hybrid_call2))]$barcode), opts$output.barcodes)

}
