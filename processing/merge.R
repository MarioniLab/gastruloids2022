#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-S", "--seurat_files"), type="character", default=NULL, 
                help="list of seurat files", metavar="character"),
    make_option(c("-F", "--final_calls_files"), type="character", default=NULL, 
                help="list of final calls files", metavar="character"),
    make_option(c("-G", "--gene_meta_files"), type="character", default=NULL, 
                help="list of gene metadata files", metavar="character"),
    make_option(c("-O", "--output.seurat"), type="character", default=NULL, 
                help="path to output seurat file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-g", "--output.genemetadata"), type="character", default=NULL, 
                help="path to output gene metadata file", metavar="character"),
    make_option(c("-s", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (is.null(opts$seurat_files)){
    print_help(opt_parser)
    stop("An list of seurat files must be supplied.n", call.=FALSE)
} else if (is.null(opts$final_calls_files)) {
    print_help(opt_parser)
    stop("A list of final calls files must be supplied.n", call.=FALSE)
} else if (is.null(opts$gene_meta_files)) {
    print_help(opt_parser)
    stop("A list of gene metadata files must be supplied.n", call.=FALSE)
} else if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
}


##########################
## Source settings file ##
##########################

source(opts$settings)


if (is.string(opts$seurat_files)) opts$seurat_files <- str_trim(strsplit(opts$seurat_files, ",")[[1]])
if (is.string(opts$final_calls_files)) opts$final_calls_files <- str_trim(strsplit(opts$final_calls_files, ",")[[1]])
if (is.string(opts$gene_meta_files)) opts$gene_meta_files <- str_trim(strsplit(opts$gene_meta_files, ",")[[1]])

srat_list <- list()
mtx_list <- list()
cell.info_list <- list()
gene_list <- list()
for (f in opts$seurat_files) {
    
    srat_list[[f]] <- readRDS(f)
    
    cell.info_list[[f]] <- as.data.table(srat_list[[f]]@meta.data)
    mtx_list[[f]] <- srat_list[[f]]@assays$RNA@counts
    gene_list[[f]] <- rownames(mtx_list[[f]])
    
}

genes <- Reduce(intersect, gene_list)

for (f in opts$seurat_files) {
    
    mtx_list[[f]] <- mtx_list[[f]][genes,]
    
}

gene.meta_list <- list()
for (f in opts$gene_meta_files) {
    
    gene.meta_list[[f]] <- fread(f)[ens_id %in% genes]
    
}
gene.meta <- Reduce(union, gene.meta_list)
gene.meta <- unique(gene.meta, by = c('ens_id'))

cell.info <- rbindlist(cell.info_list, fill=TRUE)
rownames(cell.info) <- cell.info$cell


############################################
## Add MULTI-seq Calls to Seurat metadata ##
############################################

myfunc1 <- function(x) {
    if (x[["barcode"]] %in% names(final.calls.tmp)) {
        bar <- final.calls.tmp[[x[["barcode"]]]]
        if (bar %in% c("Doublet", "Negative")) {
            return(bar)
        } else {
            return(paste0(final.calls.tmp[[x[["barcode"]]]], "_", x[["experiment"]]))
        }
        
    } else {
        return(NA)
    }
}

myfunc2 <- function(x) {
    if (x[["experiment"]] %in% io$mixedexps) {
        if (grepl(paste(paste0(io$earlierbars[[x[["experiment"]]]], "_", x[["experiment"]]), collapse="|"), x[["MULTI_class"]])) {
            return(strsplit(x[["timepoint"]], '_')[[1]][1])
        } else if (grepl(paste(paste0(io$laterbars[[x[["experiment"]]]], "_", x[["experiment"]]), collapse="|"), x[["MULTI_class"]])) {
            return(strsplit(x[["timepoint"]], '_')[[1]][2])
        } else {
            return(x[["timepoint"]])
        }
    } else {
        return(x[["timepoint"]])
    }
}

final.calls <- list()
final.calls.tmp <- list()
exps <- unique(cell.info$experiment)
for (e in exps) {
    fs <- opts$final_calls_files[grepl(e, opts$final_calls_files)]
    sfs <- opts$seurat_files[grepl(e, opts$seurat_files)]
    
    samps <- unlist(lapply(unique(cell.info[experiment==e]$batch), function(b) strsplit(b, "_")[[1]] %>% .[length(.)]))
    for (s in samps) {
        f <- fs[grepl(s, fs)]
        if (length(f) == 0) {
            next
        }
        sf <- sfs[grepl(s, sfs)]
    
        final.calls[[f]] <- readRDS(f)
        final.calls.tmp <- final.calls[[f]]
        names(final.calls.tmp) <- paste0(names(final.calls.tmp),"-1")

        cell.info_list[[sf]]$MULTI_class <- apply(cell.info_list[[sf]], 1, myfunc1)
        
    }
    
}
                           
cell.info <- rbindlist(cell.info_list, fill=TRUE)
rownames(cell.info) <- cell.info$cell
cell.info[is.na(MULTI_class), MULTI_class := "unidentifiable"]

cell.info[,timepoint.base:=timepoint]
cell.info$timepoint <- apply(cell.info, 1, myfunc2)

mtx <- do.call(cbind, mtx_list)
srat <- CreateSeuratObject(mtx, meta.data = cell.info)


##########
## Save ##
##########

saveRDS(srat, opts$output.seurat)
fwrite(cell.info, opts$output.metadata, quote=F, na="NA", sep="\t")
fwrite(gene.meta, opts$output.genemetadata, quote=F, na="NA", sep="\t")