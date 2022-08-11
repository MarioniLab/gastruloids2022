#!/usr/bin/env Rscript

# run in MULTIseq conda env
# run e.g. with: bsub -M 100000 ./barcode_align.R -i /nfs/research1/marioni/Leah/scNMTseq/ModelSystems/10x_old_EBsvGastruloids/data/MULTI/exp{??}_d* -S /nfs/research1/marioni/Leah/scNMTseq/ModelSystems/gastruloids_scRNAseq/pipeline/processing/settings.R -O /nfs/research1/marioni/Leah/scNMTseq/ModelSystems/10x_old_EBsvGastruloids/processed_files_Leah/MULTI/exp{??}_d* -e exp{??}_d* -s sample{?}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(deMULTIplex))
suppressPackageStartupMessages(library(optparse))

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-i", "--inputdir"), type="character", default=NULL, 
                help="directory of the input files", metavar="character"),
    make_option(c("-b", "--barcodes"), type="character", default=NULL,
                help="path to a file with a list of cell barcodes", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-O", "--output.base"), type="character", default=NULL, 
                help="output directory path", metavar="character"),
    make_option(c("-e", "--experiment"), type="character", default="MULTI", 
                help="specify experiment (currently MULTI)", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default=NULL, 
                help="which sample to process", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (!is.null(opts$sample)) { opts$sample <-substr(opts$sample,1,nchar(opts$sample)-1) }

# source settings file
source(opts$settings)

if (is.null(opts$inputdir)){
    print_help(opt_parser)
    stop("An input directory must be supplied.n", call.=FALSE)
} else if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.base)) {
    print_help(opt_parser)
    stop("A base output directory must be supplied.n", call.=FALSE)
} else if (is.null(opts$experiment)) {
    print_help(opt_parser)
    stop("The experiment must be specified.n", call.=FALSE)
}

if(is.null(opts$barcodes)) {
    cell.barcode.loc <- list.files(path = paste0(opts$inputdir, "/", opts$sample), pattern="*_barcodes.tsv.gz", full.names=TRUE)[1]
} else {
    cell.barcode.loc <- opts$barcodes
}

cell.id.vec <- read.delim(cell.barcode.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")$V1
cell.id.vec <- as.vector(t(data.frame(strsplit(cell.id.vec, "-"))[1,]))

gastr.barcode.loc <- sprintf("%s/BClist.tsv",opts$inputdir)
bar.ref <- read.delim(gastr.barcode.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")$V1

#if (((opts$experiment == "MULTI/exp4_d3") & (opts$sample == "sampleA")) | ((opts$experiment == "MULTI/exp5_d3.5_d4") & (opts$sample == "sampleA"))) {
if (1 == 2) {
    R1 <- paste0(opts$inputdir, "/MiSeqNano_barcodes/", io$barcode_files[[opts$experiment]][[opts$sample]], "_R1.fastq.gz")
    R2 <- paste0(opts$inputdir, "/MiSeqNano_barcodes/", io$barcode_files[[opts$experiment]][[opts$sample]], "_R2.fastq.gz")
} else if (opts$experiment == "MULTI/exp4_d3") {
    R1 <- paste0(opts$inputdir, "/DeepSeq_barcodes/", io$barcode_files[[opts$experiment]][[opts$sample]], "_R1.fastq.gz")
    R2 <- paste0(opts$inputdir, "/DeepSeq_barcodes/", io$barcode_files[[opts$experiment]][[opts$sample]], "_R2.fastq.gz")
} else {
    R1 <- paste0(opts$inputdir, "/barcodes/", io$barcode_files[[opts$experiment]][[opts$sample]], "_R1.fastq.gz")
    
    if (opts$experiment %in% c("MULTI/exp4_d3", "MULTI/exp5_d3.5_d4", "MULTI/exp6_d4.5_d5")) {
        R2 <- paste0(opts$inputdir, "/barcodes/", io$barcode_files[[opts$experiment]][[opts$sample]], "_R2.fastq.gz")
    } else {
        R2 <- paste0(opts$inputdir, "/barcodes/", io$barcode_files[[opts$experiment]][[opts$sample]], "_R3.fastq.gz")
    }
}



readTable <- MULTIseq.preProcess(R1 = R1,
                                 R2 = R2,
                                 cellIDs = cell.id.vec,
                                 cell=io$cell_pos,
                                 umi=io$umi_pos,
                                 tag=io$tag_pos
                                )

bar.table <- MULTIseq.align(readTable, cell.id.vec, bar.ref)

saveRDS(bar.table, paste0(opts$output.base,"barTable.rds"))
saveRDS(readTable, paste0(opts$output.base,"readTable.rds"))