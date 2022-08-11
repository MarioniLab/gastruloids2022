#!/usr/bin/env Rscript

# /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/load_ct_DA_files.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/20220628/ -c "Early Spinal Cord" -d "poisson"

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-o", "--outbase"), type="character", default=NULL, 
                help="basic outpath", metavar="character"),
    make_option(c("-c", "--ct1"), type="character", default=NULL, 
                help="cell type that was perturbed", metavar="character"),
    make_option(c("-C", "--ct2"), type="character", default=NULL, 
                help="DA testing method", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

gastr_types <- c('All', 'mesodermal', 'neural', 'intermediate')
tps <- c('All', 'd3', 'd3.5', 'd4', 'd4.5', 'd5')
for (gt1 in gastr_types) {
    for (gt2 in gastr_types) {
        for (tp1 in tps) {
            for (tp2 in tps) {
                file <- paste0(opts$outbase,
                               tp1, '_',
                               gt1, '_',
                               opts$ct1, '_vs_',
                               tp2, '_',
                               gt2, '_',
                               opts$ct2, '_out.csv.gz'
                              )
                if (file.exists(file)) {

                    tmp <- fread(file)[!is.na(sig)]
                    tmp[,detection_rate_group1 := signif(detection_rate_group1, 3)]
                    tmp[,detection_rate_group2 := signif(detection_rate_group2, 3)]
                    tmp[,logFC := signif(logFC, 3)]
                    tmp[,p.value := signif(p.value, 3)]
                    tmp[,log_padj_fdr := signif(log_padj_fdr, 3)]

                    tmp[,padj_fdr := NULL]

                    tmp$ct1 <- opts$ct1
                    tmp$ct2 <- opts$ct2
                    tmp$tp1 <- tp1
                    tmp$tp2 <- tp2
                    tmp$gt1 <- gt1
                    tmp$gt2 <- gt2

                    if (exists("out")) {
                        out <- rbind(out, tmp)
                    } else {
                        out <- tmp
                    }

                }
            }
        }
    }
}

if (exists("out")) {
    fwrite(out, paste0(opts$outbase, opts$ct1, '_', opts$ct2, '_out.txt.gz'))
}