#!/usr/bin/env Rscript

# /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/heterogeneity_implications/load_ct_DA_files.R -o /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/heterogeneity_simulations/20220628/ -c "Early Spinal Cord" -d "poisson"

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tibble))


###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-o", "--outbase"), type="character", default=NULL, 
                help="basic outpath", metavar="character"),
    make_option(c("-c", "--ct_perturb"), type="character", default=NULL, 
                help="cell type that was perturbed", metavar="character"),
    make_option(c("-d", "--da_meth"), type="character", default=NULL, 
                help="DA testing method", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

ntrials <- 100
sim_types <- c('full')
timepoints <- c('d4.5')
IoDs <- c()
fc_seq <- c( 0, 0.25, 0.5, 0.75 )
nsamps <- c( 2, 3, 4, 5, 6, 7, 8, 9, 10 )
norgs <- c( 1, 5, 10, 15, 20, 25, 30 )
IoD_contrs <- c(500)

for (simtype in sim_types) {
    for (da_meth in c(opts$da_meth)) {
        for (tp in timepoints) {
            for (IoD_contr in IoD_contrs) {
                for (norg in norgs) {
                    for (nsamp in nsamps) {
                        for (t in 1:ntrials) {
                            for (fc in fc_seq) {
                                for (ct in c(opts$ct_perturb)) {
                                    if (simtype == 'gampois') {
                                        for (IoD in IoDs) {
                                            f <- paste0(opts$outbase,
                                                        'trial', t, '_',
                                                        'gastr_sim_',
                                                        simtype, '_',
                                                        tp, '_',
                                                        da_meth, '_',
                                                        ct, '_',
                                                        IoD, '_',
                                                        nsamp, '_',
                                                        norg, '_',
                                                        IoD_contr, '_out.txt.gz')
                                            if (file.exists(f)) {
                                                tmp <- fread(f)
                                                tmp$trial <- t
                                                tmp$timepoint <- tp
                                                tmp$DA_meth <- da_meth
                                                tmp$ct_perturb <- ct
                                                tmp$norg <- norg
                                                tmp$nsamp <- nsamp
                                                tmp$fc <- fc
                                                
                                                if (exists("out")) {
                                                    out <- rbind(out, tmp)
                                                } else {
                                                    out <- tmp
                                                }
                                            }
                                        }
                                    } else {
                                        f <- paste0(opts$outbase,
                                                    'trial', t, '_',
                                                    'gastr_sim_',
                                                    simtype, '_',
                                                    tp, '_',
                                                    da_meth, '_',
                                                    ct, '_',
                                                    fc, '_',
                                                    NA, '_',
                                                    nsamp, '_',
                                                    norg, '_',
                                                    IoD_contr, '_out.txt.gz')
                                        if (file.exists(f)) {
                                            tmp <- fread(f)
                                            tmp$trial <- t
                                            tmp$timepoint <- tp
                                            tmp$DA_meth <- da_meth
                                            tmp$ct_perturb <- ct
                                            tmp$norg <- norg
                                            tmp$nsamp <- nsamp
                                            tmp$fc <- fc

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
                    }
                }
            }
        }
    }
}

fwrite(out, paste0(opts$outbase, 'gastr_sim_', opts$da_meth, '_', opts$ct_perturb, '_out.txt.gz'))