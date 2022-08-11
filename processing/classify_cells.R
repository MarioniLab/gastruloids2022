#!/usr/bin/env Rscript

# run in MULTIseq conda env

suppressPackageStartupMessages(library(ggplot2))
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
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-e", "--experiment"), type="character", default="MULTI", 
                help="specify experiment (currently MULTI)", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default=NULL, 
                help="which sample to process", metavar="character"),
    make_option(c("-p", "--plot.type"), type="character", default="pdf",
                help="what device should be used in ggsave", metavar="character")
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
} else if (is.null(opts$experiment)) {
    print_help(opt_parser)
    stop("The experiment must be specified.n", call.=FALSE)
}

bar.table <- readRDS(paste0(opts$inputdir,opts$experiment,"/",opts$sample,"_barTable.rds"))

full_outbase <- paste0(opts$inputdir, "processing/MULTIseq_nUMIs/", opts$experiment, "/")
if (!dir.exists(full_outbase)) {
    dir.create(full_outbase, recursive = TRUE)
}

p <- ggplot(bar.table, aes(x=nUMI)) +
    geom_density(fill="darkgrey") +
    scale_x_continuous(trans="log2")
ggsave(paste0(full_outbase, opts$sample, "_nUMI.", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")

p <- ggplot(bar.table, aes(x=nUMI_total)) +
    geom_density(fill="darkgrey") +
    scale_x_continuous(trans="log2")
ggsave(paste0(full_outbase, opts$sample, "_nUMI_total.", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")

p <- ggplot(bar.table, aes(x=nUMI_total/nUMI)) +
    geom_density(fill="darkgrey") +
    scale_x_continuous(trans="log2")
ggsave(paste0(full_outbase, opts$sample, "_nUMI_total_logFC.", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")

## Visualize barcode space
bar.tsne <- barTSNE(bar.table[,io$bars_used[[opts$experiment]][[opts$sample]]])

full_outbase <- paste0(opts$inputdir, "processing/MULTIseq_bar_tsnes/", opts$experiment, "/")
if (!dir.exists(full_outbase)) {
    dir.create(full_outbase, recursive = TRUE)
}

for (i in 3:ncol(bar.tsne)) {
    p <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
    ggsave(paste0(full_outbase, opts$sample, "_bar",i,".", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")
}

## ASSIGN CELLS

bar.table.full <- bar.table[,1:96]
good.bars <- paste("Bar",io$bars_used[[opts$experiment]][[opts$sample]],sep="")  # Barcodes 1:24 were detected
bar.table <- bar.table.full[, good.bars]  # Remove missing bars and summary columns

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  #print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

full_outbase <- paste0(opts$inputdir, "processing/MULTIseq_qsweep/", opts$experiment, "/")
if (!dir.exists(full_outbase)) {
    dir.create(full_outbase, recursive = TRUE)
}

## Round 1 -----------------------------------------------------------------------------------------------------
## Perform Quantile Sweep

i <- 1
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds

threshold.results <- findThresh(call.list=bar.table_sweep.list)
p <- ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
ggsave(paste0(full_outbase, opts$sample, "_round",i,".", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")

## Finalize round 1 classifications, remove negative cells
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

while (length(names(round.calls)[which(round.calls == "Negative")]) > 0) {
    
    i <- i+1
    bar.table_sweep.list <- list()
    n <- 0
    for (q in seq(0.01, 0.99, by=0.02)) {
      #print(q)
      n <- n + 1
      bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
      names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
    }
    
    error1 <- FALSE
    tryCatch(threshold.results <- findThresh(call.list=bar.table_sweep.list),
             error = function(e) { error1 <<- TRUE}
            )
    if (error1==TRUE) {
        break
    }

    threshold.results <- findThresh(call.list=bar.table_sweep.list)
    p <- ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
      geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
    ggsave(paste0(full_outbase, opts$sample, "_round",i,".", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")
    round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
    neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
    bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
    
}

final.calls <- c(round.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round.calls),neg.cells)

if (opts$experiment == 'MULTI/exp4_d3') {
    
    if (opts$sample == 'sampleA') {
        message('in exp4_d3 sampleA')
        final.calls2 <- final.calls
        tmp <- bar.table.full[max.col(bar.table.full[names(final.calls2),good.bars], ties.method='first') == (good.bars == 'Bar61'),] %>% rownames
        final.calls2[names(final.calls2[tmp][final.calls2[tmp] %in% c('Bar6', 'Bar14')])] <- 'Bar61'
        for (bar in c('Bar6', 'Bar14')) {
            tmp <- bar.table.full[names(final.calls2)[final.calls2 == bar],good.bars]
            final.calls2[rownames(tmp)[good.bars[max.col(tmp, ties.method="first")] != bar]] <- 'Negative'
            tmp2 <- tmp[good.bars[max.col(tmp, ties.method="first")] == bar,bar]
            names(tmp2) <- rownames(tmp[good.bars[max.col(tmp, ties.method="first")] == bar,])
            final.calls2[names(tmp2[tmp2<mean(tmp2)])] <- 'Negative'
        }
        tmp <- bar.table.full[names(final.calls)[final.calls == 'Doublet'],good.bars]
        final.calls2[rownames(tmp)[good.bars[max.col(tmp, ties.method="first")] %in% c("Bar6", "Bar14")]] <- 'Negative'
        tmp2 <- tmp[good.bars[max.col(tmp, ties.method="first")] %in% c("Bar6", "Bar14"),good.bars]
        argmax2 <- max.col(replace(as.matrix(tmp2), cbind(1:dim(tmp2)[1], max.col(as.matrix(tmp2))), -Inf))
        final.calls2[rownames(tmp2)[good.bars[argmax2] %in% c("Bar6", "Bar14")]] <- 'Negative'
        final.calls <- final.calls2
        good.bars.sub <- good.bars[!(good.bars %in% c("Bar6", "Bar14", "Bar61"))]
    }
    
    if (opts$sample == 'sampleB') {
        message('in exp4_d3 sampleB')
        final.calls2 <- final.calls
        for (bar in c("Bar12", "Bar24", "Bar30")) {
            tmp <- bar.table.full[names(final.calls2)[final.calls2 == bar],good.bars]
            final.calls2[rownames(tmp)[good.bars[max.col(tmp, ties.method="first")] != bar]] <- 'Negative'
            tmp2 <- tmp[good.bars[max.col(tmp, ties.method="first")] == bar,bar]
            names(tmp2) <- rownames(tmp[good.bars[max.col(tmp, ties.method="first")] == bar,])
            final.calls2[names(tmp2[tmp2<mean(tmp2)])] <- 'Negative'
        }
        tmp <- bar.table.full[names(final.calls)[final.calls == 'Doublet'],good.bars]
        final.calls2[rownames(tmp)[good.bars[max.col(tmp, ties.method="first")] %in% c("Bar12", "Bar24", "Bar30")]] <- 'Negative'
        tmp2 <- tmp[good.bars[max.col(tmp, ties.method="first")] %in% c("Bar12", "Bar24", "Bar30"),good.bars]
        argmax2 <- max.col(replace(as.matrix(tmp2), cbind(1:dim(tmp2)[1], max.col(as.matrix(tmp2))), -Inf))
        final.calls2[rownames(tmp2)[good.bars[argmax2] %in% c("Bar12", "Bar24", "Bar30")]] <- 'Negative'
        final.calls <- final.calls2
        good.bars.sub <- good.bars[!(good.bars %in% c("Bar12", "Bar13", "Bar24", "Bar30"))]
    }
} else {
    good.bars.sub <- good.bars
}

full_outbase <- paste0(opts$inputdir, opts$experiment, "/")
if (!dir.exists(full_outbase)) {
    dir.create(full_outbase, recursive = TRUE)
}

saveRDS(final.calls, paste0(full_outbase, opts$sample, "_finalCalls_noreclass.rds"))
saveRDS(final.calls, paste0(full_outbase, opts$sample, "_finalCalls.rds")) # this will be overwritten if reclassification works

bar.table.doublets <- bar.table.full[names(final.calls)[which(final.calls=="Doublet")],io$bars_used[[opts$experiment]][[opts$sample]]]
bar.table.doublets <- data.table(bar.table.doublets)
rownames(bar.table.doublets) <- rownames(bar.table.doublets)

argmax <- max.col(bar.table.doublets)
bar.table.doublets.sorted <- bar.table.doublets[sort(argmax, index.return=TRUE)$ix,]
argmax2 <- max.col(replace(as.matrix(bar.table.doublets.sorted), cbind(1:dim(bar.table.doublets.sorted)[1], max.col(as.matrix(bar.table.doublets.sorted))), -Inf))

sorted_argmax <- max.col(bar.table.doublets.sorted)
for (i in unique(argmax)) {
    bar.table.doublets.sorted[sorted_argmax==i,] = bar.table.doublets.sorted[sorted_argmax==i,][sort(argmax2[sorted_argmax==i], index.return=TRUE)$ix,]
}

full_outbase <- paste0(opts$inputdir, "processing/DoubletHeatmaps/", opts$experiment, "/")
if (!dir.exists(full_outbase)) {
    dir.create(full_outbase, recursive = TRUE)
}

p <- ggplot(melt(as.matrix(bar.table.doublets.sorted)), aes(Var2, Var1, fill= value)) + 
  geom_tile() +
  #scale_fill_gradient( trans = 'log2' ) +
  #scale_fill_gradient(low="white", high="#2171b5", trans = 'log2')
  scale_fill_continuous(type = "viridis", trans = 'log2') +
  theme_classic()
ggsave(paste0(full_outbase, opts$sample, "_doublet_hm.", opts$plot.type), plot = p, device=opts$plot.type, width = 10, height = 10, units = "in")

message('experiment=', opts$experiment, ', sample=', opts$sample, ', good.bars.sub=', paste(good.bars.sub, collapse=', '))

## Perform semi-supervised negative cell reclassification
if (!error1) {
    
    error2 <- FALSE
    tryCatch(reclass.cells <- findReclassCells(bar.table.full[,good.bars.sub], names(final.calls)[which(final.calls=="Negative")]),
             error = function(e) { error2 <<- TRUE}
            )
    
    if (!error2) {
    
        reclass.cells <- findReclassCells(bar.table.full[,good.bars.sub], names(final.calls)[which(final.calls=="Negative")])
        error3 <- FALSE
        tryCatch(reclass.res <- rescueCells(bar.table.full[,good.bars.sub], final.calls, reclass.cells),
                 error = function(e) { error3 <<- TRUE}
                )
        if (!error3) {

            reclass.res <- rescueCells(bar.table.full[,good.bars.sub], final.calls, reclass.cells)

            full_outbase <- paste0(opts$inputdir, "processing/NegReclass/", opts$experiment, "/")
            if (!dir.exists(full_outbase)) {
                dir.create(full_outbase, recursive = TRUE)
            }

            ## Visualize Results
            p <- ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
                geom_point() + xlim(c(nrow(reclass.res)-1,1)) + 
                ylim(c(0,1.05)) +
                geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
                geom_vline(xintercept = io$reclass_stability[[opts$experiment]][[opts$sample]], color="red") +
                geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
                geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
                geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
                theme_classic()
            ggsave(paste0(full_outbase, opts$sample, "_negreclass.", opts$plot.type), plot = p, device=opts$plot.type, width = 5, height = 5, units = "in")

            ## Finalize negative cell rescue results
            final.calls.rescued <- final.calls
            rescue.ind <- which(reclass.cells$ClassStability >= io$reclass_stability[[opts$experiment]][[opts$sample]]) ## Note: Value will be dataset-specific
            final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

            full_outbase <- paste0(opts$inputdir, opts$experiment, "/")
            if (!dir.exists(full_outbase)) {
                dir.create(full_outbase, recursive = TRUE)
            }

            saveRDS(final.calls.rescued, paste0(full_outbase, opts$sample, "_finalCalls.rds"))

        }
    }
}












