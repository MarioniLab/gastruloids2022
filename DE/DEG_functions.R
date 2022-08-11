suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tools))

do_DE <- function(sce, gene_metadata, min_detection_rate_per_group=0.1, max_detection_rate_per_group=1, min.logFC=2, threshold_fdr=0.1) {
    
    # input should already have groups
    groups <- unique(sce$group)
    
    # calculate detection rate per gene
    cdr.dt <- data.table(
      ens_id = rownames(sce),
      detection_rate_A = rowMeans(counts(sce[,sce$group==groups[1]])>0),
      detection_rate_B = rowMeans(counts(sce[,sce$group==groups[2]])>0)
    ) %>% setnames(c("ens_id",sprintf("detection_rate_%s",groups[1]),sprintf("detection_rate_%s",groups[2])))
    # .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",io$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",io$groups[2])),with=F][[1]])] %>%

    # Filter genes
    sce <- sce[rownames(sce)%in%gene_metadata$ens_id,]

    #sce <- sce[rownames(metadata(mapping$pca)$rotation),]
    
    # Filter genes by min detection rate per group
    cdr_A <- rowMeans(counts(sce[,sce$group==groups[1]])>0) >= min_detection_rate_per_group
    cdr_B <- rowMeans(counts(sce[,sce$group==groups[2]])>0) >= min_detection_rate_per_group
    sce <- sce[cdr_B | cdr_A,]

    # Filter genes by max detection rate per group
    cdr_A <- rowMeans(counts(sce[,sce$group==groups[1]])>0) < max_detection_rate_per_group
    cdr_B <- rowMeans(counts(sce[,sce$group==groups[2]])>0) < max_detection_rate_per_group
    sce <- sce[cdr_B | cdr_A,]
    
    # Convert SCE to DGEList
    sce_edger <- scran::convertTo(sce, type="edgeR")

    # Define design matrix (with intercept)
    cdr <- colMeans(counts(sce)>0)
    if (length(unique(sce$batch)) > 1) {
        design <- model.matrix(~ cdr + sce$batch + sce$group)
    } else {
        design <- model.matrix(~ cdr + sce$group)
    }

    # Estimate dispersions
    sce_edger  <- estimateDisp(sce_edger,design)

    # Fit GLM
    #fit <- glmQLFit(sce_edger, design, coef=grep('sce\\$group', colnames(design), value=TRUE))
    fit <- glmQLFit(sce_edger, design, coef=paste0('sce$group', groups[2]))

    # Likelihood ratio test
    #lrt <- glmQLFTest(fit, coef=grep('sce\\$group', colnames(design), value=TRUE))
    lrt <- glmQLFTest(fit, coef=paste0('sce$group', groups[2]))
    
    # Construct output data.frame
    out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
        setnames(c("ens_id","logFC","logCPM","LR","p.value","padj_fdr")) %>%
        .[,c("logCPM","LR"):=NULL]
    out %>% .[,log_padj_fdr:= -log10(padj_fdr)]
    out[,c("groupA_N","groupB_N"):=list(sum(sce$group == groups[1]),sum(sce$group == groups[2]))]
    setnames(out, c("groupA_N","groupB_N"),c(sprintf("N_%s",groups[1]),sprintf("N_%s",groups[2])))
    out <- merge(cdr.dt, out, all.y=TRUE, by="ens_id") %>% # Add gene statistics
        merge(gene_metadata, all.y=TRUE, by="ens_id")
    out[, sig := (padj_fdr<=threshold_fdr & abs(logFC)>=min.logFC)] # Calculate statistical significance
    setorder(out, -sig, padj_fdr, na.last=T)
    
    return(out)
    
}

runDEG <- function(io) {
    
    # io should be a list containing:
    # io$DEG_met
    # io$DEG_srat
    # io$genemetadata
    # io$test
    # io$celltype # Define cell types
    # io$groups) # Define groups
    # io$threshold_fdr # Define FDR threshold
    # io$min.logFC # Define minimum logFC for significance
    # io$min_detection_rate_per_group # For a given gene, the minimum fraction of cells that must express it in at least one group
    # io$max_detection_rate_per_group # For a given gene, the maximum fraction of cells that can express it in one group
    
    
    ###############
    ## Load data ##
    ###############

    stopifnot(all(io$groups%in%unique(io$DEG_meta$group1)))

    # Update cell metadata
    sample_metadata <- io$DEG_meta %>%
      .[group1%in%io$groups & celltype.mapped%in%io$celltype] %>%
      setnames("group1","group") %>%
      .[,c("cell","celltype.mapped","group")]

    # Sort cells so that groupA comes before groupB
    sample_metadata[,group:=factor(group,levels=io$groups)] %>% setorder(group)
    table(sample_metadata$group)

    # Load SingleCellExperiment object
    sce <- SingleCellExperiment(list(counts=io$DEG_srat@assays$RNA@counts[,sample_metadata$cell],
                                     logcounts=io$DEG_srat@assays$RNA@data[,sample_metadata$cell]),
                                colData=sample_metadata
                                )
    sce$group <- sample_metadata$group

    # Load gene metadata
    gene_metadata <- io$genemetadata %>%
      .[ens_id%in%rownames(sce)] %>%
      .[,c("symbol","ens_id")] %>% 
      setnames("symbol","gene")

    ################
    ## Parse data ##
    ################

    # calculate detection rate per gene
    cdr.dt <- data.table(
      ens_id = rownames(sce),
      detection_rate_A = rowMeans(logcounts(sce[,sce$group==io$groups[1]])>0),
      detection_rate_B = rowMeans(logcounts(sce[,sce$group==io$groups[2]])>0)
    ) %>% setnames(c("ens_id",sprintf("detection_rate_%s",io$groups[1]),sprintf("detection_rate_%s",io$groups[2])))
    # .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",io$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",io$groups[2])),with=F][[1]])] %>%

    # Filter genes
    sce <- sce[rownames(sce)%in%gene_metadata$ens_id,]

    ################################################
    ## Differential expression testing with edgeR ##
    ################################################

    out <- doDiffExpr(sce, io$groups, io$test, io$min_detection_rate_per_group, io$max_detection_rate_per_group)
    
    out <- out %>% .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% # Add sample statistics
        setnames(c("groupA_N","groupB_N"),c(sprintf("N_%s",io$groups[1]),sprintf("N_%s",io$groups[2])))
    out <- merge(cdr.dt, out, all.y=TRUE, by="ens_id") %>% # Add gene statistics
        merge(gene_metadata, all.y=TRUE, by="ens_id")
    out <- out %>%
        .[, sig := (padj_fdr<=io$threshold_fdr & abs(logFC)>=io$min.logFC)] %>% # Calculate statistical significance
        setorder(-sig, padj_fdr, na.last=T)
    
    return(out)
}

# Function to differential expression
# - sce: SingleCellExperiment object with the column "group" in the colData
# - groups: the names of the two groups
# - test: one of "edgeR","t-test","wilcoxon".
# - min_detection_rate_per_group: minimum detection rate per group
doDiffExpr <- function(sce, groups, test=c("edgeR","t-test","wilcoxon"), min_detection_rate_per_group = 0.50, max_detection_rate_per_group=1) {
    
  # Sanity checks
  if (!is(sce, "SingleCellExperiment")) stop("'sce' has to be an instance of SingleCellExperiment")
  stopifnot(length(groups)==2)
  test <- match.arg(test)

  # Filter genes by min detection rate per group
  cdr_A <- rowMeans(logcounts(sce[,sce$group==groups[1]])>0) >= min_detection_rate_per_group
  cdr_B <- rowMeans(logcounts(sce[,sce$group==groups[2]])>0) >= min_detection_rate_per_group
  sce <- sce[cdr_B | cdr_A,]

  # Filter genes by max detection rate per group
  cdr_A <- rowMeans(logcounts(sce[,sce$group==groups[1]])>0) < max_detection_rate_per_group
  cdr_B <- rowMeans(logcounts(sce[,sce$group==groups[2]])>0) < max_detection_rate_per_group
  sce <- sce[cdr_B | cdr_A,]
  
  if (test=="edgeR") {
    print("doing edgeR")
    out <- .edgeR(sce)
  } else if (test=="t-test") {
    out <- .t_test(sce)
  } else if (test=="wilcoxon") {
    out <- .wilcoxon(sce)
  } else {
    stop("Test not recognised")
  }
  
  out %>% .[,log_padj_fdr:= -log10(padj_fdr)]
  
  return(out)
}



.t_test <- function(sce) {
  # left.result1 <- wilcox.test(host.vals, target.vals, alternative="less", mu=-lfc, exact=FALSE)
  # expect_equal(p.adjust(pval, method="BH"), curres$FDR)
}

.edgeR <- function(sce) {
  
  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce, type="edgeR")
  
  # Define design matrix (with intercept)
  cdr <- colMeans(logcounts(sce)>0)
  design <- model.matrix(~cdr+sce$group)
  
  # Estimate dispersions
  sce_edger  <- estimateDisp(sce_edger,design)
  
  # Fit GLM
  fit <- glmQLFit(sce_edger,design)
  
  # Likelihood ratio test
  lrt <- glmQLFTest(fit)
  
  # Construct output data.frame
  out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
    setnames(c("ens_id","logFC","logCPM","LR","p.value","padj_fdr")) %>%
    .[,c("logCPM","LR"):=NULL]
  
  return(out)
}

################
## Plot utils ##
################

gg_volcano_plot <- function(tmp, top_genes=10, xlim=NULL, ylim=NULL, gene_list=NULL) {
  negative_hits <- tmp[sig==TRUE & logFC<0,ens_id]
  positive_hits <- tmp[sig==TRUE & logFC>0,ens_id]
  all <- nrow(tmp)
  all_red <- FALSE
  
  if (is.null(gene_list)) {
    all_red <- TRUE
    setorder(tmp, -sig, p.value, na.last=T)
    gene_list <- head(tmp[sig==T],n=top_genes)$gene
  }
  if (is.null(xlim)) {
    xlim <- max(abs(tmp$logFC), na.rm=T)
  }
  if (is.null(ylim)) {
    ylim <- max(-log10(tmp$p.value), na.rm=T)
  }
  
  tmp <- tmp[order(-match(gene,gene_list), na.last=FALSE)]
  
  p <- ggplot(tmp, aes(x=logFC, y=-log10(p.value))) +
    labs(title="", x="Log Fold Change", y=expression(paste("-log"[10],"(p.value)"))) +
    # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange")
    #ggrastr::geom_point_rast(aes(color=sig), size=1) +
  if (all_red) {
      p <- p + geom_point(aes(color=sig), size=1)
  } else {
      p <- p + geom_point(aes(color=gene %in% gene_list)) +
          scale_size_manual(values = c(0.1,1))
  }
  p <- p + scale_color_manual(values=c("black","red")) +
    scale_x_continuous(limits=c(-xlim-2,xlim+2)) +
    scale_y_continuous(limits=c(0,ylim+1)) +
    annotate("text", x=0, y=ylim+1, size=7, label=sprintf("(%d)", all)) +
    annotate("text", x=-8, y=ylim+1, size=7, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=8, y=ylim+1, size=7, label=sprintf("%d (+)",length(positive_hits))) +
    ggrepel::geom_text_repel(data=tmp[gene %in% gene_list], aes(x=logFC, y=-log10(p.value), label=gene), size=5) +
    theme_bw() +
    theme(
      plot.title=element_text(size=rel(1.5), face='bold', margin=margin(0,0,10,0), hjust=0.5),
      axis.text=element_text(size=rel(1.75), color='black'),
      axis.title=element_text(size=rel(1.95), color='black'),
      axis.title.y = element_text(margin=margin(0,10,0,0)),
      axis.title.x = element_text(margin=margin(10,0,0,0)),
      legend.position="none",
      # panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # panel.background = element_blank()
    )
  return(p)
}


matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}