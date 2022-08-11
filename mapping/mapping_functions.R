doBatchCorrect <- function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, pc_override = NULL, BPPARAM = SerialParam()){
  require(BiocParallel)
  
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    #message('running PCA')
    pca = irlba::prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  #message('making pc_list')
  pc_list    <- lapply(unique(timepoints), function(tp){
    #message(tp)
    sub_pc   <- pca[timepoints == tp, , drop = FALSE]
    sub_samp <- samples[timepoints == tp]
    list     <- lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) <- unique(sub_samp)
    return(list)
  })
  
  names(pc_list) <- unique(timepoints)
  
  #arrange to match timepoint order
  #message('arranging pc_list')
  pc_list <- pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list <- lapply(pc_list, function(x){
    #message(names(x))
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  #message('correcting within pc_list elements')
  correct_list <- lapply(pc_list, function(x){
    #message(names(x))
    if(length(x) > 1){
        #return(do.call(scran::fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
        return(do.call(reducedMNN, c(x, BPPARAM = BPPARAM))$corrected) # edited 09.02.2020 because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02.2020 changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
      #message('running full correct')
      #correct <- do.call(scran::fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
      correct <- do.call(reducedMNN, c(correct_list, BPPARAM = BPPARAM))$corrected # edited 09.02.2020 because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02.2020 changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
  } else {
    correct <- correct_list[[1]]
  }
  
  #message('outputting')
  correct <- correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}

getHVGs <- function(sce, block, min.mean = 1e-3){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getcelltypes <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getMappingScore <- function(mapping){
    out <- list()
    celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped),
                                 nrow=length(mapping$celltypes.mapped[[1]]),
                                 ncol=length(mapping$celltypes.mapped))
    cellstages_accrossK <- matrix(unlist(mapping$cellstages.mapped),
                                  nrow=length(mapping$cellstages.mapped[[1]]),
                                  ncol=length(mapping$cellstages.mapped))
    out$celltype.score <- NULL
    for (i in 1:nrow(celltypes_accrossK)){
        p <- max(table(celltypes_accrossK[i,]))
        index <- which(table(celltypes_accrossK[i,]) == p)
        p <- p/length(mapping$celltypes.mapped)
        out$celltype.score <- c(out$celltype.score,p)
    }
    out$cellstage.score <- NULL
    for (i in 1:nrow(cellstages_accrossK)){
        p <- max(table(cellstages_accrossK[i,]))
        index <- which(table(cellstages_accrossK[i,]) == p)
        p <- p/length(mapping$cellstages.mapped)
        out$cellstage.score <- c(out$cellstage.score,p)
    }
    return(out)  
}

get_meta.original <- function(correct_atlas, atlas_meta, correct_map, map_meta, k_map = 30){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
    get.distance = FALSE)
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype_old[match(x, atlas_meta$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         stage.mapped     = stage.mapped[x],
         celltypes.mapped = celltypes[x,],
         stages.mapped    = stages[x,])
  })
  names(out) <- map_meta$cell
  return(out)  
}

get_meta.extended <- function(correct_atlas, atlas_meta, correct_map, map_meta, k_map = 30){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
    get.distance = FALSE)
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype.clustering[match(x, atlas_meta$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage.mapped[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         stage.mapped     = stage.mapped[x],
         celltypes.mapped = celltypes[x,],
         stages.mapped    = stages[x,])
  })
  names(out) <- map_meta$cell
  return(out)  
}

get_meta.descendant <- function(correct_atlas, atlas_meta, correct_map, map_meta, k_map = 30){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
    get.distance = FALSE)
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype.descendant[match(x, atlas_meta$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped.descendant     = k.mapped[x,],
         celltype.mapped.descendant  = celltype.mapped[x],
         stage.mapped.descendant     = stage.mapped[x],
         celltypes.mapped.descendant = celltypes[x,],
         stages.mapped.descendant    = stages[x,])
  })
  names(out) <- map_meta$cell
  return(out)  
}

get_meta.meso <- function(correct_atlas, atlas_meta, correct_map, map_meta, k_map = 30){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
    get.distance = FALSE)
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped.meso     = k.mapped[x,],
         celltype.mapped.meso  = celltype.mapped[x],
         stage.mapped.meso     = stage.mapped[x],
         celltypes.mapped.meso = celltypes[x,],
         stages.mapped.meso    = stages[x,])
  })
  names(out) <- map_meta$cell
  return(out)  
}
                           
getMappingScore.descendant <- function(mapping){
    out <- list()
    celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped.descendant),
                                 nrow=length(mapping$celltypes.mapped.descendant[[1]]),
                                 ncol=length(mapping$celltypes.mapped.descendant))
    cellstages_accrossK <- matrix(unlist(mapping$cellstages.mapped.descendant),
                                  nrow=length(mapping$cellstages.mapped.descendant[[1]]),
                                  ncol=length(mapping$cellstages.mapped.descendant))
    out$celltype.score.descendant <- NULL
    for (i in 1:nrow(celltypes_accrossK)){
        p <- max(table(celltypes_accrossK[i,]))
        index <- which(table(celltypes_accrossK[i,]) == p)
        p <- p/length(mapping$celltypes.mapped.descendant)
        out$celltype.score.descendant <- c(out$celltype.score.descendant,p)
    }
    out$cellstage.score.descendant <- NULL
    for (i in 1:nrow(cellstages_accrossK)){
        p <- max(table(cellstages_accrossK[i,]))
        index <- which(table(cellstages_accrossK[i,]) == p)
        p <- p/length(mapping$cellstages.mapped.descendant)
        out$cellstage.score.descendant <- c(out$cellstage.score.descendant,p)
    }
    return(out)  
}
                           
getMappingScore.meso <- function(mapping){
    out <- list()
    celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped.meso),
                                 nrow=length(mapping$celltypes.mapped.meso[[1]]),
                                 ncol=length(mapping$celltypes.mapped.meso))
    cellstages_accrossK <- matrix(unlist(mapping$cellstages.mapped.meso),
                                  nrow=length(mapping$cellstages.mapped.meso[[1]]),
                                  ncol=length(mapping$cellstages.mapped.meso))
    out$celltype.score.meso <- NULL
    for (i in 1:nrow(celltypes_accrossK)){
        p <- max(table(celltypes_accrossK[i,]))
        index <- which(table(celltypes_accrossK[i,]) == p)
        p <- p/length(mapping$celltypes.mapped.meso)
        out$celltype.score.meso <- c(out$celltype.score.meso,p)
    }
    out$cellstage.score.meso <- NULL
    for (i in 1:nrow(cellstages_accrossK)){
        p <- max(table(cellstages_accrossK[i,]))
        index <- which(table(cellstages_accrossK[i,]) == p)
        p <- p/length(mapping$cellstages.mapped.meso)
        out$cellstage.score.meso <- c(out$cellstage.score.meso,p)
    }
    return(out)  
}

mapWrap <- function(atlas_sce, atlas_meta, map_sce, map_meta, tp_order, k = 30, npcs = 50, return.list = FALSE){
    
  message("Normalizing joint dataset...")
  
  #easier to avoid directly binding sce objects as it is a lot more likely to have issues
  sce_all <- SingleCellExperiment::SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(atlas_sce),counts(map_sce)),sparse=TRUE)))
  #big_sce <- scater::normalize(sce_all)
  #big_sce <- scater::logNormCounts(sce_all) # edited 09.02.2020 because normalize deprecated in favour of logNormCounts
  big_sce <- multiBatchNorm(sce_all, batch=c(atlas_meta$sample, map_meta$batch)) # edited 17.02.2020 because now multibatchnorm exists
  message("Done\n")
  
  message("Computing highly variable genes...")
  hvgs    <- getHVGs(big_sce, block=c(atlas_meta$sample, map_meta$batch))
  message("Done\n")
  
  message("Performing PCA...")
  full_pca <- multiBatchPCA(big_sce,
                            batch=c(atlas_meta$sample, map_meta$batch),
                            subset.row = hvgs,
                            d = npcs,
                            preserve.single = TRUE,
                            assay.type = "logcounts")
  big_pca <- full_pca[[1]]
  rownames(big_pca) <- colnames(big_sce) 
  atlas_pca <- big_pca[1:ncol(atlas_sce),]
  map_pca   <- big_pca[-(1:ncol(atlas_sce)),]
  message("Done\n")
  
  message("Batch effect correction for the atlas...")  
  order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage  <- factor(order_df$stage, 
                            levels = rev(c("E9.5","E9.25","E9.0","E8.75","E8.5","E8.25","E8.0","E7.75","E7.5","E7.25","mixed_gastrulation","E7.0","E6.75","E6.5")))
  order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage <- as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected <- doBatchCorrect(counts         = logcounts(atlas_sce[hvgs,]), 
                                    timepoints      = atlas_meta$stage, 
                                    samples         = atlas_meta$sample, 
                                    timepoint_order = order_df$stage, 
                                    sample_order    = order_df$sample, 
                                    pc_override     = atlas_pca,
                                    npc             = npcs)
  message("Done\n")
  
  message("MNN mapping...")
  order_df        <- map_meta[!duplicated(map_meta$batch), c("timepoint", "batch")]
  order_df$ncells <- sapply(order_df$batch, function(x) sum(map_meta$batch == x))
  order_df$timepoint  <- factor(order_df$timepoint, 
                                levels = rev(tp_order)
                               )
  order_df       <- order_df[order(order_df$timepoint, order_df$ncells, decreasing = TRUE),]
  order_df$timepoint <- as.character(order_df$timepoint)
  message(order_df$batch)
  correct <- reducedMNN(rbind(atlas_corrected, map_pca),
                      batch=c(rep("ATLAS", dim(atlas_meta)[1]), map_meta$batch),
                      merge.order=c("ATLAS", order_df$batch))$corrected
  atlas   <- 1:nrow(atlas_pca)
  correct_atlas <- correct[atlas,]
  correct_map   <- correct[-atlas,]
  message("Done\n")
  
  message("Writing output...") 
  out <- list()
  out$correct <- correct
  out$correct_atlas <- correct_atlas
  out$correct_map <- correct_map
  out$pca <- full_pca
  message("Done\n")
  
  return(out)
  
}

mapNoAtlas <- function(sce, meta, order, k = 30, npcs = 50, return.list = FALSE){
    
  message("Normalizing joint dataset...")
  sce_all <- SingleCellExperiment::SingleCellExperiment(list(counts=counts(sce)))
  big_sce <- multiBatchNorm(sce_all, batch=c(meta$batch))
  message("Done\n")
  
  message("Computing highly variable genes...")
  hvgs    <- getHVGs(big_sce, block=c(meta$batch))
  message("Done\n")
  
  message("Performing PCA...")
  full_pca <- multiBatchPCA(big_sce,
                           batch=c(meta$batch),
                           subset.row = hvgs,
                           d = npcs,
                           preserve.single = TRUE,
                           get.variance=TRUE,
                           assay.type = "logcounts")
  big_pca <- full_pca[[1]]
  rownames(big_pca) <- colnames(big_sce)
  message("Done\n")
  
  message("Batch effect correction...")  
  order_df        <- meta[!duplicated(meta$batch), c("timepoint", "batch")]
  order_df$ncells <- sapply(order_df$batch, function(x) sum(meta$batch == x))
  order_df$timepoint  <- factor(order_df$timepoint, 
                                levels = rev(order)
                               )
  order_df       <- order_df[order(order_df$timepoint, order_df$ncells, decreasing = TRUE),]
  order_df$timepoint <- as.character(order_df$timepoint)
  message(order_df)
  
  set.seed(42)
  corrected <- doBatchCorrect(counts          = logcounts(big_sce[hvgs,]), 
                              timepoints      = meta$timepoint, 
                              samples         = meta$batch, 
                              timepoint_order = order_df$timepoint, 
                              sample_order    = order_df$batch, 
                              pc_override     = big_pca,
                              npc             = npcs)
  message("Done\n")
  out <- list()
  out$corrected <- corrected
  out$pca <- full_pca
  message("Done\n")
  
  return(out)
  
}