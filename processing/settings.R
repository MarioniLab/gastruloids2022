io <- list()

io$s.genes.ensid <- '/nfs/research/marioni/Leah/data/cell_cycle_genes/s_genes_ENS.txt'
io$g2m.genes.ensid <- '/nfs/research/marioni/Leah/data/cell_cycle_genes/g2m_genes_ENS.txt'
io$s.genes.symbol <- '/nfs/research/marioni/Leah/data/cell_cycle_genes/s_genes.txt'
io$g2m.genes.symbol <- '/nfs/research/marioni/Leah/data/cell_cycle_genes/g2m_genes.txt'

io$k <- 30

io$barcode_files[["MULTI/exp1_d5"]] <- c(
  "sample1" = "lane2_ATCACGAT_CS1_10x_BC_merged_L002",
  "sample2" = "lane2_CGATGTAT_CS2_10x_BC_merged_L002",
  "sample3" = "lane2_TTAGGCAT_CS3_10x_BC_merged_L002"
)
io$barcode_files[["MULTI/exp2A_d3_d3.5"]] <- c(
  "sample1" = "lane1_ATCACGAT_BC1_gastr_d3_A1_merged_L001",
  "sample2" = "lane1_CGATGTAT_BC2_gastr_d3_A2_merged_L001",
  "sample3" = "lane1_TTAGGCAT_BC3_gastr_d3_A3_merged_L001"
)
io$barcode_files[["MULTI/exp2C_d3_d3.5"]] <- c(
  "sample1" = "lane1_TGACCAAT_BC4_gastr_d3_C1_merged_L001",
  "sample2" = "",
  "sample3" = ""
)
io$barcode_files[["MULTI/exp2D_d3_d3.5"]] <- c(
  "sample1" = "",
  "sample2" = ""
)
io$barcode_files[["MULTI/exp3A_d4_d4.5"]] <- c(
  "sample1" = "",
  "sample2" = "",
  "sample3" = ""
)
io$barcode_files[["MULTI/exp3B_d4_d4.5"]] <- c(
  "sample1" = "",
  "sample2" = "lane1_ACAGTGAT_BC13_gastr_d4_B2_merged_L001",
  "sample3" = "lane1_GCCAATAT_BC14_gastr_d4_B3_merged_L001"
)
io$barcode_files[["MULTI/exp3C_d4_d4.5"]] <- c(
  "sample1" = "lane1_ACTTGAAT_BC16_gastr_d4_C2_merged_L001",
  "sample2" = "lane1_CAGATCAT_BC15_gastr_d4_C1_merged_L001"
)
io$barcode_files[["MULTI/exp4_d3"]] <- c(
  #"sampleA" = "SITTA12_lipidtag_gastr_d3_MULTIseq_3prime_RNA_A_L002",
  #"sampleA" = "SLX20586_SIGAA11",
  #"sampleA" = "SIGAA11_gastr_d3_MULTIseq_3prime_RNA_A",
  "sampleA" = "SIGAA11_gastr_d3_MULTIseq_3primeRNA_A_L001",
  #"sampleB" = "SITTB12_lipidtag_gastr_d3_MULTIseq_3prime_RNA_B_L002"
  #"sampleB" = "SIGAB11_gastr_d3_MULTIseq_3prime_RNA_B"
  "sampleB" = "SIGAB11_gastr_d3_MULTIseq_3primeRNA_B_L001"
)
io$barcode_files[["MULTI/exp5_d3.5_d4"]] <- c(
  #"sampleA" = "SITTC12_lipidtag_gastr_d3_5_d4_MULTIseq_3prime_RNA_A_L002",
  #"sampleA" = "SLX20586_SIGAC11",
  "sampleA" = "SIGAC11_gastr_d3.5_d4_MULTIseq_3prime_RNA_A",
  #"sampleB" = "SITTD12_lipidtag_gastr_d3_5_d4_MULTIseq_3prime_RNA_B_L002",
  "sampleB" = "SIGAD11_gastr_d3.5_d4_MULTIseq_3prime_RNA_B",
  #"sampleC" = "SITTE12_lipidtag_gastr_d3_5_d4_MULTIseq_3prime_RNA_C_L002",
  "sampleC" = "SIGAE11_gastr_d3.5_d4_MULTIseq_3prime_RNA_C",
  #"sampleD" = "SITTF12_lipidtag_gastr_d3_5_d4_MULTIseq_3prime_RNA_D_L002"
  "sampleD" = "SIGAF11_gastr_d3.5_d4_MULTIseq_3prime_RNA_D"
)
io$barcode_files[["MULTI/exp6_d4.5_d5"]] <- c(
  #"sampleA" = "SITTE10_lipidtag_gastr_d4_5_d5_MULTIseq_3prime_RNA_A_L002",
  "sampleA" = "SIGAA10_gastr_d4.5_d5_MULTIseq_3prime_RNA_A",
  #"sampleB" = "SITTF10_lipidtag_gastr_d4_5_d5_MULTIseq_3prime_RNA_B_L002",
  "sampleB" = "SIGAB10_gastr_d4.5_d5_MULTIseq_3prime_RNA_B",
  #"sampleC" = "SITTG10_lipidtag_gastr_d4_5_d5_MULTIseq_3prime_RNA_C_L002"
  "sampleC" = "SIGAC10_gastr_d4.5_d5_MULTIseq_3prime_RNA_C"
)

io$cell_pos=c(1,16)
io$umi_pos=c(17,28)
io$tag_pos=c(1,8)

io$bars_used <- list()
io$bars_used[["MULTI/exp1_d5"]] <- list(
    sample1 = c(1:24),
    sample2 = c(1:24),
    sample3 = c(1:24)
)
io$bars_used[["MULTI/exp2A_d3_d3.5"]] <- list(
    sample1 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample2 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample3 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)
io$bars_used[["MULTI/exp2C_d3_d3.5"]] <- list(
    sample1 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample2 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample3 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)
io$bars_used[["MULTI/exp2D_d3_d3.5"]] <- list(
    sample1 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample2 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)
io$bars_used[["MULTI/exp3A_d4_d4.5"]] <- list(
    sample1 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample2 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample3 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)
io$bars_used[["MULTI/exp3B_d4_d4.5"]] <- list(
    sample1 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample2 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sample3 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)
io$bars_used[["MULTI/exp3C_d4_d4.5"]] <- list(
    #sample1 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47),
    sample1 = c(2, 5, 6, 12, 13, 14, 15, 19),
    #sample2 = c(2, 5, 6, 12, 13, 14, 15, 19, 22, 25, 28, 30, 39, 41, 45, 47)
    sample2 = c(2, 5, 6, 12, 13, 14, 15, 19)
)
io$bars_used[["MULTI/exp4_d3"]] <- list(
    sampleA = c(2, 6, 12, 13, 14, 15, 18, 19, 20, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sampleB = c(2, 6, 12, 13, 14, 15, 18, 19, 20, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)
io$bars_used[["MULTI/exp5_d3.5_d4"]] <- list(
    sampleA = c(1, 2, 4, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sampleB = c(1, 2, 4, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sampleC = c(1, 2, 4, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sampleD = c(1, 2, 4, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)
io$bars_used[["MULTI/exp6_d4.5_d5"]] <- list(
    sampleA = c(1, 2, 4, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sampleB = c(1, 2, 4, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66),
    sampleC = c(1, 2, 4, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 30, 41, 45, 47, 49, 52, 55, 56, 59, 61, 65, 66)
)

io$laterbars <- list()
io$laterbars[["exp2A_d3_d3.5"]] <- c('Bar13', 'Bar14', 'Bar15', 'Bar19', 'Bar39', 'Bar41', 'Bar45', 'Bar47', 'Bar59', 'Bar61', 'Bar65', 'Bar66')
io$laterbars[["exp2C_d3_d3.5"]] <- c('Bar13', 'Bar14', 'Bar15', 'Bar19', 'Bar39', 'Bar41', 'Bar45', 'Bar47', 'Bar59', 'Bar61', 'Bar65', 'Bar66')
io$laterbars[["exp3B_d4_d4.5"]] <- c('Bar13', 'Bar14', 'Bar15', 'Bar19', 'Bar39', 'Bar41', 'Bar45', 'Bar47', 'Bar59', 'Bar61', 'Bar65', 'Bar66')
io$laterbars[["exp3C_d4_d4.5"]] <- c('Bar13', 'Bar14', 'Bar15', 'Bar19', 'Bar39', 'Bar41', 'Bar45', 'Bar47')
io$laterbars[["exp5_d3.5_d4"]] <- c('Bar1', 'Bar2', 'Bar4', 'Bar6', 'Bar7', 'Bar10', 'Bar12', 'Bar13', 'Bar14', 'Bar15', 'Bar17', 'Bar18', 'Bar19', 'Bar20', 'Bar21', 'Bar22')
io$laterbars[["exp6_d4.5_d5"]] <- c('Bar1', 'Bar2', 'Bar4', 'Bar6', 'Bar7', 'Bar10', 'Bar12', 'Bar13', 'Bar14', 'Bar15', 'Bar17', 'Bar18', 'Bar19', 'Bar20', 'Bar21', 'Bar22')
io$earlierbars <- list()
io$earlierbars[["exp2A_d3_d3.5"]] <- c('Bar2', 'Bar5', 'Bar6', 'Bar12', 'Bar22', 'Bar25', 'Bar28', 'Bar30', 'Bar49', 'Bar52', 'Bar55', 'Bar56')
io$earlierbars[["exp2C_d3_d3.5"]] <- c('Bar2', 'Bar5', 'Bar6', 'Bar12', 'Bar22', 'Bar25', 'Bar28', 'Bar30', 'Bar49', 'Bar52', 'Bar55', 'Bar56')
io$earlierbars[["exp3B_d4_d4.5"]] <- c('Bar2', 'Bar5', 'Bar6', 'Bar12', 'Bar22', 'Bar25', 'Bar28', 'Bar30', 'Bar49', 'Bar52', 'Bar55', 'Bar56')
io$earlierbars[["exp3C_d4_d4.5"]] <- c('Bar2', 'Bar5', 'Bar6', 'Bar12', 'Bar22', 'Bar25', 'Bar28', 'Bar30')
io$earlierbars[["exp5_d3.5_d4"]] <- c('Bar23', 'Bar24', 'Bar25', 'Bar28', 'Bar30', 'Bar41', 'Bar45', 'Bar47', 'Bar49', 'Bar52', 'Bar55', 'Bar56', 'Bar59', 'Bar61', 'Bar65', 'Bar66')
io$earlierbars[["exp6_d4.5_d5"]] <- c('Bar23', 'Bar24', 'Bar25', 'Bar28', 'Bar30', 'Bar41', 'Bar45', 'Bar47', 'Bar49', 'Bar52', 'Bar55', 'Bar56', 'Bar59', 'Bar61', 'Bar65', 'Bar66')
io$mixedexps <- names(io$laterbars)

io$reclass_stability <- list()
io$reclass_stability[["MULTI/exp1_d5"]] <- c(
    sample1 = NA,
    sample2 = NA,
    sample3 = NA
)
#io$reclass_stability[["MULTI/exp2A_d3_d3.5"]] <- c(
#    sample1 = 0,
#    sample2 = 0,
#    sample3 = 0
#)
#io$reclass_stability[["MULTI/exp2C_d3_d3.5"]] <- c(
#    sample1 = 0,
#    sample2 = 0,
#    sample3 = 0
#)
#io$reclass_stability[["MULTI/exp2D_d3_d3.5"]] <- c(
#    sample1 = 0,
#    sample2 = 0
#)
#io$reclass_stability[["MULTI/exp3A_d4_d4.5"]] <- c(
#    sample1 = 0,
#    sample2 = 0,
#    sample3 = 0
#)
io$reclass_stability[["MULTI/exp3B_d4_d4.5"]] <- c(
    sample1 = NA,
    sample2 = 3,
    sample3 = 6
)
#io$reclass_stability[["MULTI/exp3C_d4_d4.5"]] <- c(
#    sample1 = 0,
#    sample2 = 0
#)
io$reclass_stability[["MULTI/exp4_d3"]] <- c(
    sampleA = 6,
    sampleB = 6
)
io$reclass_stability[["MULTI/exp5_d3.5_d4"]] <- c(
    sampleA = 9,
    sampleB = 10,
    sampleC = 18,
    sampleD = 8
)
io$reclass_stability[["MULTI/exp6_d4.5_d5"]] <- c(
    sampleA = 4,
    sampleB = 4,
    sampleC = 5
)

io$subset.proteincoding <- NULL # Choose between NULL and path to proteincoding genes
io$qc <- TRUE # Choose between NULL and TRUE

io$minUMIs <- 500

io$min_nFeature_RNA <- list()
io$min_nFeature_RNA[["gastr_d3"]] <- 1500
io$min_nFeature_RNA[["gastr_d4"]] <- 1500
io$min_nFeature_RNA[["MULTI/exp1_d5"]] <- c(
    sample1 = 1500,
    sample2 = 1500,
    sample3 = 1500
)
io$min_nFeature_RNA[["MULTI/exp2A_d3_d3.5"]] <- c(
    sample1 = 2000,
    sample2 = 2000,
    sample3 = 2000
)
io$min_nFeature_RNA[["MULTI/exp2C_d3_d3.5"]] <- c(
    sample1 = 500,
    sample2 = 500,
    sample3 = 500
)
io$min_nFeature_RNA[["MULTI/exp2D_d3_d3.5"]] <- c(
    sample1 = 500,
    sample2 = 500
)
io$min_nFeature_RNA[["MULTI/exp3A_d4_d4.5"]] <- c(
    sample1 = 2500,
    sample2 = 2500,
    sample3 = 2500
)
io$min_nFeature_RNA[["MULTI/exp3B_d4_d4.5"]] <- c(
    sample1 = 2500,
    sample2 = 2500,
    sample3 = 2500
)
io$min_nFeature_RNA[["MULTI/exp3C_d4_d4.5"]] <- c(
    sample1 = 2500,
    sample2 = 2500
)
io$min_nFeature_RNA[["MULTI/exp4_d3"]] <- c(
    sampleA = 2500,
    sampleB = 2500
)
io$min_nFeature_RNA[["MULTI/exp5_d3.5_d4"]] <- c(
    sampleA = 2500,
    sampleB = 2500,
    sampleC = 2500,
    sampleD = 2500
)
io$min_nFeature_RNA[["MULTI/exp6_d4.5_d5"]] <- c(
    sampleA = 2500,
    sampleB = 2500,
    sampleC = 2500
)

io$min_nCount_RNA <- list()
io$min_nCount_RNA[["gastr_d3"]] <- 4000
io$min_nCount_RNA[["gastr_d4"]] <- 4000
io$min_nCount_RNA[["MULTI/exp1_d5"]] <- c(
    sample1 = 4000,
    sample2 = 4000,
    sample3 = 4000
)
io$min_nCount_RNA[["MULTI/exp2A_d3_d3.5"]] <- c(
    sample1 = 4000,
    sample2 = 4000,
    sample3 = 4000
)
io$min_nCount_RNA[["MULTI/exp2C_d3_d3.5"]] <- c(
    sample1 = 4000,
    sample2 = 4000,
    sample3 = 4000
)
io$min_nCount_RNA[["MULTI/exp2D_d3_d3.5"]] <- c(
    sample1 = 4000,
    sample2 = 4000
)
io$min_nCount_RNA[["MULTI/exp3A_d4_d4.5"]] <- c(
    sample1 = 8000,
    sample2 = 8000,
    sample3 = 8000
)
io$min_nCount_RNA[["MULTI/exp3B_d4_d4.5"]] <- c(
    sample1 = 8000,
    sample2 = 8000,
    sample3 = 8000
)
io$min_nCount_RNA[["MULTI/exp3C_d4_d4.5"]] <- c(
    sample1 = 8000,
    sample2 = 8000
)
io$min_nCount_RNA[["MULTI/exp4_d3"]] <- c(
    sampleA = 8000,
    sampleB = 8000
)
io$min_nCount_RNA[["MULTI/exp5_d3.5_d4"]] <- c(
    sampleA = 8000,
    sampleB = 8000,
    sampleC = 8000,
    sampleD = 8000
)
io$min_nCount_RNA[["MULTI/exp6_d4.5_d5"]] <- c(
    sampleA = 8000,
    sampleB = 8000,
    sampleC = 8000
)

io$min_percent.mt <- list()
io$min_percent.mt[["gastr_d3"]] <- 3
io$min_percent.mt[["gastr_d4"]] <- 3
io$min_percent.mt[["MULTI/exp1_d5"]] <- c(
    sample1 = 3,
    sample2 = 3,
    sample3 = 3
)
io$min_percent.mt[["MULTI/exp2A_d3_d3.5"]] <- c(
    sample1 = 3,
    sample2 = 3,
    sample3 = 3
)
io$min_percent.mt[["MULTI/exp2C_d3_d3.5"]] <- c(
    sample1 = 3,
    sample2 = 3,
    sample3 = 3
)
io$min_percent.mt[["MULTI/exp2D_d3_d3.5"]] <- c(
    sample1 = 3,
    sample2 = 3
)
io$min_percent.mt[["MULTI/exp3A_d4_d4.5"]] <- c(
    sample1 = 3,
    sample2 = 3,
    sample3 = 3
)
io$min_percent.mt[["MULTI/exp3B_d4_d4.5"]] <- c(
    sample1 = 3,
    sample2 = 3,
    sample3 = 3
)
io$min_percent.mt[["MULTI/exp3C_d4_d4.5"]] <- c(
    sample1 = 3,
    sample2 = 3
)
io$min_percent.mt[["MULTI/exp4_d3"]] <- c(
    sampleA = 1,
    sampleB = 1
)
io$min_percent.mt[["MULTI/exp5_d3.5_d4"]] <- c(
    sampleA = 1,
    sampleB = 1,
    sampleC = 1,
    sampleD = 1
)
io$min_percent.mt[["MULTI/exp6_d4.5_d5"]] <- c(
    sampleA = 1,
    sampleB = 1,
    sampleC = 1
)

io$max_percent.mt <- list()
io$max_percent.mt[["gastr_d3"]] <- 10
io$max_percent.mt[["gastr_d4"]] <- 10
io$max_percent.mt[["MULTI/exp1_d5"]] <- c(
    sample1 = 10,
    sample2 = 10,
    sample3 = 10
)
io$max_percent.mt[["MULTI/exp2A_d3_d3.5"]] <- c(
    sample1 = 10,
    sample2 = 10,
    sample3 = 10
)
io$max_percent.mt[["MULTI/exp2C_d3_d3.5"]] <- c(
    sample1 = 10,
    sample2 = 10,
    sample3 = 10
)
io$max_percent.mt[["MULTI/exp2D_d3_d3.5"]] <- c(
    sample1 = 10,
    sample2 = 10
)
io$max_percent.mt[["MULTI/exp3A_d4_d4.5"]] <- c(
    sample1 = 10,
    sample2 = 10,
    sample3 = 10
)
io$max_percent.mt[["MULTI/exp3B_d4_d4.5"]] <- c(
    sample1 = 10,
    sample2 = 10,
    sample3 = 10
)
io$max_percent.mt[["MULTI/exp3C_d4_d4.5"]] <- c(
    sample1 = 10,
    sample2 = 10
)
io$max_percent.mt[["MULTI/exp4_d3"]] <- c(
    sampleA = 10,
    sampleB = 10
)
io$max_percent.mt[["MULTI/exp5_d3.5_d4"]] <- c(
    sampleA = 10,
    sampleB = 10,
    sampleC = 10,
    sampleD = 10
)
io$max_percent.mt[["MULTI/exp6_d4.5_d5"]] <- c(
    sampleA = 10,
    sampleB = 10,
    sampleC = 10
)
