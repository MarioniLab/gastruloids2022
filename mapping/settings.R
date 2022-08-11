io <- list()

io$min.counts.per.gene <- 50
io$k <- 30
io$min.mean <- 1e-3
io$npcs <- 50

io$atlas.RDS <- NULL # Choose between NULL and path to atlas SCE RDS
io$atlas.metadata <- NULL # Choose between NULL and path to atlas metadata txt
io$corrected.atlas.RDS <- NULL # Choose between NULL and path to atlas SCE RDS
io$corrected.atlas.metadata <- NULL  # Choose between NULL and path to atlas metadata txt
io$atlas_stages <- NULL # Choose between NULL and atlas stages
io$testing <- NULL # Choose between NULL and TRUE


io$order <- list()
io$order <- c("ATLAS", "exp1_d5_sample2", "exp1_d5_sample1", "exp1_d5_sample3", "exp6_d4.5_d5_sampleB", "exp6_d4.5_d5_sampleC", "exp6_d4.5_d5_sampleA", "exp3B_d4_d4.5_sample3", "exp3B_d4_d4.5_sample2", "exp3C_d4_d4.5_sample1", "exp3C_d4_d4.5_sample2", "gastr_d4", "exp5_d3.5_d4_sampleC", "exp5_d3.5_d4_sampleD", "exp5_d3.5_d4_sampleA", "exp5_d3.5_d4_sampleB", "exp2C_d3_d3.5_sample1", "exp2A_d3_d3.5_sample2", "exp2A_d3_d3.5_sample3", "exp2A_d3_d3.5_sample1", "exp4_d3_sampleB", "exp4_d3_sampleA", "gastr_d3" )