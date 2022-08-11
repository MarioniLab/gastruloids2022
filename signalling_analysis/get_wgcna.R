get_wgcna <- function(path_wgcna, organism) {
  
  # Takes output from WGCNA analysis and produces a list of relevant modules with the gene contained within them. 
  # 
  # AUTHOR:	C.Barker; Edited by L. Rosen 12.07.2021
  # INPUT: 
  #          path_wgcna,          path containing the ENSG of WGCNA modules 
  #          is.full              logical - do you want the full list of wgna modules
  #
  # OUTPUT:   list of the relevant WGCNA modules and the gene names contained within them. 
  
  wgcna<-read.delim(file = path_wgcna,header = FALSE) #~/phenotype_networks/data/modules
  wgcna.genes<-wgcna$V1[grep("ENSG", wgcna$V1)]
  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
  wgcna.split<-splitAt(wgcna$V1, grep("ME", wgcna$V1))
  names(wgcna.split)<-wgcna$V1[grep("ME", wgcna$V1)]
  for (x in c(1:length(wgcna.split))) {
    gconvert.table<-gconvert(query = wgcna.split[[x]][-1], 
                             organism = organism, 
                             target="ENSG", 
                             mthreshold = Inf, 
                             filter_na = FALSE)
    wgcna.split[[x]]<-gconvert.table$name
    cat("\r", "WGCNA module conversion progress:", x/length(wgcna.split)*100, "%")
  }
  return(wgcna.split)
}
