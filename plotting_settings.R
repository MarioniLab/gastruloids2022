# gastruloid type colours
gastr_type_colours <- c("mesodermal" = "#8DB5CE",
                        "neural" = "#65A83E",
                        "intermediate" = "#7f7f7f"
                       )

# gastruloid timepoint colours
timepoint_colours <- c("d3" = "#cd1400",
                       "d3_d3.5" = "#DA4800",
                       "d3.5" = "#E67D00",
                       "d3.5_d4" = "#EEA416",
                       "d4" = "#f6ca2b",
                       "d4_d4.5" = "#BEBC50",
                       "d4.5" = "#85ad76",
                       "d4.5_d5" = "#699D95",
                       "d5" = "#4d8db4"
                      )



# gastruloid experiment + timepoint colours (so experiments are split by which timepoint the cell is)
tp_experiment_colours <- c("d3_10X" = "#A91200",
                           "d3_exp4" = "#F21A00",
                           "d3.5_exp5" = "#E67D00",
                           "d4_10X" = "#ffdf54",
                           "d4_exp5" = "#ffd21f",
                           #"4_exp3B_d4_d4.5" = "#e3b80e",
                           "d4_exp3B" = "#e3b10e",
                           "d4.5_exp3B" = "#a3c995",
                           "d4.5_exp6" = "#669156",
                           "d5_exp6" = "#78b7c5",
                           "d5_exp1" = "#2162A3"
                          )
experiment_tp_colours <- c("10X_d3" = "#A91200",
                           "exp4_d3" = "#F21A00",
                           "exp5_d3.5" = "#E67D00",
                           "10X_d4" = "#ffdf54",
                           "exp5_d4" = "#ffd21f",
                           #"4_exp3B_d4_d4.5" = "#e3b80e",
                           "exp3B_d4" = "#e3b10e",
                           "exp3B_d4.5" = "#a3c995",
                           "exp6_d4.5" = "#669156",
                           "exp6_d5" = "#78b7c5",
                           "exp1_d5" = "#2162A3"
                          )

gastruloid_order <- list()
gastruloid_order$d3 <- c('Bar49_exp4_d3', 'Bar24_exp4_d3', 'Bar41_exp4_d3', 'Bar28_exp4_d3', 'Bar66_exp4_d3', 'Bar45_exp4_d3', 'Bar20_exp4_d3', 'Bar25_exp4_d3', 'Bar47_exp4_d3', 'Bar59_exp4_d3', 'Bar52_exp4_d3', 'Bar56_exp4_d3', 'Bar30_exp4_d3', 'Bar55_exp4_d3', 'Bar19_exp4_d3', 'Bar13_exp4_d3', 'Bar14_exp4_d3', 'Bar2_exp4_d3', 'Bar6_exp4_d3', 'Bar18_exp4_d3', 'Bar61_exp4_d3', 'Bar65_exp4_d3', 'Bar15_exp4_d3', 'Bar12_exp4_d3')
gastruloid_order$d5 <- c('Bar16_exp1_d5', 'Bar4_exp1_d5', 'Bar13_exp1_d5', 'Bar7_exp6_d4.5_d5', 'Bar15_exp1_d5', 'Bar2_exp6_d4.5_d5', 'Bar12_exp1_d5', 'Bar10_exp6_d4.5_d5', 'Bar22_exp1_d5', 'Bar21_exp6_d4.5_d5', 'Bar3_exp1_d5', 'Bar17_exp1_d5', 'Bar14_exp6_d4.5_d5', 'Bar7_exp1_d5', 'Bar19_exp6_d4.5_d5', 'Bar18_exp1_d5', 'Bar10_exp1_d5', 'Bar19_exp1_d5', 'Bar20_exp1_d5', 'Bar5_exp1_d5', 'Bar14_exp1_d5', 'Bar24_exp1_d5', 'Bar2_exp1_d5', 'Bar1_exp6_d4.5_d5', 'Bar6_exp1_d5', 'Bar22_exp6_d4.5_d5', 'Bar9_exp1_d5', 'Bar13_exp6_d4.5_d5', 'Bar12_exp6_d4.5_d5', 'Bar23_exp1_d5', 'Bar11_exp1_d5', 'Bar15_exp6_d4.5_d5', 'Bar8_exp1_d5', 'Bar21_exp1_d5', 'Bar18_exp6_d4.5_d5', 'Bar6_exp6_d4.5_d5', 'Bar1_exp1_d5', 'Bar17_exp6_d4.5_d5', 'Bar4_exp6_d4.5_d5', 'Bar20_exp6_d4.5_d5')
gastruloid_order$d3.5 <- c('Bar30_exp5_d3.5_d4', 'Bar23_exp5_d3.5_d4', 'Bar47_exp5_d3.5_d4', 'Bar66_exp5_d3.5_d4', 'Bar45_exp5_d3.5_d4', 'Bar59_exp5_d3.5_d4', 'Bar25_exp5_d3.5_d4', 'Bar41_exp5_d3.5_d4', 'Bar24_exp5_d3.5_d4', 'Bar52_exp5_d3.5_d4', 'Bar55_exp5_d3.5_d4', 'Bar65_exp5_d3.5_d4', 'Bar56_exp5_d3.5_d4', 'Bar28_exp5_d3.5_d4', 'Bar49_exp5_d3.5_d4', 'Bar61_exp5_d3.5_d4')
gastruloid_order$d4 <- c('Bar22_exp3B_d4_d4.5', 'Bar49_exp3B_d4_d4.5', 'Bar2_exp3B_d4_d4.5', 'Bar25_exp3B_d4_d4.5', 'Bar20_exp5_d3.5_d4', 'Bar22_exp5_d3.5_d4', 'Bar6_exp3B_d4_d4.5', 'Bar5_exp3B_d4_d4.5', 'Bar30_exp3B_d4_d4.5', 'Bar13_exp5_d3.5_d4', 'Bar56_exp3B_d4_d4.5', 'Bar1_exp5_d3.5_d4', 'Bar7_exp5_d3.5_d4', 'Bar55_exp3B_d4_d4.5', 'Bar12_exp3B_d4_d4.5', 'Bar21_exp5_d3.5_d4', 'Bar4_exp5_d3.5_d4', 'Bar12_exp5_d3.5_d4', 'Bar14_exp5_d3.5_d4', 'Bar28_exp3B_d4_d4.5', 'Bar10_exp5_d3.5_d4', 'Bar18_exp5_d3.5_d4', 'Bar17_exp5_d3.5_d4', 'Bar6_exp5_d3.5_d4', 'Bar2_exp5_d3.5_d4', 'Bar52_exp3B_d4_d4.5', 'Bar19_exp5_d3.5_d4', 'Bar15_exp5_d3.5_d4')
gastruloid_order$d4.5 <- c('Bar65_exp6_d4.5_d5', 'Bar28_exp6_d4.5_d5', 'Bar15_exp3B_d4_d4.5', 'Bar45_exp6_d4.5_d5', 'Bar47_exp3B_d4_d4.5', 'Bar19_exp3B_d4_d4.5', 'Bar61_exp3B_d4_d4.5', 'Bar55_exp6_d4.5_d5', 'Bar25_exp6_d4.5_d5', 'Bar23_exp6_d4.5_d5', 'Bar13_exp3B_d4_d4.5', 'Bar56_exp6_d4.5_d5', 'Bar66_exp6_d4.5_d5', 'Bar45_exp3B_d4_d4.5', 'Bar52_exp6_d4.5_d5', 'Bar61_exp6_d4.5_d5', 'Bar49_exp6_d4.5_d5', 'Bar41_exp6_d4.5_d5', 'Bar66_exp3B_d4_d4.5', 'Bar14_exp3B_d4_d4.5', 'Bar39_exp3B_d4_d4.5', 'Bar59_exp3B_d4_d4.5', 'Bar24_exp6_d4.5_d5', 'Bar47_exp6_d4.5_d5', 'Bar65_exp3B_d4_d4.5', 'Bar59_exp6_d4.5_d5', 'Bar30_exp6_d4.5_d5', 'Bar41_exp3B_d4_d4.5')

# which gastruloids are assigned which class
mesodermal <- list('d3' = c('Bar12_exp4_d3', 'Bar6_exp4_d3', 'Bar15_exp4_d3', 'Bar65_exp4_d3', 'Bar18_exp4_d3', 'Bar2_exp4_d3', 'Bar14_exp4_d3', 'Bar13_exp4_d3', 'Bar61_exp4_d3', 'Bar19_exp4_d3', 'Bar55_exp4_d3'),
        'd3.5' = c('Bar28_exp5_d3.5_d4', 'Bar56_exp5_d3.5_d4', 'Bar65_exp5_d3.5_d4', 'Bar55_exp5_d3.5_d4', 'Bar49_exp5_d3.5_d4', 'Bar61_exp5_d3.5_d4'),
        'd4' = c('Bar21_exp5_d3.5_d4', 'Bar28_exp3B_d4_d4.5', 'Bar12_exp3B_d4_d4.5', 'Bar12_exp5_d3.5_d4', 
               'Bar4_exp5_d3.5_d4', 'Bar17_exp5_d3.5_d4', 'Bar55_exp3B_d4_d4.5', 'Bar18_exp5_d3.5_d4', 
               'Bar56_exp3B_d4_d4.5', 'Bar10_exp5_d3.5_d4', 'Bar14_exp5_d3.5_d4', 'Bar13_exp5_d3.5_d4',
               'Bar52_exp3B_d4_d4.5', 'Bar2_exp5_d3.5_d4', 'Bar6_exp5_d3.5_d4', 'Bar7_exp5_d3.5_d4',
               'Bar1_exp5_d3.5_d4', 'Bar5_exp3B_d4_d4.5', 'Bar15_exp5_d3.5_d4', 'Bar30_exp3B_d4_d4.5',
               'Bar19_exp5_d3.5_d4', 'Bar6_exp3B_d4_d4.5'),
        'd4.5' = c('Bar66_exp6_d4.5_d5', 'Bar45_exp3B_d4_d4.5', 'Bar61_exp6_d4.5_d5', 'Bar14_exp3B_d4_d4.5', 'Bar39_exp3B_d4_d4.5', 'Bar49_exp6_d4.5_d5',
                 'Bar66_exp3B_d4_d4.5', 'Bar41_exp6_d4.5_d5', 'Bar52_exp6_d4.5_d5', 'Bar47_exp6_d4.5_d5', 'Bar59_exp3B_d4_d4.5', 'Bar24_exp6_d4.5_d5',
                 'Bar59_exp6_d4.5_d5', 'Bar65_exp3B_d4_d4.5', 'Bar41_exp3B_d4_d4.5', 'Bar30_exp6_d4.5_d5'),
        'd5' = c('Bar1_exp1_d5', 'Bar17_exp6_d4.5_d5', 'Bar21_exp1_d5', 'Bar20_exp6_d4.5_d5', 'Bar8_exp1_d5', 'Bar6_exp6_d4.5_d5', 'Bar4_exp6_d4.5_d5',
               'Bar23_exp1_d5', 'Bar18_exp6_d4.5_d5', 'Bar15_exp6_d4.5_d5', 'Bar11_exp1_d5', 'Bar9_exp1_d5', 'Bar6_exp1_d5', 'Bar22_exp6_d4.5_d5', 'Bar13_exp6_d4.5_d5')
       )
neural <- list('d3' = c('Bar56_exp4_d3', 'Bar30_exp4_d3', 'Bar52_exp4_d3', 'Bar59_exp4_d3', 'Bar25_exp4_d3', 'Bar47_exp4_d3', 'Bar45_exp4_d3', 'Bar20_exp4_d3', 'Bar66_exp4_d3', 'Bar28_exp4_d3', 'Bar41_exp4_d3', 'Bar24_exp4_d3', 'Bar49_exp4_d3'),
          'd3.5' = c('Bar30_exp5_d3.5_d4', 'Bar47_exp5_d3.5_d4', 'Bar59_exp5_d3.5_d4', 'Bar24_exp5_d3.5_d4', 'Bar25_exp5_d3.5_d4', 'Bar41_exp5_d3.5_d4',
                   'Bar52_exp5_d3.5_d4', 'Bar45_exp5_d3.5_d4', 'Bar23_exp5_d3.5_d4', 'Bar66_exp5_d3.5_d4'),
          'd4' = c('Bar22_exp3B_d4_d4.5', 'Bar22_exp5_d3.5_d4', 'Bar25_exp3B_d4_d4.5', 'Bar20_exp5_d3.5_d4',
                 'Bar2_exp3B_d4_d4.5', 'Bar49_exp3B_d4_d4.5'),
          'd4.5' = c('Bar65_exp6_d4.5_d5', 'Bar28_exp6_d4.5_d5', 'Bar45_exp6_d4.5_d5', 'Bar47_exp3B_d4_d4.5', 'Bar19_exp3B_d4_d4.5', 'Bar15_exp3B_d4_d4.5',
                   'Bar13_exp3B_d4_d4.5', 'Bar23_exp6_d4.5_d5', 'Bar55_exp6_d4.5_d5', 'Bar61_exp3B_d4_d4.5', 'Bar25_exp6_d4.5_d5', 'Bar56_exp6_d4.5_d5'),
          'd5' = c('Bar17_exp1_d5', 'Bar3_exp1_d5', 'Bar10_exp6_d4.5_d5', 'Bar21_exp6_d4.5_d5', 'Bar2_exp6_d4.5_d5', 'Bar22_exp1_d5', 'Bar12_exp1_d5', 'Bar15_exp1_d5',
                 'Bar13_exp1_d5', 'Bar7_exp6_d4.5_d5', 'Bar4_exp1_d5', 'Bar16_exp1_d5')
         )

intermediate <- list('d5' = c('Bar10_exp1_d5', 'Bar12_exp6_d4.5_d5', 'Bar14_exp1_d5', 'Bar14_exp6_d4.5_d5',
                              'Bar18_exp1_d5', 'Bar19_exp1_d5', 'Bar19_exp6_d4.5_d5', 'Bar1_exp6_d4.5_d5', 'Bar20_exp1_d5',
                              'Bar24_exp1_d5', 'Bar2_exp1_d5', 'Bar5_exp1_d5', 'Bar7_exp1_d5')
                    )

# OG 116,312-cell atlas colours
celltype_colours = c(
 "Epiblast" = "#635547",
 "Primitive Streak" = "#DABE99",
 "Caudal epiblast" = "#9e6762",
 "PGC" = "#FACB12",
 "Anterior Primitive Streak" = "#c19f70",
 "Notochord" = "#0F4A9C",
 "Def. endoderm" = "#F397C0",
 "Gut" = "#EF5A9D",
 "Nascent mesoderm" = "#C594BF",
 "Mixed mesoderm" = "#DFCDE4",
 "Intermediate mesoderm" = "#139992",
 "Caudal Mesoderm" = "#3F84AA",
 "Paraxial mesoderm" = "#8DB5CE",
 "Somitic mesoderm" = "#005579",
 "Pharyngeal mesoderm" = "#C9EBFB",
 "Cardiomyocytes" = "#B51D8D",
 "Allantois" = "#532C8A",
 "ExE mesoderm" = "#8870ad",
 "Mesenchyme" = "#cc7818",
 "Haematoendothelial progenitors" = "#FBBE92",
 "Endothelium" = "#ff891c",
 "Blood progenitors 1" = "#f9decf",
 "Blood progenitors 2" = "#c9a997",
 "Erythroid1" = "#C72228",
 "Erythroid2" = "#f79083",
 "Erythroid3" = "#EF4E22",
 "NMP" = "#8EC792",
 "Rostral neurectoderm" = "#65A83E",
 "Caudal neurectoderm" = "#354E23",
 "Neural crest" = "#C3C388",
 "Forebrain/Midbrain/Hindbrain" = "#647a4f",
 "Spinal cord" = "#CDE088",
 "Surface ectoderm" = "#f7f79e",
 "Visceral endoderm" = "#F6BFCB",
 "ExE endoderm" = "#7F6874",
 "ExE ectoderm" = "#989898",
 "Parietal endoderm" = "#1A1A1A",
 "Anterior-most_somites" = "red",
 "Head_mesoderm" = "darkgreen",
 "Dermomyotome" = "brown",
 "Presomitic_mesoderm" = "blue",
 "Sclerotome" = "yellow",
 "Posterior-most_somites" = "turquoise"
)

# extended atlas colours
celltype_colours_final = c(
    "Epiblast" = "#635547",
    "Primitive Streak" = "#DABE99",
    "Caudal epiblast" = "#9E6762",
    "PGC" = "#FACB12",
    "Anterior Primitive Streak" = "#C19F70",
    "Node" = "#153B3D",
    "Notochord" = "#0F4A9C",
    "Gut tube" = "#EF5A9D",
    "Hindgut" = "#F397C0",
    "Midgut" = "#FF00B2",
    "Foregut" = "#FFB7FF",
    "Pharyngeal endoderm" = "#95E1FF",
    "Thyroid primordium" = "#97BAD3",
    "Nascent mesoderm" = "#C594BF",
    "Intermediate mesoderm" = "#139992",
    "Caudal mesoderm" = "#3F84AA",
    "Lateral plate mesoderm" = "#F9DFE6",
    "Limb mesoderm" = "#E35F82",
    "Forelimb" = "#D02D75",
    "Kidney primordium" = "#E85639",
    "Presomitic mesoderm" = "#5581CA", #"0000ff", #"blue",
    "Somitic mesoderm" = "#005579",
    "Posterior somitic tissues" = "#5ADBE4", #"40e0d0",#"turquoise",
    "Paraxial mesoderm" = "#8DB5CE",
    "Cranial mesoderm" = "#456722",#“#006400”,#darkgreen",
    "Anterior somitic tissues" = "#D5E839",
    "Sclerotome" = "#E3CB3A", #"ffff00", #"yellow",
    "Dermomyotome" = "#00BFC4", #"a52a2a", #"brown",
    "Pharyngeal mesoderm" = "#C9EBFB",
    "Cardiopharyngeal progenitors" = "#556789",
    "Anterior cardiopharyngeal progenitors" = "#683ED8",
    "Allantois" = "#532C8A",
    "Mesenchyme" = "#CC7818",
    "YS mesothelium" = "#FF7F9C",
    "Epicardium" = "#F79083",
    "Embryo proper mesothelium" = "#FF487D",
    "Cardiopharyngeal progenitors FHF" = "#D780B0",
    "Cardiomyocytes FHF 1" = "#A64D7E",
    "Cardiomyocytes FHF 2" = "#B51D8D",
    "Cardiopharyngeal progenitors SHF" = "#4B7193",
    "Cardiomyocytes SHF 1" = "#5D70DC",
    "Cardiomyocytes SHF 2" = "#332C6C",
    "Haematoendothelial progenitors" = "#FBBE92",
    "Blood progenitors" = "#6C4B4C",
    "Erythroid" = "#C72228",
    "Chorioallantoic-derived erythroid progenitors" = "#E50000",
    "Megakaryocyte progenitors" = "#E3CB3A",
    "MEP" = "#EF4E22",
    "EMP" = "#7C2A47",
    "YS endothelium" = "#FF891C",
    "YS mesothelium-derived endothelial progenitors" = "#AE3F3F",
    "Allantois endothelium" = "#2F4A60",
    "Embryo proper endothelium" = "#90E3BF",
    "Venous endothelium" = "#BD3400",
    "Endocardium" = "#9D0049",
    "NMPs/Mesoderm-biased" = "#89C1F5",
    "NMPs" = "#8EC792",
    "Ectoderm" = "#FF675C",
    "Optic vesicle" = "#BD7300",
    "Ventral forebrain progenitors" = "#A0B689",
    "Early dorsal forebrain progenitors" = "#0F8073",
    "Late dorsal forebrain progenitors" = "#7A9941",
    "Midbrain/Hindbrain boundary" = "#8AB3B5",
    "Midbrain progenitors" = "#9BF981",
    "Dorsal midbrain neurons" = "#12ED4C",
    "Ventral hindbrain progenitors" = "#7E907A",
    "Dorsal hindbrain progenitors" = "#2C6521",
    "Hindbrain floor plate" = "#BF9DA8",
    "Hindbrain neural progenitors" = "#59B545",
    "Neural tube" = "#233629",
    "Migratory neural crest" = "#4A6798",
    "Branchial arch neural crest" = "#BD84B0",
    "Frontonasal mesenchyme" = "#D3B1B1",
    "Spinal cord progenitors" = "#6B2035",
    "Dorsal spinal cord progenitors" = "#E273D6",
    "Non-neural ectoderm" = "#F7F79E",
    "Surface ectoderm" = "#FCFF00",
    "Epidermis" = "#FFF335",
    "Limb ectoderm" = "#FFD731",
    "Amniotic ectoderm" = "#DBB400",
    "Placodal ectoderm" = "#FF5C00",
    "Otic placode" = "#F1A262",
    "Otic neural progenitors" = "#00B000",
    "Visceral endoderm" = "#F6BFCB",
    "ExE endoderm" = "#7F6874",
    "ExE ectoderm" = "#989898",
    "Parietal endoderm" = "#1A1A1A"
)

cluster_colours = c(
    "Epiblast" = "#635547",
    "Primitive Streak" = "#DABE99",
    "Anterior Primitive Streak" = "#a68250",
    "Early Nascent Mesoderm" = "#DFCDE4",
    "Late Nascent Mesoderm" = "#C594BF",
    "Head Mesoderm" = "#7E508D",
    "Cardiopharyngeal Mesoderm" = "#B51D8D",
    "Endothelium" = "#ff891c",
    "Early Posterior PSM" = "#6eb8b5",
    "Late Posterior PSM" = "#479b98",
    "Early Anterior PSM" = "#8ac5e3",
    "Late Anterior PSM" = "#459AC0",
    "Somites" = "#006f9e",
    #"Somites" = "#00405c",
    "Mature Somites" = "#00405c",
    #"Mature Somites" = "#00BFC4",
    "Caudal Epiblast" = "#9E6762",
    "NMPs" = "#8EC792",
    "Caudal Neurectoderm" = "#65A83E",
    "Early Spinal Cord" = "#99C463",
    "Late Spinal Cord" = "#CDE088",
    "Early Neurectoderm" = "#649146",
    "Late Neurectoderm" = "#647a4f",
    "Neurons" = "#414f33",
    "Notochord" = "#0F4A9C",
    "Mature Endoderm" = "#EF5A9D",
    "PGCs" = "#FACB12"
)

# the order in which I plot the clusters
cluster_order <- c("PGCs", "Epiblast", "Primitive Streak", "Anterior Primitive Streak",
                   "Mature Endoderm", "Notochord",
                   "Early Nascent Mesoderm", "Late Nascent Mesoderm", "Head Mesoderm", "Cardiopharyngeal Mesoderm", "Endothelium",
                   "Early Posterior PSM", "Late Posterior PSM", "Early Anterior PSM", "Late Anterior PSM", "Somites", "Mature Somites",
                   "Caudal Epiblast", "NMPs", "Caudal Neurectoderm", "Early Neurectoderm", "Late Neurectoderm", "Early Spinal Cord", "Late Spinal Cord", "Neurons"
                  )

celltype_order <- c("PGC",
                    "Anterior somitic tissues", "Limb mesoderm",
                    "Presomitic mesoderm", "Posterior somitic tissues", "Somitic mesoderm", "Sclerotome", "Dermomyotome",
                    "Pharyngeal mesoderm", "Embryo proper endothelium", "Venous endothelium", "Anterior cardiopharyngeal progenitors", "Cardiopharyngeal progenitors", "Cardiopharyngeal progenitors SHF", "Cardiopharyngeal progenitors FHF", "Cardiomyocytes FHF 1", "Epicardium",
                    "Caudal mesoderm", "NMPs/Mesoderm-biased", "NMPs", 
                    "Dorsal spinal cord progenitors", "Spinal cord progenitors", "Hindbrain neural progenitors", "Dorsal hindbrain progenitors", "Ventral hindbrain progenitors", "Midbrain/Hindbrain boundary", "Midbrain progenitors", "Dorsal midbrain neurons", "Optic vesicle", "Thyroid primordium",
                    "Non-neural ectoderm 5", 
                    "Notochord", "Pharyngeal endoderm", "Hindgut"
                   )

gastr_type_order <- c("mesodermal", "intermediate", "neural")

timepoint_order <- c("d3", "d3_d3.5", "d3.5", "d3.5_d4", "d4", "d4_d4.5", "d4.5", "d4.5_d5", "d5")

# I think these are the OG atlas stage colours but idk why they go to E9.5...
stage_colours <- c("E6.5" = "#D53E4F",
                   "E6.75" = "#F46D43",
                   "E7.0" = "#FDAE61",
                   "E7.25" = "#FEE08B",
                   "E7.5" = "#FFFFBF",
                   "E7.75" = "#E6F598",
                   "E8.0" = "#ABDDA4",
                   "E8.25" = "#66C2A5",
                   "E8.5" = "#3288BD",
                   "E8.75" = "#3250BD",
                   "E9.0" = "#5532bd",
                   "E9.25" = "#9a32bd",
                   "E9.5" = "#bd329a",
                   "mixed_gastrulation" = "#A9A9A9"
                  )

# Ivan's extended atlas stage colours
library(wesanderson)
# For original stage annotation
stage_colours_extension  <- c(rev(wes_palette("Zissou1", 13, type = "continuous")),"#A9A9A9")
stage_colours_extension[11:13] <- c("#3399FF", "#297ACC", "#2162A3")
#stage_colours_extension <- c(colorRampPalette(c("red", "orange", "yellow", "green", "blue"), space = "Lab")(13),"#A9A9A9")
names(stage_colours_extension) <- c("E6.5", "E6.75", "E7.0", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "mixed_gastrulation")

barplot.pub <- function(to.plot, x="cluster", colors=NULL, xlim.max=NULL, y="N") {
  p <- ggplot(to.plot, aes_string(x=x, y=y)) +
    scale_x_discrete(drop=FALSE) + 
    # coord_flip() +
    labs(y="Number of cells") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.3)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
    
    #if (is.null(colors)) {
    #    p <- p + geom_bar(stat="identity", color="black")
    #}
    if (is.null(colors)) {
        p <- p + geom_bar(stat="identity", aes_string(fill=x))
    } else {
        p <- p + geom_bar(aes_string(fill=x), stat="identity", color="black") + 
            scale_fill_manual(values=colors, drop=F)
    }

    if (!is.null(xlim.max)) {
      p <- p + coord_flip(ylim=c(0,xlim.max))
    } else {
      p <- p + coord_flip()
    }
  
    return(p)
}

plotExpression_UMAP <- function(srat,
                                gene_ensid=NULL, gene_symbol=NULL,
                                smoothed=FALSE, type='normalised',
                                assay='RNA', graph_name='RNA_nn', dim1_name='X1', dim2_name='X2',
                                cluster_subset = NULL, tp_subset = NULL, gastr_type_subset = NULL, dataset_subset = NULL,
                                facet_by_r=NULL, facet_by_c=NULL
                               ) {
    
    if (is.null(gene_ensid) & is.null(gene_symbol)) {
        stop("either ensembl id or symbol of gene must be given")
    }
    if (!(type %in% c('normalised', 'counts', 'log2'))) {
        stop("type must be one of normalised, counts, or log2")
    }
    if ((smoothed) & is.null(srat@graphs[[graph_name]])) {
        stop("graph_name not among srat graphs (did you remember to first FindNeighbors?")
    }
    
    if (is.null(gene_symbol)) {
        if (!(gene_ensid %in% srat@assays[[assay]]@meta.features$ens_id)) {
            stop("can't find gene_ensid in meta.features ens_id column")
        }
        gene_symbol <- srat@assays[[assay]]@meta.features$symbol[srat@assays[[assay]]@meta.features$ens_id == gene_ensid][1]
    }
    if (is.null(gene_ensid)) {
        if (!(gene_symbol %in% srat@assays[[assay]]@meta.features$symbol)) {
            stop("can't find gene_symbol in meta.features symbol column")
        }
        gene_ensid <- srat@assays[[assay]]@meta.features$ens_id[srat@assays[[assay]]@meta.features$symbol == gene_symbol][1]
    }
    
    if (type == 'normalised') {
        if (smoothed) {
            to.plot <- as.matrix(srat@assays[[assay]]@data[c(gene_ensid),])
            idx <- split(as(t(srat@graphs[[graph_name]]), "dgTMatrix")@i+1, rownames(srat@graphs[[graph_name]])[as(t(srat@graphs[[graph_name]]), "dgTMatrix")@j+1])
            gene_ave <- lapply(idx, function(x) mean(to.plot[x]))
            to.plot <- do.call(rbind, gene_ave) %>% as.data.frame
        } else {
            to.plot <- as.data.frame(srat@assays[[assay]]@data[c(gene_ensid),])
        }
        
    } else if (type %in% c('counts', 'log2')) {
        if (smoothed) {
            to.plot <- as.matrix(srat@assays[[assay]]@counts[c(gene_ensid),])
            idx <- split(as(t(srat@graphs[[graph_name]]), "dgTMatrix")@i+1, rownames(srat@graphs[[graph_name]])[as(t(srat@graphs[[graph_name]]), "dgTMatrix")@j+1])
            gene_ave <- lapply(idx, function(x) mean(to.plot[x]))
            to.plot <- do.call(rbind, gene_ave) %>% as.data.frame
        } else {
            to.plot <- as.data.frame(srat@assays[[assay]]@counts[c(gene_ensid),])
        }
        
    }
    
    colnames(to.plot) <- c(str_replace(gene_symbol, '-', '_'))
    to.plot$cell <- rownames(to.plot)
    to.plot <- merge(to.plot, srat@meta.data, by="cell") %>% as.data.table

    if (!(is.null(cluster_subset))) {
        to.plot <- to.plot[cluster %in% cluster_subset]
    }
    if (!(is.null(tp_subset))) {
        to.plot <- to.plot[timepoint %in% tp_subset]
    }
    if (!(is.null(gastr_type_subset))) {
        to.plot <- to.plot[gastr_type %in% gastr_type_subset]
    }
    if (!(is.null(dataset_subset))) {
        to.plot <- to.plot[dataset %in% dataset_subset]
    }
    
    #to.plot <- to.plot[cluster %in% cluster_order]
    #to.plot$cluster <- factor(to.plot$cluster, levels = cluster_order)
    
    p <- ggplot(data=to.plot, mapping = aes_string(x=dim1_name, y=dim2_name, colour=str_replace(gene_symbol, '-', '_')))
    
    if (is.null(facet_by_r) & is.null(facet_by_c)) {
        tmp <- NULL
    } else if ((!is.null(facet_by_r)) & (is.null(facet_by_c))) {
        if (facet_by_r == "gastr_type") {
            p <- p + geom_point(data=select(to.plot,-gastr_type), colour="grey", size=0.1, alpha=1)
        } else if (facet_by_r == "timepoint") {
            p <- p + geom_point(data=select(to.plot,-timepoint), colour="grey", size=0.1, alpha=1)
        } else if (facet_by_r == "dataset") {
            p <- p + geom_point(data=select(to.plot,-dataset), colour="grey", size=0.1, alpha=1)
        }
    } else if ((!is.null(facet_by_c)) & (is.null(facet_by_r))) {
        if (facet_by_c == "gastr_type") {
            p <- p + geom_point(data=select(to.plot,-gastr_type), colour="grey", size=0.1, alpha=1)
        } else if (facet_by_c == "timepoint") {
            p <- p + geom_point(data=select(to.plot,-timepoint), colour="grey", size=0.1, alpha=1)
        } else if (facet_by_c == "dataset") {
            p <- p + geom_point(data=select(to.plot,-dataset), colour="grey", size=0.1, alpha=1)
        }
    } else if (((facet_by_r == "timepoint") & (facet_by_c == "gastr_type")) | ((facet_by_c == "timepoint") & (facet_by_r == "gastr_type"))) {
        p <- p + geom_point(data=select(to.plot,-c(timepoint, gastr_type)), colour="grey", size=0.1, alpha=1)
    } else if (((facet_by_r == "timepoint") & (facet_by_c == "dataset")) | ((facet_by_c == "timepoint") & (facet_by_r == "dataset"))) {
        p <- p + geom_point(data=select(to.plot,-c(timepoint, dataset)), colour="grey", size=0.1, alpha=1)
    } else if (((facet_by_r == "gastr_type") & (facet_by_c == "dataset")) | ((facet_by_c == "gastr_type") & (facet_by_r == "dataset"))) {
        p <- p + geom_point(data=select(to.plot,-c(dataset, gastr_type)), colour="grey", size=0.1, alpha=1)
    } else if ((!(is.null(facet_by_c))) | (!(is.null(facet_by_r)))) {
        stop('unknown facet_by_c or unknown facet_by_r')
    }
    
    p <- p + geom_point(size=0.1, alpha=1) +
      # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
      labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
      theme_classic() +
      theme(axis.text=element_blank(),
            axis.ticks = element_blank()
           )
    
    if (type == 'log2') {
        p <- p + scale_colour_viridis(trans='log2')
    } else {
        p <- p + scale_colour_viridis()
    }
    
    if (is.null(facet_by_r) & is.null(facet_by_c)) {
        tmp <- NULL
    } else if ((!is.null(facet_by_r)) & (is.null(facet_by_c))) {
        if (facet_by_r == "gastr_type") {
            p <- p + facet_wrap(~ gastr_type, ncol=1)
        } else if (facet_by_r == "timepoint") {
            p <- p + facet_wrap(~ timepoint, ncol=1)
        } else if (facet_by_r == "dataset") {
            p <- p + facet_wrap(~ dataset, ncol=1)
        }
    } else if ((!is.null(facet_by_c)) & (is.null(facet_by_r))) {
        if (facet_by_c == "gastr_type") {
            p <- p + facet_wrap(~ gastr_type, nrow=1)
        } else if (facet_by_c == "timepoint") {
            p <- p + facet_wrap(~ timepoint, nrow=1)
        } else if (facet_by_c == "dataset") {
            p <- p + facet_wrap(~ dataset, nrow=1)
        }
    } else if ((facet_by_r == "timepoint") & (facet_by_c == "gastr_type")) {
        p <- p + facet_grid(timepoint ~ gastr_type)
    } else if ((facet_by_r == "gastr_type") & (facet_by_c == "timepoint")) {
        p <- p + facet_grid(gastr_type ~ timepoint)
    } else if ((facet_by_r == "timepoint") & (facet_by_c == "dataset")) {
        p <- p + facet_grid(timepoint ~ dataset)
    } else if ((facet_by_r == "dataset") & (facet_by_c == "timepoint")) {
        p <- p + facet_grid(dataset ~ timepoint)
    } else if ((facet_by_r == "gastr_type") & (facet_by_c == "dataset")) {
        p <- p + facet_grid(gastr_type ~ dataset)
    } else if ((facet_by_r == "dataset") & (facet_by_c == "gastr_type")) {
        p <- p + facet_grid(dataset ~ gastr_type)
    } else if ((!(is.null(facet_by_c))) | (!(is.null(facet_by_r)))) {
        stop('unknown facet_by_c or unknown facet_by_r')
    }
    
    return(p)
}

plotGeneCorrelation <- function(srat, gene1_ensid=NULL, gene1_symbol=NULL, gene2_ensid=NULL, gene2_symbol=NULL, smoothed=TRUE, type='normalised', assay='RNA', graph_name='RNA_nn', cluster_subset = NULL, tp_subset = NULL, gastr_type_subset = NULL) {
    
    if (is.null(gene1_ensid) & is.null(gene1_symbol)) {
        stop("either ensembl id or symbol of gene 1 must be given")
    }
    if (is.null(gene2_ensid) & is.null(gene2_symbol)) {
        stop("either ensembl id or symbol of gene 2 must be given")
    }
    if (!(type %in% c('normalised', 'counts', 'log2'))) {
        stop("type must be one of normalised, counts, or log2")
    }
    if ((smoothed) & is.null(srat@graphs[[graph_name]])) {
        stop("graph_name not among srat graphs (did you remember to first FindNeighbors?")
    }
    
    if (is.null(gene1_symbol)) {
        if (!(gene1_ensid %in% srat@assays[[assay]]@meta.features$ens_id)) {
            stop("can't find gene1_ensid in meta.features ens_id column")
        }
        gene1_symbol <- srat@assays[[assay]]@meta.features$symbol[srat@assays[[assay]]@meta.features$ens_id == gene1_ensid][1]
    }
    if (is.null(gene1_ensid)) {
        if (!(gene1_symbol %in% srat@assays[[assay]]@meta.features$symbol)) {
            stop("can't find gene1_symbol in meta.features symbol column")
        }
        gene1_ensid <- srat@assays[[assay]]@meta.features$ens_id[srat@assays[[assay]]@meta.features$symbol == gene1_symbol][1]
    }
    if (is.null(gene2_symbol)) {
        if (!(gene2_ensid %in% srat@assays[[assay]]@meta.features$ens_id)) {
            stop("can't find gene2_ensid in meta.features ens_id column")
        }
        gene2_symbol <- srat@assays[[assay]]@meta.features$symbol[srat@assays[[assay]]@meta.features$ens_id == gene2_ensid][1]
    }
    if (is.null(gene2_ensid)) {
        if (!(gene2_symbol %in% srat@assays[[assay]]@meta.features$symbol)) {
            stop("can't find gene2_symbol in meta.features symbol column")
        }
        gene2_ensid <- srat@assays[[assay]]@meta.features$ens_id[srat@assays[[assay]]@meta.features$symbol == gene2_symbol][1]
    }
    
    if (type == 'normalised') {
        if (smoothed) {
            to.plot <- as.matrix(srat@assays[[assay]]@data[c(gene1_ensid, gene2_ensid),])
            idx <- split(as(t(srat@graphs[[graph_name]]), "dgTMatrix")@i+1, rownames(srat@graphs[[graph_name]])[as(t(srat@graphs[[graph_name]]), "dgTMatrix")@j+1])
            gene_ave <- lapply(idx, function(x) rowMeans(to.plot[,x]))
            to.plot <- do.call(rbind, gene_ave) %>% as.data.frame
        } else {
            to.plot <- as.data.frame(t(srat@assays[[assay]]@data[c(gene1_ensid, gene2_ensid),]))
        }
        colnames(to.plot) <- c(str_replace(gene1_symbol, '-', '_'), str_replace(gene2_symbol, '-', '_'))
        to.plot$cell <- rownames(to.plot)
        to.plot <- merge(to.plot, srat@meta.data, by="cell") %>% as.data.table
        
        if (!(is.null(cluster_subset))) {
            to.plot <- to.plot[cluster %in% cluster_subset]
        }
        if (!(is.null(tp_subset))) {
            to.plot <- to.plot[timepoint %in% tp_subset]
        }
        if (!(is.null(gastr_type_subset))) {
            to.plot <- to.plot[gastr_type %in% gastr_type_subset]
        }
        
        p <- ggplot(data=to.plot, mapping = aes_string(x=str_replace(gene1_symbol, '-', '_'), y=str_replace(gene2_symbol, '-', '_'), colour="cluster")) +
          geom_point(size=0.1, alpha=1) +
          # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
          scale_colour_manual(values = cluster_colours[names(cluster_colours) %in% unique(to.plot$cluster)], name = "Celltype") +
          guides(colour = guide_legend(override.aes = list(size=6))) +
          theme_classic()
    } else if (type == 'counts') {
        if (smoothed) {
            to.plot <- as.matrix(srat@assays[[assay]]@counts[c(gene1_ensid, gene2_ensid),])
            idx <- split(as(t(srat@graphs[[graph_name]]), "dgTMatrix")@i+1, rownames(srat@graphs[[graph_name]])[as(t(srat@graphs[[graph_name]]), "dgTMatrix")@j+1])
            gene_ave <- lapply(idx, function(x) rowMeans(to.plot[,x]))
            to.plot <- do.call(rbind, gene_ave) %>% as.data.frame
        } else {
            to.plot <- as.data.frame(t(srat@assays[[assay]]@counts[c(gene1_ensid, gene2_ensid),]))
        }
        colnames(to.plot) <- c(str_replace(gene1_symbol, '-', '_'), str_replace(gene2_symbol, '-', '_'))
        to.plot$cell <- rownames(to.plot)
        to.plot <- merge(to.plot, srat@meta.data, by="cell") %>% as.data.table
                               
        
        if (!(is.null(cluster_subset))) {
            to.plot <- to.plot[cluster %in% cluster_subset]
        }
        if (!(is.null(tp_subset))) {
            to.plot <- to.plot[timepoint %in% tp_subset]
        }
        if (!(is.null(gastr_type_subset))) {
            to.plot <- to.plot[gastr_type %in% gastr_type_subset]
        }
        
        p <- ggplot(data=to.plot, mapping = aes_string(x=str_replace(gene1_symbol, '-', '_'), y=str_replace(gene2_symbol, '-', '_'), colour="cluster")) +
          geom_point(size=0.1, alpha=1) +
          # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
          scale_colour_manual(values = cluster_colours[names(cluster_colours) %in% unique(to.plot$cluster)], name = "Celltype") +
          guides(colour = guide_legend(override.aes = list(size=6))) +
          theme_classic()
    } else if (type == 'log2') {
        if (smoothed) {
            to.plot <- as.matrix(srat@assays[[assay]]@counts[c(gene1_ensid, gene2_ensid),])
            idx <- split(as(t(srat@graphs[[graph_name]]), "dgTMatrix")@i+1, rownames(srat@graphs[[graph_name]])[as(t(srat@graphs[[graph_name]]), "dgTMatrix")@j+1])
            gene_ave <- lapply(idx, function(x) rowMeans(to.plot[,x]))
            to.plot <- do.call(rbind, gene_ave) %>% as.data.frame
        } else {
            to.plot <- as.data.frame(t(srat@assays[[assay]]@counts[c(gene1_ensid, gene2_ensid),]))
        }
        colnames(to.plot) <- c(str_replace(gene1_symbol, '-', '_'), str_replace(gene2_symbol, '-', '_'))
        to.plot$cell <- rownames(to.plot)
        to.plot <- merge(to.plot, srat@meta.data, by="cell") %>% as.data.table
        
        if (!(is.null(cluster_subset))) {
            to.plot <- to.plot[cluster %in% cluster_subset]
        }
        if (!(is.null(tp_subset))) {
            to.plot <- to.plot[timepoint %in% tp_subset]
        }
        if (!(is.null(gastr_type_subset))) {
            to.plot <- to.plot[gastr_type %in% gastr_type_subset]
        }
        
        p <- ggplot(data=to.plot, mapping = aes_string(x=str_replace(gene1_symbol, '-', '_'), y=str_replace(gene2_symbol, '-', '_'), colour="cluster")) +
          geom_point(size=0.1, alpha=1) +
          # ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
          scale_colour_manual(values = cluster_colours[names(cluster_colours) %in% unique(to.plot$cluster)], name = "Celltype") +
          guides(colour = guide_legend(override.aes = list(size=6))) +
          scale_x_continuous(trans='log2') +
          scale_y_continuous(trans='log2') +
          theme_classic()
    }
    
    return(p)
}

# NOTE: FOR THIS PLOT DEFAULT TYPE IS 'log2' BECAUSE I FIND IT WORKS BETTER
#fill_by: should be NULL or gastr_type, or timepoint, other columns may work but also could throw errors
#plot_type: one of violin, boxplot, or violin_w_boxplot (default)
plotGeneViolins <- function(srat, gene_ensid=NULL, gene_symbol=NULL, smoothed=TRUE, type='log2', assay='RNA', graph_name='RNA_nn', cluster_subset = NULL, tp_subset = NULL, gastr_type_subset = NULL, fill_by=NULL, plot_type='violin_w_boxplot', xaxis='cluster') {
    
    if (is.null(gene_ensid) & is.null(gene_symbol)) {
        stop("either ensembl id or symbol of gene must be given")
    }
    if (!(type %in% c('normalised', 'counts', 'log2'))) {
        stop("type must be one of normalised, counts, or log2")
    }
    if ((smoothed) & is.null(srat@graphs[[graph_name]])) {
        stop("graph_name not among srat graphs (did you remember to first FindNeighbors?")
    }
    
    if (is.null(gene_symbol)) {
        if (!(gene_ensid %in% srat@assays[[assay]]@meta.features$ens_id)) {
            stop("can't find gene_ensid in meta.features ens_id column")
        }
        gene_symbol <- srat@assays[[assay]]@meta.features$symbol[srat@assays[[assay]]@meta.features$ens_id == gene_ensid][1]
    }
    if (is.null(gene_ensid)) {
        if (!(gene_symbol %in% srat@assays[[assay]]@meta.features$symbol)) {
            stop("can't find gene_symbol in meta.features symbol column")
        }
        gene_ensid <- srat@assays[[assay]]@meta.features$ens_id[srat@assays[[assay]]@meta.features$symbol == gene_symbol][1]
    }
    
    if (type == 'normalised') {
        if (smoothed) {
            to.plot <- as.matrix(srat@assays[[assay]]@data[c(gene_ensid),])
            idx <- split(as(t(srat@graphs[[graph_name]]), "dgTMatrix")@i+1, rownames(srat@graphs[[graph_name]])[as(t(srat@graphs[[graph_name]]), "dgTMatrix")@j+1])
            gene_ave <- lapply(idx, function(x) mean(to.plot[x]))
            to.plot <- do.call(rbind, gene_ave) %>% as.data.frame
        } else {
            to.plot <- as.data.frame(srat@assays[[assay]]@data[c(gene_ensid),])
            colnames(to.plot) <- c(str_replace(gene_symbol, '-', '_'))
        }
        
    } else if (type %in% c('counts', 'log2')) {
        if (smoothed) {
            to.plot <- as.matrix(srat@assays[[assay]]@counts[c(gene_ensid),])
            idx <- split(as(t(srat@graphs[[graph_name]]), "dgTMatrix")@i+1, rownames(srat@graphs[[graph_name]])[as(t(srat@graphs[[graph_name]]), "dgTMatrix")@j+1])
            gene_ave <- lapply(idx, function(x) mean(to.plot[x]))
            to.plot <- do.call(rbind, gene_ave) %>% as.data.frame
        } else {
            to.plot <- as.data.frame(srat@assays[[assay]]@counts[c(gene_ensid),])
            colnames(to.plot) <- c(str_replace(gene_symbol, '-', '_'))
        }
        
    }
    
    colnames(to.plot) <- c(str_replace(gene_symbol, '-', '_'))
    to.plot$cell <- rownames(to.plot)
    to.plot <- merge(to.plot, srat@meta.data, by="cell") %>% as.data.table

    if (!(is.null(cluster_subset))) {
        to.plot <- to.plot[cluster %in% cluster_subset]
    }
    if (!(is.null(tp_subset))) {
        to.plot <- to.plot[timepoint %in% tp_subset]
    }
    if (!(is.null(gastr_type_subset))) {
        to.plot <- to.plot[gastr_type %in% gastr_type_subset]
    }
    
    to.plot <- to.plot[cluster %in% cluster_order]
    to.plot$cluster <- factor(to.plot$cluster, levels = cluster_order)
    
    if (is.null(fill_by)) {
        p <- ggplot(data=to.plot[cluster %in% cluster_order], mapping = aes_string(x=xaxis, y=str_replace(gene_symbol, '-', '_'), fill="cluster"))
        if (plot_type == 'violin') {
            p <- p + geom_violin(trim=FALSE)
        } else if (plot_type == 'violin_w_boxplot') {
            p <- p +
              geom_violin(trim=FALSE) +
              geom_boxplot(width=0.1, outlier.size=0.2, fill="white")
        } else if (plot_type == 'boxplot') {
            p <- p + geom_boxplot(outlier.size=0.2)
        } else {
            stop("unknown plot_type")
        }
        p <- p + scale_fill_manual(values = cluster_colours[names(cluster_colours) %in% unique(to.plot$cluster)], name = "Celltype")
    } else if (fill_by == 'gastr_type') {
        to.plot$gastr_type <- factor(to.plot$gastr_type, levels = gastr_type_order)
        p <- ggplot(data=to.plot[gastr_type %in% gastr_type_order], mapping = aes_string(x=xaxis, y=str_replace(gene_symbol, '-', '_'), fill="gastr_type"))
        if (plot_type == 'violin') {
            p <- p + geom_violin(trim=FALSE)
        } else if (plot_type == 'violin_w_boxplot') {
            p <- p +
              geom_violin(trim=FALSE) +
              geom_boxplot(width=0.1, outlier.size=0.2, fill="white")
        } else if (plot_type == 'boxplot') {
            p <- p + geom_boxplot(outlier.size=0.2)
        } else {
            stop("unknown plot_type")
        }
        p <- p + scale_fill_manual(values = gastr_type_colours[names(gastr_type_colours) %in% unique(to.plot$gastr_type)], name = "Gastruloid Type") +
          guides(colour = guide_legend(override.aes = list(size=6)))
    } else if (fill_by == 'timepoint') {
        to.plot$timepoint <- factor(to.plot$timepoint, levels = timepoint_order)
        p <- ggplot(data=to.plot[timepoint %in% timepoint_order], mapping = aes_string(x=xaxis, y=str_replace(gene_symbol, '-', '_'), fill="timepoint"))
        if (plot_type == 'violin') {
            p <- p + geom_violin(trim=FALSE)
        } else if (plot_type == 'violin_w_boxplot') {
            p <- p +
              geom_violin(trim=FALSE) +
              geom_boxplot(width=0.1, outlier.size=0.2, fill="white")
        } else if (plot_type == 'boxplot') {
            p <- p + geom_boxplot(outlier.size=0.2)
        } else {
            stop("unknown plot_type")
        }
        p <- p + scale_fill_manual(values = timepoint_colours[names(timepoint_colours) %in% unique(to.plot$timepoint)], name = "Timepoint") +
          guides(colour = guide_legend(override.aes = list(size=6)))
    } else {
        p <- ggplot(data=to.plot, mapping = aes_string(x=xaxis, y=str_replace(gene_symbol, '-', '_'), fill=fill_by))
        if (plot_type == 'violin') {
            p <- p + geom_violin(trim=FALSE)
        } else if (plot_type == 'violin_w_boxplot') {
            p <- p +
              geom_violin(trim=FALSE) +
              geom_boxplot(width=0.1, outlier.size=0.2, fill="white")
        } else if (plot_type == 'boxplot') {
            p <- p + geom_boxplot(outlier.size=0.2)
        } else {
            stop("unknown plot_type")
        }
        p <- p + guides(colour = guide_legend(override.aes = list(size=6)))
    }
    
    if (type == 'log2') {
        p <- p + scale_y_continuous(trans='log2')
    }
    
    p <- p +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    if (is.null(fill_by)) {
        p <- p + theme(legend.position = "None")
    }
    
    return(p)
}