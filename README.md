# gastruloids2022
Code to reproduce the analyses from the 2022 gastruloid timecourse paper. For the intial processing there is a snakemake pipeline. Beyond that, there is a lot of code in the make_plots.ipynb notebook.

In general, please be aware that the code was designed to run on EBI's codon cluster, so some file directories/conda environments may be hard coded into the scripts, so please do check all code before running.

In general, most scripts call ./plotting_settings.R which involves colours used for plots, the categorisation of each gastruloid into mesodermal, neural, or intermediate, as well as the order in which the gastruloids, cell types, timepoints, and gastruloid types are plot. Additionally, it includes some functions for exploratory plotting of the data.

To explore the data, I recommend downloading `seurat_final.rds` from GEO and using the `explore_data.ipynb`.

## Snakemake Pipeline:
There is a snakemake pipeline to load the data, run the per-sample QC and create Seurat files, align the MULTI-seq barcodes, do MULTI-seq demultiplexing, do the joint QC, and map to both the original and extended atlases. The files involved in this are:
- ./Snakefile
- ./config.yaml
- ./processing
    - ./processing/merge.R this merges all the Seurat objects and MULTI-seq classifications for all samples
    - ./processing/create_seurat.R this creates a separate Seurat for each sample
    - ./processing/classify_cells.R this is part of the MULTI-seq demultiplexing pipeline, and is run separately for each sample
    - ./processing/barcode_align.R this is part of the MULTI-seq demultiplexing pipeline, and is run separately for each sample
    - ./processing/settings.R
- ./mapping
    - ./mapping/mapping_mnn_extendedAtlas.R
    - ./mapping/mapping_mnn_originalAtlas.R
    - ./mapping/settings.R
    - ./mapping/mapping_functions.R

Furthermore, the following scripts are not part of the pipeline, but relevant to the initial data processing:
- ./processing/exp1_d5_DoubletBarcodes.ipynb this notebook explores the doublets in the first MULTI-seq experiment
- ./processing/select_MULTIseq_barcodes.ipynb this notebook selects the optimal barcodes that are certain distance metrics apart in sequence and reverse complement
- ./processing/check_MULTIseq_assignment is a directory of R jupyter notebooks where I chose the optimal reclassification stability for each sample, to recover as many MULTI-seq assignments as possible.

Please note that the results of the reference-based celltype annotation are merged and added to files in the ./annotate_cells.ipynb as well as the embryonic day assignment.

## QC Plots:
The QC plots that are in the paper are in the big make_plots.ipynb along with code to normalise across samples. The calculation is performed in ./annotate_cells.ipynb

## Annotating cell types:
Cell types were annotated in the following jupyter notebook: ./annotate_cells.ipynb. This notebook also takes the reference-based celltype annotations, merges them, and adds them to the files, as well as the trajectory inference annotations (this is a bit circular when first running the code), and does the embryonic day assignment. This notebook also creates Figure 1D,E, and Supplementary Figure 3A.

## Assigning equivalent embryonic days:
The equivalent embryonic days were annotated and explored in the following jupyter notebook: ./annotate_cells.ipynb

## Comparing with other gastruloid data:
Here you will find comparisons with the [van den Brink et al. 2020](https://www.nature.com/articles/s41586-020-2024-3) data (in paper), but also (not in paper) integration with the [Veenvliet et al. 2020](https://www.science.org/doi/full/10.1126/science.aba4937#) and [Rossi et al. 2020](https://doi.org/10.1016/j.stem.2020.10.013) datasets.

## Trajectory Inference
First, velocyto was run to obtain spliced and unspliced count matrices. All scripts for this are in ./TrajectoryInference/run_velocyto and are run in the following order:
1. ./TrajectoryInference/run_velocyto/run_all_split_by_chrom.sh this script calls ./TrajectoryInference/run_velocyto/run_split_by_chrom.sh which in turn calls ./TrajectoryInference/run_velocyto/split_by_chrom.sh For each experiment, this takes the outs/possorted_genome_bam.bam and splits it by chromosome.
2. ./TrajectoryInference/run_velocyto/run_all_velocyto.sh which calls run_velocyto.sh which runs the actual velocyto for each experiment and each chromosome separately, and does require the gene annotation (references/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf) and masked areas of the genome (mm10_rmsk.gtf)
3. run_merge_chroms.sh calls merge_chroms.py which merges the separate chromosomes of the velocyto outputs for each experiment
This produces a loom for each experiment. These are then combined into a single loom with the following code:
```
files = [y for x in os.walk(indir) for y in glob(os.path.join(x[0], 'all.loom')) if not 'human' in y]
ldata_sub = {}
for f in files:
    
    ldata_sub[f] = scv.read(f, cache=True)

ldata = anndata.concat([ldata_sub[f] for f in files], axis=0, join='outer', merge="first")

ldata.write_loom(outfile)
```

Trajectory inference was done using CellRank, for mesodermal and neural gastruloids separately.
Notebooks running CellRank
./TrajectoryInference/CellRank_mesogastr.ipynb
./TrajectoryInference/CellRank_neurogastr.ipynb

For the vector plot in figure 1, the following notebook merges the 2 gastruloid types
./TrajectoryInference/merge_meso_neural.ipynb

The following 2 notebooks generate the figure 2 WOT flows.
./TrajectoryInference/plot_WOT_meso.ipynb
./TrajectoryInference/plot_WOT_neuro.ipynb

## PCA on gastruloid cell type proportions
This is in the big make_plots.ipynb

## Gastruloid linear order
This was determined in the big make_plots.ipynb based on fitting a principle curve

## Combining cell types into broad lineages
This is in the big make_plots.ipynb

## Calculating the AICc
This is in the big make_plots.ipynb

## Differential abundance testing between gastruloid classes
This is in the big make_plots.ipynb

## Inferring signalling in day 3 gastruloids
First, the gastruloids are pseudobulked in ./signalling_analysis/pseudobulk.R including creating random pseudobulks.

This approach is based on [Barker et al. 2020](https://genome.cshlp.org/content/32/4/750). To run my modification of his method, use ./signalling_analysis/run_Charlie_method.R. This runs scWGCNA, selects significant modules, and his signalling inference. However, I only use the WGCNA part of the output ("[day]_[experiment]_WGCNA_data_[date].rds").

This is then input to ./signalling_analysis/d3_improved_analysis.ipynb where the GO enrichment analysis is performed.

## Comparing to [Gouti 2014](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001937) bulk, 2D differentiation data
- ./signalling_analysis/Gouti2014/run_tophat_persample.sh runs the tophat alignment
- ./signalling_analysis/Gouti2014/run_htseq_persample.sh runs the htseq alignment
- ./signalling_analysis/Gouti2014/merge_counts.ipynb merges the resulting count tables
- ./signalling_analysis/Gouti2014/RunComparison.ipynb runs DEseq2 and Limma to compare the datasets.

## Comparing day 4.5 and day 5 gastruloid PCA distances
This is in the big make_plots.ipynb

## Running perturbation simulations
The simulations are run from ./perturbation_simulations/run_DA_simulations.sh and subsequently loaded using ./perturbation_simulations/load_DA_simulation_files.sh which calls ./perturbation_simulations/load_ct_DA_files.R please be warned that running the full simulations requires 1.9 million 5GB jobs. This will make any other user of your cluster sad.

All simulations were done on d4.5 gastruloids using the "full" simulation type, an contribution IoD of 500 and 100 trials.

Each individual simulation is done by ./perturbation_simulations/DA_simulations_gastruloid.R which loads the ./perturbation_simulations/DA_functions.py (the deconvolution python functions) and ./perturbation_simulations/DA_simulations_functions.R (all other functions, in R).
This takes as input an already-calculated PCA space signal (it isn't recalculated each time), which is calculated using: ./perturbation_simulations/Run_PCA.ipynb

## Scripts not directly contributing to figures in the paper
- ./DE These scripts were used for data exploration, and there outputs are explorable in the shinyapp. `run_all_DEG_gastrtype_tp.sh` does the calling (though submits many, many jobs), by calling `run_DEG_gastrtype_tp.R` for each cell type, gastruloid type, and time point combination. It gives it as input the `DEG_functions.R`. The results are then loaded and merged using `load_DE_results.sh` which goes cell type by cell type and calling `load_DE_results.R` to load the cell type - cell type combination results.
- ./explore_data.ipynb this notebook contains example code to use the functions in plotting_settings.R to explore the data (especially gene expression). I hope it's useful!

## Where to find all paper figures:
All figures in make_plots.ipynb unless otherwise stated

Figure 1
- D: `./annotate_cells.ipynb`
- E: `./annotate_cells.ipynb`
- F: `./TrajectoryInference/merge_meso_neural.ipynb`

Supplementary Figure 2
- A: `./other_gastruloid_scRNAseq/vdbrink.ipynb`
- B: `./other_gastruloid_scRNAseq/vdbrink.ipynb`
- C: `./other_gastruloid_scRNAseq/vdbrink.ipynb`
- D: `./other_gastruloid_scRNAseq/vdbrink.ipynb`

Supplementary Figure 3
- A: `./annotate_cells.ipynb`

Figure 2
- F: `./TrajectoryInference/plot_WOT_meso.ipynb`
- G: `./TrajectoryInference/plot_WOT_neuro.ipynb`

Supplementary Figure 4
- D: `./processing/exp1_d5_DoubletBarcodes.ipynb`
- F: `./processing/exp1_d5_DoubletBarcodes.ipynb`

Supplementary Figure 8
- A: `./TrajectoryInference/CellRank_mesogastr.ipynb`
- B: `./TrajectoryInference/CellRank_neurogastr.ipynb`
- C: `./TrajectoryInference/plot_WOT_meso.ipynb`
- D: `./TrajectoryInference/plot_WOT_neuro.ipynb`

Figure 3
- A: `./signalling_analysis/d3_improved_analysis.ipynb`
- B: `./signalling_analysis/d3_improved_analysis.ipynb`
- C: `./signalling_analysis/Gouti2014/RunComparison.ipynb`

Figure 4
- C: `./perturbation_simulations/makeplots.ipynb`
- D: `./perturbation_simulations/makeplots.ipynb`
- E: `./perturbation_simulations/makeplots.ipynb`

Supplementary Figure 8
- B: `./perturbation_simulations/makeplots.ipynb`














