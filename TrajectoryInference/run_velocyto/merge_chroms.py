#!/usr/bin/env python3

# run in scvelo environment using bsub -M 100000 /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/mapping_eval/scvelo/merge_chroms.py -i /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/cellranger_outs/MULTI/exp4_d3/SITTA12_gastr_d3_MULTIseq_3prime_RNA_A/velocyto/ -m /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/MULTI/exp4_d3/sampleA_sample_metadata.txt.gz

import argparse
import os
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import copy
import anndata

parser = argparse.ArgumentParser(description='Merge individual loom files that were output by velocity. Also saves modified loom files, whose obs and var have been modified to align with my naming conventions.')
parser.add_argument('--indir', '-i', type=str, required=True,
                    help='the directory into which velocyto wrote its output')
parser.add_argument('--metadata', '-m', type=str, required=True,
                    help='the metadata of the sample')

args = parser.parse_args()

md = pd.read_csv(args.metadata,
                 sep='\t'
                )

files = [ os.path.join(args.indir, name, name2) for name in os.listdir(args.indir) if os.path.isdir(os.path.join(args.indir, name)) for name2 in os.listdir(os.path.join(args.indir, name))]
files_working = copy.deepcopy(files)
ldata_sub = {}
for f in files_working:
    
    ldata_sub[f] = scv.read(f, cache=True)
    
    if (sum(np.array(ldata_sub[f].X.sum(axis=0))[0]) == 0):
        files.remove(f)
        continue
    
    s = ldata_sub[f].obs_names[0]
    start = s.find(":") + len(":")
    end = s.find("x")
    ldata_sub[f].obs_names = pd.Index([on[start:end] + '-1' for on in ldata_sub[f].obs_names], name='CellID')
    if len(ldata_sub[f].obs_names[0]) != 18:
        print('error: barcode length not 18')
        break
    
    md_sub = md.loc[md.barcode.isin(ldata_sub[f].obs_names),]
    ldata_sub[f] = ldata_sub[f][md_sub.barcode]
    ldata_sub[f].obs['barcode'] = ldata_sub[f].obs_names
    ldata_sub[f].obs['cell'] = [ md_sub.cell[md_sub.barcode == bc].values[0] for bc in ldata_sub[f].obs.barcode ]
    ldata_sub[f].obs_names = pd.Index(ldata_sub[f].obs.cell, name='CellID')
    ldata_sub[f] = ldata_sub[f][:,np.array(ldata_sub[f].X.sum(axis=0))[0] > 0]
    
    ldata_sub[f].var["symbol"] = ldata_sub[f].var_names.values
    ldata_sub[f].var_names = pd.Index(ldata_sub[f].var["Accession"].values, name='Gene')
    
    chrom = os.path.dirname(f).split('/')[-1]
    ldata_sub[f].write_loom(args.indir + chrom + '.loom')

ldata = anndata.concat([ldata_sub[f] for f in files], axis=1, join='outer')
ldata.write_loom(args.indir + 'all.loom')