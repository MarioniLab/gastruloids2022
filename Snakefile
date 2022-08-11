# bsub -M 500000 -q bigmem -n 5 snakemake --cores 5 --latency-wait 600 -k --resources mem_mb=500000

configfile: "config.yaml"

rule all:
    input:
        config["directories"]["base"] + config["directories"]["output"] + "mapping/all_extended.rds",
        config["directories"]["base"] + config["directories"]["output"] + "mapping/all_original.rds",
        expand(config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_extended.rds", batch=config["samples"].keys()),
        expand(config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_original.rds", batch=config["samples"].keys()),
        #config["directories"]["base"] + config["directories"]["output"] + "mapping/all_umap.rds",
        #expand(config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_umap.rds", batch=config["batches"])
        #config["directories"]["base"] + config["directories"]["output"] + "seurat_QCed.rds"
        

rule CreateSeurat:
    input:
        script=config["scripts"]["CreateSeurat"],
        sets=config["settings"]["Processing"],
        dir=config["directories"]["base"] + config["directories"]["extension"] + "{GEXexperiment}"
    output:
        seurat=config["directories"]["base"] + config["directories"]["output"] + "{GEXexperiment}/{GEXsample}seurat.rds",
        metadata=config["directories"]["base"] + config["directories"]["output"] + "{GEXexperiment}/{GEXsample}sample_metadata.txt.gz",
        featmetadata=config["directories"]["base"] + config["directories"]["output"] + "{GEXexperiment}/{GEXsample}feature_metadata.txt.gz",
        barcodes=config["directories"]["base"] + config["directories"]["output"] + "{GEXexperiment}/{GEXsample}cleaned_barcodes.tsv.gz"
    wildcard_constraints:
        GEXsample=".*"
    params:
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    resources:
        mem_mb=50000
    run:
        if wildcards.GEXsample:
            sample = "-s {wildcards.GEXsample}"
        else:
            sample = ""
        shell("{params.cmd} {input.script} -i {input.dir} -S {input.sets} -O {output.seurat} -o {output.metadata} -g {output.featmetadata} -b {output.barcodes} -e {wildcards.GEXexperiment} "+sample)


rule MULTIBarAlign:
    input:
        script=config["scripts"]["MULTIBarAlign"],
        sets=config["settings"]["Processing"],
        dir=config["directories"]["base"] + config["directories"]["extension"] + "{experiment}",
        #barcodes=rules.CreateSeurat.output.barcodes
        barcodes=config["directories"]["base"] + config["directories"]["output"] + "{experiment}/{sample}cleaned_barcodes.tsv.gz"
    output:
        barTable=config["directories"]["base"] + config["directories"]["output"] + "{experiment}/{sample}barTable.rds",
        readTable=config["directories"]["base"] + config["directories"]["output"] + "{experiment}/{sample}readTable.rds"
    params:
        outbase = config["directories"]["base"] + config["directories"]["output"] + "{experiment}/{sample}",
        cmd       = f"conda run -n " + config["environments"]["MULTIseq"] + " Rscript"
    resources:
        mem_mb=100000
    shell:
        "{params.cmd} {input.script} -i {input.dir} -b {input.barcodes} -S {input.sets} -O {params.outbase} -e {wildcards.experiment} -s {wildcards.sample}"


rule MULTIAssignCells:
    input:
        script=config["scripts"]["MULTIAssignCells"],
        sets=config["settings"]["Processing"],
        barTable=rules.MULTIBarAlign.output.barTable
    output:
        finalCalls=config["directories"]["base"] + config["directories"]["output"] + "{experiment}/{sample}finalCalls.rds"
    params:
        dir = config["directories"]["base"] + config["directories"]["output"],
        cmd       = f"conda run -n " + config["environments"]["MULTIseq"] + " Rscript"
    resources:
        mem_mb=25000
    shell:
        "{params.cmd} {input.script} -i {params.dir} -S {input.sets} -e {wildcards.experiment} -s {wildcards.sample}"

rule MergePreprocessed:
    input:
        script=config["scripts"]["MergePreprocessed"],
        seurat_files=[config["directories"]["base"] + config["directories"]["output"] + "{GEXexperiment}/{GEXsample}seurat.rds".format(GEXsample=GEXsample, GEXexperiment=GEXexperiment) for GEXexperiment in config["samples"].keys() for GEXsample in config["samples"][GEXexperiment]],
        genemeta_files=[config["directories"]["base"] + config["directories"]["output"] + "{GEXexperiment}/{GEXsample}feature_metadata.txt.gz".format(GEXsample=GEXsample, GEXexperiment=GEXexperiment) for GEXexperiment in config["samples"].keys() for GEXsample in config["samples"][GEXexperiment]],
        final_calls_files=[config["directories"]["base"] + config["directories"]["output"] + "{experiment}/{sample}finalCalls.rds".format(sample=sample, experiment=experiment) for experiment in config["samples_MULTI"].keys() for sample in config["samples_MULTI"][experiment]],
        sets=config["settings"]["Processing"]
    output:
        seurat=config["directories"]["base"] + config["directories"]["output"] + "seurat.rds",
        genemeta=config["directories"]["base"] + config["directories"]["output"] + "feature_metadata.txt.gz",
        metadata=config["directories"]["base"] + config["directories"]["output"] + "sample_metadata.txt.gz"
    params:
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    resources:
        mem_mb=100000
    run:
        Sfiles = ','.join(input.seurat_files)
        Ffiles = ','.join(input.final_calls_files)
        Gfiles = ','.join(input.genemeta_files)
        shell("{params.cmd} {input.script} -S "+Sfiles+" -F "+Ffiles+" -G "+Gfiles+" -s {input.sets} -O {output.seurat} -o {output.metadata} -g {output.genemeta}")


rule jointQC:
    input:
        script=config["scripts"]["JointQC"],
        seurat=rules.MergePreprocessed.output.seurat,
        meta=rules.MergePreprocessed.output.metadata,
        genemeta=rules.MergePreprocessed.output.genemeta,
        sets=config["settings"]["Mapping"],
        fns=config["settings"]["MappingFunctions"]
    output:
        seurat=config["directories"]["base"] + config["directories"]["output"] + "seurat_QCed.rds",
        #anndata=config["directories"]["base"] + config["directories"]["output"] + "seurat_QCed.h5ad",
        genemeta=config["directories"]["base"] + config["directories"]["output"] + "feature_metadata_QCed.txt.gz",
        metadata=config["directories"]["base"] + config["directories"]["output"] + "sample_metadata_QCed.txt.gz"
    params:
        dir = config["directories"]["base"] + config["directories"]["output"] + "processing/jointQC/",
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    resources:
        mem_mb=100000
    shell:
        #"mkdir -p {params.dir} && {params.cmd} {input.script} -I {input.seurat} -i {input.meta} -G {input.genemeta} -S {input.sets} -f {input.fns} -O {output.seurat} -o {output.metadata} -g {output.genemeta} -p {params.dir} -h {output.anndata}"
        "mkdir -p {params.dir} && {params.cmd} {input.script} -I {input.seurat} -i {input.meta} -G {input.genemeta} -S {input.sets} -f {input.fns} -O {output.seurat} -o {output.metadata} -g {output.genemeta} -p {params.dir}"


rule MappingAll_ExtendedAtlas:
    input:
        script=config["scripts"]["Mapping"]["Extended"],
        query=rules.jointQC.output.seurat,
        meta=rules.jointQC.output.metadata,
        genemeta=rules.jointQC.output.genemeta,
        sets=config["settings"]["Mapping"],
        fns=config["settings"]["MappingFunctions"],
        atlasrds=config["directories"]["atlas"]["Extended"]+"sce_atlas.rds",
        atlasmd=config["directories"]["atlas"]["Extended"]+"atlas_metadata.txt.gz"
    output:
        rds=config["directories"]["base"] + config["directories"]["output"] + "mapping/all_extended.rds",
        txtgz=config["directories"]["base"] + config["directories"]["output"] + "mapping/all_extended.txt.gz",
        genemeta=config["directories"]["base"] + config["directories"]["output"] + "mapping/all_extended_feature_metadata.txt.gz"
    params:
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    resources:
        mem_mb=500000
    shell:
        "{params.cmd} {input.script} -A {input.atlasrds} -a {input.atlasmd} -Q {input.query} -q {input.meta} -G {input.genemeta} -S {input.sets} -f {input.fns} -O {output.rds} -o {output.txtgz} -g {output.genemeta}"


rule MappingBatch_ExtendedAtlas:
    input:
        script=config["scripts"]["Mapping"]["Extended"],
        query=rules.jointQC.output.seurat,
        genemeta=rules.jointQC.output.genemeta,
        meta=rules.jointQC.output.metadata,
        sets=config["settings"]["Mapping"],
        fns=config["settings"]["MappingFunctions"],
        atlasrds=config["directories"]["atlas"]["Extended"]+"sce_atlas.rds",
        atlasmd=config["directories"]["atlas"]["Extended"]+"atlas_metadata.txt.gz"
    output:
        rds=config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_extended.rds",
        txtgz=config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_extended.txt.gz",
        genemeta=config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_extended_feature_metadata.txt.gz"
    params:
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    resources:
        mem_mb=500000
    shell:
        "{params.cmd} {input.script} -A {input.atlasrds} -a {input.atlasmd} -Q {input.query} -q {input.meta} -G {input.genemeta} -S {input.sets} -f {input.fns} -b {wildcards.batch} -O {output.rds} -o {output.txtgz} -g {output.genemeta}"


rule MappingAll_OriginalAtlas:
    input:
        script=config["scripts"]["Mapping"]["Original"],
        query=rules.jointQC.output.seurat,
        meta=rules.jointQC.output.metadata,
        genemeta=rules.jointQC.output.genemeta,
        sets=config["settings"]["Mapping"],
        fns=config["settings"]["MappingFunctions"],
        atlasrds=config["directories"]["atlas"]["Original"]+"sce_atlas.rds",
        atlasmd=config["directories"]["atlas"]["Original"]+"atlas_metadata.txt.gz"
    output:
        rds=config["directories"]["base"] + config["directories"]["output"] + "mapping/all_original.rds",
        txtgz=config["directories"]["base"] + config["directories"]["output"] + "mapping/all_original.txt.gz",
        genemeta=config["directories"]["base"] + config["directories"]["output"] + "mapping/all_original_feature_metadata.txt.gz"
    params:
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    resources:
        mem_mb=100000
    shell:
        "{params.cmd} {input.script} -A {input.atlasrds} -a {input.atlasmd} -Q {input.query} -q {input.meta} -G {input.genemeta} -S {input.sets} -f {input.fns} -O {output.rds} -o {output.txtgz} -g {output.genemeta}"


rule MappingBatch_OriginalAtlas:
    input:
        script=config["scripts"]["Mapping"]["Original"],
        query=rules.jointQC.output.seurat,
        meta=rules.jointQC.output.metadata,
        genemeta=rules.jointQC.output.genemeta,
        sets=config["settings"]["Mapping"],
        fns=config["settings"]["MappingFunctions"],
        atlasrds=config["directories"]["atlas"]["Original"]+"sce_atlas.rds",
        atlasmd=config["directories"]["atlas"]["Original"]+"atlas_metadata.txt.gz"
    output:
        rds=config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_original.rds",
        txtgz=config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_original.txt.gz",
        genemeta=config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_original_feature_metadata.txt.gz"
    params:
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    resources:
        mem_mb=100000
    shell:
        "{params.cmd} {input.script} -A {input.atlasrds} -a {input.atlasmd} -Q {input.query} -q {input.meta} -G {input.genemeta} -S {input.sets} -f {input.fns} -b {wildcards.batch} -O {output.rds} -o {output.txtgz} -g {output.genemeta}"


'''
rule PlottingBatch:
    input:
        script=config["scripts"]["Plotting"],
        rds=rules.MappingBatch.output.rds,
        meta=rules.MappingBatch.output.txtgz,
        sets=config["settings"]["Plotting"],
        atlasrds=config["directories"]["atlas"]+"sce_atlas.rds",
        atlasmd=config["directories"]["atlas"]+"atlas_metadata.txt.gz"
    output:
        umap=config["directories"]["base"] + config["directories"]["output"] + "mapping/{batch}_umap.rds"
    params:
        base=config["directories"]["base"] + config["directories"]["output"] + "plots/mapping/{batch}\/",
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    shell:
        "mkdir -p {params.base} && {params.cmd} {input.script} -R {input.rds} -m {input.meta} -a {input.atlasmd} -A {input.atlasrds} -S {input.sets} -u {output.umap} -O {params.base}"


rule PlottingAll:
    input:
        script=config["scripts"]["Plotting"],
        rds=rules.MappingAll.output.rds,
        meta=rules.MergePreprocessed.output.metadata,
        sets=config["settings"]["Plotting"],
        atlasrds=config["directories"]["atlas"]+"sce_atlas.rds",
        atlasmd=config["directories"]["atlas"]+"atlas_metadata.txt.gz"
    output:
        umap=config["directories"]["base"] + config["directories"]["output"] + "mapping/all_umap.rds"
    params:
        base=config["directories"]["base"] + config["directories"]["output"] + "plots/mapping/all\/",
        cmd       = f"conda run -n " + config["environments"]["CreateSeurat"] + " Rscript"
    shell:
        "mkdir -p {params.base} && {params.cmd} {input.script} -R {input.rds} -m {input.meta} -a {input.atlasmd} -A {input.atlasrds} -S {input.sets} -u {output.umap} -O {params.base}"
'''






























