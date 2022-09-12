# I'm testing if there's a huge difference if i aggregate at gene level using 
# -g flag in `salmon quant -g` or if i aggregate them in R using tximport as suggested.
# There'd be multiple rules for that untill i decided which one to choose

rule create_tables_g:
    input:
        expand("results/02_salmon_g/{sample}/quant.genes.sf", sample = SAMPLES)
    output:
        tpm         = "results/04deseq2_g/tpm.tsv",
        fpkm        = "results/04deseq2_g/fpkm.tsv",
        raw_counts  = "results/04deseq2_g/Raw_counts.tsv"
    params:
        sample_names  = expand("{sample}", sample = SAMPLES),
        exclude       = config["diffexp"].get("exclude", None)
    log:
        "results/00log/deseq2_g/create_tables.log"
    script:
        "../scripts/createTables_count_rpkm_tpm.R"


rule deseq2_g:
    input:
        expand("results/02_salmon_g/{sample}/quant.sf", sample = SAMPLES)
    output:
        rds         	= "results/04deseq2_g/all.rds",
        norm_counts 	= "results/04deseq2_g/Normalized_counts.tsv",
    params:
        sample_names  	= expand("{sample}", sample = SAMPLES),
        samples  		= config["samples"],
        exclude  		= config["diffexp"].get("exclude", None),
        tx2gene  		= config["ref"]["annotation"],
        annot_col 		= config["ref"]["annot_type"]
    log:
        "results/00log/deseq2_g/init.log"
    script:
        "../scripts/DESeq2.R"

rule get_contrasts_g:
    input:
        rds     = rules.deseq2_g.output.rds,
        fpkm    = rules.create_tables_g.output.fpkm
    output:
        table     = "results/04deseq2_g/{contrast}/{contrast}_diffexp.tsv",
        ma_plot   = "results/04deseq2_g/{contrast}/{contrast}_ma-plot.pdf",
        pval_hist = "results/04deseq2_g/{contrast}/{contrast}_pval-hist.pdf",
    params:
        contrast        = lambda w: config["diffexp"]["contrasts"][w.contrast],
        lfcShrink       = config["lfcShrink"],
        samples         = config["samples"],
        exclude         = config["diffexp"].get("exclude", None),
        annot           = config["ref"]["geneInfo"].get("file", None),
        column_used     = config["ref"]["geneInfo"]["column_used"],
        column_toAdd    = config["ref"]["geneInfo"]["column_toAdd"],
        name_annotation = config["ref"]["geneInfo"]["name_annotation"],
    log:
        "results/00log/deseq2_g/{contrast}.diffexp.log"
    script:
        "../scripts/get_DESeq2_contrasts.R"

rule pca:
    input:
        rules.deseq2_g.output.rds
    output:
        "results/04deseq2_g/pca_elipse_names_top{ntop}.pdf",
        "results/04deseq2_g/pca_top{ntop}.pdf",
        "results/04deseq2_g/pca_names_top{ntop}.pdf"
    params:
        pca_labels = config["pca"]["labels"],
    log:
        "results/00log/pca_{ntop}.log"
    script:
        "../scripts/plot-PCA.R"

rule filter_deg_g:
    input:
        diffExp = rules.get_contrasts_g.output.table,
    output:
        "results/04deseq2_g/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_diffexp_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.tsv"
    params:
        pval      = lambda w: w.pvalue,
        log2fc    = lambda w: w.log2fc,
        fpkm_filt = lambda w: w.fpkm
    log:
        "results/00log/deseq2_g/{contrast}.{log2fc}.{pvalue}_fpkm{fpkm}.filter_deg.log"
    script:
        "../scripts/filter_deg.R"

# rule enrichments:
#     input:
#         rules.filter_deg.output
#     output:
#         enrichments = "results/04deseq2/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_enrichments_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.xlsx",
#     params:
#         genome       = config["ref"]["genome"],
#         pvalue       = config["enrichments"]["pval"],
#         qvalue       = config["enrichments"]["qval"],
#         set_universe = config["enrichments"]["set_universe"],
#     log:
#         "results/00log/deseq2/{contrast}.log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.enrichments.log"
#     script:
#         "../scripts/enrichments.R"   


rule volcano_g:
    input:
        rules.filter_deg_g.output
    output:
        volcano_pdf	 = "results/04deseq2_g/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_volcano_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.pdf",
        volcano_png	 = "results/04deseq2_g/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_volcano_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.png",
    params:
        pval     = lambda w: w.pvalue,
        log2fc   = lambda w: w.log2fc,
        contrast = lambda w: w.contrast,
    log:
        "results/00log/deseq2_g/{contrast}.log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.volcano.log"
    script:
        "../scripts/volcano.R"  


############################################################
############# RULES FOR SALMON + TXIMPORT ##################
############################################################

# rule create_tables_tx:
#     input:
#          expand("results/02_salmon_tx/{sample}/quant.sf", sample = SAMPLES)
#     output:
#         tpm         = "results/04deseq2_tx/tpm.tsv",
#         fpkm        = "results/04deseq2_tx/fpkm.tsv",
#         raw_counts  = "results/04deseq2_tx/Raw_counts.tsv"
#     params:
#         exclude = config["diffexp"].get("exclude", None)
#     log:
#         "results/00log/deseq2_tx/create_tables.log"
#     script:
#         "../scripts/createTables_count_rpkm_tpm.R"
 
# rule deseq2_tx:
#     input:
#         rules.create_tables.output.raw_counts
#     output:
#         rds         = "results/04deseq2_tx/all.rds",
#         norm_counts = "results/04deseq2_tx/Normalized_counts.tsv",
#     params:
#         samples  = config["samples"],
#         exclude  = config["diffexp"].get("exclude", None),
#     log:
#         "results/00log/deseq2_tx/init.log"
#     script:
#         "../scripts/DESeq2.R"

# rule get_contrasts_tx:
#     input:
#         rds     = rules.deseq2.output.rds,
#         fpkm    = rules.create_tables.output.fpkm
#     output:
#         table     = "results/04deseq2_tx/{contrast}/{contrast}_diffexp.tsv",
#         ma_plot   = "results/04deseq2_tx/{contrast}/{contrast}_ma-plot.pdf",
#         pval_hist = "results/04deseq2_tx/{contrast}/{contrast}_pval-hist.pdf",
#     params:
#         contrast        = lambda w: config["diffexp"]["contrasts"][w.contrast],
#         lfcShrink       = config["lfcShrink"],
#         samples         = config["samples"],
#         exclude         = config["diffexp"].get("exclude", None),
#         annot           = config["ref"]["geneInfo"].get("file", None),
#         column_used     = config["ref"]["geneInfo"]["column_used"],
#         column_toAdd    = config["ref"]["geneInfo"]["column_toAdd"],
#         name_annotation = config["ref"]["geneInfo"]["name_annotation"],
#     log:
#         "results/00log/deseq2_tx/{contrast}.diffexp.log"
#     script:
#         "../scripts/get_DESeq2_contrasts.R"


# rule filter_deg_g:
#     input:
#         diffExp = rules.get_contrasts.output.table,
#     output:
#         "results/04deseq2_tx/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_diffexp_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.tsv"
#     params:
#         pval      = lambda w: w.pvalue,
#         log2fc    = lambda w: w.log2fc,
#         fpkm_filt = lambda w: w.fpkm
#     log:
#         "results/00log/deseq2_tx/{contrast}.{log2fc}.{pvalue}_fpkm{fpkm}.filter_deg.log"
#     script:
#         "../scripts/filter_deg.R"

# rule volcano_g:
#     input:
#         rules.filter_deg.output
#     output:
#         volcano_pdf	 = "results/04deseq2_tx/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_volcano_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.pdf",
#         volcano_png	 = "results/04deseq2_tx/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_volcano_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.png",
#     params:
#         pval     = lambda w: w.pvalue,
#         log2fc   = lambda w: w.log2fc,
#         contrast = lambda w: w.contrast,
#     log:
#         "results/00log/deseq2_tx/{contrast}.log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.volcano.log"
#     script:
#         "../scripts/volcano.R"  


rule star:
    input:
        get_fq
    output:
        bam   = "results/02alignments/{sample}/{sample}.bam",
        log   = "results/02alignments/{sample}/Log.final.out"
    log:
        align   = "results/00log/alignments/{sample}.log",
        rm_dups = "results/00log/alignments/rm_dup/{sample}.log",
    params:
        out_dir      = "results/02alignments/{sample}/",
        star_params  = config["params"]["star"],
        # path to STAR reference genome index
        index        = config["ref"]["index"],
        samtools_mem = config["params"]["samtools_mem"]
    threads:
        CLUSTER["star"]["cpu"]
    shadow: 
        "minimal"
    shell: 
        """
        STAR --genomeDir {params.index} \
        --runThreadN {threads} \
        --readFilesIn {input} \
        --outFileNamePrefix {params.out_dir} \
        --outSAMtype None \
        --outStd Log \
        --quantMode TranscriptomeSAM \
        --outSAMunmapped Within \
        --quantTranscriptomeBan Singleend \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --outFilterMultimapNmax 20 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        {params.star_params} > {log.align}
        samtools view -h {params.out_dir}Aligned.toTranscriptome.out.bam \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb - > {output.bam} 2>> {log.align}
        """
