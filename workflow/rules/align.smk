########################################################
# This file contains rules to:
## 1) Align reads with STAR (remove duplicates)
## 2) Produce counts with featureCounts
## 3) Produce bigWig files (STAR + FC)
## 4) Create BAm files marked for QC
########################################################

## Alignment
# Mapping fastq file to reference genome with STAR. In this branch of the pipeline
# this rule is here only to output BW files.
rule star:
    input:
        get_fq
    output:
        bam   = "results/02alignments/{sample}/{sample}.bam",
        index = temp("results/02alignments/{sample}/{sample}.bam.bai"),
        log   = "results/02alignments/{sample}/Log.final.out"
    log:
        align   = "results/00log/alignments/{sample}.log",
        rm_dups = "results/00log/alignments/rm_dup/{sample}.log",
    params:
        out_dir      = directory("results/02alignments/{sample}/"),
        star_params  = config["params"]["star_noSalmon"],
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
        --outSAMtype SAM \
        --outStd SAM \
        {params.star_params} 2> {log.align} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam} 2>> {log.align}
        """

## CountTables
# producing count tables with featureCounts. In this branch of the pipeline, this rule is here only to
# output an annotated sam/bam to be used w/ deeptols
rule featureCounts:
    input:
        rules.star.output.bam
    output: 
        annot_sam     = temp("results/03featureCounts/{sample}/{sample}.bam.featureCounts.sam"),
        featureCounts = temp("results/03featureCounts/{sample}/{sample}.featureCounts"),
        summary       = temp("results/03featureCounts/{sample}/{sample}.featureCounts.summary")
    log:
        "results/00log/featureCounts/{sample}.log"
    params:
        tmp      = "results/03featureCounts/{sample}/{sample}.seqDepth",
        gtf      = config["ref"]["annotation"],
        options  = config["params"]["featureCounts"],
        # Add -p option for pair-end data if it's the case
        pair_end = lambda w: "-p" if not is_single_end(w.sample) else str()
    threads: 
        CLUSTER["featureCounts"]["cpu"]
    shell:
        """
        featureCounts -T {threads} \
        -a {params.gtf} \
        -o {output.featureCounts} \
        {params.pair_end} {params.options} -R SAM {input} > {log} 2>&1
        """

## BigWig file
# This rule output a normlaized bigwig file
rule bam2bigwig:
    input:
        rules.featureCounts.output.annot_sam
    output:
        "results/05bigwig/{sample}.bw"
    params:
        tmp = "results/05bigwig/{sample}_tmp"
    log:
        "results/00log/bam2bigwig/{sample}.log"
    threads:
        CLUSTER["bam2bigwig"]["cpu"]
    shell:
        """
        fgrep -v "Unassigned_" {input} | samtools view -@ {threads} -Sb - | samtools sort -@ {threads} -o {params.tmp}
    samtools index {params.tmp}
    bamCoverage --normalizeUsing CPM -p {threads} -bs 1 -b {params.tmp} -o {output} 2> {log}
        rm {params.tmp} {params.tmp}.bai
        """


# rule to remove all the previous outputs to clean up the work dir.

## Alignment for QC
# This rule output a bam file that has duplicates marked insted of reemoved (like the first rule). This is necessary
# because QC rule library_complexity (qc.smk) requires a bam file with duplicates marked.
rule star_dupRadar:
    input:
        get_fq
    output:
        bam           = "results/02alignments/{sample}/{sample}_marked.bam"
    log:
        align_dup     = "results/00log/alignments/{sample}_marked.log",
        dup_mark      = "results/00log/alignments/dup_mark/{sample}_marked.log"
    params:
        out_dir       = directory("results/02alignments/marked/{sample}/"),
        star_params   = config["params"]["star_noSalmon"],
        # path to STAR reference genome index
        index         = config["ref"]["index"],
        samtools_mem  = config["params"]["samtools_mem"],
    threads:
        CLUSTER["star_dupRadar"]["cpu"]
    shadow: 
        "minimal"
    shell: 
        """
        STAR --genomeDir {params.index} \
        --runThreadN {threads} \
        --readFilesIn {input} \
        --outFileNamePrefix {params.out_dir} \
        --outSAMtype SAM \
        --outStd SAM \
        {params.star_params} 2> {log.align_dup} \
        | samblaster 2> {log.dup_mark} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align_dup}
        samtools index {output.bam} 2>> {log.align_dup}
        """