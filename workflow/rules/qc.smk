# ------- FASTQC ------- #
rule fastqc:
    """This rule is a little bit tricky to fit PE fastq in multiqc.
    Basically the probelm is that fastqc needs to run R1 and R2 separatelly,
    which means 2 fastqc_zip files with different names. This will be recognized
    by multiqc as different samples, so the report will be a mess.
    My workaround has been to use just the forward fastq to create the report.
    For this I need to change the fastq file name (because it has .1.) to fit
    what multiqc is expecting as name. If multiqc reads A.1.fastq it won't know
    that that file must match A.bam in the report, and they will be in different
    rows. I know it's a pain of workaround but it's the only solution I found.
    For this I need to create a symlink to the fastq to be able to change it's name
    (without duplicating the file which would be less efficient). To make sure that 
    the symlink will work it needs to be created from the folder where it's going to be,
    that's why the cd command of the rule it's imporant. Since the fastq folder can change
    this step needs to work always, it's the only solution I came up with.
    """
    input:  
        get_fq_forward
    output: 
        "results/01qc/fqc/{sample}_fastqc.zip"
    log:    
        "results/00log/fqc/{sample}.log"
    params:
        folder_name = "results/01qc/fqc/",
        tmp = "{sample}.fastq.gz"
    threads: 
        CLUSTER["fastqc"]["cpu"]
    message: 
        "Running fastqc for {input}"
    shadow: 
        "minimal"
    shell:
        """
        cd {params.folder_name} # Move to folder where symlink is going to be created
        ln -s {input} {params.tmp} # Create symlink to fastq file. Imporant to set the desired file name.
        cd - # Go back to workdir
        fastqc -o {params.folder_name} -f fastq -t {threads} --noextract {params.folder_name}/{params.tmp} 2> {log}
        """


##------- RSEQC -------##
rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"]
    output:
        bed  = "results/01qc/rseqc/annotation.bed",
        pred = temp("results/01qc/rseqc/annotation.pred")
    log:
        "results/00log/rseqc_gtf2bed.log"
    shell:
        "resources/gtfToGenePred {input} {output.pred} && resources/genePredToBed {output.pred} {output.bed}"

rule rseqc_stat:
    input:
        expand("results/02alignments/{sample}/{sample}.bam", sample = SAMPLES)
    output:
        "results/01qc/rseqc/{sample}.stats.txt"
    priority: 1
    log:
        "results/00log/rseqc/rseqc_stat/{sample}.log"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

rule rseqc_innerdis:
    input:
        bam =  expand("results/02alignments/{sample}/{sample}.bam", sample = SAMPLES),
        bed = "results/01qc/rseqc/annotation.bed"
    output:
        "results/01qc/rseqc/{sample}.inner_distance_freq.inner_distance_freq.txt"
    priority: 1
    log:
        "results/00log/rseqc/rseqc_innerdis/{sample}.log"
    params:
        prefix="results/01qc/rseqc/{sample}.inner_distance_freq"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam =  expand("results/02alignments/{sample}/{sample}.bam", sample = SAMPLES),
        bed = "results/01qc/rseqc/annotation.bed"
    output:
        "results/01qc/rseqc/{sample}.read_distribution.txt"
    priority: 1
    log:
        "results/00log/rseqc/rseqc_readdis/{sample}.log"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule rseqc_geneCoverage:
    input:
        bam   = expand("results/02alignments/{sample}/{sample}.bam", sample = SAMPLES),
        index = expand("results/02alignments/{sample}/{sample}.bam.bai", sample = SAMPLES),
        bed   = "results/01qc/rseqc/annotation.bed"
    output:
        "results/01qc/rseqc/{sample}.geneBodyCoverage.geneBodyCoverage.txt"
    params:
        prefix="results/01qc/rseqc/{sample}.geneBodyCoverage"
    priority: 1
    shadow:
        "minimal"
    log:
        "results/00log/rseqc/rseqc_geneCoverage/{sample}.log"
    shell:
        "geneBody_coverage.py -r {input.bed} -i {input.bam}  -o {params.prefix} 2> {log}"


## Complexity Curve
# This rule uses the Preseq package () to compute and plot an observed/expected genomic library complexity curve
# Can be used to assess the sequencing depth of an experiment and to evaluate further depth
# Takes in a BAM file with duplicates
rule library_complexity:
    input:
        "results/02alignments/{sample}/{sample}_marked.bam"
    output:
        c_curve    = "results/01qc/preseq/curve/{sample}_c_curve.txt",
        lc_extrap  = "results/01qc/preseq/extrap/{sample}_lc_extrap.txt",
    message:
        "creating creating complexity tables for current and future genomic libraries"
    params:
        curve      = config["params"]["curve"],
        extrap     = config["params"]["extrap"]
    log:
        c_curve    = "results/00log/preseq/{sample}_c_curve.log",
        lc_extrap  = "results/00log/preseq/{sample}_lc_extrap.log",
    shell:
        """
       preseq c_curve {params.curve} \
       -o {output.c_curve} {input} 2> {log.c_curve}
       preseq lc_extrap {params.extrap} \
       -o {output.lc_extrap} {input} 2> {log.lc_extrap}
        """

rule plot_lib_complex:
    input:
       c_curve   = expand("results/01qc/preseq/curve/{sample}_c_curve.txt", sample = SAMPLES),
       lc_extrap = expand("results/01qc/preseq/extrap/{sample}_lc_extrap.txt", sample = SAMPLES)
    output:
        c_curve_plot    = "results/01qc/preseq/curve/c_curve.pdf",
        lc_extrap_plot  = "results/01qc/preseq/extrap/lc_extrap.pdf",
        obs_exp_plot    = "results/01qc/preseq/obs_exp_plot.pdf"
    params:
        sample_names = expand("{sample}", sample = SAMPLES)
    log:
        "results/00log/preseq/plot_complexity_curve.log"
    script:
        "../scripts/plot_library_complexity.R"


## Duplication Plot
# This rule uses and R package to plot and assess duplication problems. Read duplication has a strong
# biological source that's more than the usual technical PCR duplication. Gives insights into the duplication problem
# by relating gene expression level and duplication rate present on it.
rule dupRadar:
    input:
        marked_bam = "results/02alignments/{sample}/{sample}_marked.bam"
    output:
        plot = "results/01qc/dupradar/{sample}.pdf"
    params:
        # GTF to counts the read falling on a feature, use a different GTF cause gencode one seems to give a less clear plot
        gtf = config["ref"]["annot_refseq"],
        strand = 1,
        #'0' (unstranded), '1' (stranded) and '2' (reversely stranded)
        threads = 8
    log:
        "results/00log/dupradar/{sample}.log"
    script:
        "../scripts/dupRadar.R"

    

# ---------------- MultiQC report ----------------- #
rule multiQC_inputs:
    input:
        expand("results/01qc/fqc/{sample}_fastqc.zip", sample = SAMPLES),
        expand("results/02alignments/{sample}/Log.final.out", sample = SAMPLES),
        expand("results/03featureCounts/{sample}/{sample}.featureCounts.summary", sample = SAMPLES),
        expand("results/01qc/rseqc/{sample}.stats.txt", sample = SAMPLES),
        expand("results/01qc/rseqc/{sample}.inner_distance_freq.inner_distance_freq.txt", sample = SAMPLES),
        expand("results/01qc/rseqc/{sample}.read_distribution.txt", sample = SAMPLES),
        expand("results/00log/alignments/rm_dup/{sample}.log", sample = SAMPLES),
        expand("results/01qc/rseqc/{sample}.geneBodyCoverage.geneBodyCoverage.txt", sample = SAMPLES),
        expand("results/01qc/dupradar/{sample}.pdf", sample = SAMPLES),
        "results/01qc/preseq/curve/c_curve.pdf"
    output: 
        file = "results/01qc/multiqc/multiqc_inputs.txt"
    message:
        "create file containing all multiqc input files"
    run:
        with open(output.file, 'w') as outfile:
            for fname in input:
                    outfile.write(fname + "\n")



rule multiQC:
    input:
        "results/01qc/multiqc/multiqc_inputs.txt"
    output: 
        "results/01qc/multiqc/multiqc_report.html"
    params:
        log_name = "multiqc_report",
        folder   = "results/01qc/multiqc"
    log:
        "results/00log/multiqc/multiqc.log"
    message:
        "multiqc for all logs"
    shell:
        """
        multiqc -o {params.folder} -l {input} -f -v -n {params.log_name} 2> {log}
        """
