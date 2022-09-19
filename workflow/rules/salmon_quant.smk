########################################################
# This file contains rules to:
## 1) Index the transcriptome, it gets stored in DATA/DP, move it where you want
##		and update the config file @ ref-salmon_index
## 2) Quantify transcript abundance. Output both quant.sf and quant.genes.sf
## 		Quant.sf is used to then aggregate to gene level abundacnies w/ tximport
##		and downstream with DEseq2
########################################################

# RUN THIS ONLY IF THE INDEX IS DIFFERENT, salmon_quant will crash because it'll not found
# the index. I use this part of the pipeline to avoid running this step in a different script
# Find a way to make it as a conditional execution based on the existence of the index

# The commands of the index rule were obtained from here: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

rule salmon_index:
    input:
        primary_assembly = config["ref"]["assembly"],
        transcripts      = config["ref"]["transcriptome"],
    output:
        index = directory("results/02salmon/salmon_index"),
    log:
        "results/00log/salmonIndex/log"
    params:
        destination = config["ref"]["salmon_index"],
    threads:
        CLUSTER["salmon_index"]["cpu"]
    shell:
        """
        grep "^>" <(gunzip -c {input.primary_assembly}) \
        | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt
        cat {input.transcripts} {input.primary_assembly} > gentrome.fa.gz
        salmon index -t gentrome.fa.gz -d decoys.txt -p {threads} -i {output} --gencode
        rm decoys.txt* gentrome.fa.gz
        """


def set_reads(wildcards, input):
        n = len(input.fastq)
        if n == 1:
            reads = "-r {}".format(*input.fastq)
            return reads
        else:
            reads = "-1 {} -2 {}".format(*input.fastq)
            return reads

### Quasi-Mapping Mode
# Quantification using "salmon quant -g" to aggregate to gene level (to have the file ready)
# Output also quant.sf and use it to aggregate at gene level w/ tximport in R as suggested by the creator Rob Paltro
# https://crazyhottommy.blogspot.com/2016/07/comparing-salmon-kalliso-and-star-htseq.html

# FIX SAM OUTPUT
rule salmon_quant:
    input:
        index = check_index,
        #index = config["salmon_index"],
        fastq = get_fq
    output:
        quant       = "results/02salmon/{sample}/quant.sf",
        quant_genes = "results/02salmon/{sample}/quant.genes.sf",
        sam         = "results/02salmon/{sample}/{sample}.sam",
    log: 
        "results/00log/salmonQuant/{sample}.log"
    params:
        gtf           = config["ref"]["annotation"],
        options       = config["params"]["salmon"],
        library       = config["params"]["salmon_library"],
        out_fold      = "results/02salmon/{sample}",
        reads  	      = set_reads,
    threads:    
        CLUSTER["salmon_quant"]["cpu"]
    shell:
        """
        salmon quant \
        -i {input.index} \
        -p {threads} \
        -g {params.gtf} \
        -l {params.library}\
        {params.reads} \
        -o {params.out_fold} \
        {params.options} \
        2> {log} \
        > {output.sam}
        """

# ### SALMON BigWig
# rule salmon_bw:
#     input:
#         rules.salmon_quant.output.sam
#     output:
#         "results/05bigwig_salmon/{sample}.bw"
#     params:
#         tmp = "results/05bigwig_salmon/{sample}_tmp"
#     log:
#         "results/00log/bam2bigwig_salmon/{sample}.log"
#     threads:
#         CLUSTER["bam2bigwig"]["cpu"]
#     run:
#     """
#     take correct lines of the SAM file {input} | samtools view -@ {threads} -Sb - | \
#     samtools sort -@ {threads} -o {params.tpm}
#     samtools index {params.tmp}
#     bamCoverage --normalizeUsing CPM -p {threads} -bs 1 -b {params.tmp} -o {output} 2> {log}
#     rm {params.tmp} {params.tmp}.bai
#     """





### Salmon Alignment mode (i'm not planning on using it)
# Fall back to Alignment based quantification with STAR but use Salmon in align mode instead of featureCounts
# # Most star parameters taken from https://www.biorxiv.org/content/biorxiv/early/2019/10/31/657874.full.pdf

# rule star:
#     input:
#         get_fq
#     output:
#         bam   = "results/02alignments/{sample}/{sample}.bam",
#         log   = "results/02alignments/{sample}/Log.final.out"
#     log:
#         align   = "results/00log/alignments/{sample}.log",
#         rm_dups = "results/00log/alignments/rm_dup/{sample}.log",
#     params:
#         out_dir      = "results/02alignments/{sample}/",
#         star_params  = config["params"]["star"],
#         # path to STAR reference genome index
#         index        = config["ref"]["index"],
#         samtools_mem = config["params"]["samtools_mem"]
#     threads:
#         CLUSTER["star"]["cpu"]
#     shadow: 
#         "minimal"
#     shell: 
#         """
#         STAR --genomeDir {params.index} \
#         --runThreadN {threads} \
#         --readFilesIn {input} \
#         --outFileNamePrefix {params.out_dir} \
#         --outSAMtype None \
#         --outStd Log \
#         --quantMode TranscriptomeSAM \
#         --outSAMunmapped Within \
#         --quantTranscriptomeBan Singleend \
#         --outFilterType BySJout \
#         --alignSJoverhangMin 8 \
#         --outFilterMultimapNmax 20 \
#         --alignSJDBoverhangMin 1 \
#         --outFilterMismatchNmax 999 \
#         --outFilterMismatchNoverReadLmax 0.04 \
#         --alignIntronMin 20 \
#         --alignIntronMax 1000000 \
#         --alignMatesGapMax 1000000 \
#         {params.star_params} > {log.align}
#         samtools view -h {params.out_dir} Aligned.toTranscriptome.out.bam \
#         | samblaster --removeDups 2> {log.rm_dups} \
#         | samtools view -Sb - > {output.bam} 2>> {log.align}
#         """

