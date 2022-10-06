log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

#------------------------------------------------------------------------------------------
# Load in Libraries + custom themes
#------------------------------------------------------------------------------------------

# Lib
library("ggplot2")
library("dplyr")
library("tibble")
library("purrr")
library("tidyr")
library("dupRadar")

#------------------------------------------------------------------------------------------
# Read in marked bam file and GTF annotation
#------------------------------------------------------------------------------------------

# Reading in params
marked_bam <- snakemake@input$marked_bam
gtf <- snakemake@params[["gtf"]]
stranded <- snakemake@params[["strand"]]
paired <- TRUE
threads <- snakemake@params[["threads"]]


# Run Duplicate rate Analysis
duplication_matrix <- analyzeDuprates(marked_bam, gtf, stranded, paired, threads)


#------------------------------------------------------------------------------------------
# Plot result
#------------------------------------------------------------------------------------------
pdf(file = snakemake@output[["plot"]], width = 10, height = 10)

duprateExpDensPlot(DupMat=duplication_matrix)
title("Duplication Rate Plot")

dev.off()

