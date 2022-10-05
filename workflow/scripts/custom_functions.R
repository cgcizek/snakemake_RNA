#----------------------------------------------------------------------------------------------
# Annotation of differential expression table
#----------------------------------------------------------------------------------------------
Annot_DE <- function(df, log2FC = 2, padjust = 0.05, tpm = 0) {
  require(dplyr)
  require(purrr)
    
    df <- mutate(df, max_tpm = reduce(select(df, contains("_TPM")), pmax)) # Add column with max tpm value between conditions
    
    df$DEG <- "NS"
    df$DEG[which(df$log2FoldChange >= log2FC & df$padj <= padjust & df$max_tpm >= tpm)] <- "Upregulated" #Annotate UPregulated genes
    df$DEG[which(df$log2FoldChange <= -log2FC & df$padj <= padjust & df$max_tpm >= tpm)] <- "Downregulated" #Annotate DOWNregulated genes
    
    df <- select(df, -max_tpm) # Remove temporary max tpm column
  
  return(df)
}


#----------------------------------------------------------------------------------------------
# Functions to transform counts into fpkm and tpm
#----------------------------------------------------------------------------------------------
do_fpkm = function (counts, effective_lengths) {
  exp(log(counts) - log(effective_lengths) - log(sum(counts)) + log(1E9))
}

do_tpm = function (counts, effective_lengths) {
  rate = log(counts) - log(effective_lengths)
  exp(rate - log(sum(exp(rate))) + log(1E6))
}


#----------------------------------------------------------------------------------------------
# Function to downsample count matrix
#----------------------------------------------------------------------------------------------
Down_Sample_Matrix <-
function (expr_mat) 
{
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
        prob <- min_lib_size/sum(x)
        return(unlist(lapply(x, function(y) {
            rbinom(1, y, prob)
        })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
}


#----------------------------------------------------------------------------------------------
# Enrichment analyses functions
#----------------------------------------------------------------------------------------------
goEnrichment <- function(
  df, 
  ont      = "BP", 
  db       = org.Mm.eg.db, 
  pvalue   = 0.05, 
  qvalue   = 0.1,
  universe = NULL
  ) {
  require(clusterProfiler)
    ego <- enrichGO(gene        = df,
                  OrgDb         = db,
                  keyType       = 'ENTREZID',
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = pvalue,
                  qvalueCutoff  = qvalue,
                  universe      = universe,
                  readable      = TRUE)
  return(ego)
}

KEGGenrichment <- function(
  df, 
  org      = "mmu", 
  db       = org.Mm.eg.db, 
  pvalue   = 0.05, 
  qvalue   = 0.1,
  universe = NULL
  ) {
  require(clusterProfiler)
  ekgg <- enrichKEGG(gene          = df,
                     organism      = org,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = pvalue,
                     qvalueCutoff  = qvalue,
                     universe      = universe)
  ekgg <- DOSE::setReadable(ekgg, OrgDb = db, keyType="ENTREZID")
  return(ekgg)
}

PAenrichment <- function(
  df, 
  org      = "mouse", 
  pvalue   = 0.05, 
  qvalue   = 0.1,
  universe = NULL
  ) {
  require(clusterProfiler)
  require(ReactomePA)
  ePA <- enrichPathway(gene         = df,
                       pvalueCutoff = pvalue,
                       qvalueCutoff = qvalue,
                       organism     = org,
                       universe     = universe,
                       readable     = TRUE)
  return(ePA)
}

msig_db_enrichment <- function(
  df, 
  db        = org.Mm.eg.db, 
  pvalue    = 0.05, 
  qvalue    = 0.1,
  universe  = NULL,
  term2gene = term2gene
) {
  require(clusterProfiler)
  GSEA_hyp  <- enricher(gene        = df,
                       pvalueCutoff = pvalue,
                       qvalueCutoff = qvalue,
                       TERM2GENE    = term2gene)
  if(is.null(GSEA_hyp)){
    return(data.frame())
  }
  else{
    GSEA_hyp <- DOSE::setReadable(GSEA_hyp, OrgDb = db, keyType="ENTREZID")
    return(GSEA_hyp)
  }
}

GSEA_enrichment <- function(df, pathways.gmt) {
  require(clusterProfiler)
  pathways <- fgsea::gmtPathways(pathways.gmt)
  df <- df %>% filter(!is.na(log2FoldChange))
  # Create a list containing a named vector (with genenames) of log2fc
  geneList <- df$log2FoldChange
  names(geneList) <- toupper(df$SYMBOL)
  
  # Run GSEA algorithm
  fgseaRes <- fgsea::fgsea(pathways = pathways, 
                           stats   = geneList,
                           minSize = 15,
                           maxSize = 2000)
                           #nperm   = 10000)
  return(fgseaRes)
}


#----------------------------------------------------------------------------------------------
# Volcano plot code
#----------------------------------------------------------------------------------------------
VolcanoPlot <- function(df, xlim=NULL, ylim=NULL, main = NULL, labelSize = 8, pval = 0.05, log2FC = 1) {
  require(ggplot2)
  require(dplyr)
  # require(ggrastr)

  df <- mutate(df, shape = "circle")

  df <- mutate(df, shape = ifelse(-log10(padj) >  ylim[2], "triangle", shape))
  df <- mutate(df, padj = ifelse(-log10(padj) >  ylim[2], 10^-ylim[2], padj))

  df <- mutate(df, shape = ifelse(log2FoldChange > xlim[2], "triangle", shape))
  df <- mutate(df, shape = ifelse(log2FoldChange < -xlim[2], "triangle", shape))

  df <- mutate(df, log2FoldChange = ifelse(log2FoldChange > xlim[2], xlim[2], log2FoldChange))
  df <- mutate(df, log2FoldChange = ifelse(log2FoldChange < -xlim[2], -xlim[2], log2FoldChange))


  p <-  ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape) ) +

    # geom_point_rast(alpha=0.7, size=1.7, raster.height = 5.15, raster.width = 6, raster.dpi = 400) +
    geom_point(alpha=0.7, size=1.7) +

    annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = xlim[2],
             vjust="inward",hjust="inward", size = labelSize) +
    annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = xlim[1],
             vjust="inward",hjust="inward", size = labelSize) +

    theme_classic(base_size = 20) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "top") +

    ggtitle(main) +
    theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

    xlim(xlim) + ylim(ylim) +

    geom_hline(yintercept = -log10(pval), linetype = 2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

    xlab(bquote(~Log[2]~ "fold change")) + ylab(bquote(~-Log[10]~italic(P))) +

    scale_colour_manual(values=c("Downregulated" = "darkgreen", "NS" = "gray", "Upregulated" = "red"),
                        labels = c("Downregulated" = "Downregulated", "NS" = "NS", "Upregulated" = "Upregulated"),
                        drop = FALSE) + #Force legend to show always

    guides(shape=FALSE) # Remove legend for shapes

  return(p)
}


# Themes from ggtech, code from Ricardo Bion, FONT NECESSARY!!
# https://github.com/ricardo-bion/ggtech

theme_tech <- function(theme="airbnb", tech_key = list(
                  airbnb = list(
                    family_title="Circular Air Bold"
                    , family_text = "Circular Air Medium"
                    , colour_title = "#F14000"
                    , colour_text = "#535353"),
                  facebook = list(
                    family_title="Facebook Letter Faces"
                    , family_text = "Facebook Letter Faces"
                    , colour_title = "#3D579D"
                    , colour_text = "#535353"),
                  google = list(
                    family_title="Product Sans"
                    , family_text = "Product Sans"
                    , colour_title = "#dd4b39"
                    , colour_text = "black"),
                  etsy = list(
                    family_title="Georgia"
                    , family_text = "Georgia"
                    , colour_title = "#F14000"
                    , colour_text = "#535353"),
                  twitter = list(
                    family_title="PicoBlackAl"
                    , family_text = "[z] Arista Light"
                    , colour_title = "#5380E4"
                    , colour_text = "black"),
                  X23andme = list(
                    family_title = "Helvetica Neue Bold"
                    , family_text = "Helvetica"
                    , colour_title = "#7BC143"
                    , colour_text = "#DA1963"
                  )
                    )) {

  theme_classic() + 
    theme(text=element_text(size=18, family=tech_key[[theme]]$family_text)) +
    theme(legend.title=element_blank()) + 
    theme(plot.title = element_text(size = 25, colour = tech_key[[theme]]$colour_title, family=tech_key[[theme]]$family_title)) + 
    theme(plot.subtitle = element_text(size = 15, colour = tech_key[[theme]]$colour_title, family=tech_key[[theme]]$family_title)) + 
    theme(axis.text.x=element_text(color=tech_key[[theme]]$colour_text)) +
    theme(axis.text.y=element_text(color=tech_key[[theme]]$colour_text)) +
    theme(axis.title.x=element_text(color=tech_key[[theme]]$colour_text, vjust=0)) +
    theme(axis.title.y=element_text(color=tech_key[[theme]]$colour_text, vjust=1.25)) + 
    theme(plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
    theme(panel.border= element_blank())+
    theme(axis.line.x = element_line(color=tech_key[[theme]]$colour_text, size = 0.5),
        axis.line.y = element_line(color=tech_key[[theme]]$colour_text, size = 0.5)) +
  theme(line = element_line(color=tech_key[[theme]]$colour_text)) +
  theme(rect = element_rect(color=tech_key[[theme]]$colour_text)) +
  theme(axis.ticks.x = element_line(color=tech_key[[theme]]$colour_text)) +
  theme(axis.ticks.y = element_line(color=tech_key[[theme]]$colour_text))

}

# VolcanoPlot <- function(df, xlim=NULL, ylim=NULL, main = NULL, labelSize = 8, pval = 0.05, log2FC = 1) {
#   require(ggplot2)
#   require(ggrastr)

#   df <- mutate(df, shape = "circle")
#   df <- mutate(df, shape = ifelse(-log10(padj) >  ylim[2], "triangle_up", shape))
#   df <- mutate(df, padj = ifelse(-log10(padj) >  ylim[2], 10^-ylim[2], padj))

#   df <- mutate(df, shape = ifelse(log2FoldChange > xlim[2], "triangle_right", shape))
#   df <- mutate(df, shape = ifelse(log2FoldChange < -xlim[2], "triangle_left", shape))

#   df <- mutate(df, log2FoldChange = ifelse(log2FoldChange > xlim[2], xlim[2], log2FoldChange))
#   df <- mutate(df, log2FoldChange = ifelse(log2FoldChange < -xlim[2], -xlim[2], log2FoldChange))


#   p <-  ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape) ) +

#     # geom_point_rast(alpha=0.7, size=1.7, raster.height = 5.15, raster.width = 6, raster.dpi = 400) +
#     geom_point(alpha=0.7, size=1.7) +
    
#     annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = xlim[2],
#              vjust="inward",hjust="inward", size = labelSize) +
#     annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = xlim[1],
#              vjust="inward",hjust="inward", size = labelSize) +

#     theme_classic(base_size = 20) +
#     theme(legend.title = element_blank()) +
#     theme(legend.position = "top") +


#     ggtitle(main) +
#     theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

#     xlim(xlim) + ylim(ylim) +

#     geom_hline(yintercept = -log10(pval), linetype = 2) +
#     geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

#     xlab(bquote(~Log[2]~ "fold change")) + ylab(bquote(~-Log[10]~italic(P))) +

#     scale_colour_manual(values=c("Downregulated" = "darkgreen", "NotDE" = "gray", "Upregulated" = "red"),
#                         labels = c("Downregulated" = "Downregulated", "NS" = "NS", "Upregulated" = "Upregulated"),
#                         drop = FALSE) + #Force legend to show always

#     scale_shape_manual(values=c("triangle_up" = "\u25B2", "triangle_right" = "\u25BA", "triangle_left" = "\u25C4", "circle" = "\u25CF",
#                                 drop = FALSE)) +

#     guides(shape=FALSE) # Remove legend for shapes

#   return(p)
# }