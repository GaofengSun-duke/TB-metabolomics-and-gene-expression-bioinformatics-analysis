library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggfun)
library(tidyverse)
volcano_plot_enhanced <- function(df, num_symbol = 10, logFC="logFC", FDR="adj.P.Val", logFC_Value=1.5,
                                  Symbol="Symbol", y_increased=20, labs="Genes") {
  ## Define the up and down genes/symbols
  df <- df %>% mutate(DEG = case_when(
    !!sym(logFC) > logFC_Value & !!sym(FDR) < 0.05 ~ "Up",
    abs(!!sym(logFC)) < logFC_Value | !!sym(FDR) > 0.05 ~ "None",
    !!sym(logFC) < -logFC_Value & !!sym(FDR) < 0.05 ~ "Down"
  )) %>%
    mutate(!!sym(Symbol) := rownames(.)) 
  
  ## Check that belongs to the column name is correct
  if (!logFC %in% colnames(df)) {
    cat(paste0("The '", logFC, "' column is not found in the dataset.\n"))
    logFC <- readline("Please enter the name of the 'logFC' column: ")
  }
  
  if (!FDR %in% colnames(df)) {
    cat(paste0("The '", FDR, "' column is not found in the dataset.\n"))
    FDR <- readline("Please enter the name of the 'FDR' column: ")
  }
  
  if (!Symbol %in% colnames(df)) {
    cat(paste0("The '", Symbol, "' column is not found in the dataset.\n"))
    Symbol <- readline("Please enter the name of the 'Symbol' column: ")
  }
  
  print(paste0("Column name of logFC is: ", logFC))
  print(paste0("Column name of FDR is: ", FDR))
  print(paste0("Column name of Symbol is: ", Symbol))
  
  ## Calculate the number of upregulated and downregulated genes
  upregulated_genes <- sum(df[[logFC]] > logFC_Value & df[[FDR]] < 0.05)
  downregulated_genes <- sum(df[[logFC]] < -logFC_Value & df[[FDR]] < 0.05)
  
  y_upper <- ceiling(max(-log10(df[[sym(FDR)]])))
  x_abs <- max(abs(c(floor(min(df[[sym(logFC)]])), ceiling(max(df[[sym(logFC)]])))))
  
  plot <- ggplot(data = df) + 
    geom_point(aes(x = !!sym(logFC), y = -log10(!!sym(FDR)), color = !!sym(logFC), size = -log10(!!sym(FDR)))) + 
    #up
    geom_text_repel(data = df %>%
                      tidyr::drop_na() %>%
                      dplyr::filter(!!sym(logFC) > 0 & !!sym(FDR) < 0.05) %>%
                      dplyr::arrange(desc(-log10(!!sym(FDR)))) %>%
                      dplyr::slice(1:num_symbol),
                    aes(x = !!sym(logFC), y = -log10(!!sym(FDR)), label = !!sym(Symbol)), 
                    nudge_x = 0.5, nudge_y = 0.2, segment.curvature = -0.1, segment.ncp = 3,
                    direction = "y", hjust = "left",
                    max.overlaps = 200) +
    #down
    geom_text_repel(data = df %>%
                      tidyr::drop_na() %>%
                      dplyr::filter(!!sym(logFC) < 0 & !!sym(FDR) < 0.05) %>%
                      dplyr::arrange(desc(-log10(!!sym(FDR)))) %>%
                      dplyr::slice(1:num_symbol),
                    aes(x = !!sym(logFC), y = -log10(!!sym(FDR)), label = !!sym(Symbol)), 
                    box.padding = 0.5, nudge_x = -0.2, nudge_y = 0.2, segment.curvature = -0.1, segment.ncp = 3,
                    segment.angle = 20,
                    max.overlaps = 200) +
    scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"), values = seq(0, 1, 0.2)) +
    scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"), values = seq(0, 1, 0.2)) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 4) + 
    scale_size(range = c(1,6)) + 
    ggtitle(label = "Volcano Plot") + 
    xlim(-x_abs, x_abs) + 
    ylim(c(-1, y_upper + y_increased)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.background = element_roundrect(color = "#808080", linetype = 1),
          axis.text = element_text(size = 13, color = "#000000"),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5)) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#3288bd")), 
      xmin = -x_abs, xmax = -1, ymin = y_upper + y_increased, ymax = y_upper + y_increased) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Down",
        gp = grid::gpar(col = "#3288bd")),
      xmin = -x_abs, xmax = -1, ymin = y_upper + y_increased, ymax = y_upper + y_increased) +
    annotation_custom(
      grob = grid::segmentsGrob(
        y0 = unit(-10, "pt"),
        y1 = unit(-10, "pt"),
        arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
        gp = grid::gpar(lwd = 3, col = "#d73027")), 
      xmin = x_abs, xmax = 1, ymin = y_upper + y_increased, ymax = y_upper + y_increased) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Up",
        gp = grid::gpar(col = "#d73027")),
      xmin = x_abs, xmax = 1, ymin = y_upper + y_increased, ymax = y_upper + y_increased) +
    annotation_custom(
      grob = grid::textGrob(
        label = paste(downregulated_genes, labs),
        gp = grid::gpar(col = "#3288bd")),
      xmin = -x_abs, xmax = -1, ymin = y_upper + y_increased - 2, ymax = y_upper + y_increased - 4) +
    annotation_custom(
      grob = grid::textGrob(
        label = paste(upregulated_genes, labs),
        gp = grid::gpar(col = "#d73027")),
      xmin = x_abs, xmax = 1, ymin = y_upper + y_increased - 2, ymax = y_upper + y_increased - 4) 
  
  return(plot)
}

dat <- metadata
group_list <- as.factor(group$group)
YM_data <- metadata[,1:74]
NC_data <- metadata[,75:180]
pvals <- apply(metadata,1,function(x){
  t.test(as.numeric(x)~group_list)$p.value
})
p.adj = p.adjust(pvals,method = 'BH')
YM_mean <- rowMeans(YM_data)
NC_mean <- rowMeans(NC_data)
FC <- YM_mean/NC_mean
log2FC <- log2(FC)
result <- cbind(YM_mean,NC_mean,FC,log2FC,pvals,p.adj,vip)
write.table(result,'TB.vs.NC.neg.Result.xls',sep = '\t')

df <- read.table('YM.vs.NC.Volcano.xls',sep = '\t',header = T,row.names = 1,check.names = F)
colnames(df) <- c("Type","ID","YMmean","NCmean","FC","logFC","P.Value","P.adj","VIP")
volcano_plot_enhanced(df,y_increased = 20,labs = "Metabolite",logFC_Value = 1,num_symbol = 0,FDR="P.Value")
