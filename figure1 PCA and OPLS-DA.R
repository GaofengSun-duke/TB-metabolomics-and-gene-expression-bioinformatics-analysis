  library(tidyselect)
  library(dplyr)
  library(ropls)
  library(Rserve)
  library(ggplot2)
  metadata <- read.table('TB.vs.NC.txt',sep = '\t',header = T,row.names = 1,check.names = F)
  exprset <- t(metadata)
  pData <- data.frame(c(rep('TB',74),rep('NC',106)))
  pData$sample <- row.names(exprset)
  colnames(pData) <- c('group','sample')
  group <- pData
  
  plsda <- opls(exprset,group$group,predI=1,orthoI=1,log10L=F,crossvalI=nrow(exprset),scaleC="pareto",permI=200)
  
  vip <- plsda@vipVn%>% as.data.frame() %>% rename("VIP"=".") %>% tibble::rownames_to_column(var = 'variables')
  
  plsdaScore <- data.frame(
    t1 =plsda@scoreMN,
    to1 =plsda@orthoScoreMN 
  ) %>% scale(center = T,scale = T) %>% 
    as.data.frame() %>% 
    rename(
      "t1"="p1",
      "to1"="o1"
    ) %>%
    tibble::rownames_to_column(var = 'sample') %>%
    merge(group,by="sample")
  
  
  t1Weight=sprintf("%.1f%%", plsda@modelDF[1,1]*100);t1Weight
  to1Weight=sprintf("%.1f%%", plsda@modelDF[2,1]*100);to1Weight
  
  R2X=plsda@modelDF[1,1]+plsda@modelDF[2,1]
  R2Y=plsda@modelDF[1,3]+plsda@modelDF[2,3]
  Q2Y=plsda@modelDF[1,6]+plsda@modelDF[2,6]
  
  subTitle <- paste0("R2X=",R2X,"R2Y=",R2Y,"Q2Y=",Q2Y)
  ##opls-da
  oplsdaFig <- ggplot(plsdaScore,aes(x=t1,y=to1,color=group))+
    geom_point(aes(shape=group),size=3)+
    stat_ellipse(aes(fill=group),alpha=0.2,geom = "polygon")+
    theme_classic()+
    scale_fill_manual(values=c("#2B8BC34E","#4E0C664E"))+
    scale_color_manual(values=c("#2B8BC34E","#4E0C664E"))+
    labs(title = "OPLS-DA", 
         subtitle = subTitle,
         x = paste0("t1(",t1Weight,")"),
         y = paste0("to1(",to1Weight,")"))+
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1.2)), 
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)) 
    )
  #saving
  ggsave("outputfile/oplsda.png",width = 4,height = 4)
  print(oplsdaFig)
  ggsave("outputfile/oplsda.pdf",width = 4,height = 4,onefile=F)
  print(oplsdaFig)
  dev.off()
  
  permutation <- plsda@suppLs[["permMN"]] %>% as.data.frame() %>% 
    select(`R2Y(cum)`,`Q2(cum)`) %>% 
    mutate("number"=1:nrow(.)) %>% 
    tidyr::pivot_longer(!number,names_to = "Type",values_to = "value") %>% 
    mutate(Type = stringr::str_replace(Type, "\\(cum\\)", ""))
  permutation %>% write.csv("outputfile/200permutation.csv")
  
  R2Xlabel <- paste0("R2X:",plsda@summaryDF[1,1])
  R2Ylabel <- paste0("R2Y:",plsda@summaryDF[1,2],"\np=",
                     plsda@summaryDF[1,7],"(",plsda@summaryDF[1,7]*200,"/200)")
  Q2lablel <- paste0("Q2:",plsda@summaryDF[1,3],"\np=",plsda@summaryDF[1,8],
                     "(",plsda@summaryDF[1,8]*200,"/200)")
  
  calculate_bin_counts <- function(data,#Long-format data of permutation
                                   column, #Data column in the long format of permutation
                                   bins=30#Consistent with the 'bins' parameter of geom_histogram
  ) {
    if (!column %in% names(data)) {
      stop("zcp:Specified column does not exist in the dataframe!!!")
    }
    data_vector <- data[[column]]
    range_data <- range(data_vector, na.rm = TRUE)
    bin_width <- (range_data[2] - range_data[1]) / bins
    bin_breaks <- seq(from = range_data[1], to = range_data[2], by = bin_width)
    binned_data <- cut(data_vector, breaks = bin_breaks, include.lowest = TRUE, right = FALSE)
    
    bin_counts <- table(binned_data)
    
    return(bin_counts)
  }
  
  bin_counts <- calculate_bin_counts(data=permutation,column="value", bins=30)
  max_counts <- max(bin_counts)
  
  labels_data <- data.frame(
    x = c(plsda@summaryDF[1,1], plsda@summaryDF[1,2], plsda@summaryDF[1,3]),
    y = c(max_counts*0.4, max_counts*0.8,max_counts*1.2), 
    label = c(R2Xlabel, R2Ylabel,Q2lablel)
  )
  #Plotting Permutation Histogram
  PermutationFig <- ggplot()+
    geom_histogram(data = permutation,aes(x=value,fill=Type),
                   bins=30,# Keep consistent with the bins parameter in calculate_bin_counts function!
                   alpha=0.6,
                   size=0.3,
                   color="black")+
    scale_fill_manual(values=c("#2B8BC34E","#4E0C664E"))+
    xlim(c(NA, max(labels_data$x) * 1.2))+
    ylim(c(NA, max(labels_data$y) * 1.2))+
    geom_label(data = labels_data, 
               aes(x = x, y = y, label = label), 
               size = 2.5, color = "darkred",
               fill = "white", label.padding = unit(0.1, "lines"),
               label.size = 0.25)+
    geom_segment(data = labels_data, aes(x = x, xend = x, y = y, yend = 0,color=label),
                 arrow = arrow(type = "closed", length = unit(0.15, "cm")),
                 alpha=0.5,
                 show.legend = FALSE)+
    labs(title = "Permutation of OPLS-DA", 
         x = "Permutations",
         y ="Frequency")+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1.2)))
  
  #saving
  ggsave("outputfile/Permutation.png", width = 6, height = 4)
  print(PermutationFig)
  ggsave("outputfile/Permutation.pdf", width = 6, height = 4,onefile=F)
  print(PermutationFig)
  dev.off()    
  
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
  write.table(result,'outputfile/TB.vs.NC.Result.xls',sep = '\t')
  
  library(factoextra)
 library(FactoMineR)
 group <- as.factor(pData$group)
 df <- as.data.frame(t(metadata))
 dat.pca <- PCA(df,graph = F)
 fviz_pca_ind(dat.pca,
             geom.ind = 'point',
             col.ind = group,
             addEllipses = T,
             legend.title = 'Groups')
