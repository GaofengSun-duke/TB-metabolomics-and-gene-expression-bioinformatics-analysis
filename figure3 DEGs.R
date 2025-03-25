###GSE19491 DEGs
if(T){
exprset <- read.table('GSE19491_Matrix.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pdata <- read.table('GSE19491_ATBvsCon_Group.txt',sep = '\t',header = T,check.names = F)
group <- pdata$group
library(limma)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
row.names(design) = colnames(exprset)
design
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels = design)
contrast.matrix
#step1
fit <- lmFit(exprset,design)
#step2
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput =topTable(fit2,coef = 1,n = Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)
write.table(nrDEG,'GSE19491_DEGs.xls',sep = '\t')

###GSE107994 DEGs
if(T){
  exprset <- read.table('GSE107994_ATBvsCon_nomolization.txt',sep = '\t',header = T,row.names = 1,check.names = F)
  pdata <- read.table('GSE107994_Group.txt',sep = '\t',header = T,check.names = F)
  group <- pdata$Group
  library(limma)
  design <- model.matrix(~0+factor(group))
  colnames(design) = levels(factor(group))
  row.names(design) = colnames(exprset)
  design
  contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels = design)
  contrast.matrix
  #step1
  fit <- lmFit(exprset,design)
  #step2
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput =topTable(fit2,coef = 1,n = Inf)
  nrDEG = na.omit(tempOutput)
  head(nrDEG)
  write.table(nrDEG,'GSE107994_DEGs.xls',sep = '\t')
  
