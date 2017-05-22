## Monocle for HSC
rm(list = ls())
source("~/Documents/Project/Network.comparison/functions.R")
setwd("~/Documents/Project/Jan_2017/")
library(monocle)
## data ################################################################################
meso.ct <- read.csv("./Data/mesoderm.1.csv")
rownames(meso.ct) <- meso.ct[,1]
meso.ct <- meso.ct[,c(-1)]
dim(meso.ct) # 3934 * 46
gene.list <- colnames(read.csv("./Data/gene33.meso.csv"))[c(-1)]
meso.model <- meso.ct[colnames(meso.ct) %in% gene.list]
colnames(meso.model) <- Transform.standrdize.names(colnames(meso.model)); rm(meso.ct)
meso.model <- meso.model[,order(colnames(meso.model))]
meso.model[meso.model == -14] <- min(meso.model[meso.model != -14]) - 0.5; meso.model <- meso.model - min(meso.model);range(meso.model) # -11.020677   5.346537
data <- t(meso.model)
cell.name <- c(NA, ncol(data))
for(i in 1:ncol(data))
{
  cell.name.full <- colnames(data)[i];
  cell.name[i] <- strsplit(cell.name.full,"_")[[1]][1]
}

length(grep("PS",cell.name)) + length(grep("SG",cell.name)) + 
  length(grep("SFG",cell.name)) + length(grep("HF",cell.name)) + length(grep("NP",cell.name))
# should be 3934
##phenoData, an AnnotatedDataFrame object, where rows are cells, 
##and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
temp <- matrix(data = NA, ncol = ncol(data), nrow = 3)
rownames(temp)<-c("cellID", "celltype", "stage"); colnames(temp) <- colnames(data); 
temp <- as.data.frame(temp); temp[1, ] <- cell.name
temp[2:3,grep("PS",cell.name)] <- c("PS",1);temp[2:3,grep("NP",cell.name)] <- c("NP",2);
temp[2:3,grep("HF",cell.name)] <- c("HF",3);temp[2:3,grep("SG",cell.name)] <- c("SG",4);temp[2:3,grep("SFG",cell.name)] <- c("SFG",5);
phenoData <- new("AnnotatedDataFrame", data = as.data.frame(t(temp)))

##featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes),
##and columns are gene attributes, such as biotype, gc content, etc.

temp <- matrix(data = NA, ncol =1, nrow = nrow(data))
colnames(temp) <- "GeneID"; rownames(temp) <- temp[,1] <- rownames(data); 
temp <- as.data.frame(temp); 
featureData <- new("AnnotatedDataFrame", data = temp)

HSC_monole <- newCellDataSet(data, phenoData = phenoData, featureData = featureData, expressionFamily = gaussianff())
print(head(pData(HSC_monole)))

HSC_monole <- estimateSizeFactors(HSC_monole)
##HSC_monole <- estimateDispersions(HSC_monole) 
pData(HSC_monole)$Total_mRNAs <- Matrix::colSums(exprs(HSC_monole))
upper_bound <- 10^(mean(log10(pData(HSC_monole)$Total_mRNAs)) + 2*sd(log10(pData(HSC_monole)$Total_mRNAs)));lower_bound <- 10^(mean(log10(pData(HSC_monole)$Total_mRNAs)) - 2*sd(log10(pData(HSC_monole)$Total_mRNAs)))
qplot(Total_mRNAs, data=pData(HSC_monole), color=stage, geom="density") +
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound)


# this will be very slow, almost an hour??? VERY SLOW (prbrbly bc sample size)
HSC_monole <- reduceDimension(HSC_monole, pseudo_expr = 0, norm_method = "none")
HSC_monole <- orderCells(HSC_monole)
plot_cell_trajectory(HSC_monole)

