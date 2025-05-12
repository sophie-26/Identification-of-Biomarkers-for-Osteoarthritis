# import libraries
library(readxl)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(openxlsx)
library(dplyr) 
library(gplots)
library(readr)

# set working directory
setwd("/Users/sophiezhuang/Documents")

# read data 
ctrl_pat = read_excel("OA summer research project data GSE114007_raw_counts.xlsx", "Normal")
OA_pat = read_excel("OA summer research project data GSE114007_raw_counts.xlsx", "OA")
data = bind_cols(OA_pat, ctrl_pat[,-c(1)])
# take out outlier from data
data = data %>% select(-18, -27)

# export the combined data into excel file
write.xlsx(data, "combined research project data.xlsx")

# data without gene name
count = (data[,-c(1)])
dim(count)

# gene info 
genes = data[,1]
dim(genes)

# plot expression value vs the genes
boxplot(log2(count+1),las=2)
ggFormat = melt(count)
ggplot(ggFormat, aes(x=variable, y =log2(value+1))) + geom_boxplot() 


# sample info 
samples = colnames(data) # column name is the name of samples 
group = c(rep("OA", times = 19), rep("HC", times = 17)) # vector of healthy as HC and disease as OA
sampleinfo = data.frame(samples[-c(1)],group) # get correct data and combine into data frame
sampleinfo
dim(sampleinfo)

# create DGEList
all_cds <- DGEList(count, group =sampleinfo$group,genes=genes )
names( all_cds )
dim(all_cds)

# filter
all_keep <- rowSums(cpm(all_cds)>30) >= 6 
all_count_filtered <- all_cds[all_keep,]
dim(all_cds) # before 
dim(all_count_filtered) # after

all_count_filtered = all_cds
# visualizing expresion profile
cpmRNA<-cpm(all_count_filtered)
cpmRNA<-as.data.frame(cpmRNA)
# plot 1: log2counts
boxplot(log2(all_count_filtered$counts+1), las=2, main="log2 counts") # expand plots UI to load
# plot 2: log2CPM
boxplot(log2(cpmRNA+1),las=2,main="log2CPM") # expand plots UI to load

# TMM Normalization
cf_all <- calcNormFactors(all_count_filtered, method="TMM")

# export normalized data
normalized_data = bind_cols(data.frame(cf_all$genes), data.frame(cf_all$counts))
write_csv(normalized_data, file = file.path("Normalized OA data.csv"))



# NEW
# read data
correct_normal = read.table("GSE114007_normal_normalized.counts.txt", header = TRUE, sep = "\t")
correct_OA = read.table("GSE114007_OA_normalized.counts.txt",  header = TRUE, sep = "\t")
correct_data = bind_cols(correct_OA, correct_normal)
export_correct_data = correct_data[order(correct_data$symbol...1),]
write_csv(export_correct_data, file = file.path("OA vs HC correct combined normalized counts.csv"))








# Multidimensional scaling
plotMDS(cf_all, method="bcv", col=as.numeric(cf_all$samples$group))
#we start to see samples that are similar to each other

# PCA dimensional reduction
par(bty = 'n')
FCpca <- prcomp(t(cpm(cf_all, log=TRUE))) # compute principal components
head(FCpca$x)
# plot PCA of RNA-seq
plot(FCpca$x[,1], FCpca$x[,2])  ## plot  PC1 and PC2
text(FCpca$x[,"PC1"], FCpca$x[,"PC2"], rownames(FCpca$x), cex=0.7, pos=1, col=rep(c("red","blue"),each = 19))
title("PCA of RNA-seq", adj = 0, line = 0.1)

# make group column
pcaValues<-data.frame(FCpca$x[,1:2])
pcaValues$group<-rownames(pcaValues)
pcaValues$group[grep("^OA",pcaValues$group)]<-"OA"
pcaValues$group[grep("^.ormal",pcaValues$group)]<-"HC"



# FINISHED PLOT: sample relations based on multidimensional scaling
ggplot(pcaValues, aes(PC1, PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1")) +
  ylab(paste0("PC2")) + 
  coord_fixed() +
  scale_color_manual(values = c("blue", "red"))

# convert edgeR to DESeq2
edgeR.DDS <- DESeqDataSetFromMatrix(countData = round(cf_all$counts), 
                                    colData = cf_all$samples, design = ~group)
# rlog Apply a 'regularized log' transformation
transform.edgeR.DDS <- rlog(edgeR.DDS, blind = TRUE)
# plot PCA
pcaData <- plotPCA(transform.edgeR.DDS, intgroup = c("group"), returnData = TRUE, ntop = 1000)  
percentVar <- round(100 * attr(pcaData, "percentVar"))



# FINISHED PLOT: alternative to getting sample relations based on multidimensional scaling
ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  scale_color_manual(values = c("blue", "red"))

# Estimate Dispersion
cf_all <-estimateDisp(cf_all) 
cf_all <- estimateCommonDisp(cf_all)
cf_all <- estimateTagwiseDisp(cf_all)   



# FINISHED PLOT: dispersion of data --> Larger dispersions = higher variance between replicates
plotBCV(cf_all)  
# coefficient of biological variance (TOO HIGH)
sqrt(cf_all$common.dispersion) # true abundance for each gene can vary up or down by 56.8 % between replicates

# Testing for DE genes

cf_et <- exactTest(cf_all, pair = c("HC","OA"))
DE<- topTags( cf_et , n = nrow( cf_et$table ) )$table
head(DE)
de1 <- decideTestsDGE(cf_et, adjust.method="BH", p.value=0.01) # lower p value to get more significant
summary(de1)
# plot
de1tags12 <- rownames(cf_et)[as.logical(de1)] 



# FINISHED PLOT: average expression of genes and their full change --> shows up reg or down reg
plotSmear(cf_et, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")

PTH_genes = DE[c(grep("PTH", DE$symbol)),]
write_csv(PTH_genes, file = file.path("PTH_genes.csv"))
sig_PTH_genes = subset(PTH_genes, PValue<0.05)





# looking at expression value of one gene
cpmTMM<-cpm(cf_all)
cpmTMM<-data.frame(cpmTMM)
cpmTMM<-cbind(cpmTMM,cf_all$genes) # add symbol to the normalized expression matrix
cpmTMM<-data.frame(cpmTMM)
head(cpmTMM);tail(cpmTMM)


boxplot_data = cpmTMM[c(16866,16867,16870),]
dim(boxplot_data)
boxplot_data = data.frame(boxplot_data)
rownames(boxplot_data) = c("PTH","PTH1R","PTHLH")

# Heatmap of DE genes
# up reg top 5
colnames(cpmTMM)
mypal <- colorRampPalette(c("blue","black", "yellow"))
genes1<-PTH_genes$symbol[1:5]

# FINISHED PLOT: blue = lowly expressed, yellow = highly expressed
heatmap.2(log2(as.matrix(cpmTMM[cpmTMM$symbol%in%c(genes1),c(1:36)])+1),
          labRow=genes1, key=T, trace="none",symkey = T,symbreaks=T,
          Colv=FALSE, 
          breaks=seq(-2,+2,0.2),
          col=mypal(20), Rowv=T,scale="row",dendrogram ="row",
          cexCol=1,cexRow=1,
          margin=c(10, 10), main="OA vs HC PTH Genes")


