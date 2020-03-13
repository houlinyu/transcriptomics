
library(Rsubread)
file_vector <- list.files(path="/project/uma_lijun_ma/Hyu/RNASeq/TopHatOut/single_end", full.names=TRUE, recursive=FALSE)
list_out <- featureCounts(files = file_vector, 
            annot.ext = '/project/uma_lijun_ma/Hyu/RNASeq/TopHatOut/Magnaporthe_oryzae.MG8.46.gtf',
            isGTFAnnotationFile = T)


library(DESeq2)
countdata <- read.csv("matrix.csv", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)


coldata <- as.matrix(read.csv("coldata.csv", header = F))
condition <- as.factor(coldata[,2])
  
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ V2)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts[as.vector(rnaid$Ref),], 'normalized_counts_interests.csv')

colnames(normalized_counts)
DESeq(dds)
vsd <- vst(dds, blind=FALSE)
library(tidyverse)



pcaData <- plotPCA(vsd, intgroup="V2", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=V2)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
plotPCA(DESeqTransform(se),intgroup="V2")


cds <- DESeq(dds)
vst <- getVarianceStabilizedData(cds)
df <- vst
df <- df - rowMeans(df)
df <- as.data.frame(df)

colnames(df)
res_km <- kmeans(counts(dds, normalized=TRUE), 7, iter.max = 50)
res_km$size


as.numeric(res_km$cluster[as.vector(rnaid$Ref)]) -> cluster
rnaid <- read.csv('rnaseqid2.csv')
mutate(rnaid, cluster = cluster)



data_plot <- data.table(melt(data.table(class = as.factor(res_km$cluster), df)))
data_plot[, Time := rep(1:ncol(df), each = nrow(df))]
data_plot[, ID := rep(1:nrow(df), ncol(df))]
head(data_plot)
# prepare centroids
centers <- data.table(melt(res_km$centers))
setnames(centers, c("Var1", "Var2"), c("class", "Time"))
centers[, ID := class]
centers[, gr := as.numeric(as.factor(Time))]
head(centers)
head(data_plot)
# plot the results
ggplot(data_plot, aes(variable, value, group = ID)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.65) +
  
  geom_line(data = centers, aes(gr, value),
            color = "firebrick1", alpha = 0.80, size = 1.2) +
  labs(x = "Time", y = "Load (normalised)") +
  theme_bw()




library(DESeq2)
countdata <- read.csv("normalized_counts_interests.csv", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)

kmeans(countdata, 7) -> cluster_1
as.numeric(cluster_1$cluster[as.vector(rnaid$Ref)]) -> cluster
rnaid <- read.csv('rnaseqid2.csv')
mutate(rnaid, cluster = cluster) %>% group_by(tps) %>% summarise()




##                  tps6      tpss
##    in             8        17
##    not in         11       63 


matrix(c(8,17,11,63), nrow = 2)
fisher.test(matrix(c(8,17,11,63), nrow = 2))
cluster_1$size


coldata <- as.matrix(read.csv("coldata.csv", header = F))
condition <- as.factor(coldata[,2])
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ V2)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts[as.vector(rnaid$Ref),], 'normalized_counts_interests.csv')

colnames(normalized_counts)

