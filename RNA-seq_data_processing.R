#reads_count
library(Rsubread)
file_vector <- list.files(path="/project/uma_lijun_ma/Hyu/RNASeq/TopHatOut/single_end", 
                          full.names=TRUE, recursive=FALSE)
list_out <- featureCounts(files = file_vector, 
            annot.ext = '/project/uma_lijun_ma/Hyu/RNASeq/TopHatOut/Magnaporthe_oryzae.MG8.46.gtf',
            isGTFAnnotationFile = T)
#write.csv(list_out$count)



####DESeq2_pipeline
# readin_data
library(DESeq2)
countdata <- as.matrix(read.csv("matrix.csv", header=TRUE, row.names=1))
coldata <- as.matrix(read.csv("coldata.csv"))

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ Condition)
# output_normalized_counts
ddsf <- estimateSizeFactors(dds)
normalized_counts <- counts(ddsf, normalized=TRUE)
write.csv(normalized_counts, 'normalized_counts_all.csv')

# pre-filtering
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]


###DESeq
deseq <- DESeq(dds)
resultsNames(deseq) ## what is this?

##1.bacteria
Condition_bacteria_C33h_vs_bacteria_control3h <- 
  results(deseq, contrast = c('Condition','bacteria_C33h','bacteria_control3h'))
Condition_bacteria_DCA3h_vs_bacteria_control3h <- 
  results(deseq, contrast = c('Condition','bacteria_DCA3h','bacteria_control3h'))
Condition_bacteria_C39h_vs_bacteria_control9h <- 
  results(deseq, contrast = c('Condition','bacteria_C39h','bacteria_control9h'))
Condition_bacteria_DCA9h_vs_bacteria_control9h <- 
  results(deseq, contrast = c('Condition','bacteria_DCA9h','bacteria_control9h'))

##2.conidia
Condition_conidia_macroconidia_vs_conidia_vegetativehyphae <- 
  results(deseq, contrast = c('Condition','conidia_macroconidia','conidia_vegetativehyphae'))
Condition_conidia_microconidia_vs_conidia_vegetativehyphae <- 
  results(deseq, contrast = c('Condition','conidia_microconidia','conidia_vegetativehyphae'))
Condition_conidia_microconidia_vs_conidia_macroconidia <- 
  results(deseq, contrast = c('Condition','conidia_microconidia','conidia_macroconidia'))

##3.Z2C6AraR1
Condition_Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h <- 
  results(deseq, contrast = c('Condition','Z2C6AraR1_ara2h','Z2C6AraR1_fru2h'))
Condition_Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h <- 
  results(deseq, contrast = c('Condition','Z2C6AraR1_ara4h','Z2C6AraR1_fru2h'))
Condition_Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h <- 
  results(deseq, contrast = c('Condition','Z2C6AraR1_ara8h','Z2C6AraR1_fru2h'))

##4.Z2C6Gpf1Cnf2
Condition_Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT <- 
  results(deseq, contrast = c('Condition','Z2C6Gpf1Cnf2_Gpf1','Z2C6Gpf1Cnf2_WT'))
Condition_Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT <- 
  results(deseq, contrast = c('Condition','Z2C6Gpf1Cnf2_Cnf2','Z2C6Gpf1Cnf2_WT'))

##5.bHTHCrf1
Condition_bHTHCrf1_Crf1_vs_bHTHCrf1_WT <- 
  results(deseq, contrast = c('Condition','bHTHCrf1_Crf1','bHTHCrf1_WT'))

##6.bHTHCrf1
Condition_MoRfx1_MoRfx1_vs_MoRfx1_WT <- 
  results(deseq, contrast = c('Condition','MoRfx1_MoRfx1','MoRfx1_WT'))

test_all <- list(Condition_bacteria_C33h_vs_bacteria_control3h,
              Condition_bacteria_DCA3h_vs_bacteria_control3h,
              Condition_bacteria_C39h_vs_bacteria_control9h,
              Condition_bacteria_DCA9h_vs_bacteria_control9h,
              Condition_conidia_macroconidia_vs_conidia_vegetativehyphae,
              Condition_conidia_microconidia_vs_conidia_vegetativehyphae,
              Condition_conidia_microconidia_vs_conidia_macroconidia,
              Condition_Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h,
              Condition_Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h,
              Condition_Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h,
              Condition_Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT,
              Condition_Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT,
              Condition_bHTHCrf1_Crf1_vs_bHTHCrf1_WT,
              Condition_MoRfx1_MoRfx1_vs_MoRfx1_WT)

test_all_id <- c('bacteria_C33h_vs_bacteria_control3h',
                 'bacteria_DCA3h_vs_bacteria_control3h',
                 'bacteria_C39h_vs_bacteria_control9h',
                 'bacteria_DCA9h_vs_bacteria_control9h',
                 'conidia_macroconidia_vs_conidia_vegetativehyphae',
                 'conidia_microconidia_vs_conidia_vegetativehyphae',
                 'conidia_microconidia_vs_conidia_macroconidia',
                 'Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h',
                 'Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h',
                 'Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h',
                 'Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT',
                 'Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT',
                 'bHTHCrf1_Crf1_vs_bHTHCrf1_WT',
                 'MoRfx1_MoRfx1_vs_MoRfx1_WT')


#tips_cluster <- read.csv('tps_cluster_id.csv')

library(tidyverse)
#tip_cluster_withFCtest <- tips_cluster
for (i in 1:14){
assign(paste(test_all_id[i], '_log2FC', sep = ''),
as_tibble(test_all[[i]][,c('log2FoldChange','pvalue')])$log2FoldChange)
assign(paste((test_all_id[i]), '_pvalue', sep = ''),
    as_tibble(test_all[[i]][,c('log2FoldChange','pvalue')])$pvalue)                              
}

testandp <- data.frame(test_all,
           bacteria_C33h_vs_bacteria_control3h_log2FC, bacteria_C33h_vs_bacteria_control3h_pvalue,
           bacteria_C39h_vs_bacteria_control9h_log2FC, bacteria_C39h_vs_bacteria_control9h_pvalue,
           bacteria_DCA3h_vs_bacteria_control3h_log2FC, bacteria_DCA3h_vs_bacteria_control3h_pvalue,
           bacteria_DCA9h_vs_bacteria_control9h_log2FC, bacteria_DCA9h_vs_bacteria_control9h_pvalue,
           conidia_macroconidia_vs_conidia_vegetativehyphae_log2FC, conidia_macroconidia_vs_conidia_vegetativehyphae_pvalue,
           conidia_microconidia_vs_conidia_macroconidia_log2FC, conidia_microconidia_vs_conidia_macroconidia_pvalue,
           conidia_microconidia_vs_conidia_vegetativehyphae_log2FC, conidia_microconidia_vs_conidia_vegetativehyphae_pvalue,
           Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_log2FC, Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_pvalue,
           Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_log2FC, Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_pvalue,
           Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_log2FC, Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_pvalue,
           Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_log2FC, Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_pvalue,
           Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_log2FC, Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_pvalue,
           bHTHCrf1_Crf1_vs_bHTHCrf1_WT_log2FC, bHTHCrf1_Crf1_vs_bHTHCrf1_WT_pvalue,
           MoRfx1_MoRfx1_vs_MoRfx1_WT_log2FC, MoRfx1_MoRfx1_vs_MoRfx1_WT_pvalue)
           
write.csv(testandp, 'testandp.csv')
worktab <- read.csv('testandp.csv')

worktab$bacteria_C33h_vs_bacteria_control3h_pvalue[!(worktab$bacteria_C33h_vs_bacteria_control3h_pvalue < 0.05)] = 0
worktab$bacteria_C39h_vs_bacteria_control9h_pvalue[!(worktab$bacteria_C39h_vs_bacteria_control9h_pvalue < 0.05)] = 0
worktab$bacteria_DCA3h_vs_bacteria_control3h_pvalue[!(worktab$bacteria_DCA3h_vs_bacteria_control3h_pvalue < 0.05)] = 0
worktab$bacteria_DCA9h_vs_bacteria_control9h_pvalue[!(worktab$bacteria_DCA9h_vs_bacteria_control9h_pvalue < 0.05)] = 0
worktab$conidia_macroconidia_vs_conidia_vegetativehyphae_pvalue[!(worktab$conidia_macroconidia_vs_conidia_vegetativehyphae_pvalue < 0.05)] = 0
worktab$conidia_microconidia_vs_conidia_macroconidia_pvalue[!(worktab$conidia_microconidia_vs_conidia_vegetativehyphae_pvalue < 0.05)] = 0
worktab$conidia_microconidia_vs_conidia_vegetativehyphae_pvalue[!(worktab$conidia_microconidia_vs_conidia_vegetativehyphae_pvalue < 0.05)] = 0
worktab$Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_pvalue[!(worktab$Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_pvalue < 0.05)] = 0
worktab$Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_pvalue[!(worktab$Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_pvalue < 0.05)] = 0
worktab$Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_pvalue[!(worktab$Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_pvalue < 0.05)] = 0
worktab$Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_pvalue[!(worktab$Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_pvalue < 0.05)] = 0
worktab$Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_pvalue[!(worktab$Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_pvalue < 0.05)] = 0
worktab$bHTHCrf1_Crf1_vs_bHTHCrf1_WT_pvalue[!(worktab$bHTHCrf1_Crf1_vs_bHTHCrf1_WT_pvalue < 0.05)] = 0
worktab$MoRfx1_MoRfx1_vs_MoRfx1_WT_pvalue[!(worktab$MoRfx1_MoRfx1_vs_MoRfx1_WT_pvalue < 0.05)] = 0

write.csv(worktab, 'worktab2.csv')

worktab <- read.csv('testandp.csv')
worktab1 <- read.csv('worktab2.csv')

worktab1$bacteria_C33h_vs_bacteria_control3h_log2FC[!(worktab$bacteria_C33h_vs_bacteria_control3h_pvalue < 0.05)] = 0
worktab1$bacteria_C39h_vs_bacteria_control9h_log2FC[!(worktab$bacteria_C39h_vs_bacteria_control9h_pvalue < 0.05)] = 0
worktab1$bacteria_DCA3h_vs_bacteria_control3h_log2FC[!(worktab$bacteria_DCA3h_vs_bacteria_control3h_pvalue < 0.05)] = 0
worktab1$bacteria_DCA9h_vs_bacteria_control9h_log2FC[!(worktab$bacteria_DCA9h_vs_bacteria_control9h_pvalue < 0.05)] = 0
worktab1$conidia_macroconidia_vs_conidia_vegetativehyphae_log2FC[!(worktab$conidia_macroconidia_vs_conidia_vegetativehyphae_pvalue < 0.05)] = 0
worktab1$conidia_microconidia_vs_conidia_macroconidia_log2FC[!(worktab$conidia_microconidia_vs_conidia_vegetativehyphae_pvalue < 0.05)] = 0
worktab1$conidia_microconidia_vs_conidia_vegetativehyphae_log2FC[!(worktab$conidia_microconidia_vs_conidia_vegetativehyphae_pvalue < 0.05)] = 0
worktab1$Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_log2FC[!(worktab$Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_pvalue < 0.05)] = 0
worktab1$Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_log2FC[!(worktab$Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_pvalue < 0.05)] = 0
worktab1$Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_log2FC[!(worktab$Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_pvalue < 0.05)] = 0
worktab1$Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_log2FC[!(worktab$Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_pvalue < 0.05)] = 0
worktab1$Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_log2FC[!(worktab$Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_pvalue < 0.05)] = 0
worktab1$bHTHCrf1_Crf1_vs_bHTHCrf1_WT_log2FC[!(worktab$bHTHCrf1_Crf1_vs_bHTHCrf1_WT_pvalue < 0.05)] = 0
worktab1$MoRfx1_MoRfx1_vs_MoRfx1_WT_log2FC[!(worktab$MoRfx1_MoRfx1_vs_MoRfx1_WT_pvalue < 0.05)] = 0

write.csv(worktab1, 'workingtab2.csv')




library(readxl)
tab3 <- read.csv('workingtab2.csv')
tab3$bacteria_C33h_vs_bacteria_control3h_log2FC[!(tab3$bacteria_C33h_vs_bacteria_control3h_log2FC == 0)] = 1
tab3$bacteria_C39h_vs_bacteria_control9h_log2FC[!(tab3$bacteria_C39h_vs_bacteria_control9h_log2FC == 0)] = 1
tab3$bacteria_DCA3h_vs_bacteria_control3h_log2FC[!(tab3$bacteria_DCA3h_vs_bacteria_control3h_log2FC == 0)] = 1
tab3$bacteria_DCA9h_vs_bacteria_control9h_log2FC[!(tab3$bacteria_DCA9h_vs_bacteria_control9h_log2FC == 0)] = 1
tab3$conidia_macroconidia_vs_conidia_vegetativehyphae_log2FC[!(tab3$conidia_macroconidia_vs_conidia_vegetativehyphae_log2FC == 0)] = 1
tab3$conidia_microconidia_vs_conidia_macroconidia_log2FC[!(tab3$conidia_microconidia_vs_conidia_vegetativehyphae_log2FC == 0)] = 1
tab3$conidia_microconidia_vs_conidia_vegetativehyphae_log2FC[!(tab3$conidia_microconidia_vs_conidia_vegetativehyphae_log2FC == 0)] = 1
tab3$Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_log2FC[!(tab3$Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_log2FC == 0)] = 1
tab3$Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_log2FC[!(tab3$Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_log2FC == 0)] = 1
tab3$Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_log2FC[!(tab3$Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_log2FC == 0)] = 1
tab3$Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_log2FC[!(tab3$Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_log2FC == 0)] = 1
tab3$Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_log2FC[!(tab3$Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_log2FC == 0)] = 1
tab3$bHTHCrf1_Crf1_vs_bHTHCrf1_WT_log2FC[!(tab3$bHTHCrf1_Crf1_vs_bHTHCrf1_WT_log2FC == 0)] = 1
tab3$MoRfx1_MoRfx1_vs_MoRfx1_WT_log2FC[!(tab3$MoRfx1_MoRfx1_vs_MoRfx1_WT_log2FC == 0)] = 1

mutate(tab3, biotic = bacteria_C33h_vs_bacteria_control3h_log2FC + 
         bacteria_C33h_vs_bacteria_control3h_log2FC + 
         bacteria_DCA3h_vs_bacteria_control3h_log2FC +
         bacteria_DCA9h_vs_bacteria_control9h_log2FC,
       lifecycle = conidia_macroconidia_vs_conidia_vegetativehyphae_log2FC +
         conidia_microconidia_vs_conidia_macroconidia_log2FC +
         conidia_microconidia_vs_conidia_vegetativehyphae_log2FC,
       nutrition = Z2C6AraR1_ara2h_vs_Z2C6AraR1_fru2h_log2FC +
         Z2C6AraR1_ara4h_vs_Z2C6AraR1_fru2h_log2FC +
         Z2C6AraR1_ara8h_vs_Z2C6AraR1_fru2h_log2FC,
       TF = Z2C6Gpf1Cnf2_Gpf1_vs_Z2C6Gpf1Cnf2_WT_log2FC +
         Z2C6Gpf1Cnf2_Cnf2_vs_Z2C6Gpf1Cnf2_WT_log2FC +
         bHTHCrf1_Crf1_vs_bHTHCrf1_WT_log2FC +
         MoRfx1_MoRfx1_vs_MoRfx1_WT_log2FC
       ) -> tab4

tab4 %>% select(X, biotic, lifecycle, nutrition,TF) -> tab5




#13470 gene   ### with 3664 no expression.
#12617 has transcript support

write.csv(tab5, 'workingtab5.csv')

library(tidyverse)
tab6 <- read.csv('workingtab5.csv', header = T)
tab6 %>% filter(TF != 1000) -> tab7
mean(tab7$biotic)
##0.6713957/4
mean(tab7$lifecycle)
##1.466038
mean(tab7$nutrition)
##0.7583419
mean(tab7$TF)
##1.272648

hist(tab7$biotic)
##set >= 3/4

hist(tab7$lifecycle)
##set >= 3

hist(tab7$nutrition)
##set >= 3

hist(tab7$TF)
##set >= 4




tab7 <- read_excel('workingtab7.xls')
filter(tab7, `Sum/14` != 0) -> tab8
write_csv(tab8, 'workingtab8.csv')


tab8 <- read.csv('workingtab8.csv')
tab8 %>% filter(`Sum.14` > 1) -> tab9

write_csv(tab9, 'workingtab9.csv')


sigratio <- read_csv('significant_ratio2.csv')
sigratio <- sigratio[-8,]

sigratio_long <- pivot_longer(sigratio, cols = 2:3, names_to = 'count')




sigratio_long$tps <- factor(sigratio_long$tps, levels = c('tps6','tps9','dtps1','tps7','tps10',

                                                                                                    'dtps2','tps8'))

ggplot(data = sigratio_long, aes(x = tps, y=value, fill = count)) +
  geom_bar(stat="identity") + 
  scale_y_continuous(limits = c(0,20), breaks = c(0,5,10,15,20), 
                     expand = c(0,0)) +
  xlab(label = 'terpine synthase clusters') +
  ylab(label = 'count') +
  scale_fill_manual(values=c('darkgreen','red'))
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))






#library(reshape)
#x <- melt(sigratio)


sigratio %>% 
ggplot(aes(x=tps)) + 
  geom_bar(aes(fill=factor(c('total','POI'))),position="stack") 


# position默认，堆叠








###other_codes
vsd <- vst(dds, blind=FALSE)
cds <- DESeq(dds)
vst <- getVarianceStabilizedData(cds)
df <- vst
df <- df - rowMeans(df)
df <- as.data.frame(df)

###PCA_plot
# by_package_builin
dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
plotPCA(DESeqTransform(se),intgroup="V2")
# by_ggplot2
library(tidyverse)
pcaData <- plotPCA(vsd, intgroup="V2", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=V2)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

###k-means clustering
colnames(df)
res_km <- kmeans(counts(dds, normalized=TRUE), 7, iter.max = 50)
res_km$size
as.numeric(res_km$cluster[as.vector(rnaid$Ref)]) -> cluster
rnaid <- read.csv('rnaseqid2.csv')
mutate(rnaid, cluster = cluster)



# k-means clustering and plot (not tested)
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













