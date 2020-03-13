library(tidyverse)
tab6 <- read.csv('workingtab5.csv', header = T)
tab6 %>% filter(TF != 1000) -> tab7
head(tab7)
tab7 <- filter(tab7, X %in% c('MGG_01674',
                              'MGG_01677',
                              'MGG_01687',
                              'MGG_01689',
                              'MGG_01690',
                              'MGG_01700',
                              'MGG_01701',
                              'MGG_01702',
                              'MGG_01704',
                              'MGG_01708',
                              'MGG_01711',
                              'MGG_14813',
                              'MGG_14814',
                              'MGG_01720',
                              'MGG_01721',
                              'MGG_01722',
                              'MGG_01723',
                              'MGG_01724',
                              'MGG_01728'))


mean(tab7$biotic)
##0.6713957/4
mean(tab7$lifecycle)
##1.466038
mean(tab7$nutrition)
##0.7583419
mean(tab7$TF)
##1.272648








hist(tab7$biotic, col = 'black', main = 'Biotic', xlab = '', ylab = '')
##set >= 3/4
hist(tab7$lifecycle, col = 'black', main = 'Lifecycle', xlab = '', ylab = '')
##set >= 3
hist(tab7$nutrition, col = 'black', main = 'Nutrition', xlab = '', ylab = '')
##set >= 3
hist(tab7$TF, col = 'black', main = 'TF', xlab = '', ylab = '')
##set >= 4

nrow(tab7)
nrow(tab7 %>% mutate(sum = biotic + lifecycle + nutrition + TF) %>% 
  filter(biotic > 2, lifecycle > 2, nutrition > 2, TF > 3))

tab7 %>% mutate(sum = biotic + lifecycle + nutrition + TF) %>% 
  filter(biotic > 1, lifecycle > 1, nutrition > 1, TF > 1) -> tab8


print(tab8$X)
write_lines(tab8$X, '1112.txt')


#9546/12617





#caculate reads count
library(readxl)
stat_sum <- read_xls('stat_summary.xls', sheet = 'Sheet1')
stat_sum %>% 
  ggplot(aes(x = tps, y = Normalized)) +
  #geom_boxplot() +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center',outlier.shape = NA,
               binwidth = 0.5, ) +
  theme(panel.grid.major =element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))






tpssum <- rep(NA, 7)
tpsid <- unique(stat_sum$tps)


tpssum <- list(rep(NA, 7))
for (i in 1:7){
stat_sum %>% filter(tps == tpsid[i]) -> sum1
tpssum[[i]] <- sum1$Sum/sum1$`Transcript Length`
}






tpslen <- rep(NA, 7)
for (i in 1:7){
  stat_sum %>% filter(tps == tpsid[i]) -> sum1
  tpslen[i] <- sum(sum1$`Transcript Length`)
}


tibble(TPS = tpsid, MedianExpression = tpssum) -> tpstib
tpstib$TPS <- c('TPS2', 'TPS6','TPS1','TPS7','TPS9','TPS8','TPS10')
tpstib$TPS <- factor(tpstib$TPS,
                     levels = c('TPS1','TPS2',
                                'TPS6','TPS7',
                                'TPS8','TPS9','TPS10'))
tpstib %>% 
  ggplot(aes(TPS, MedianExpression)) +
  geom_col() +
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))




tpssum/tpslen
tpssum


#write_lines(stat_sum$X1, 'tps_gene.txt')

















