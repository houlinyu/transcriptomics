library(Biostrings)
library(tidyverse)

all <- read_csv('normalized_counts_all.csv')
id <- read_csv('tps_cluster_id.csv')
all %>% filter(X1 %in% id$Ref) -> tps

write_csv(tps, 'tps_normalized_counts.csv')
write_lines(id$Ref, 'tps_id_all.txt')

