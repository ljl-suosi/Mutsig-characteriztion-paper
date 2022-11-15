
library(BSgenome.Hsapiens.UCSC.hg19)

hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_names <- str_c('chr', c(1:22, 'X', 'Y'))

# 64 trinucleotides

nucleotides <- c('A', 'T', 'C', 'G')
trinucleotides <- c()

for(i in nucleotides){
  for(j in nucleotides){
    for(k in nucleotides){
      trinucleotides <- c(trinucleotides, str_c(i, j, k))
    }
  }
}

chrom_trinuc_count <- list()

for(i in chrom_names){
  
  # i <- 'chr2'
  
  sin_chr <- hg19_genome[[i]]
  
  sin_count <- tibble(trinuc = trinucleotides) %>% 
    mutate(count = map_dbl(trinuc, ~countPattern(., sin_chr))) %>% 
    dplyr::rename(!!i := count)
  
  chrom_trinuc_count[[i]] <- sin_count
  
  print(str_c(i, 'finished!', sep = ' '))
}


hg19_trinucleotide_count <- chrom_trinuc_count %>% 
  purrr::reduce(inner_join, by = 'trinuc') %>% 
  pivot_longer(-trinuc, names_to = 'chrom', values_to = 'count') %>% 
  group_by(trinuc) %>% 
  summarise(count = sum(count))
# total count: 2,861,326,455

write_tsv(hg19_trinucleotide_count, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/hg19_trinucleotide_count.tsv')

