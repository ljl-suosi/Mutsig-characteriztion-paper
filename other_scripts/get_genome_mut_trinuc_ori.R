
library(tidyverse)

hg19_trinuc_count <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/hg19_trinucleotide_count.tsv')

hg19_mut_trinuc_count <- hg19_trinuc_count %>% 
  mutate(trinuc = map_chr(trinuc, ~str_flatten(unlist(str_split(., '')), ','))) %>% 
  separate(trinuc, into = c('first', 'mid', 'last'), sep = ',') %>% 
  mutate(alt_base = pmap_chr(list(first, mid, last), function(x, y, z) str_flatten(str_c(x, '[', str_c(y, setdiff(c('A', 'T', 'C', 'G'), y), sep = '>'), ']', z), collapse = ','))) %>% 
  select(mut_trinuc = alt_base, count) %>% 
  separate_rows(mut_trinuc, sep = ',')

write_tsv(hg19_mut_trinuc_count, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/hg19_mut_trinuc_count.tsv')

