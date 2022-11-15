
library(tidyverse)

driver_gene_nonsyno_trinuc <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/driver_gene_nonsyno_trinuc.tsv')

hg19_trinuc_count <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/hg19_trinucleotide_count.tsv')
hg19_mut_trinuc_count <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/hg19_mut_trinuc_count.tsv')

comb_mutsig_mut_dat <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_mutsig_mut_dat.tsv')
comb_driver_gene <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')

mutsig_order <- comb_mutsig_mut_dat %>% distinct(Attri_Signature) %>% pull(Attri_Signature) %>% str_sort(numeric = TRUE)

pcawg_mutsig_mut_dat <- comb_mutsig_mut_dat %>% filter(Project %in% 'PCAWG')
# 20,970,058
pcawg_sam_ttype <- pcawg_mutsig_mut_dat %>% distinct(ICGC_abbr_top, Sample)
# 1,935 samples

driver_gene_nonsyno_trinuc_count <- driver_gene_nonsyno_trinuc %>% count(tricontext, name = 'count') %>% rename(mut_trinuc = tricontext)

# limited to PCAWG

# select signature-cancer_type pair with cancer gene mutations

signature_ttype_with_dri <- pcawg_mutsig_mut_dat %>% 
  semi_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top')) %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site'))

mutsig_ttype_dri_num <- signature_ttype_with_dri %>% count(Attri_Signature, ICGC_abbr_top, name = 'dri_count')

# number of signature mutations in each cancer type

mutsig_ttype_mutnum <- pcawg_mutsig_mut_dat %>% 
  count(Attri_Signature, ICGC_abbr_top, name = 'total_count') %>% 
  inner_join(mutsig_ttype_dri_num, by = c('Attri_Signature', 'ICGC_abbr_top'))
# 108 paired signature-ttype


pcawg_mutsig_trinuc <- pcawg_mutsig_mut_dat %>% 
  select(mut_id, Attri_Signature, tricontext) %>% 
  separate(mut_id, into = c('Sample', 'chr', 'pos', 'ref', 'alt'), sep = ':') %>% 
  mutate(trinuc_1 = str_sub(tricontext, 1, 1), 
         trinuc_3 = str_sub(tricontext, 7, 7), 
         trinuc = str_c(trinuc_1, ref, trinuc_3), 
         mut_trinuc = str_c(trinuc_1, '[', ref, '>', alt, ']', trinuc_3)) %>% 
  select(-(tricontext:trinuc_3)) %>% 
  left_join(pcawg_sam_ttype, by = 'Sample') %>% 
  select(ICGC_abbr_top, everything())

vroom::vroom_write(pcawg_mutsig_trinuc, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/pcawg_mutsig_trinuc.tsv')

sig_dri_test <- function(sig, ttype){
  
  # e.g., SBS5 in Panc-AdenoCA
  
  # sig <- 'SBS5'
  # ttype <- 'Panc-AdenoCA'
  
  # get mut probability accounting for base change and 5', 3' context
  
  mut_dat <- pcawg_mutsig_trinuc %>% filter(ICGC_abbr_top %in% all_of(ttype), Attri_Signature %in% all_of(sig))
  mut_count_dat <- mutsig_ttype_mutnum %>% filter(ICGC_abbr_top %in% all_of(ttype), Attri_Signature %in% all_of(sig))
  
  mut_trinuc_prob_dat <- mut_dat %>% 
    count(trinuc, mut_trinuc) %>% 
    left_join(hg19_trinuc_count, by = 'trinuc') %>% 
    mutate(prob = n/count, 
           prob = prob/sum(prob)) %>% # normalize
    select(mut_trinuc, prob)
  
  # whole genome mut prob
  
  whole_prob <- hg19_mut_trinuc_count %>% 
    left_join(mut_trinuc_prob_dat, by = 'mut_trinuc') %>% 
    drop_na(prob) %>% 
    mutate(pseudo_prob = count*prob)
  
  # cancer gene non-synonymous mut prob
  
  dri_ny_prob <- driver_gene_nonsyno_trinuc_count %>% 
    left_join(mut_trinuc_prob_dat, by = 'mut_trinuc') %>% 
    drop_na(prob) %>% 
    mutate(pseudo_prob = count*prob)
  
  # probability of a mutation located in cancer gene non-synonymous mutations
  
  ny_prob <- sum(dri_ny_prob$pseudo_prob)/(sum(whole_prob$pseudo_prob) - sum(dri_ny_prob$pseudo_prob))
  
  # binom.test
  
  # test_p <- binom.test(mut_count_dat$dri_count, mut_count_dat$total_count, p = ny_prob)$p.value
  test_p <- binom.test(mut_count_dat$dri_count, mut_count_dat$total_count, p = ny_prob, alternative = 'greater')$p.value
  
  mut_count_dat_output <- mut_count_dat %>% 
    mutate(exp_dri_count = round(total_count * ny_prob), p_val = test_p) %>% 
    dplyr::rename(obs_dri_count = dri_count) %>% 
    select(-Attri_Signature, -ICGC_abbr_top)
  
}

# perform test

mutsig_dri_test_res <- mutsig_ttype_mutnum %>% 
  select(Attri_Signature, ICGC_abbr_top) %>% 
  mutate(data = map2(Attri_Signature, ICGC_abbr_top, ~sig_dri_test(.x, .y))) %>% 
  unnest(data) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH'))

# write_tsv(mutsig_dri_test_res, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_driver_contribution_analysis/mutsig_dri_test_res_two_sided.tsv')
write_tsv(mutsig_dri_test_res, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_driver_contribution_analysis/mutsig_dri_test_res.tsv')
# 8 significant
