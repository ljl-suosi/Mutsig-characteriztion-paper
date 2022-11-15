
library(tidyverse)

comb_mutsig_we <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')
load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/driver_gs_status.Rdata')
comb_driver_gene <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')
comb_TMB <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_TMB.tsv')


comb_mutsig_we_filter <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% 
  ungroup()

# for each signature, compute the number of cancer type where the signature presents in over 5% samples

mutsig_freq <- comb_mutsig_we_filter %>% distinct(Signature, ICGC_abbr_top) %>% count(Signature)

# for each cancer gene, compute the number of cancer type where the gene presents in over 5% samples

dri_freq <- driver_gs_status %>% 
  pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status') %>% 
  inner_join(distinct(comb_mutsig_we_filter, ICGC_abbr_top, Sample), by = 'Sample') %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(mean(Status) >= 0.05) %>% 
  ungroup() %>% 
  distinct(ICGC_abbr_top, Hugo_Symbol) %>% 
  count(Hugo_Symbol)

# for each pathway gene, compute the number of cancer type where the gene presents in over 5% samples

path_gene_freq <- key_pathway_gene_sam_status %>% 
  pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status') %>% 
  inner_join(distinct(comb_mutsig_we_filter, ICGC_abbr_top, Sample), by = 'Sample') %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(mean(Status) >= 0.05) %>% 
  ungroup() %>% 
  distinct(ICGC_abbr_top, Hugo_Symbol) %>% 
  count(Hugo_Symbol)

# age ---------------------------------------------------------------------


comb_mutsig_we_age <- comb_mutsig_we_filter %>% 
  inner_join(select(comb_cli_data, Sample, age), by = 'Sample') %>% 
  mutate(id = str_c(ICGC_abbr_top, Signature, sep = ':'))

sig_fct_lev <- comb_mutsig_we_age %>% distinct(id) %>% pull(id)

age_mutsig_lmRob <- comb_mutsig_we_age %>% 
  mutate(id = factor(id, levels = sig_fct_lev)) %>%
  group_by(id) %>% 
  group_map(possibly(~summary(robust::lmRob(Exposure_TMB ~ age, data = .))$coefficients[2, ], set_names(c(NA, NA, NA, NA), c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')))) %>% 
  bind_rows() %>% 
  rename(estimate = Estimate, std.error = `Std. Error`, statistic = `t value`, p.value = `Pr(>|t|)`) %>% 
  mutate(id = sig_fct_lev) %>% 
  separate(id, into = c('ICGC_abbr_top', 'Signature'), sep = ':') %>% 
  select(ICGC_abbr_top, Signature, everything()) %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(q_val = p.adjust(p.value, method = 'BH')) %>% 
  ungroup()

age_mutsig_lmRob_signif <- age_mutsig_lmRob %>% filter(q_val < 0.05, q_val != 0)

# write_tsv(age_mutsig_lmRob_signif, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/age_mutsig_lmRob_signif.tsv')


# Signature     n
# <chr>     <int>
# SBS5         20
# SBS1         19
# SBS40         4
# SBS13         1
# SBS2          1
# SBS4          1
# SBS7a         1


# other clinical characteristics ------------------------------------------

source('/pub6/Temp/Liaojl/Code/common_R_functions/var_relate_als.R')


comb_cli_data_t <- comb_cli_data %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(TMB = ifelse(TMB > median(TMB, na.rm = TRUE), 'High_TMB', 'Low_TMB')) %>% 
  ungroup() %>% 
  semi_join(comb_mutsig_we_filter, by = 'Sample') %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(histological_type = fct_lump_min(histological_type, 15)) %>% # lump histological type whose factor level is less than 15 to other
  ungroup()
# 8,662

comb_mutsig_we_wide <- comb_mutsig_we %>% 
  select(-Exposure, -Exposure_TMB) %>% 
  pivot_wider(ICGC_abbr_top:Sample, names_from = Signature, values_from = Weight)


mutsig_cli_relate_result <- list()

for(i in unique(comb_mutsig_we_wide$ICGC_abbr_top)){
  
  # i <- 'Cervix'
  # i <- 'Panc-AdenoCA'
  
  sin_mutsig <- comb_mutsig_we_wide %>% 
    filter(ICGC_abbr_top %in% i) %>% 
    select_if(~mean(. > 0) >= 0.05)
  
  sin_mutsig_group <- comb_mutsig_we_filter %>% 
    filter(ICGC_abbr_top %in% i) %>% 
    group_by(Signature) %>% 
    mutate(group = ifelse(Weight > median(Weight), 'High', 'Low')) %>% 
    ungroup() %>% 
    mutate(group = factor(group, levels = c('Low', 'High'))) %>% 
    select(ICGC_abbr_top:Signature, group) %>% 
    pivot_wider(names_from = Signature, values_from = group)
  
  ana_sigs <- sin_mutsig %>% select(starts_with('SBS')) %>% colnames()
  
  # gender, histological_type, 'stage', 'grade', wgd_status, 'first_therapy_response', 'TMB'
  
  sin_result1 <- var_relate_als(sin_mutsig, comb_cli_data_t, ana_sigs, c('gender', 'histological_type', 'stage', 'grade', 'wgd_status', 'first_therapy_response', 'TMB'), 'Sample', 'Sample')
  
  # purity, ploidy
  
  sin_result2 <- var_relate_als(sin_mutsig_group, comb_cli_data_t, ana_sigs, c('purity', 'ploidy'), 'Sample', 'Sample')
  sin_result2_trans <- sin_result2 %>% bind_rows() %>% column_to_rownames('Variable') %>% t() %>% as.data.frame() %>% rownames_to_column('Variable') %>% as_tibble()
  
  sin_result_wrap <- sin_result1 %>% 
    bind_rows() %>% 
    left_join(sin_result2_trans, by = 'Variable') %>% 
    mutate(ICGC_abbr_top = i) %>% 
    select(ICGC_abbr_top, Signature = Variable, everything()) %>% 
    pivot_longer(-(ICGC_abbr_top:Signature), names_to = 'Clinical_Variable', values_to = 'p_val') %>% 
    group_by(Clinical_Variable) %>% 
    mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
    ungroup()
  
  mutsig_cli_relate_result[[i]] <- sin_result_wrap
  
  # sin_result_wrap %>% filter(Clinical_Variable %in% 'gender')
  
  print(i)
}


mutsig_cli_relate_signif <- mutsig_cli_relate_result %>% bind_rows() %>% filter(q_val < 0.05)
# 375

# write_tsv(mutsig_cli_relate_signif, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_cli_relate_signif.tsv')

# get mean weight of each test group for significantly clinical characteristics related mutational signatures

# gender, histological_type, 'stage', 'grade', wgd_status, 'first_therapy_response', 'TMB'


comb_mutsig_cli_dat <- comb_mutsig_we_filter %>% 
  select(ICGC_abbr_top:Signature, Weight) %>% 
  inner_join(select(comb_cli_data_t, Sample:first_therapy_response, wgd_status:TMB), by = 'Sample') %>% 
  pivot_longer(gender:TMB, names_to = 'Clinical_Variable', values_to = 'variable_value')

mutsig_cli_relate_signif_cat <- mutsig_cli_relate_signif %>% filter(Clinical_Variable %in% c('gender', 'histological_type', 'stage', 'grade', 'wgd_status', 'first_therapy_response', 'TMB'))

cli_signif_cat_stat <- comb_mutsig_cli_dat %>% 
  semi_join(mutsig_cli_relate_signif_cat, by = c('ICGC_abbr_top', 'Signature', 'Clinical_Variable')) %>% 
  filter(!is.na(variable_value)) %>% 
  group_by(ICGC_abbr_top, Signature, Clinical_Variable, variable_value) %>% 
  summarise(size = n(), size_non_zero = sum(Weight > 0), mean_Weight = mean(Weight), median_Weight = median(Weight)) %>% 
  ungroup() %>% 
  left_join(mutsig_cli_relate_signif_cat, by = c('ICGC_abbr_top', 'Signature', 'Clinical_Variable'))
# 931

write_tsv(filter(cli_signif_cat_stat, !Clinical_Variable %in% c('stage', 'grade')), '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/cli_signif_cat_stat.tsv')
write_tsv(filter(cli_signif_cat_stat, Clinical_Variable %in% c('stage', 'grade')), '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/cli_signif_stage_stat.tsv')


# cli_signif_cat_stat %>% filter(Clinical_Variable %in% 'first_therapy_response')


# purity, ploidy

comb_mutsig_pp_dat <- comb_mutsig_we_filter %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  mutate(group = ifelse(Weight > median(Weight), 'High', 'Low')) %>% 
  ungroup() %>% 
  select(ICGC_abbr_top:Signature, group) %>% 
  inner_join(select(comb_cli_data_t, Sample, purity, ploidy), by = 'Sample') %>% 
  filter_all(~!is.na(.)) %>% 
  group_by(ICGC_abbr_top, Signature, group) %>% 
  summarise(purity = median(purity), ploidy = median(ploidy)) %>% 
  ungroup() %>% 
  pivot_longer(purity:ploidy, names_to = 'Clinical_Variable', values_to = 'median_value')

mutsig_cli_relate_signif_con <- mutsig_cli_relate_signif %>% filter(Clinical_Variable %in% c('purity', 'ploidy')) # 69

cli_signif_pp_stat <- comb_mutsig_pp_dat %>% 
  inner_join(mutsig_cli_relate_signif_con, by = c('ICGC_abbr_top', 'Signature', 'Clinical_Variable')) %>% 
  arrange(ICGC_abbr_top, Signature, Clinical_Variable)

# write_tsv(cli_signif_pp_stat, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/cli_signif_pp_stat.tsv')


# Cancer Gene Mutation Status ----------------------------------------------------

comb_driver_gene_cnv_status <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene_cnv_status.tsv')

mutsig_dri_relate_result <- list()

for(i in unique(comb_mutsig_we_wide$ICGC_abbr_top)){
  
  # i <- 'Biliary-AdenoCA'
  
  sin_mutsig <- comb_mutsig_we_filter %>% 
    filter(ICGC_abbr_top %in% i) %>% 
    select(-Exposure, -Exposure_TMB)
  
  sin_driver <- comb_driver_gene %>% filter(ICGC_abbr_top %in% i) %>% pull(Hugo_Symbol) %>% intersect(colnames(driver_gs_status))
  
  sin_gs <- driver_gs_status[, c('Sample', sin_driver)] %>% 
    semi_join(sin_mutsig, by = 'Sample') %>% 
    column_to_rownames('Sample') %>% 
    select_if(~mean(.) >= 0.05 & sum(.) >= 5) # only focus on cancer genes of which mutation count >= 5 & mutation frequency >= 0.05
  
  if(nrow(sin_gs) < 30 | ncol(sin_gs) < 1){ # delete cancer type of which sample count < 30
    
    next
    
  }else{
    
    sin_gs_t <- sin_gs %>% 
      as.data.frame() %>% 
      rownames_to_column('Sample') %>% 
      as_tibble() %>% 
      pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status')
    
    sin_mutsig_gs_comb <- sin_mutsig %>% 
      inner_join(sin_gs_t, by = 'Sample') %>% 
      left_join(comb_driver_gene_cnv_status, by = c('Sample', 'Hugo_Symbol'))
    
    sin_res <- sin_mutsig_gs_comb %>% 
      group_by(ICGC_abbr_top, Signature, Hugo_Symbol) %>% 
      filter(mean(Weight > 0) >= 0.05 & sum(Weight > 0) >= 5) %>% # only retain signatures of which frequency over 0.05 and number more than 5
      nest() %>% 
      ungroup() %>% 
      mutate(p_val_wilcox = map_dbl(data, ~wilcox.test(Weight~Status, data = ., alternative = 'less')$p.value), 
             p_val_glm_uni = map_dbl(data, ~pnorm(summary(glm(Status ~ Weight, family = binomial, data = .))$coefficients[2, 3], lower.tail = FALSE)), 
             p_val_glm_multi = map_dbl(data, ~pnorm(summary(glm(Status ~ Weight + cnv_status, family = binomial, data = .))$coefficients[2, 3], lower.tail = FALSE)), 
             fc = map_dbl(data, ~group_by(., Status) %>% summarise(mean_weight = mean(Weight)) %>% mutate(fc = mean_weight[2]/mean_weight[1]) %>% pull(fc) %>% .[1])) %>% 
      select(-data) %>%
      group_by(Signature) %>%
      mutate(q_val_wilcox = p.adjust(p_val_wilcox, method = 'BH'), 
             q_val_glm_uni = p.adjust(p_val_glm_uni, method = 'BH'), 
             q_val_glm_multi = p.adjust(p_val_glm_multi, method = 'BH')) %>% 
      ungroup()
    
    mutsig_dri_relate_result[[i]] <- sin_res
  }
  
  print(i)
}

mutsig_dri_relate_res <- mutsig_dri_relate_result %>% bind_rows()
mutsig_dri_relate_signif <- mutsig_dri_relate_res %>% filter(q_val_glm_uni < 0.05, q_val_glm_multi < 0.05)


write_tsv(mutsig_dri_relate_res, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_res.tsv')
write_tsv(mutsig_dri_relate_signif, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_signif.tsv')

# mutsig_dri_relate_signif <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_signif.tsv')
# write_tsv(mutsig_dri_relate_res, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_res.tsv')


# onco pathway aberration --------------------------------------------

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/onco_pathway_sam_status.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/onco_pathway_gene_sam_status.Rdata')

mutsig_onco_pathway_gene_relate_result <- list()
mutsig_onco_pathway_relate_result <- list()

for(i in unique(comb_mutsig_we_wide$ICGC_abbr_top)){
  
  # i <- 'Biliary-AdenoCA'
  # i <- 'Liver-HCC'

  sin_mutsig <- comb_mutsig_we_filter %>% 
    filter(ICGC_abbr_top %in% i) %>% 
    select(-Exposure, -Exposure_TMB)
  
  sin_path_gene_status <- onco_pathway_gene_sam_status %>% 
    semi_join(sin_mutsig, by = 'Sample') %>% 
    column_to_rownames('Sample') %>% 
    select_if(~mean(.) >= 0.05 & sum(.) >= 5) %>% # only focus on pathways of which mutation count >= 5 & mutation frequency >= 0.05
    rownames_to_column('Sample') %>% 
    as_tibble()
  
  sin_path_status <- onco_pathway_sam_status %>% 
    semi_join(sin_mutsig, by = 'Sample') %>% 
    column_to_rownames('Sample') %>% 
    select_if(~mean(.) >= 0.05 & sum(.) >= 5) %>% # only focus on pathways of which mutation count >= 5 & mutation frequency >= 0.05
    rownames_to_column('Sample') %>% 
    as_tibble()
  
  
  if(nrow(sin_path_gene_status) < 30 | ncol(sin_path_gene_status) == 1){ # delete cancer type of which sample count < 30
    
    next
    
  }else{
    
    sin_path_gene_status_t <- sin_path_gene_status %>% pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status')
    sin_path_status_t <- sin_path_status %>% pivot_longer(-Sample, names_to = 'pathway', values_to = 'Status')
    
    sin_mutsig_path_gene_comb <- sin_mutsig %>% inner_join(sin_path_gene_status_t, by = 'Sample')
    sin_mutsig_path_comb <- sin_mutsig %>% inner_join(sin_path_status_t, by = 'Sample')
    
    sin_path_gene_res <- sin_mutsig_path_gene_comb %>% 
      group_by(ICGC_abbr_top, Signature, Hugo_Symbol) %>% 
      nest() %>% 
      ungroup() %>% 
      mutate(p_val = map_dbl(data, ~wilcox.test(Weight~Status, data = ., alternative = 'less')$p.value), 
             fc = map_dbl(data, ~group_by(., Status) %>% summarise(mean_weight = mean(Weight)) %>% mutate(fc = mean_weight[2]/mean_weight[1]) %>% pull(fc) %>% .[1])) %>% 
      select(-data) %>%
      group_by(Signature) %>% 
      mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
      ungroup()
    
    sin_path_res <- sin_mutsig_path_comb %>% 
      group_by(ICGC_abbr_top, Signature, pathway) %>% 
      nest() %>% 
      ungroup() %>% 
      mutate(p_val = map_dbl(data, ~wilcox.test(Weight~Status, data = ., alternative = 'less')$p.value), 
             fc = map_dbl(data, ~group_by(., Status) %>% summarise(mean_weight = mean(Weight)) %>% mutate(fc = mean_weight[2]/mean_weight[1]) %>% pull(fc) %>% .[1])) %>% 
      select(-data) %>%
      group_by(Signature) %>% 
      mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
      ungroup()
    
    mutsig_onco_pathway_gene_relate_result[[i]] <- sin_path_gene_res
    mutsig_onco_pathway_relate_result[[i]] <- sin_path_res
    
  }
  
  print(i)
}

mutsig_onco_pathway_gene_relate <- mutsig_onco_pathway_gene_relate_result %>% bind_rows()
mutsig_onco_pathway_relate <- mutsig_onco_pathway_relate_result %>% bind_rows()

mutsig_onco_pathway_gene_relate_signif <- mutsig_onco_pathway_gene_relate %>% filter(q_val < 0.05)
mutsig_onco_pathway_relate_signif <- mutsig_onco_pathway_relate %>% filter(q_val < 0.05)

# write_tsv(mutsig_onco_pathway_gene_relate_signif, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_onco_pathway_gene_relate_signif.tsv')
# write_tsv(mutsig_onco_pathway_relate_signif, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_onco_pathway_relate_signif.tsv')
# write_tsv(mutsig_onco_pathway_gene_relate, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_onco_pathway_gene_relate.tsv')
# write_tsv(mutsig_onco_pathway_relate, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_onco_pathway_relate.tsv')


# DDR pathway aberration --------------------------------------------

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_pathway_sam_status.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_pathway_gene_sam_status.Rdata')

mutsig_DDR_pathway_gene_relate_result <- list()
mutsig_DDR_pathway_relate_result <- list()

for(i in unique(comb_mutsig_we_wide$ICGC_abbr_top)){
  
  # i <- 'Bladder-TCC'
  
  sin_mutsig <- comb_mutsig_we_filter %>% 
    filter(ICGC_abbr_top %in% i) %>% 
    select(-Exposure, -Exposure_TMB)
  
  sin_path_gene_status <- DDR_pathway_gene_sam_status %>% 
    semi_join(sin_mutsig, by = 'Sample') %>% 
    column_to_rownames('Sample') %>% 
    select_if(~mean(.) >= 0.05 & sum(.) >= 5) %>% # only focus on pathways of which mutation count >= 5 & mutation frequency >= 0.05
    rownames_to_column('Sample') %>% 
    as_tibble()
  
  sin_path_status <- DDR_pathway_sam_status %>% 
    semi_join(sin_mutsig, by = 'Sample') %>% 
    column_to_rownames('Sample') %>% 
    select_if(~mean(.) >= 0.05 & sum(.) >= 5) %>% # only focus on pathways of which mutation count >= 5 & mutation frequency >= 0.05
    rownames_to_column('Sample') %>% 
    as_tibble()
  
  
  if(nrow(sin_path_gene_status) < 30 | ncol(sin_path_gene_status) == 1){ # delete cancer type of which sample count < 30
    
    next
    
  }else{
    
    sin_path_gene_status_t <- sin_path_gene_status %>% pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status')
    sin_path_status_t <- sin_path_status %>% pivot_longer(-Sample, names_to = 'pathway', values_to = 'Status')
    
    sin_mutsig_path_gene_comb <- sin_mutsig %>% inner_join(sin_path_gene_status_t, by = 'Sample')
    sin_mutsig_path_comb <- sin_mutsig %>% inner_join(sin_path_status_t, by = 'Sample')
    
    sin_path_gene_res <- sin_mutsig_path_gene_comb %>% 
      group_by(ICGC_abbr_top, Signature, Hugo_Symbol) %>% 
      nest() %>% 
      ungroup() %>% 
      mutate(p_val = map_dbl(data, ~wilcox.test(Weight~Status, data = ., alternative = 'less')$p.value), 
             fc = map_dbl(data, ~group_by(., Status) %>% summarise(mean_weight = mean(Weight)) %>% mutate(fc = mean_weight[2]/mean_weight[1]) %>% pull(fc) %>% .[1])) %>% 
      select(-data) %>%
      group_by(Signature) %>% 
      mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
      ungroup()
    
    sin_path_res <- sin_mutsig_path_comb %>% 
      group_by(ICGC_abbr_top, Signature, pathway) %>% 
      nest() %>% 
      ungroup() %>% 
      mutate(p_val = map_dbl(data, ~wilcox.test(Weight~Status, data = ., alternative = 'less')$p.value), 
             fc = map_dbl(data, ~group_by(., Status) %>% summarise(mean_weight = mean(Weight)) %>% mutate(fc = mean_weight[2]/mean_weight[1]) %>% pull(fc) %>% .[1])) %>% 
      select(-data) %>%
      group_by(Signature) %>% 
      mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
      ungroup()
    
    mutsig_DDR_pathway_gene_relate_result[[i]] <- sin_path_gene_res
    mutsig_DDR_pathway_relate_result[[i]] <- sin_path_res
    
  }
  
  print(i)
}

mutsig_DDR_pathway_gene_relate <- mutsig_DDR_pathway_gene_relate_result %>% bind_rows()
mutsig_DDR_pathway_relate <- mutsig_DDR_pathway_relate_result %>% bind_rows()

mutsig_DDR_pathway_gene_relate_signif <- mutsig_DDR_pathway_gene_relate %>% filter(q_val < 0.05)
mutsig_DDR_pathway_relate_signif <- mutsig_DDR_pathway_relate %>% filter(q_val < 0.05)

# HRD status ---------------------------------------------------------------

# HR deficiency: HRD score >42 and/or tumor BRCA1/2 mutation

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/scoreHRD/comb_HRDscore_res.Rdata')

comb_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mut_data.tsv')

BRCA_ns_mut_sam <- comb_mut_data %>% 
  filter(Hugo_Symbol %in% c('BRCA1', 'BRCA2'), Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  distinct(Sample) %>% 
  mutate(BRCA_status = 'MUT')

HRDscore_mutsig_we <- comb_mutsig_we_filter %>% 
  inner_join(comb_HRDscore_res, by = 'Sample') %>% 
  left_join(BRCA_ns_mut_sam, by = 'Sample') %>% 
  # mutate(HRD_status = ifelse(`HRD score` >= 42, TRUE, FALSE)) %>% 
  mutate(HRD_status = ifelse(`HRD score` > 42 | BRCA_status %in% 'MUT', TRUE, FALSE))

HRDscore_test_res <- HRDscore_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~wilcox.test(Weight~HRD_status, data = ., alternative = 'less')$p.value, NA))) %>%
  select(-data) %>% 
  unnest(p_val) %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
  ungroup()

HRDscore_test_res_signif <- HRDscore_test_res %>% filter(q_val < 0.05)

# save(HRDscore_test_res_signif, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/HRDscore_test_res_signif.Rdata')

HRDscore_mutsig_we_signif <- HRDscore_mutsig_we %>% semi_join(HRDscore_test_res_signif, by = c('ICGC_abbr_top', 'Signature'))

# save(HRDscore_mutsig_we_signif, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/HRDscore_mutsig_we_signif.Rdata')


# Figure 1 plots ---------------------------------------------------------------

library(ggh4x)
library(ggrepel)

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')
load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/onco_pathway_gene.Rdata')
load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_gene.Rdata')
load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/driver_gs_status.Rdata')
load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/onco_pathway_gene_sam_status.Rdata')
load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_pathway_gene_sam_status.Rdata')

comb_mutsig_we <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

comb_mutsig_we_filter <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% 
  ungroup()

comb_cli_data_t <- comb_cli_data %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(TMB = ifelse(TMB > median(TMB, na.rm = TRUE), 'High_TMB', 'Low_TMB')) %>% 
  ungroup() %>% 
  semi_join(comb_mutsig_we_filter, by = 'Sample') %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(histological_type = fct_lump_min(histological_type, 15)) %>% # lump histological type whose factor level is less than 15 to other
  ungroup()
# 8,662

comb_mutsig_we_wide <- comb_mutsig_we %>% 
  select(-Exposure, -Exposure_TMB) %>% 
  pivot_wider(ICGC_abbr_top:Sample, names_from = Signature, values_from = Weight)



# for each signature, compute the number of cancer type where the signature presents in over 5% samples

mutsig_freq <- comb_mutsig_we_filter %>% distinct(Signature, ICGC_abbr_top) %>% count(Signature)

# for each cancer gene, compute the number of cancer type where the gene presents in over 5% samples

dri_freq <- driver_gs_status %>% 
  pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status') %>% 
  inner_join(distinct(comb_mutsig_we_filter, ICGC_abbr_top, Sample), by = 'Sample') %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(mean(Status) >= 0.05) %>% 
  ungroup() %>% 
  distinct(ICGC_abbr_top, Hugo_Symbol) %>% 
  count(Hugo_Symbol)

# for each onco pathway gene, compute the number of cancer type where the gene presents in over 5% samples

onco_path_gene_freq <- onco_pathway_gene_sam_status %>% 
  pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status') %>% 
  inner_join(distinct(comb_mutsig_we_filter, ICGC_abbr_top, Sample), by = 'Sample') %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(mean(Status) >= 0.05) %>% 
  ungroup() %>% 
  distinct(ICGC_abbr_top, Hugo_Symbol) %>% 
  count(Hugo_Symbol)


# for each DDR pathway gene, compute the number of cancer type where the gene presents in over 5% samples

DDR_path_gene_freq <- DDR_pathway_gene_sam_status %>% 
  pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'Status') %>% 
  inner_join(distinct(comb_mutsig_we_filter, ICGC_abbr_top, Sample), by = 'Sample') %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(mean(Status) >= 0.05) %>% 
  ungroup() %>% 
  distinct(ICGC_abbr_top, Hugo_Symbol) %>% 
  count(Hugo_Symbol)


my_legend_theme <- theme(legend.background = element_blank(), 
                         legend.title.align = 0.5,   
                         legend.key.size = unit(0.45, 'cm'), 
                         legend.key.height = unit(0.45, 'cm'), 
                         legend.key.width = unit(0.45, 'cm'), 
                         legend.title = element_text(size = 10), 
                         legend.text = element_text(size = 8))


# Figure 1A data

age_mutsig_lmRob_signif <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/age_mutsig_lmRob_signif.tsv')

# Figure 1D data

mutsig_dri_relate_res <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_res.tsv')
mutsig_dri_relate_signif <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_signif.tsv')

# Figure 1A and 1D cancer type color

age_dri_cancer_type <- union(age_mutsig_lmRob_signif$ICGC_abbr_top, setdiff(mutsig_dri_relate_signif$ICGC_abbr_top, c('ColoRect-AdenoCA', 'Stomach-AdenoCA', 'Uterus-AdenoCA')))
# 'ColoRect-AdenoCA', 'Stomach-AdenoCA', 'Uterus-AdenoCA', cancer driver for Figure S1


all_cancer_colors <- tribble(
  ~cancer_type, ~color, 
  'Adrenal-neoplasm', '#0C487A', 
  'Biliary-AdenoCA', '#BAA133', 
  'Bladder-TCC', '#F8D0D6', 
  'Bone-Osteosarc', '#EE9028', 
  'Bone-Other', '#0C8D46', 
  'Breast', '#D92F86', 
  'Cervix', '#ECB065', 
  'CNS-GBM', '#A74C93', 
  'CNS-LGG', '#DD2929', 
  'CNS-Medullo', '#374F99', 
  'CNS-Oligo', '#90AD1C', 
  'CNS-PiloAstro', '#CC98BF', 
  'ColoRect-AdenoCA', '#9DD6F0', 
  'Eso-AdenoCA', '#0E78AA', 
  'Eye-Melanoma', '#DC6B72', 
  'Head-SCC', '#93C9A3', 
  'Kidney-ChRCC', '#7A378B', 
  'Kidney-Papillary', '#B2212D', 
  'Kidney-RCC', '#EEAAAE', 
  'Liver-HCC', '#C7C9D7', 
  'Lung-AdenoCA', '#CEC0DC', 
  'Lung-SCC', '#997FB5', 
  'Lymph-BNHL', '#724C29', 
  'Lymph-CLL', '#502D80', 
  'Mesothelium-Mesothelioma', '#EDE33E', 
  'Myeloid-AML', '#1F8C43', 
  'Myeloid-MDS/MPN', '#0FA094', 
  'Ovary-AdenoCA', '#CF7A29', 
  'Panc-AdenoCA', '#6B789B', 
  'Panc-Endocrine', '#DCBE29',
  'Pheochromocytoma', '#CD6600',  
  'Prost-AdenoCA', '#7A1A1D', 
  'Sarcoma', '#CAA98D', 
  'Skin-Melanoma', '#B3CB47', 
  'Stomach-AdenoCA', '#2FA0D2', 
  'Testis-Ca', '#008B8B', 
  'Thy-AdenoCA', '#9370DB', 
  'Thymoma', '#E066FF', 
  'UCS', '#FF8C69', 
  'Uterus-AdenoCA', '#F7DEC5'
)

age_dri_cancer_type_color <- all_cancer_colors %>% filter(cancer_type %in% age_dri_cancer_type)
# 32 cancer types


# Figure 1A ---------------------------------------------------------------
# Scatterplot
# Age-related mutational signatures 


age_signature_sort <- age_mutsig_lmRob_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)
age_cancer_sort <- age_mutsig_lmRob_signif %>% distinct(ICGC_abbr_top) %>% pull(ICGC_abbr_top) %>% str_sort(numeric = TRUE)


pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure 1A/age_mutsig_point_plot.pdf', width = 6, height = 6)

age_mutsig_lmRob_signif %>% 
  filter(estimate > 0) %>%
  mutate(log_p_val = -log10(p.value), 
         Signature = factor(Signature, levels = age_signature_sort)) %>% 
  ggplot(aes(estimate, abs(statistic))) +
  geom_point(aes(col = Signature, size = log_p_val)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_text_repel(aes(label = ICGC_abbr_top), size = 2) +
  ggsci::scale_color_npg() +
  labs(x = 'Robust regression coefficient', y = 'Robust regression statistic', col = NULL) +
  theme_classic() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank())
  
dev.off()


# Figure 1B ---------------------------------------------------------------
# Heatmap
# Association between clinical characteristics and mutational signatures

cli_signif_cat_stat <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/cli_signif_cat_stat.tsv')
cli_signif_pp_stat <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/cli_signif_pp_stat.tsv')

## TMB

tmb_mutsig_bi_count <- cli_signif_cat_stat %>% 
  filter(Clinical_Variable %in% 'TMB') %>% 
  select(ICGC_abbr_top, Signature, variable_value, mean_Weight, median_Weight) %>% 
  pivot_wider(names_from = variable_value, values_from = c(mean_Weight, median_Weight)) %>% 
  mutate(type = ifelse(mean_Weight_High_TMB >= mean_Weight_Low_TMB & median_Weight_High_TMB >= median_Weight_Low_TMB, 'positive', 'negative')) %>% 
  count(Signature, type, name = 'TMB')

## wgd_status

wgd_mutsig_bi_count <- cli_signif_cat_stat %>% 
  filter(Clinical_Variable %in% 'wgd_status') %>% 
  select(ICGC_abbr_top, Signature, variable_value, mean_Weight, median_Weight) %>% 
  pivot_wider(names_from = variable_value, values_from = c(mean_Weight, median_Weight)) %>% 
  mutate(type = ifelse(mean_Weight_wgd >= mean_Weight_no_wgd & median_Weight_wgd >= median_Weight_no_wgd, 'positive', 'negative')) %>% 
  count(Signature, type, name = 'wgd_status')

## purity

purity_mutsig_bi_count <- cli_signif_pp_stat %>% 
  filter(Clinical_Variable %in% 'purity') %>% 
  select(-p_val, -q_val) %>% 
  pivot_wider(names_from = group, values_from = median_value) %>% 
  mutate(type = ifelse(High > Low, 'positive', 'negative')) %>% 
  count(Signature, type, name = 'purity')

## ploidy

ploidy_mutsig_bi_count <- cli_signif_pp_stat %>% 
  filter(Clinical_Variable %in% 'ploidy') %>% 
  select(-p_val, -q_val) %>% 
  pivot_wider(names_from = group, values_from = median_value) %>% 
  mutate(type = ifelse(High > Low, 'positive', 'negative')) %>% 
  count(Signature, type, name = 'ploidy')


## combine 4 variables

comb_bi_count <- list(tmb_mutsig_bi_count, wgd_mutsig_bi_count, purity_mutsig_bi_count, ploidy_mutsig_bi_count) %>% 
  reduce(full_join, by = c('Signature', 'type')) %>% 
  replace(is.na(.), 0)


### total data

mutsig_cli_relate_signif <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_cli_relate_signif.tsv')

mutsig_cli_relate_signif_rmsg <- mutsig_cli_relate_signif %>% filter(!Clinical_Variable %in% c('stage', 'grade'))
cli_heatmap_sig_sort <- mutsig_cli_relate_signif_rmsg %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

cli_heatmap_dat <- mutsig_cli_relate_signif_rmsg %>%
  count(Signature, Clinical_Variable) %>%
  mutate(Signature = factor(Signature, levels = cli_heatmap_sig_sort)) %>% 
  pivot_wider(names_from = Clinical_Variable, values_from = n, values_fill = 0) %>% 
  pivot_longer(-Signature, names_to = 'Clinical_Variable', values_to = 'n') %>% 
  mutate(label = ifelse(n > 0, n, ''))

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure 1B/mutsig_cli_relate_heatmap_total.pdf', width = 8, height = 1.6)

cli_heatmap_dat %>%
  ggplot(aes(Signature, Clinical_Variable, fill = n)) +
  geom_tile(col = 'grey60') +
  geom_text(aes(label = label), size = 3) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_gradientn(colours = c('white', '#5082AF')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

### positive association

pos_heatmap_dat <- comb_bi_count %>% 
  filter(type %in% 'positive') %>% 
  select(-type) %>% 
  pivot_longer(-Signature, names_to = 'Clinical_Variable', values_to = 'pos_n') %>% 
  right_join(select(cli_heatmap_dat, Signature, Clinical_Variable), by = c('Signature', 'Clinical_Variable')) %>% 
  mutate(new_n = ifelse(is.na(pos_n), 0, pos_n), 
         new_label = ifelse(new_n > 0, new_n, ''), 
         Signature = factor(Signature, levels = cli_heatmap_sig_sort)) %>% 
  select(Signature, Clinical_Variable, n = new_n, label = new_label)


pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure 1B/mutsig_cli_relate_heatmap_pos.pdf', width = 8, height = 1.6)

pos_heatmap_dat %>%
  ggplot(aes(Signature, Clinical_Variable, fill = n)) +
  geom_tile(col = 'grey60') +
  geom_text(aes(label = label), size = 3) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_gradientn(colours = c('white', '#C9372E')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


### negative association

neg_heatmap_dat <- comb_bi_count %>% 
  filter(type %in% 'negative') %>% 
  select(-type) %>% 
  pivot_longer(-Signature, names_to = 'Clinical_Variable', values_to = 'neg_n') %>% 
  right_join(select(cli_heatmap_dat, Signature, Clinical_Variable), by = c('Signature', 'Clinical_Variable')) %>% 
  mutate(new_n = ifelse(is.na(neg_n), 0, neg_n), 
         new_label = ifelse(new_n > 0, new_n, ''), 
         Signature = factor(Signature, levels = cli_heatmap_sig_sort)) %>% 
  select(Signature, Clinical_Variable, n = new_n, label = new_label)

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure 1B/mutsig_cli_relate_heatmap_neg.pdf', width = 8, height = 1.6)

neg_heatmap_dat %>%
  ggplot(aes(Signature, Clinical_Variable, fill = n)) +
  geom_tile(col = 'grey60') +
  geom_text(aes(label = label), size = 3) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_gradientn(colours = c('white', '#8151A0')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


# Figure 1C ---------------------------------------------------------------
# Boxplot
# Examples for clinical features-signatures association

comb_cli_mutsig_data <- comb_cli_data_t %>% inner_join(comb_mutsig_we_wide, by = c('ICGC_abbr_top', 'Sample'))
# save(comb_cli_mutsig_data, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_cli_mutsig_data.Rdata')


forBoxplot <- function(comb_data, Cancer_type, which_signature, group_var, is_var_gene_or_pathway = FALSE){
  
  # comb_data <- comb_cli_mutsig_data
  # Cancer_type <- 'Panc-AdenoCA'
  # which_signature <- 'SBS18'
  # group_var <- 'first_therapy_response'
  
  sin_data <- comb_data %>% filter(ICGC_abbr_top %in% Cancer_type) %>% 
    select(all_of(c(which_signature, group_var))) %>% 
    filter_all(~!is.na(.))
  
  colnames(sin_data) <- c('Signature', 'Variable')
  
  if(is_var_gene_or_pathway){
#    sin_data <- sin_data %>% mutate(Variable = ifelse(Variable == TRUE, str_c(group_var, ' mutant'), str_c(group_var, ' wildtype')))
    sin_data <- sin_data %>% mutate(Variable = ifelse(Variable == TRUE, 'MUT', 'WT'), 
                                    Variable = factor(Variable, levels = c('WT', 'MUT')))
  }
  
  sin_plot <- sin_data %>% 
    ggplot(aes(Variable, Signature, fill = Variable)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = 'white', outlier.color = NA) +
    labs(x = NULL, y = 'Signature contribution', title = str_c(which_signature, Cancer_type, sep = ', ')) +
    scale_fill_manual(values = c('#5873A8', '#C98035')) +
    theme_bw() +
    theme(legend.position = 'none', 
          axis.ticks.x = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  if(is_var_gene_or_pathway){
    sin_plot + annotate('text', x = 1.5, y = -Inf, label = group_var, vjust = -0.5, hjust = 0.5)
  }else{
    sin_plot
  }
  
}


Stomach_SBS1_gender <- forBoxplot(comb_cli_mutsig_data, 'Stomach-AdenoCA', 'SBS1', 'gender') + annotate('text', x = 1.5, y = Inf, label = 'p=8.64e-3', vjust = 1.5)
Medullo_SBS1_histype <- forBoxplot(comb_cli_mutsig_data, 'CNS-Medullo', 'SBS1', 'histological_type') + 
  annotate('text', x = 2, y = Inf, label = 'p=4.42e-4', vjust = 1.5) + 
  scale_fill_manual(values = c('#00A1D5', '#5873A8', '#C98035')) + 
  scale_x_discrete(labels = c('Desmoplastic', 'Large cell', 'NOS'))


# purity and ploidy

Prost_SBS18_purity_dat <- comb_cli_mutsig_data %>% 
  filter(ICGC_abbr_top %in% 'Prost-AdenoCA') %>% 
  mutate(SBS18 = ifelse(SBS18 > median(SBS18), 'High', 'Low'), 
         SBS18 = factor(SBS18, levels = c('Low', 'High')))

Ovary_SBS40_ploidy_dat <- comb_cli_mutsig_data %>%
  filter(ICGC_abbr_top %in% 'Ovary-AdenoCA') %>%
  mutate(SBS40 = ifelse(SBS40 > median(SBS40), 'High', 'Low'), 
         SBS40 = factor(SBS40, levels = c('Low', 'High')))

Prost_SBS18_purity <- forBoxplot(Prost_SBS18_purity_dat, 'Prost-AdenoCA', 'purity', 'SBS18') + annotate('text', x = 1.5, y = Inf, label = 'p=3.22e-4', vjust = 1.5) + labs(y = 'Tumor purity', title = 'SBS18, Prost-AdenoCA')
Uterus_SBS40_ploidy <- forBoxplot(Ovary_SBS40_ploidy_dat, 'Ovary-AdenoCA', 'ploidy', 'SBS40') + annotate('text', x = 1.5, y = Inf, label = 'p=5.80e-3', vjust = 1.5) + labs(y = 'Tumor ploidy', title = 'SBS40, Ovary-AdenoCA')

# combine boxplots

cli_sig_example_plots <- list(Stomach_SBS1_gender, Prost_SBS18_purity, Medullo_SBS1_histype, Uterus_SBS40_ploidy)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_clinical_features_boxplot.pdf', width = 6, height = 6)

wrap_plots(cli_sig_example_plots) + plot_layout(nrow = 2)

dev.off()


# Figure 1D --------------------------------------------------------------- 

# dot plots of each cancer type


mutsig_dri_dot_dat <- mutsig_dri_relate_res %>% 
  semi_join(mutsig_dri_relate_signif, by = 'ICGC_abbr_top') %>% 
  mutate(dot_size = case_when(
    q_val_glm_multi > 0.05 ~ 0.5, 
    q_val_glm_multi < 0.05 & q_val_glm_multi > 0.01 ~ 1, 
    q_val_glm_multi < 0.01 & q_val_glm_multi > 0.001 ~ 1.25, 
    q_val_glm_multi < 0.001 & q_val_glm_multi > 0.0001 ~ 1.5, 
    q_val_glm_multi < 0.0001 ~ 1.75
  ), fc = case_when(
    is.na(fc) ~ 0, 
    fc %in% Inf ~ max(fc[is.finite(fc)]), 
    TRUE ~ fc
  ), 
  log2fc = ifelse(fc < 1, 0, log2(fc)))

dri_dot_sig_sort <- mutsig_dri_dot_dat %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)


mutsig_dri_dot_dat_re <- mutsig_dri_dot_dat %>% 
  group_by(Signature) %>% 
  filter(any(q_val_glm_multi < 0.05)) %>% # filter signatures without significant results
  ungroup() %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(any(q_val_glm_multi < 0.05)) %>% # filter genes without significant results
  ungroup()


mutsig_dri_dot_dat_rep1 <- mutsig_dri_dot_dat_re %>% 
  filter(!ICGC_abbr_top %in% c('ColoRect-AdenoCA', 'Stomach-AdenoCA', 'Uterus-AdenoCA')) %>% 
  mutate(Signature = factor(Signature, levels = dri_dot_sig_sort), 
         ttype_gene = str_c(ICGC_abbr_top, Hugo_Symbol, sep = '_')) %>% 
  left_join(age_dri_cancer_type_color, by = c('ICGC_abbr_top' = 'cancer_type'))


ttype_gene_idx <- mutsig_dri_dot_dat_rep1 %>% 
  distinct(ttype_gene, Hugo_Symbol, color) %>% 
  arrange(ttype_gene)


pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_dri_relate_analysis/mutsig_dri_dot_part1.pdf', width = 15, height = 4.5)

mutsig_dri_dot_dat_rep1 %>% 
  ggplot(aes(ttype_gene, Signature, fill = log2fc, size = dot_size)) +
  geom_point(shape = 21, color = 'grey60') +
  labs(x = NULL, y = NULL, size = NULL) +
  scale_fill_gradientn(colors = c('white', '#C9372E')) +
  scale_x_discrete(breaks = ttype_gene_idx$ttype_gene, labels = ttype_gene_idx$Hugo_Symbol) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = ttype_gene_idx$color), 
        axis.ticks = element_line(color = 'grey60'), 
        legend.position = 'top')

dev.off()

# Figure 1E ---------------------------------------------------------------
# Boxplot
# Examples for cancer gene-signatures association

comb_mutsig_dri_data <- comb_mutsig_we_wide %>% inner_join(driver_gs_status, by = 'Sample')

Bladder_SBS5_ERCC2 <- forBoxplot(comb_mutsig_dri_data, 'Bladder-TCC', 'SBS5', 'ERCC2', is_var_gene_or_pathway = TRUE) + annotate('text', x = 1.5, y = Inf, label = 'p=1.25e-3', vjust = 1.5)
Melanoma_SBS7a_COL5A1 <- forBoxplot(comb_mutsig_dri_data, 'Skin-Melanoma', 'SBS7a', 'COL5A1', is_var_gene_or_pathway = TRUE) + annotate('text', x = 1.5, y = Inf, label = 'p=3.87e-10', vjust = 1.5)

cli_sig_dri_example_plots <- list(Bladder_SBS5_ERCC2, Melanoma_SBS7a_COL5A1)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_cancer_gene_boxplot.pdf', width = 3.3, height = 6)

wrap_plots(cli_sig_dri_example_plots) + plot_layout(nrow = 2)

dev.off()


# Figure 1F ---------------------------------------------------------------
# Heatmap
# Association between Oncogenic and DDR pathways and mutational signatures

# Oncogenic pathways

mutsig_onco_pathway_relate_signif <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_onco_pathway_relate_signif.tsv')

mutsig_onco_sig_sort <- mutsig_onco_pathway_relate_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

mutsig_onco_pathway_plot_dat <- mutsig_onco_pathway_relate_signif %>% 
  count(pathway, Signature) %>% 
  mutate(Signature = factor(Signature, levels = mutsig_onco_sig_sort)) %>% 
  pivot_wider(names_from = pathway, values_from = n, values_fill = 0) %>% 
  pivot_longer(-Signature, names_to = 'pathway', values_to = 'n') %>% 
  mutate(label = ifelse(n > 0, n, ''))

onco_pathway_heatmap <- mutsig_onco_pathway_plot_dat %>% 
  ggplot(aes(Signature, pathway, fill = n)) +
  geom_tile(col = 'grey60') +
  geom_text(aes(label = label), size = 3) +
  labs(x = NULL, y = NULL, title = 'Curated Oncogenic Pathways', fill = NULL) +
  scale_fill_gradientn(colours = c('white', '#1F8C43')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))


# DDR pathways

mutsig_DDR_pathway_relate_signif <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_DDR_pathway_relate_signif.tsv')

mutsig_DDR_sig_sort <- mutsig_DDR_pathway_relate_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

mutsig_DDR_pathway_plot_dat <- mutsig_DDR_pathway_relate_signif %>% 
  count(pathway, Signature) %>% 
  mutate(Signature = factor(Signature, levels = mutsig_DDR_sig_sort)) %>% 
  pivot_wider(names_from = pathway, values_from = n, values_fill = 0) %>% 
  pivot_longer(-Signature, names_to = 'pathway', values_to = 'n') %>% 
  mutate(label = ifelse(n > 0, n, ''))

DDR_pathway_heatmap <- mutsig_DDR_pathway_plot_dat %>% 
  ggplot(aes(Signature, pathway, fill = n)) +
  geom_tile(col = 'grey60') +
  geom_text(aes(label = label), size = 3) +
  labs(x = NULL, y = NULL, title = 'DNA Damage Repair Pathways', fill = NULL) +
  scale_fill_gradientn(colours = c('white', '#8151A0')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1))


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_pathway_relate_heatmap.pdf', width = 16, height = 3)

onco_pathway_heatmap + DDR_pathway_heatmap

dev.off()


# Figure S1 plots ---------------------------------------------------------------

library(patchwork)

# dot plots of each cancer type

mutsig_dri_relate_res <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_res.tsv')
mutsig_dri_relate_signif <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_relate_signif.tsv')

mutsig_dri_dot_dat <- mutsig_dri_relate_res %>% 
  semi_join(mutsig_dri_relate_signif, by = 'ICGC_abbr_top') %>% 
  mutate(dot_size = case_when(
    q_val_glm_multi > 0.05 ~ 0.5, 
    q_val_glm_multi < 0.05 & q_val_glm_multi > 0.01 ~ 1, 
    q_val_glm_multi < 0.01 & q_val_glm_multi > 0.001 ~ 1.25, 
    q_val_glm_multi < 0.001 & q_val_glm_multi > 0.0001 ~ 1.5, 
    q_val_glm_multi < 0.0001 ~ 1.75
  ), fc = case_when(
    is.na(fc) ~ 0, 
    fc %in% Inf ~ max(fc[is.finite(fc)]), 
    TRUE ~ fc
  ), 
  log2fc = ifelse(fc < 1, 0, log2(fc)))

dri_dot_sig_sort <- mutsig_dri_dot_dat %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)


mutsig_dri_dot_dat_re <- mutsig_dri_dot_dat %>% 
  group_by(Signature) %>% 
  filter(any(q_val_glm_multi < 0.05)) %>% # filter signatures without significant results
  ungroup() %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(any(q_val_glm_multi < 0.05)) %>% # filter genes without significant results
  ungroup()


mutsig_dri_dot_dat_rep1 <- mutsig_dri_dot_dat_re %>% filter(!ICGC_abbr_top %in% c('ColoRect-AdenoCA', 'Stomach-AdenoCA', 'Uterus-AdenoCA'))


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_dot_part1.pdf', width = 15, height = 6)

dot_plot_part1 <- mutsig_dri_dot_dat_rep1 %>% 
  mutate(Signature = factor(Signature, levels = dri_dot_sig_sort)) %>% 
  ggplot(aes(Hugo_Symbol, Signature, fill = log2fc, size = dot_size)) +
  geom_point(shape = 21, color = 'grey60') +
  labs(x = NULL, y = NULL, size = NULL) +
  scale_fill_gradientn(colors = c('white', '#C9372E')) +
  facet_grid2(. ~ ICGC_abbr_top, scales = 'free', space = 'free', strip = strip_vanilla(clip = 'off')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks = element_line(color = 'grey60'), 
        legend.position = 'bottom', 
        panel.border = element_rect(size = 0.5, color = 'grey60'),
        strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3),
        strip.background = element_blank(), 
        panel.spacing = unit(0, 'cm'))

dot_plot_part1

dev.off()


dot_plot_part2 <- list()

for(i in c('ColoRect-AdenoCA', 'Stomach-AdenoCA', 'Uterus-AdenoCA')){
  
  # i <- 'Uterus-AdenoCA'
  # 
  # pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_dri_dot_uterus.pdf', width = 12, height = 4)
  
  sin_dot_plot_dat <- mutsig_dri_dot_dat_re %>% 
    mutate(Signature = factor(Signature, levels = dri_dot_sig_sort)) %>% 
    filter(ICGC_abbr_top %in% i) %>% 
    group_by(Signature) %>% 
    filter(any(q_val_glm_multi < 0.05)) %>% # filter signatures without significant results
    ungroup() %>% 
    group_by(Hugo_Symbol) %>% 
    filter(any(q_val_glm_multi < 0.05)) %>% # filter genes without significant results
    ungroup()
  
  sin_dot_plot <- sin_dot_plot_dat %>% 
    ggplot(aes(Hugo_Symbol, Signature, fill = log2fc, size = dot_size)) +
    geom_point(shape = 21, color = 'grey60') +
    labs(x = NULL, y = NULL, size = NULL, title = i) +
    scale_fill_gradientn(colors = c('white', '#C9372E')) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5)) +
    my_legend_theme
  
  dot_plot_part2[[i]] <- sin_dot_plot
  # sin_dot_plot
  
  # dev.off()
  
}

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_dri_relate_analysis/mutsig_dri_dot_part2.pdf', width = 14, height = 4)

dot_plot_part2[[3]] + (dot_plot_part2[[1]] / dot_plot_part2[[2]]) + plot_layout(widths = c(3, 1))

dev.off()


# heatmap (pathway gene)

# DDR

DDR_pathway_color <- tribble(
  ~pathway, ~color,
  'BER', '#A5CCE1',
  'DR', '#2177B2',
  'DS', '#B2DE8B',
  'FA', '#379E2B',
  'HR', '#FA9897',
  'MMR', '#E21A19',
  'NER', '#FDBD6F',
  'NHEJ', '#FD7F00', 
  'TLS', '#CAB0D5'
)

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_gene.Rdata')

mutsig_DDR_pathway_gene_relate_signif <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_DDR_pathway_gene_relate_signif.tsv')
mutsig_DDR_pathway_gene_relate <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_DDR_pathway_gene_relate.tsv')

mutsig_DDR_gene_dot_dat <- mutsig_DDR_pathway_gene_relate %>% 
  semi_join(mutsig_dri_relate_signif, by = 'ICGC_abbr_top') %>% 
  mutate(dot_size = case_when(
    q_val > 0.05 ~ 0.5, 
    q_val < 0.05 & q_val > 0.01 ~ 1, 
    q_val < 0.01 & q_val > 0.001 ~ 1.25, 
    q_val < 0.001 & q_val > 0.0001 ~ 1.5, 
    q_val < 0.0001 ~ 1.75
  ), fc = case_when(
    is.na(fc) ~ 0, 
    fc %in% Inf ~ max(fc[is.finite(fc)]), 
    TRUE ~ fc
  ), 
  log2fc = ifelse(fc < 1, 0, log2(fc))) %>% 
  left_join(select(DDR_gene, -type), by = 'Hugo_Symbol')


DDR_dot_sig_sort <- mutsig_DDR_gene_dot_dat %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)


mutsig_DDR_gene_dot_dat_re <- mutsig_DDR_gene_dot_dat %>% 
  group_by(Signature) %>% 
  filter(any(q_val < 0.05)) %>% # filter signatures without significant results
  ungroup() %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  filter(any(q_val < 0.05)) %>% # filter genes without significant results
  ungroup()


mutsig_DDR_tile_dat <- mutsig_DDR_gene_dot_dat_re %>% mutate(ttype_gene_id = str_c(ICGC_abbr_top, Hugo_Symbol, sep = ':')) %>% distinct(ICGC_abbr_top, ttype_gene_id, pathway)


mutsig_DDR_dot_plot <- mutsig_DDR_gene_dot_dat_re %>% 
  mutate(Signature = factor(Signature, levels = DDR_dot_sig_sort)) %>% 
  ggplot(aes(Hugo_Symbol, Signature, fill = log2fc, size = dot_size)) +
  geom_point(shape = 21, color = 'grey60') +
  labs(x = NULL, y = NULL, size = NULL) +
  scale_fill_gradientn(colors = c('white', '#C9372E')) +
  facet_grid2(. ~ ICGC_abbr_top, scales = 'free', space = 'free', strip = strip_vanilla(clip = 'off')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.border = element_rect(size = 0.5),
        # strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
        strip.background = element_blank(), 
        panel.spacing = unit(0, 'cm'))

mutsig_DDR_anno_plot <- mutsig_DDR_tile_dat %>% 
  ggplot(aes(ttype_gene_id, 1, fill = pathway)) +
  geom_bar(stat = 'identity') +
  labs(x = NULL, y = NULL, fill = NULL) +
  # scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(breaks = DDR_pathway_color$pathway, values = DDR_pathway_color$color) +
  ggthemes::theme_few() +
  facet_grid( ~ ICGC_abbr_top, scales = 'free', space = 'free') +
  guides(fill = guide_legend(nrow = 2)) +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing = unit(0, 'cm'), 
        strip.background = element_blank(), 
        panel.border = element_rect(size = 0.5), 
        strip.text = element_blank())


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_DDR_path_gene_dotplot.pdf', width = 8, height = 5)

mutsig_DDR_dot_plot / mutsig_DDR_anno_plot + plot_layout(heights = c(20, 1), guides = 'collect') & theme(legend.position = 'bottom')

dev.off()


# Cancer gene interaction in UCEC(Uterus-AdenoCA), ColoRect-AdenoCA, Stomach-AdenoCA -------------------------

library(maftools)

comb_mut_data <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_mut_data.tsv')
tcga_mut_data <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv')
comb_driver_gene <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')
comb_mutsig_we <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')

comb_mut_data_ttype <- comb_mut_data %>% inner_join(select(comb_cli_data, Sample, ICGC_abbr_top), by = 'Sample')

comb_mutsig_we_filter <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% 
  ungroup()

UCEC_dri_gene <- comb_driver_gene %>% filter(ICGC_abbr_top %in% 'Uterus-AdenoCA') %>% pull(Hugo_Symbol)
ColoRect_dri_gene <- comb_driver_gene %>% filter(ICGC_abbr_top %in% 'ColoRect-AdenoCA') %>% pull(Hugo_Symbol)
Stomach_dri_gene <- comb_driver_gene %>% filter(ICGC_abbr_top %in% 'Stomach-AdenoCA') %>% pull(Hugo_Symbol)

UCEC_mut_data <- tcga_mut_data %>%
  filter(ttype %in% 'UCEC') %>% # only Uterus-AdenoCA from tcga has mutation data
  semi_join(comb_mutsig_we_filter, by = c('patient_id' = 'Sample')) %>%
  select(Tumor_Sample_Barcode = patient_id, Hugo_Symbol:Tumor_Seq_Allele2)

ColoRect_mut_data <- comb_mut_data_ttype %>%
  filter(ICGC_abbr_top %in% 'ColoRect-AdenoCA') %>% 
  semi_join(comb_mutsig_we_filter, by = c('Sample')) %>%
  select(Tumor_Sample_Barcode = Sample, Hugo_Symbol:Tumor_Seq_Allele2)

Stomach_mut_data <- comb_mut_data_ttype %>%
  filter(ICGC_abbr_top %in% 'Stomach-AdenoCA') %>% 
  semi_join(comb_mutsig_we_filter, by = c('Sample')) %>%
  select(Tumor_Sample_Barcode = Sample, Hugo_Symbol:Tumor_Seq_Allele2)


# vroom::vroom_write(UCEC_mut_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/UCEC_mut_data.maf')
# vroom::vroom_write(ColoRect_mut_data, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/ColoRect_mut_data.maf')
# vroom::vroom_write(Stomach_mut_data, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/Stomach_mut_data.maf')

# UCEC

UCEC_maf <- read.maf(maf = '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/UCEC_mut_data.maf')

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/UCEC_dri_interaction.pdf')

UCEC_dri_interation_res <- somaticInteractions(maf = UCEC_maf, top = 55, genes = UCEC_dri_gene, pvalue = c(0.05, 0.1))

dev.off()

# ColoRect-AdenoCA

ColoRect_maf <- read.maf(maf = '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/ColoRect_mut_data.maf')

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_dri_relate_analysis/ColoRect_dri_interaction.pdf')

ColoRect_dri_interation_res <- somaticInteractions(maf = ColoRect_maf, top = 55, genes = ColoRect_dri_gene, pvalue = c(0.05, 0.1))

dev.off()

# Stomach-AdenoCA

Stomach_maf <- read.maf(maf = '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/Stomach_mut_data.maf')

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_dri_relate_analysis/Stomach_dri_interaction.pdf')

Stomach_dri_interation_res <- somaticInteractions(maf = Stomach_maf, top = 55, genes = Stomach_dri_gene, pvalue = c(0.05, 0.1))

dev.off()

# HRD status ---------------------------------------------------------------

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/HRDscore_test_res_signif.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/HRDscore_mutsig_we_signif.Rdata')

HRD_p_anno_dat <- HRDscore_test_res_signif %>% 
  mutate(start = 'HR proficiency', end = 'HR deficiency', y = 1.1, 
         label = case_when(
           p_val < 0.0001 ~ '****', 
           p_val >= 0.0001 & p_val < 0.001 ~ '***', 
           p_val >= 0.001 & p_val < 0.01 ~ '**', 
           p_val >= 0.01 & p_val < 0.05 ~ '*', 
           p_val >= 0.05 ~ 'ns'
         ))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/HRDscore_mutsig_boxplot.pdf', width = 16, height = 4)

HRDscore_mutsig_we_signif %>% 
  mutate(HRD_status = ifelse(HRD_status == TRUE, 'HR deficiency', 'HR proficiency'), 
         HRD_status = factor(HRD_status, levels = c('HR proficiency', 'HR deficiency')), 
         Signature = factor(Signature, levels = str_sort(unique(HRD_p_anno_dat$Signature), numeric = TRUE))) %>% 
  arrange(Signature) %>% 
  ggplot(aes(HRD_status, Weight)) +
  geom_point(aes(col = HRD_status), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1), shape = 16, alpha = 7/10, size = 2) + 
  geom_boxplot(aes(col = HRD_status), fill = NA, position = position_dodge(width = 1), outlier.color = NA) +
  labs(x = NULL, y = 'Signature contribution', col = NULL) +
  scale_color_manual(breaks = c('HR proficiency', 'HR deficiency'), values = c('#5082AF', '#C9372E')) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggsignif::geom_signif(data = HRD_p_anno_dat, aes(xmin = start, xmax = end, annotations = label, y_position = y), manual = TRUE) +
  facet_nested(~ ICGC_abbr_top + Signature) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = 'top', 
        strip.text = element_text(size = 10))

dev.off()

# SBS18 reactive oxygen species ------------------------------------

# KEGG

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/SBS18_ros_score_group_data.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/SBS18_ros_score_test.Rdata')

SBS18_p_anno_dat <- SBS18_ros_score_test %>% 
  mutate(start = 'Low', end = 'High', y = 3.75, 
         label = case_when(
           p_val < 0.0001 ~ '****', 
           p_val >= 0.0001 & p_val < 0.001 ~ '***', 
           p_val >= 0.001 & p_val < 0.01 ~ '**', 
           p_val >= 0.01 & p_val < 0.05 ~ '*', 
           p_val >= 0.05 ~ 'ns'
         ))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/SBS18_ROS_ssgsea_boxplot.pdf', width = 16, height = 4)

SBS18_ros_score_group_data %>% 
  ggplot(aes(group, `Reactive Oxygen Species`)) +
  geom_point(aes(col = group), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1), shape = 16, alpha = 7/10, size = 2) + 
  geom_boxplot(fill = NA, position = position_dodge(width = 1), outlier.color = NA, show.legend = FALSE) +
  labs(x = NULL, y = 'ROS ssGSEA score', col = 'SBS18') +
  scale_color_manual(values = c("#C98035", "#5873A8")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggsignif::geom_signif(data = SBS18_p_anno_dat, aes(xmin = start, xmax = end, annotations = label, y_position = y), manual = TRUE) +
  facet_grid(.~ICGC_abbr_top) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))

dev.off()

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/SBS18_ros_score_group_data_alt.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/SBS18_ros_score_test_alt.Rdata')

SBS18_p_anno_dat_alt <- SBS18_ros_score_test_alt %>% 
  mutate(start = 'Low', end = 'High', y = 3.75, 
         label = case_when(
           p_val < 0.0001 ~ '****', 
           p_val >= 0.0001 & p_val < 0.001 ~ '***', 
           p_val >= 0.001 & p_val < 0.01 ~ '**', 
           p_val >= 0.01 & p_val < 0.05 ~ '*', 
           p_val >= 0.05 ~ 'ns'
         ))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/SBS18_ROS_ssgsea_boxplot_alt.pdf', width = 16, height = 4)

SBS18_ros_score_group_data_alt %>% 
  ggplot(aes(group, `Reactive Oxygen Species`)) +
  geom_point(aes(col = group), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1), shape = 16, alpha = 7/10, size = 2) + 
  geom_boxplot(fill = NA, position = position_dodge(width = 1), outlier.color = NA, show.legend = FALSE) +
  labs(x = NULL, y = 'ROS ssGSEA score', col = 'SBS18') +
  scale_color_manual(values = c("#C98035", "#5873A8")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggsignif::geom_signif(data = SBS18_p_anno_dat_alt, aes(xmin = start, xmax = end, annotations = label, y_position = y), manual = TRUE) +
  facet_grid(.~ICGC_abbr_top) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))

dev.off()

