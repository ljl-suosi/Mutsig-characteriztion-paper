library(tidyverse)

Cancer_type_conv_tab <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/ICGC_TCGA_Cancer_Type_Convert.csv')

# PanCanAtlas

tcga_mutsig_ori <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/TCGA_WES_sigProfiler_SBS_signatures_in_samples.tsv')
tcga_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv')
tcga_cli_data <- read_tsv('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_cli.tsv', guess_max = 5000)

sam_ttype <- tcga_mut_data %>% distinct(ttype, Tumor_Sample_Barcode) # 10,103 samples/patients
tcga_sam_mutcount <- tcga_mut_data %>% count(patient_id, name = 'Mut_count')

tcga_mutsig_wt <- tcga_mutsig_ori %>% 
  inner_join(sam_ttype, by = c('Sample Names' = 'Tumor_Sample_Barcode')) %>% 
  mutate(`Sample Names` = str_sub(`Sample Names`, 9, 12)) %>% 
  left_join(Cancer_type_conv_tab, by = c('ttype' = 'TCGA_abbr')) %>% 
  select(ICGC_abbr_top, Sample = `Sample Names`, Accuracy:SBS60) %>% 
  pivot_longer(SBS1:SBS60, names_to = 'Signature', values_to = 'Exposure') %>% 
  group_by(Sample) %>% 
  mutate(Weight = Exposure/sum(Exposure)) %>% 
  ungroup()

tcga_mutsig_we <- tcga_mutsig_wt %>% 
  left_join(tcga_sam_mutcount, by = c('Sample' = 'patient_id')) %>% 
  mutate(Exposure = Weight * Mut_count, Exposure_TMB = Exposure/38) %>% 
  select(-Mut_count)

vroom::vroom_write(tcga_mutsig_we, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_mutsig_we.tsv')

# tcga_mutsig_we <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_mutsig_we.tsv', guess_max = 10000)

# PCAWG

pcawg_mutsig_ori <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_Driver/Downloaded_Data/PCAWG/pcawg_sigProfiler_SBS_signatures_in_samples.tsv') # 2,780
load('/pub6/Temp/Liaojl/Data/SBS18_pancancer_research/Data/Processed/pcawg_wgs_cli.Rdata') # pcawg_wgs_cli
pcawg_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_mut_data.tsv')
# only 1,950 samples, but there are 2,780 samples in pcawg_mutsig_ori

########## 2022.3.7
########## Attention!!! 8 samples' (such as SP99289) cancer type are consistent between pcawg_mutsig_ori and pcawg_wgs_cli. Our study is based on pcawg clinical data and modify the cancer type information of specific samples.
##########

pcawg_sam_mutcount <- pcawg_mut_data %>% count(icgc_specimen_id, name = 'Mut_count')

pcawg_sam_ttype <- pcawg_wgs_cli %>% distinct(histology_abbreviation, icgc_specimen_id)
ICGC_tab <- Cancer_type_conv_tab %>% distinct(ICGC_abbr_top, ICGC_abbr) %>% mutate(Ctype_lowercase = str_to_lower(ICGC_abbr))

pcawg_mutsig_we <- pcawg_mutsig_ori %>% 
  inner_join(pcawg_sam_ttype, by = c('Sample_Names' = 'icgc_specimen_id')) %>% 
  select(-Cancer_Types) %>% 
  select(Cancer_Types = histology_abbreviation, Sample_Names:SBS60) %>% 
  mutate(Ctype_lowercase = str_to_lower(Cancer_Types)) %>% 
  left_join(ICGC_tab, by = 'Ctype_lowercase') %>% 
  select(ICGC_abbr_top, Sample = Sample_Names, Accuracy:SBS60) %>% 
  pivot_longer(SBS1:SBS60, names_to = 'Signature', values_to = 'Exposure') %>% 
  group_by(Sample) %>% 
  mutate(Weight = Exposure/sum(Exposure)) %>% 
  ungroup() %>% 
  mutate(Exposure_TMB = Exposure/2800)

vroom::vroom_write(pcawg_mutsig_we, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/pcawg_mutsig_we.tsv')

# combine PanCanAtlas and PCAWG mutsig data

comb_mutsig_we_all <- pcawg_mutsig_we %>% bind_rows(tcga_mutsig_we)

vroom::vroom_write(comb_mutsig_we_all, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we_all.tsv')

comb_mutsig_we <- pcawg_mutsig_we %>% 
  bind_rows(tcga_mutsig_we) %>% 
  group_by(Signature) %>% 
  mutate(max_sigexp = max(Exposure)) %>% 
  ungroup() %>% 
  filter(Accuracy >= 0.8, max_sigexp >= 1) %>% 
  select(-Accuracy, -max_sigexp)

vroom::vroom_write(comb_mutsig_we, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

# final samples for downstream analsyis

final_focus_sample <- comb_mutsig_we %>% distinct(ICGC_abbr_top, Sample)
# 42 cancer types, 8836 samples, each cancer contain 3 samples to 790 samples

save(final_focus_sample, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/final_focus_sample.Rdata')

# prevalent signatures of each cancer type

cancer_mutsig <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Exposure >= 1) >= 0.05) %>% 
  ungroup() %>% 
  distinct(ICGC_abbr_top, Signature)

write_tsv(cancer_mutsig, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/cancer_mutsig.tsv')


