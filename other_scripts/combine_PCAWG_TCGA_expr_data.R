
library(tidyverse)
library(sva)

load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/Transcriptomic_Data/pcawg_gene_fpkm.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/Transcriptomic_Data/GRCh37_protein_coding_gene.Rdata')
tcga_mutsig_ori <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/TCGA_WES_sigProfiler_SBS_signatures_in_samples.tsv')
tcga_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv')

# TCGA samples with both mutation and mutsig data

tcga_sam_ttype <- tcga_mut_data %>% distinct(ttype, Tumor_Sample_Barcode) # 10,103 samples/patients

tcga_mutsig_sam <- tcga_mutsig_ori %>% 
  filter(Accuracy >= 0.8) %>% 
  semi_join(tcga_sam_ttype, by = c('Sample Names' = 'Tumor_Sample_Barcode')) %>% 
  mutate(`Sample Names` = str_sub(`Sample Names`, 1, 19)) %>% 
  pull(`Sample Names`)
# 6074 samples

# TCGA tpm data

tcga_cancer_type <- tcga_mut_data %>% distinct(ttype) %>% pull(ttype)
tcga_tpm_list <- list()

for(i in tcga_cancer_type){
  
  tcga_tpm_list[[i]] <- get(load(str_c('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/TCGA_expression/', i, '_tpm_expr.Rdata')))
  
  # print(dim(tcga_tpm_list[[i]]))
  
}

tcga_tpm_comb_all <- tcga_tpm_list %>% 
  reduce(left_join, by = 'Hugo_Symbol') %>% 
  as_tibble() %>% 
  rename_all(~str_sub(., 1, 19))

tcga_inter_sam <- intersect(tcga_mutsig_sam, colnames(tcga_tpm_comb_all))
# 5235 samples

tcga_tpm_comb <- tcga_tpm_comb_all %>% 
  select(Hugo_Symbol, all_of(tcga_inter_sam)) %>% 
  rename_all(~str_sub(., 9, 12)) %>% 
  rename(Hugo_Symbol = bol)

# funtion: FPKM value to TPM

fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

# PCAWG tpm data

pcawg_tpm_expr <- pcawg_gene_fpkm %>% 
  mutate_if(is.numeric, ~fpkmToTpm(.)) %>% 
  left_join(GRCh37_protein_coding_gene, by = 'Gene_ID') %>% 
  select(Hugo_Symbol = Gene_Symbol, starts_with('SP')) %>% 
  distinct(Hugo_Symbol, .keep_all = TRUE) # select first for duplicate gene symbols

pcawg_sample_donor <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_Driver/Downloaded_Data/PCAWG/pcawg_sample_sheet.tsv') %>% 
  distinct(icgc_specimen_id, icgc_donor_id)

# convert sample to donor

pcawg_tpm_exprd <- pcawg_tpm_expr %>% 
  column_to_rownames('Hugo_Symbol') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  as_tibble() %>% 
  inner_join(pcawg_sample_donor, by = c('Sample' = 'icgc_specimen_id')) %>% 
  select(-Sample) %>% 
  distinct(icgc_donor_id, .keep_all = TRUE) %>% 
  column_to_rownames('icgc_donor_id') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Hugo_Symbol') %>% 
  as_tibble()


# combine PCAWG and TCGA expr data

comb_tpm_expr_data <- pcawg_tpm_exprd %>% inner_join(tcga_tpm_comb, by = 'Hugo_Symbol')
# 17,512 genes, 6520 patients


# Remove batch effects between TCGA and PCAWG -----------------------------


# patients with mutsig info

comb_mutsig_we <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

mutsig_expr_group <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>%
  mutate(group = ifelse(Weight > median(Weight), 'High', 'Low'), group = factor(group, levels = c('Low', 'High'))) %>% 
  ungroup() %>% 
  left_join(pcawg_sample_donor, by = c('Sample' = 'icgc_specimen_id')) %>% 
  mutate(Patient = ifelse(is.na(icgc_donor_id), Sample, icgc_donor_id)) %>% 
  filter(Patient %in% colnames(comb_tpm_expr_data)) %>% # intersection between mutsig patients and expr patients
  mutate(Project = ifelse(is.na(icgc_donor_id), 'TCGA', 'PCAWG')) %>% 
  select(ICGC_abbr_top, Patient, Signature, group, Project) %>% 
  distinct(ICGC_abbr_top, Patient, Signature, .keep_all = TRUE)

# save(mutsig_expr_group, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mutsig_expr_group.Rdata')

mutsig_cancer_type_patient <- mutsig_expr_group %>% 
  distinct(ICGC_abbr_top, Patient, Project)
# 6453 patients, 1218 PCAWG patients, 5235 TCGA patients


# cancer types only have one project

cancer_type1pro <- mutsig_cancer_type_patient %>% distinct(ICGC_abbr_top, Project) %>% count(ICGC_abbr_top) %>% filter(n < 2)
# 15 cancer types
cancer_type1pro_patients <- mutsig_cancer_type_patient %>% semi_join(cancer_type1pro, by = 'ICGC_abbr_top') %>% pull(Patient)
# 652 patients

comb_tpm_expr_type1 <- comb_tpm_expr_data %>% select(Hugo_Symbol, all_of(cancer_type1pro_patients))


# cancer types have two projects

cancer_type2pro <- mutsig_cancer_type_patient %>% distinct(ICGC_abbr_top, Project) %>% count(ICGC_abbr_top) %>% filter(n == 2) %>% pull(ICGC_abbr_top)
# 21 cancer types

# combat

combat_expr_list <- list()

for(i in cancer_type2pro){
  
  print(str_c('Processing ', i))
  
  sin_can_data <- mutsig_cancer_type_patient %>% filter(ICGC_abbr_top %in% i)
  
  sin_can_expr <- comb_tpm_expr_data %>% select(Hugo_Symbol, all_of(sin_can_data$Patient)) %>% column_to_rownames('Hugo_Symbol')
  
  combat_expr_list[[i]] <- ComBat(dat = sin_can_expr, batch = sin_can_data$Project, prior.plots = FALSE, par.prior = TRUE) %>% as.data.frame() %>% rownames_to_column('Hugo_Symbol') %>% as_tibble()
  
}

# ComBat_seq(count_matrix, batch = batch)

comb_tpm_expr_type2 <- combat_expr_list %>% reduce(left_join, by = 'Hugo_Symbol')
comb_tpm_expr_final <- comb_tpm_expr_type1 %>% left_join(comb_tpm_expr_type2, by = 'Hugo_Symbol')

save(comb_tpm_expr_final, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_tpm_expr_final.Rdata')

