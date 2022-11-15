library(tidyverse)

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')
tcga_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv')
# 3,202,399
pcawg_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_mut_data.tsv')
# 23,159,591

# samples with info about project source

tcga_sam <- tcga_mut_data %>% distinct(patient_id) %>% mutate(Project = 'TCGA') %>% rename(Sample = patient_id)
pcawg_sam <- pcawg_mut_data %>% distinct(icgc_specimen_id) %>% mutate(Project = 'PCAWG') %>% rename(Sample = icgc_specimen_id)
mut_sam_project <- tcga_sam %>% bind_rows(pcawg_sam)

save(mut_sam_project, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mut_sam_project.Rdata')

# primary TCGA and PCAWG mutation data

comb_mut_data <- tcga_mut_data %>% 
  select(Hugo_Symbol:Tumor_Seq_Allele2, Sample = patient_id) %>% 
  bind_rows(select(pcawg_mut_data, Hugo_Symbol, Chromosome, Start_Position = Start_position, End_Position = End_position, Variant_Classification:Tumor_Seq_Allele2, Sample = icgc_specimen_id)) %>% 
  semi_join(comb_cli_data, by = 'Sample')

vroom::vroom_write(comb_mut_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mut_data.tsv')


