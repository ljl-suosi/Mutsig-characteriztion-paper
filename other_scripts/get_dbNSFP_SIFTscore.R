
library(tidyverse)

# mutation preparation for ANNOVAR

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')
tcga_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv') # 3,202,399
pcawg_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_mut_data.tsv') # 23,159,591


# primary TCGA and PCAWG mutation data
# only focus on non-silent mutations

comb_mut_datav <- tcga_mut_data %>% 
  select(Chromosome:End_Position, Reference_Allele:Tumor_Seq_Allele2, Variant_Classification, Sample = patient_id) %>% 
  bind_rows(select(pcawg_mut_data, Chromosome, Start_Position = Start_position, End_Position = End_position, Reference_Allele:Tumor_Seq_Allele2, Variant_Classification, Sample = icgc_specimen_id)) %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  select(-Variant_Classification)
# 2,129,611

comb_mut_datav_id <- comb_mut_datav %>% 
  mutate(mut_id = str_c(Sample, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ':')) %>% 
  select(mut_id)

missense_mut_id <- tcga_mut_data %>% 
  select(Chromosome:End_Position, Reference_Allele:Tumor_Seq_Allele2, Variant_Classification, Sample = patient_id) %>% 
  bind_rows(select(pcawg_mut_data, Chromosome, Start_Position = Start_position, End_Position = End_position, Reference_Allele:Tumor_Seq_Allele2, Variant_Classification, Sample = icgc_specimen_id)) %>% 
  filter(Variant_Classification %in% 'Missense_Mutation') %>% 
  mutate(mut_id = str_c(Sample, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ':')) %>% 
  select(mut_id)

# vroom::vroom_write(select(comb_mut_datav, -Sample), '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/comb_mut_datav.tsv', col_names = FALSE)

# run ANNOVAR

# perl table_annovar.pl /pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/comb_mut_datav.tsv humandb/ -build hg19 -out /pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/comb_mut_anno -remove -protocol dbnsfp30a -operation f -nastring .

comb_mut_anno <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/comb_mut_anno.hg19_multianno.txt')

missense_mut_anno <- comb_mut_datav_id %>% 
  cbind(select(comb_mut_anno, SIFT_pred, Polyphen2_HDIV_pred, LRT_pred, MutationTaster_pred, MutationAssessor_pred)) %>% 
  as_tibble() %>% 
  semi_join(missense_mut_id, by = 'mut_id')

# deleterious missense mutations predicted by various methods

SIFT_deleterious_mut <- missense_mut_anno %>% filter(SIFT_pred %in% 'D') %>% pull(mut_id)
PolyPhen2_deleterious_mut <- missense_mut_anno %>% filter(Polyphen2_HDIV_pred %in% c('D', 'P')) %>% pull(mut_id)
LRT_deleterious_mut <- missense_mut_anno %>% filter(LRT_pred %in% 'D') %>% pull(mut_id)
MutationTaster_deleterious_mut <- missense_mut_anno %>% filter(MutationTaster_pred %in% c('A', 'D')) %>% pull(mut_id)
MutationAssessor_deleterious_mut <- missense_mut_anno %>% filter(MutationAssessor_pred %in% c('H', 'M')) %>% pull(mut_id)

# consensus deleterious mutations (remove mutations predicted to be 'benign' by 3 or more of the 5 prediction algorithms)

missense_deleterious_mut <- missense_mut_anno %>% 
  rowwise() %>% 
  mutate(na_sum = sum(c(SIFT_pred, Polyphen2_HDIV_pred, LRT_pred, MutationTaster_pred, MutationAssessor_pred) %in% '.')) %>% 
  ungroup() %>% 
  filter(na_sum < 3) %>% 
  mutate(SIFT_pred = ifelse(SIFT_pred %in% 'T', 1, 0), 
         Polyphen2_HDIV_pred = ifelse(Polyphen2_HDIV_pred %in% 'B', 1, 0), 
         LRT_pred = ifelse(LRT_pred %in% 'N', 1, 0), 
         MutationTaster_pred = ifelse(MutationTaster_pred %in% c('N', 'P'), 1, 0), 
         MutationAssessor_pred = ifelse(MutationAssessor_pred %in% c('L', 'N'), 1, 0), 
         sum_benign = SIFT_pred + Polyphen2_HDIV_pred + LRT_pred + MutationTaster_pred + MutationAssessor_pred) %>% 
  filter(sum_benign < 3)
# 1,106,127

vroom::vroom_write(missense_deleterious_mut, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/missense_deleterious_mut.tsv')

