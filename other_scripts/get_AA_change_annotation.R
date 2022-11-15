
library(tidyverse)

pcawg_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_mut_data.tsv') # 23,159,591
comb_mutsig_clo_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_alt_data.tsv')

SP116478_MGAT4C_mutsig <- comb_mutsig_clo_data %>% 
  filter(Sample %in% 'SP116478', Hugo_Symbol %in% 'MGAT4C', !is.na(Attri_Signature)) %>% 
  select(mut_id, Signature = Attri_Signature) %>% 
  separate(mut_id, into = c('Sample', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'))
  
SP116478_MGAT4C_mut_dat <- pcawg_mut_data %>% 
  filter(icgc_specimen_id %in% 'SP116478', Hugo_Symbol %in% 'MGAT4C') %>% 
  select(Chromosome, Start_Position = Start_position, End_Position = End_position, Reference_Allele:Tumor_Seq_Allele2)

# write_tsv(SP116478_MGAT4C_mut_dat, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/SP116478_MGAT4C_mut_dat.tsv', col_names = FALSE)


MGAT4C_HGSV_anno <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/output/SP116478_MGAT4C_anno.exonic_variant_function.tsv', col_names = c('id', 'type', 'HGVS_anno', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')) %>% 
  select(-(id:type)) %>% 
  separate_rows(HGVS_anno, sep = ',') %>% 
  filter(HGVS_anno != '') %>% 
  separate(HGVS_anno, into = c('Hugo_Symbol', 'Refseq_ID', NA, NA, 'HGVS_short'), sep = ':') %>% 
  mutate(Chromosome = as.character(Chromosome), Start_Position = as.character(Start_Position)) %>% 
  left_join(SP116478_MGAT4C_mutsig, by = c('Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'))

save(MGAT4C_HGSV_anno, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/MGAT4C_HGSV_anno.Rdata')


