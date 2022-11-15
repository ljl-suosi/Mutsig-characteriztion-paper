
library(tidyverse)

exonic_vep_conseq_dat <- vroom::vroom('/boot3/bio_liaojl/Common_data/consequence_ranking.tsv', col_names = c('CHROMOSOME', 'POSITION', 'REF', 'ALT', 'ENSEMBL_NOVER', 'CONSEQUENCE'))
comb_driver_gene <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')
genecode_v19_protein_coding_gene <- read_tsv('/boot3/bio_liaojl/Common_data/processed_data/genecode_v19_protein_coding_gene.tsv')
ensembl_variation_severity <- read_tsv('/boot3/bio_liaojl/Common_data/ensembl_variation_severity.tsv')

ensembl_var_syno <- ensembl_variation_severity %>% 
  mutate(conseq_type = ifelse(RANK < 15, 'non_synonymous', 'synonymous')) %>% 
  select(CONSEQUENCE, conseq_type)

genecode_driver_gene <- genecode_v19_protein_coding_gene %>% 
  semi_join(comb_driver_gene, by = 'Hugo_Symbol') %>% 
  separate(Gene_ID, into = c('Gene_ID_rmv', NA), sep = '\\.', remove = FALSE)


driver_gene_conseq <- exonic_vep_conseq_dat %>% 
  inner_join(select(genecode_driver_gene, Gene_ID_rmv, Hugo_Symbol), by = c('ENSEMBL_NOVER' = 'Gene_ID_rmv')) %>% 
  inner_join(ensembl_var_syno, by = c('CONSEQUENCE'))
# 281 ensembl IDs, 280 symbols

retained_gene <- driver_gene_conseq %>% distinct(Hugo_Symbol, .keep_all = TRUE) %>% select(Hugo_Symbol, ENSEMBL_NOVER)
# 280 symbols

driver_gene_nonsyno_trinuc <- driver_gene_conseq %>% 
  filter(conseq_type %in% 'non_synonymous') %>% 
  semi_join(retained_gene, by = 'ENSEMBL_NOVER') %>% # 280 symbols
  select(chr = CHROMOSOME, pos = POSITION, ref = REF, alt = ALT, Hugo_Symbol) %>% 
  mutate(chr = str_c('chr', chr), 
         context = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chr, pos - 1, pos + 1, as.character = T), 
         tricontext = str_c(str_sub(context, 1, 1), '[', ref, '>', alt, ']', str_sub(context, 3, 3)))
# 2,096,781 substitutions


vroom::vroom_write(driver_gene_nonsyno_trinuc, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/driver_gene_nonsyno_trinuc.tsv')



