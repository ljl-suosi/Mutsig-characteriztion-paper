
library(tidyverse)
library(MutationalPatterns)
library(BSgenome)

comb_mutsig_mut_dat <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_mutsig_mut_dat.tsv')
# mutations in comb_mutsig_mut_dat is larger than comb_mut_data (because of clinical data limitation)

comb_snv_sim <- comb_mutsig_mut_dat %>% 
  select(mut_id) %>% 
  separate(mut_id, into = c(NA, 'chrom', 'start', 'ref', NA), sep = ':', remove = FALSE) %>% 
  mutate(start = as.numeric(start))

comb_snv_sim_gr <- GRanges(
  seqnames = Rle(str_c('chr', comb_snv_sim$chrom)),
  ranges = IRanges(comb_snv_sim$start, end = comb_snv_sim$start, names = comb_snv_sim$mut_id),
  strand = Rle(strand(rep('*', nrow(comb_snv_sim)))), # Attention!!! must be *, otherwise fail to map some mutations to genomic intervals
  REF = DNAStringSet(comb_snv_sim$ref))
# 23,422,116 snvs

genecode_v19_protein_coding_gene <- read_tsv('/boot3/bio_liaojl/Common_data/processed_data/genecode_v19_protein_coding_gene.tsv')

genecode_v19_gr <- GRanges(
  seqnames = Rle(genecode_v19_protein_coding_gene$chrom),
  ranges = IRanges(genecode_v19_protein_coding_gene$start, end = genecode_v19_protein_coding_gene$end, names = genecode_v19_protein_coding_gene$Gene_ID),
  strand = Rle(genecode_v19_protein_coding_gene$strand))

snv_strand <- mut_strand(comb_snv_sim_gr, genecode_v19_gr)
# the strand of mutations in comb_snv_sim_gr must be *, otherwise insufficient result be returned
# e.g., 'ACAP3, SP101724:1:1230448:G:A, +' doesn't match with 'chr1 1227756-1244989      -' due to unmatched strand
# the strand info will be added in mut_strand function, ref C or T will be assigned '+', otherwise '-', in a view of C>X or T>X

comb_snv_strand <- comb_snv_sim %>% 
  mutate(strand = as.character(snv_strand)) %>% 
  select(mut_id, strand)

vroom::vroom_write(comb_snv_strand, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_snv_strand.tsv')
