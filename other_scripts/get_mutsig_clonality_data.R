
# TCGA single mutation signature attribution ---------------------------------------------------

source('/pub6/Temp/Liaojl/Code/common_R_functions/assignTriContextSig.R')
source('/pub6/Temp/Liaojl/Code/common_R_functions/getMutTriContext.R')
load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_mutsig_w96_data.Rdata')
tcga_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv')
tcga_mutsig_we <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_mutsig_we.tsv')
# Attention!!! tcga_mutsig_we has't filtered samples whose signature accuracy were less than 0.8

tcga_tri_sig <- assignTriContextSig(tcga_mutsig_we, pcawg_mutsig_w96_data)

# 1. add mut trinucleotides info

tcga_sub_context <- getMutTriContext(tcga_mut_data, 'patient_id', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')

# 2. add signature assignment info (note: produce trinucleotide signature assignment first)

tcga_sub_data <- tcga_sub_context %>% 
  left_join(select(tcga_tri_sig, -prob), by = c('sample.id' = 'Sample', 'tricontext' = 'tri_nuc'))
# 3,039,840

# vroom::vroom_write(tcga_sub_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_sub_data.tsv')

# tcga_sub_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_sub_data.tsv')

# infer TCGA clonality  ---------------------------------------------------------------

tcga_cli_data <- read_tsv('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_cli.tsv', guess_max = 5000)

tcga_abs_cn_ori <- vroom::vroom('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Original/TCGA_mastercalls.abs_segtabs.fixed.txt')
tcga_abs_pp_ori <- read_tsv('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Original/TCGA_mastercalls.abs_tables_JSedit.fixed.txt')

tcga_abs_cn <- tcga_abs_cn_ori %>% select(Sample:End, minor_cn = Modal_HSCN_1, major_cn = Modal_HSCN_2, total_cn = Modal_Total_CN)
tcga_abs_pur <- tcga_abs_pp_ori %>% select(Sample = array, purity)

tcga_cnp <- tcga_abs_cn %>% 
  inner_join(tcga_abs_pur, by = 'Sample') %>% 
  filter_all(~!is.na(.)) %>% 
  mutate(Chromosome = ifelse(Chromosome == 23, 'X', Chromosome))
# no Y chromosome, only 1:23 chromosomes

tcga_mut_cnp <- tcga_mut_data %>% 
  filter(Variant_Type %in% 'SNP') %>% # 3,039,840 substitutions
  mutate(Sample = str_sub(Tumor_Sample_Barcode, 1, 15)) %>% 
  inner_join(tcga_cnp, by = c('Sample', 'Chromosome')) %>% 
  filter(Start_Position >= Start & Start_Position <= End) %>% 
  select(-End_Position, -(Start:End)) # 2,742,142 substitutions


############################### method from science translational medicine, based on ccf

tcga_somatics_m <- tcga_mut_cnp %>% 
  filter(total_cn > 0, Chromosome != 'X') %>% # 86 total_cn > 0 mutations and remove X chromosome mutations
  mutate(sampleId = Sample, 
         totalReadCount = t_alt_count + t_ref_count, 
         normal_cn = 2, 
         adjustedVaf = t_alt_count / (purity * totalReadCount), 
         adjustedCopyNumber = normal_cn * (1 - purity) + total_cn * purity, 
         somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber)) %>%  # mutation copy number equal to VAF/p * (Nnormal*(1-p) + Ntumor*p)
  rename(chromosome = Chromosome, alleleReadCount = t_alt_count)

absolute_ccf_estimation <- function(n.alt, depth, purity, local.copy.number){
  
  # n.alt <- unlist(mut.table$var_counts)
  # depth <- depth.t
  # purity <- cellularity
  # local.copy.number <- abs.cn
  
  f.function <- function (c,purity,local.copy.number){
    
    return((purity*c) / (2*(1-purity) + purity*local.copy.number))
  }
  
  x <- dbinom(n.alt, depth, prob = sapply(seq(0.01, 1, length.out = 100), f.function, purity, local.copy.number))
  if(min(x)==0){
    x[length(x)] <- 1
  }
  
  names(x) <- seq(0.01,1,length.out=100)
  
  # ccf distribution
  
  xnorm   <- x/sum(x)
  xsort   <- sort(xnorm, decreasing = TRUE)
  xcumLik <- cumsum(xsort)
  n = sum(xcumLik < 0.95) + 1
  LikThresh <- xsort[n]
  cint  <- x[xnorm >= LikThresh]
  cellu <- as.numeric(names(cint))
  l.t   <- cellu[1]
  r.t   <- cellu[length(cellu)]
  m     <- cellu[which.max(cint)]

  return(m)
  
}

n <- 10
nr <- nrow(tcga_somatics_m)
tcga_somatics_m_parts <- split(tcga_somatics_m, rep(1:ceiling(n), each = nr/n, length.out = nr))

save(tcga_somatics_m_parts, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/tcga_somatics_m_parts.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/tcga_somatics_m_parts.Rdata')

tcga_somatics_m_p <- tcga_somatics_m_parts[[10]]
all_muts_ccf_list10 <- list()

for(i in 1:nrow(tcga_somatics_m_p)){
  
  print(i)
  print(c(tcga_somatics_m_p$alleleReadCount[i], tcga_somatics_m_p$totalReadCount[i], tcga_somatics_m_p$purity[i], tcga_somatics_m_p$total_cn[i]))
  all_muts_ccf_list10[[i]] <- absolute_ccf_estimation(tcga_somatics_m_p$alleleReadCount[i], tcga_somatics_m_p$totalReadCount[i], tcga_somatics_m_p$purity[i], tcga_somatics_m_p$total_cn[i])
  
}

all_muts_ccf10 <- all_muts_ccf_list10 %>% bind_rows()
save(all_muts_ccf10, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf10.Rdata')


load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf1.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf2.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf3.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf4.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf5.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf6.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf7.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf8.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf9.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_muts_ccf10.Rdata')


tcga_stm_ccf_data10 <- tcga_somatics_m_parts[[10]] %>% 
  select(patient_id:Sample) %>% 
  cbind(all_muts_ccf10) %>% 
  as_tibble() %>% 
  mutate(ccf_ci95_clonality = ifelse(absolute_ccf_0.95 >= 1, 'Clonal', 'Subclonal'), 
         ccf_prob_clonality = ifelse(prob_clonal > 0.5, 'Clonal', 'Subclonal'))

tcga_stm_ccf_data <- list(tcga_stm_ccf_data1, tcga_stm_ccf_data2, tcga_stm_ccf_data3, tcga_stm_ccf_data4, tcga_stm_ccf_data5, tcga_stm_ccf_data6, tcga_stm_ccf_data7, tcga_stm_ccf_data8, tcga_stm_ccf_data9, tcga_stm_ccf_data10) %>% bind_rows()

# vroom::vroom_write(tcga_stm_ccf_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_stm_ccf_data.tsv')

# combin TCGA mutation attribution and clonality

tcga_sub_sim <- tcga_sub_data %>% 
  mutate(chr = str_remove(chr, 'chr'), mut_id = str_c(sample.id, chr, pos, ref, alt, sep = ':')) %>% 
  select(Hugo_Symbol, Chromosome = chr, Start_position = pos, Variant_Classification, Sample = sample.id, Ori_Cancer_Type = ttype, tricontext, Attri_Signature, mut_id)

tcga_stm_ccf_data_sim <- tcga_stm_ccf_data %>% 
  mutate(mut_id = str_c(patient_id, chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ':')) %>% 
  select(mut_id, clonality = ccf_ci95_clonality)
# 2,741,798

tcga_mutsig_clo_alt_data <- tcga_sub_sim %>% 
  left_join(tcga_stm_ccf_data_sim, by = 'mut_id') %>% 
  mutate(Chromosome = str_c('chr', Chromosome))
# 3,039,840

# vroom::vroom_write(tcga_mutsig_clo_alt_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_mutsig_clo_alt_data.tsv')
# tcga_mutsig_clo_alt_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/tcga_mutsig_clo_alt_data.tsv')


# PCAWG substitution data(SNV context, attributed signature and clonality) --------------------------------------------------------------

pcawg_sub_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_sub_data.tsv')
# 21,606,842

pcawg_sub_data_sim <- pcawg_sub_data %>% 
  mutate(mut_id = str_c(icgc_specimen_id, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2, sep = ':')) %>% 
  select(Hugo_Symbol:Start_position, Variant_Classification, icgc_specimen_id, histology_abbreviation:Attri_Signature, mut_id)

############################### method from science translational medicine, based on ccf

load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_sub_cnpp_data.Rdata') # 21,490,922

pcawg_tcga_somatics <- pcawg_sub_cnpp_data %>% 
  filter(total_cn > 0, !Chromosome %in% c('X', 'Y'), !is.na(t_alt_count), !is.na(t_ref_count)) %>% 
  mutate(sampleId = icgc_specimen_id, 
         totalReadCount = t_alt_count + t_ref_count, 
         normal_cn = 2, 
         adjustedVaf = t_alt_count / (purity * totalReadCount), 
         adjustedCopyNumber = normal_cn * (1 - purity) + total_cn * purity, 
         somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber)) %>%
  rename(chromosome = Chromosome, alleleReadCount = t_alt_count)
# 20,373,807

n_p <- 10
nr_p <- nrow(pcawg_tcga_somatics)
pcawg_tcga_somatics_parts <- split(pcawg_tcga_somatics, rep(1:ceiling(n_p), each = nr_p/n_p, length.out = nr_p))

for(i in 1:n_p){
  
  vroom::vroom_write(pcawg_tcga_somatics_parts[[i]], str_c('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_', i, '.tsv'))
  
}


pcawg_tcga_somatics_p <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_2.tsv')
pcawg_all_muts_ccf_list <- list()

for(i in 1:nrow(pcawg_tcga_somatics_p)){
  
  print(i)
  print(c(pcawg_tcga_somatics_p$alleleReadCount[i], pcawg_tcga_somatics_p$totalReadCount[i], pcawg_tcga_somatics_p$purity[i], pcawg_tcga_somatics_p$total_cn[i]))
  pcawg_all_muts_ccf_list[[i]] <- absolute_ccf_estimation(pcawg_tcga_somatics_p$alleleReadCount[i], pcawg_tcga_somatics_p$totalReadCount[i], pcawg_tcga_somatics_p$purity[i], pcawg_tcga_somatics_p$total_cn[i])
  
}

pcawg_all_muts_ccf2 <- pcawg_all_muts_ccf_list %>% bind_rows()
vroom::vroom_write(pcawg_all_muts_ccf2, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf2.tsv')


pcawg_tcga_somatics_p1 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_1.tsv')
pcawg_tcga_somatics_p2 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_2.tsv')
pcawg_tcga_somatics_p3 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_3.tsv')
pcawg_tcga_somatics_p4 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_4.tsv')
pcawg_tcga_somatics_p5 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_5.tsv')
pcawg_tcga_somatics_p6 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_6.tsv')
pcawg_tcga_somatics_p7 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_7.tsv')
pcawg_tcga_somatics_p8 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_8.tsv')
pcawg_tcga_somatics_p9 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_9.tsv')
pcawg_tcga_somatics_p10 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_tcga_somatics_part_10.tsv')


pcawg_all_muts_ccf1 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf1.tsv')
pcawg_all_muts_ccf2 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf2.tsv')
pcawg_all_muts_ccf3 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf3.tsv')
pcawg_all_muts_ccf4 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf4.tsv')
pcawg_all_muts_ccf5 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf5.tsv')
pcawg_all_muts_ccf6 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf6.tsv')
pcawg_all_muts_ccf7 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf7.tsv')
pcawg_all_muts_ccf8 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf8.tsv')
pcawg_all_muts_ccf9 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf9.tsv')
pcawg_all_muts_ccf10 <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/pcawg_all_muts_ccf10.tsv')


pcawg_stm_ccf_data10 <- pcawg_tcga_somatics_p10 %>% 
  select(Hugo_Symbol:histology_abbreviation) %>% 
  cbind(pcawg_all_muts_ccf10) %>% 
  as_tibble() %>% 
  mutate(ccf_ci95_clonality = ifelse(absolute_ccf_0.95 >= 1, 'Clonal', 'Subclonal'), 
         ccf_prob_clonality = ifelse(prob_clonal > 0.5, 'Clonal', 'Subclonal'))

pcawg_stm_ccf_data <- list(pcawg_stm_ccf_data1, pcawg_stm_ccf_data2, pcawg_stm_ccf_data3, pcawg_stm_ccf_data4, pcawg_stm_ccf_data5, pcawg_stm_ccf_data6, pcawg_stm_ccf_data7, pcawg_stm_ccf_data8, pcawg_stm_ccf_data9, pcawg_stm_ccf_data10) %>% bind_rows()

# vroom::vroom_write(pcawg_stm_ccf_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/pcawg_stm_ccf_data.tsv')
# pcawg_stm_ccf_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/pcawg_stm_ccf_data.tsv')


# combin PCAWG mutation attribution and clonality

pcawg_stm_ccf_data_sim <- pcawg_stm_ccf_data %>% 
  mutate(mut_id = str_c(icgc_specimen_id, chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2, sep = ':')) %>% 
  select(mut_id, clonality = ccf_ci95_clonality)
# 20,373,807

pcawg_mutsig_clo_alt_data <- pcawg_sub_data_sim %>% 
  left_join(pcawg_stm_ccf_data_sim, by = 'mut_id') %>% 
  mutate(Chromosome = str_c('chr', Chromosome))
# 21,606,842

# vroom::vroom_write(pcawg_mutsig_clo_alt_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/pcawg_mutsig_clo_alt_data.tsv')
# pcawg_mutsig_clo_alt_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/pcawg_mutsig_clo_alt_data.tsv')


# combine PanCanAtlas and PCAWG mutsig clonality data ---------------------

Cancer_type_conv_tab <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/ICGC_TCGA_Cancer_Type_Convert.csv')


############################## method from science translational medicine, based on ccf

tcga_mutsig_clo_alt_data_t <- tcga_mutsig_clo_alt_data %>% 
  fuzzyjoin::regex_left_join(Cancer_type_conv_tab %>% select(ICGC_abbr_top, TCGA_abbr), by = c('Ori_Cancer_Type' = 'TCGA_abbr'), ignore_case = TRUE) %>% 
  select(-TCGA_abbr)

comb_mutsig_clo_alt_data <- pcawg_mutsig_clo_alt_data %>% 
  rename(Sample = icgc_specimen_id, Ori_Cancer_Type = histology_abbreviation) %>% 
  fuzzyjoin::regex_left_join(Cancer_type_conv_tab %>% select(ICGC_abbr_top, ICGC_abbr), by = c('Ori_Cancer_Type' = 'ICGC_abbr'), ignore_case = TRUE) %>% 
  select(-ICGC_abbr) %>% 
  bind_rows(tcga_mutsig_clo_alt_data_t)
# 24,646,682

vroom::vroom_write(comb_mutsig_clo_alt_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_alt_data.tsv')

# Signature fitting -------------------------------------------------------

library(deconstructSigs)
source('/pub6/Temp/Liaojl/Code/common_R_functions/whichSignatures.R')
source('/pub6/Temp/Liaojl/Code/common_R_functions/mutsig_fit.R')
load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_mutsig_w96_data.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mut_sam_project.Rdata')
comb_mutsig_clo_alt_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_alt_data.tsv')
comb_mutsig_we <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')
cancer_mutsig <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/cancer_mutsig.tsv') # mutsig freq >= 0.05

sample_mutsig_names <- comb_mutsig_we %>% distinct(ICGC_abbr_top, Sample, Signature) %>% semi_join(cancer_mutsig, by = c('ICGC_abbr_top', 'Signature'))
# save(sample_mutsig_names, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/sample_mutsig_names.Rdata')

comb_mutsig_clo_alt_data_rmna <- comb_mutsig_clo_alt_data %>% filter(!is.na(clonality))

suff_sam <- comb_mutsig_clo_alt_data_rmna %>% 
  count(Sample, clonality) %>% # both > 50: 3,467 samples; both > 30: 4,685
  pivot_wider(names_from = clonality, values_from = n, values_fill = 0) %>% 
  filter(clonal > 30, subclonal > 30) %>% 
  inner_join(mut_sam_project, by = 'Sample') %>% 
  semi_join(sample_mutsig_names, by = 'Sample') %>% 
  select(Sample, Project)
# 4,276
# save(suff_sam, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/suff_sam.Rdata')

suff_sam_strict <- comb_mutsig_clo_alt_data_rmna %>% 
  count(Sample, clonality) %>% # both > 50: 3,467 samples; both > 30: 4,685
  pivot_wider(names_from = clonality, values_from = n, values_fill = 0) %>% 
  filter(clonal > 50, subclonal > 50) %>% 
  inner_join(mut_sam_project, by = 'Sample') %>% 
  semi_join(sample_mutsig_names, by = 'Sample') %>% 
  select(Sample, Project)
# 3,304
# save(suff_sam_strict, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/suff_sam_strict.Rdata')


clo_mut_data <- comb_mutsig_clo_alt_data_rmna %>% 
  filter(clonality %in% 'clonal') %>% 
  mutate(Reference_Allele = str_sub(mut_id, -3, -3), Tumor_Seq_Allele2 = str_sub(mut_id, -1, -1)) %>% 
  select(Sample, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2)
# 19,962,550
vroom::vroom_write(clo_mut_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/clo_mut_data.tsv')

sub_mut_data <- comb_mutsig_clo_alt_data_rmna %>% 
  filter(clonality %in% 'subclonal') %>% 
  mutate(Reference_Allele = str_sub(mut_id, -3, -3), Tumor_Seq_Allele2 = str_sub(mut_id, -1, -1)) %>% 
  select(Sample, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2)
# 3,153,055
vroom::vroom_write(sub_mut_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/sub_mut_data.tsv')

ssam_sig_fit <- function(sample, project = 'PCAWG', mut_data, mutsig_w_data){
  
  # sample <- '0145'
  # project <- 'TCGA'
  # mut_data <- clo_mut_data
  # mutsig_w_data <- pcawg_mutsig_w96_data
  # project: TCGA or PCAWG, it affect the input of tri_counts_method
  
  tri_counts_method_input <- ifelse(project == 'PCAWG', 'default', 'exome2genome')
  
  sam_sig <- sample_mutsig_names %>% filter(Sample %in% sample) %>% pull(Signature)
  sam_signatures <- mutsig_w_data[sam_sig, ]
  
  ssam_mut <- mut_data %>% filter(Sample %in% sample)
  
  ssam_pre_mut <- PreMutData(ssam_mut, 'Sample', 'Chromosome', 'Start_position', 'Reference_Allele', 'Tumor_Seq_Allele2')
  ssam_fit <- RunDeSigs(ssam_pre_mut, sam_signatures, tri_counts_method = tri_counts_method_input)
  ssam_mfit <- MerSamSig(ssam_pre_mut, ssam_fit)
  
  return(ssam_mfit$mutSig_we_ori)
  
}

# clonal signature fitting


clo_sig_fit <- suff_sam %>% mutate(fit_res = map2(Sample, Project, ~ssam_sig_fit(.x, .y, clo_mut_data, pcawg_mutsig_w96_data))) %>% select(fit_res) %>% unnest(fit_res)

# subclonal signature fitting

sub_sig_fit <- suff_sam %>% mutate(fit_res = map2(Sample, Project, ~ssam_sig_fit(.x, .y, sub_mut_data, pcawg_mutsig_w96_data))) %>% select(fit_res) %>% unnest(fit_res)

# merge

clonality_sig_fit <- clo_sig_fit %>% 
  mutate(clonality = 'clonal') %>% 
  bind_rows(mutate(sub_sig_fit, clonality = 'subclonal')) %>% 
  inner_join(distinct(sample_mutsig_names, ICGC_abbr_top, Sample), by = 'Sample')
# 4,276 samples

save(clonality_sig_fit, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/clonality_sig_fit.Rdata')

clonality_sig_fit_strict <- clonality_sig_fit %>% semi_join(suff_sam_strict, by = 'Sample')
# 3,304 samples

save(clonality_sig_fit_strict, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/clonality_sig_fit_strict.Rdata')
