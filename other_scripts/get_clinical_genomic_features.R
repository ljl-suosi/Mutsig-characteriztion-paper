
library(tidyverse)

Cancer_type_conv_tab <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/ICGC_TCGA_Cancer_Type_Convert.csv')

tcga_cli_data <- read_tsv('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_cli.tsv', guess_max = 5000)
tcga_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv')
tcga_mut_sam <- tcga_mut_data %>% distinct(Tumor_Sample_Barcode) %>% transmute(Sample = str_sub(Tumor_Sample_Barcode, 1, 15)) # 10,103
load('/pub6/Temp/Liaojl/Data/SBS18_pancancer_research/Data/Processed/pcawg_wgs_cli.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_pp_data.Rdata')

tcga_pp_data <- read_tsv('/pub6/Temp/Liaojl/Data/Mutation_burden_research/Original/TCGA_mastercalls.abs_tables_JSedit.fixed.txt') %>% 
  semi_join(tcga_mut_sam, by = c('array' = 'Sample')) %>% 
  mutate(patient = str_sub(array, 9, 12), 
         wgd_status = case_when(
           `Genome doublings` == 0 ~ 'no_wgd', 
           `Genome doublings` > 0 ~ 'wgd', 
           TRUE ~ NA_character_
         )) %>%
  select(Sample = patient, purity, ploidy, wgd_status)
# 9859

comb_pp_data <- pcawg_pp_data %>% rename(Sample = icgc_specimen_id) %>% bind_rows(tcga_pp_data)


# clinical variables ------------------------------------------------------


# PCAWG: only extract primary samples except for Skin-Melanoma

pcawg_pri_cli <- pcawg_wgs_cli %>% 
  filter(str_detect(dcc_specimen_type, 'Primary')) %>% 
  mutate(OS = ifelse(donor_vital_status %in% 'deceased', 1, 0)) %>% 
  select(Ori_Cancer_Type = histology_abbreviation, Sample = icgc_specimen_id, age = donor_age_at_diagnosis, gender = donor_sex, 
         histological_type = tumour_histological_type, stage = tumour_stage, grade = tumour_grade, first_therapy_response, OS, OS.time = donor_survival_time) %>% 
  fuzzyjoin::regex_left_join(Cancer_type_conv_tab %>% select(ICGC_abbr_top, ICGC_abbr), by = c('Ori_Cancer_Type' = 'ICGC_abbr'), ignore_case = TRUE) %>% 
  select(-ICGC_abbr) %>% 
  distinct_all()
# 2,748

# TMB ---------------------------------------------------------------------

comb_mut_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mut_data.tsv')

comb_TMB <- comb_mut_data %>% 
  filter(Variant_Type %in% 'SNP', Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation')) %>% 
  count(Sample, name = 'Mut_num') %>% 
  mutate(TMB = Mut_num/38) %>% 
  select(Sample, TMB)

# write_tsv(comb_TMB, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_TMB.tsv')

# combine PanCanAtlas and PCAWG clinical data

comb_cli_data <- tcga_cli_data %>% 
  mutate(gender = str_to_lower(gender)) %>% 
  select(Ori_Cancer_Type = type, Sample = patient_id, age = age_at_initial_pathologic_diagnosis, gender, histological_type, stage = ajcc_pathologic_tumor_stage, 
         grade = histological_grade, first_therapy_response = treatment_outcome_first_course, OS, OS.time) %>% 
  left_join(Cancer_type_conv_tab %>% select(ICGC_abbr_top, TCGA_abbr), by = c('Ori_Cancer_Type' = 'TCGA_abbr')) %>% 
  bind_rows(pcawg_pri_cli) %>% 
  left_join(comb_pp_data, by = 'Sample') %>% 
  mutate(stage = str_remove(stage, 'Stage '),
         stage = case_when(
    str_detect(stage, 'T1') ~ 'I', 
    str_detect(stage, 'T2') ~ 'II', 
    str_detect(stage, 'T3') ~ 'III', 
    str_detect(stage, 'T4') ~ 'IV', 
    str_detect(stage, 'TX') ~ NA_character_, 
    str_detect(stage, 'Tx') ~ NA_character_, 
    stage %in% c('1', '1a', '1b', 'IA', 'IB') ~ 'I', 
    stage %in% c('2', '2a', '2b', 'IIA', 'IIB') ~ 'II', 
    stage %in% c('3', '3a', '3b', '3c', 'IIIA', 'IIIB', 'IIIC') ~ 'III', 
    stage %in% c('4', 'IVA', 'IVB', 'IVC') ~ 'IV', 
    stage %in% c('I/II NOS', 'N/A', 'unknown') ~ NA_character_, 
    TRUE ~ stage
  ), grade = case_when(
    grade %in% c('1', 'I', 'NET-G1', '1 - Well differentiated', 'Well differentiated', 'I-II', 'I-II', 'I-III') ~ 'G1', 
    grade %in% c('2', 'II', 'NET-G2', '2 - Moderately differentiated', 'Moderate to Poor', 'Moderately differentiated', 'II-I', 'II-II', 'II-III') ~ 'G2',
    grade %in% c('III', '3', '3 - Poorly differentiated', 'NET-G3', 'Poorly differentiated') ~ 'G3', 
    grade %in% c('4', 'IV', '4 - Undifferentiated', 'Undifferentiated') ~ 'G4', 
    grade %in% c('PD', 'WD', 'needs cpath clarification', 'unknown', 'X - Cannot be assessed', 'Status Post Therapy') ~ NA_character_, 
    TRUE ~ grade
  ), first_therapy_response = case_when(
    first_therapy_response %in% c('Complete Remission/Response', 'complete response') ~ 'complete/partial response', 
    first_therapy_response %in% c('Partial Remission/Response', 'partial response') ~ 'complete/partial response', 
    first_therapy_response %in% c('Progressive Disease', 'disease progression') ~ 'stable/progressive disease', 
    first_therapy_response %in% c('Stable Disease', 'Persistent Disease', 'stable disease') ~ 'stable/progressive disease', 
    first_therapy_response %in% c('unknown', 'No Measureable Tumor or Tumor Markers', 'Normalization of Tumor Markers, but Residual Tumor Mass') ~ NA_character_, 
    TRUE ~ first_therapy_response
  )) %>% 
  left_join(comb_TMB, by = 'Sample')

save(comb_cli_data, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')



# Cancer Gene Mutation Status ----------------------------------------------------

source('/pub6/Temp/Liaojl/Code/common_R_functions/getGeneSamMat.R')

# use Cancer Genes from TCGA, as for Cancer types not found in TCGA, use PCAWG driver gene instead

# PCAWG driver gene

load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_driver.Rdata') # PCAWG driver

pcawg_driver_gene <- pcawg_driver %>% 
  filter(top_category %in% 'mutational', ttype %in% c('Bone-Benign', 'Bone-Epith', 'Bone-Osteosarc', 'CNS-Medullo', 'CNS-Oligo', 'CNS-PiloAstro', 'Lymph-CLL', 'Myeloid-MDS', 'Myeloid-MPN', 'Panc-Endocrine', 'SoftTissue-Leiomyo', 'SoftTissue-Liposarc')) %>% 
  select(Ori_Cancer_Type = ttype, Hugo_Symbol = gene) %>% 
  fuzzyjoin::regex_left_join(Cancer_type_conv_tab %>% select(ICGC_abbr_top, ICGC_abbr), by = c('Ori_Cancer_Type' = 'ICGC_abbr'), ignore_case = TRUE) %>% 
  distinct(ICGC_abbr_top, Hugo_Symbol)

# TCGA driver gene

tcga_driver_gene <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/tcga_panatlas_cancer_gene.txt') %>% 
  mutate(Cancer = ifelse(Cancer %in% 'COADREAD', 'COAD', Cancer)) %>% 
  filter(Cancer != 'PANCAN') %>% 
  inner_join(Cancer_type_conv_tab, by = c('Cancer' = 'TCGA_abbr')) %>% 
  distinct(ICGC_abbr_top, Gene)

comb_driver_gene <- tcga_driver_gene %>% rename(Hugo_Symbol = Gene) %>% bind_rows(pcawg_driver_gene)

# write_tsv(comb_driver_gene, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')

driver_gene_mut_data <- comb_mut_data %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  semi_join(comb_driver_gene, by = 'Hugo_Symbol')

total_mut_sam <- comb_mut_data %>% distinct(Sample) %>% pull(Sample)

driver_gs_mat <- getGeneSamMat(driver_gene_mut_data, 'Hugo_Symbol', 'Sample', total_mut_sam = total_mut_sam, status = 'ToF')
driver_gs_status <- driver_gs_mat %>% t() %>% as_tibble(rownames = 'Sample')

# save(driver_gs_status, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/driver_gs_status.Rdata')
# load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/driver_gs_status.Rdata')

# add driver gene copy number alteration status

## PCAWG gene copy number level data

pcawg_gene_cn_level <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/all_samples.consensus_level_calls.by_gene.170214.txt')
pcawg_sample_sheet <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/pcawg_sample_sheet.tsv')


pcawg_driver_gene_cnl <- pcawg_gene_cn_level %>% 
  filter(`Gene Symbol` %in% colnames(driver_gs_status)[-1]) %>% 
  select(-(`Locus ID`:Cytoband)) %>% 
  column_to_rownames('Gene Symbol') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('aliquot_id') %>% 
  as_tibble() %>% 
  inner_join(distinct(pcawg_sample_sheet, aliquot_id, icgc_specimen_id), by = 'aliquot_id') %>% 
  select(-aliquot_id) %>% 
  select(Sample = icgc_specimen_id, everything())


## TCGA gene copy number level data

tcga_gene_cn_level <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/all_thresholded.by_genes_whitelisted.tsv')

tcga_driver_gene_cnl <- tcga_gene_cn_level %>% 
  filter(`Gene Symbol` %in% colnames(driver_gs_status)[-1]) %>% 
  select(-(`Locus ID`:Cytoband)) %>% 
  column_to_rownames('Gene Symbol') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  as_tibble() %>% 
  mutate(Sample = str_sub(Sample, 9, 12)) %>% 
  distinct(Sample, .keep_all = TRUE)

comb_driver_gene_cnv_status <- pcawg_driver_gene_cnl %>% 
  bind_rows(tcga_driver_gene_cnl) %>% 
  semi_join(driver_gs_status, by = 'Sample') %>% 
  mutate_if(is.numeric, as.logical) %>% 
  pivot_longer(-Sample, names_to = 'Hugo_Symbol', values_to = 'cnv_status')

vroom::vroom_write(comb_driver_gene_cnv_status, '/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene_cnv_status.tsv')


# key pathway aberration --------------------------------------------

# 10 Oncogenic pathways

onco_pathway_gene <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/curated_pathways.txt') %>% 
  select(pathway = Pathway, Hugo_Symbol = Gene) %>% 
  mutate(type = 'oncogenic_pathway')
# 334 uique genes, SCRIB exists both in HIPPO and RTK-RAS pathways

# save(onco_pathway_gene, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/onco_pathway_gene.Rdata')

onco_pathway_mut_data <- comb_mut_data %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  inner_join(onco_pathway_gene, by = 'Hugo_Symbol')
# 333 genes left

total_mut_sam <- comb_mut_data %>% distinct(Sample) %>% pull(Sample)

onco_pathway_sam_mat <- getGeneSamMat(onco_pathway_mut_data, 'pathway', 'Sample', total_mut_sam = total_mut_sam, status = 'ToF')
onco_pathway_sam_status <- onco_pathway_sam_mat %>% t() %>% as_tibble(rownames = 'Sample')

save(onco_pathway_sam_status, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/onco_pathway_sam_status.Rdata')

onco_pathway_gene_sam_mat <- getGeneSamMat(onco_pathway_mut_data, 'Hugo_Symbol', 'Sample', total_mut_sam = total_mut_sam, status = 'ToF')
onco_pathway_gene_sam_status <- onco_pathway_gene_sam_mat %>% t() %>% as_tibble(rownames = 'Sample')

save(onco_pathway_gene_sam_status, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/onco_pathway_gene_sam_status.Rdata')


# 9 DDR pathways core genes

missense_deleterious_mut <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/missense_deleterious_mut.tsv')

DDR_gene <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/DDR_core_gene.tsv') %>% 
  pivot_longer(BER:DS, names_to = 'pathway', values_to = 'Hugo_Symbol') %>% 
  filter(!is.na(Hugo_Symbol)) %>% 
  arrange(pathway) %>% 
  mutate(type = 'DDR_pathway')
# 80 genes

# key_pathway_gene <- onco_pathway_gene %>% bind_rows(DDR_gene)
# 412 unique genes, except SCRIB mentioned above, ATM and CHEK2 exist both in TP53 and DS pathways
# save(key_pathway_gene, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/key_pathway_gene.Rdata')

# save(DDR_gene, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_gene.Rdata')

DDR_pathway_trunc_mut_data <- comb_mut_data %>% 
  filter(Variant_Classification %in% c('Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  inner_join(DDR_gene, by = 'Hugo_Symbol')

DDR_pathway_deleterious_missense_data <- comb_mut_data %>% 
  inner_join(DDR_gene, by = 'Hugo_Symbol') %>% 
  mutate(mut_id = str_c(Sample, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ':')) %>% 
  semi_join(missense_deleterious_mut, by = 'mut_id')

DDR_pathway_mut_data <- DDR_pathway_deleterious_missense_data %>% select(-mut_id) %>% bind_rows(DDR_pathway_trunc_mut_data)
# 8,795 mutations
# 80 genes left

total_mut_sam <- comb_mut_data %>% distinct(Sample) %>% pull(Sample)

DDR_pathway_sam_mat <- getGeneSamMat(DDR_pathway_mut_data, 'pathway', 'Sample', total_mut_sam = total_mut_sam, status = 'ToF')
DDR_pathway_sam_status <- DDR_pathway_sam_mat %>% t() %>% as_tibble(rownames = 'Sample')

save(DDR_pathway_sam_status, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_pathway_sam_status.Rdata')


DDR_pathway_gene_sam_mat <- getGeneSamMat(DDR_pathway_mut_data, 'Hugo_Symbol', 'Sample', total_mut_sam = total_mut_sam, status = 'ToF')
DDR_pathway_gene_sam_status <- DDR_pathway_gene_sam_mat %>% t() %>% as_tibble(rownames = 'Sample')

save(DDR_pathway_gene_sam_status, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_pathway_gene_sam_status.Rdata')



# ROS pathway -------------------------------------------------------------


load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Original/kegg_ros_gene.Rdata') # 223 genes

ROS_pathway_mut_data <- comb_mut_data %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins'), 
         Hugo_Symbol %in% kegg_ros_gene) %>% 
  mutate(pathway = 'Reactive Oxygen Species')

# 196 genes left

ROS_pathway_sam_mat <- getGeneSamMat(ROS_pathway_mut_data, 'pathway', 'Sample', total_mut_sam = total_mut_sam, status = 'ToF')
ROS_pathway_sam_status <- ROS_pathway_sam_mat %>% t() %>% as_tibble(rownames = 'Sample')

ROS_gene_sam_mat <- getGeneSamMat(ROS_pathway_mut_data, 'Hugo_Symbol', 'Sample', total_mut_sam = total_mut_sam, status = 'ToF')
ROS_gene_sam_status <- ROS_gene_sam_mat %>% t() %>% as_tibble(rownames = 'Sample')

save(ROS_pathway_sam_status, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/ROS_pathway_sam_status.Rdata')
save(ROS_gene_sam_status, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/ROS_gene_sam_status.Rdata')


