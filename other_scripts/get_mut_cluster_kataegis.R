
library(tidyverse)
library(data.table) # necessary

source('/pub6/Temp/Liaojl/Code/common_R_functions/detect_mut_cluster.R')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mut_sam_project.Rdata')
comb_mutsig_clo_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_alt_data.tsv')
comb_mutsig_we <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

comb_mutsig_we_filter <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% 
  ungroup()

comb_mutsig_mut_dat <- comb_mutsig_clo_data %>% 
  semi_join(comb_mutsig_we_filter, by = c('Sample', 'Attri_Signature' = 'Signature')) %>%  # Samples whose signature accuary >= 0.8
  inner_join(mut_sam_project, by = 'Sample')
# vroom::vroom_write(comb_mutsig_mut_dat, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_mutsig_mut_dat.tsv')
# comb_mutsig_mut_dat <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_mutsig_mut_dat.tsv')

cancer_sam <- comb_mutsig_mut_dat %>% distinct(Sample, ICGC_abbr_top)

# clustered mutations (>=2 mutations in less than 1000bp)

all_mutsig_sam <- comb_mutsig_mut_dat %>% distinct(Sample) %>% pull(Sample)

all_clustermut_res_list1 <- list()
all_clustermut_res_list2 <- list()
all_clustermut_res_list3 <- list()
all_clustermut_res_list4 <- list()
all_clustermut_res_list5 <- list()

all_mutsig_sam_sub1 <- all_mutsig_sam[1:1600] # this one was time-consuming and should be divided into another several parts
all_mutsig_sam_sub2 <- all_mutsig_sam[1601:3200]
all_mutsig_sam_sub3 <- all_mutsig_sam[3201:4800]
all_mutsig_sam_sub4 <- all_mutsig_sam[4801:6400]
all_mutsig_sam_sub5 <- all_mutsig_sam[6401:7993]

for(i in all_mutsig_sam_sub5){
  
  print(str_c('processing ', i))
  
  sin_mut <- comb_mutsig_mut_dat %>% 
    filter(Sample %in% i) %>% 
    select(Chromosome, Hugo_Symbol, Start_Position = Start_position, End_Position = Start_position, Tumor_Sample_Barcode = Sample, con.class = Attri_Signature) %>% 
    mutate(Chromosome = str_remove(Chromosome, 'chr'), 
           Chromosome = case_when(
             Chromosome %in% 'X' ~ '23', 
             Chromosome %in% 'Y' ~ '24', 
             TRUE ~ Chromosome
           ), 
           Chromosome = factor(Chromosome, levels = 1:24, labels = 1:24)) %>% 
    arrange(Chromosome, Start_Position) %>% 
    as.data.table() # it's a must because of many data.table operations following
  
  sin_res <- detect_mut_clusters(sin_mut, cluster_mut_num = 2)
  
  all_clustermut_res_list5[[i]] <- sin_res
  
}

# save(all_clustermut_res_list5, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_clustermut_res_list5.Rdata')

# merge results

# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_clustermut_res_list1.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_clustermut_res_list2.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_clustermut_res_list3.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_clustermut_res_list4.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_clustermut_res_list5.Rdata')

all_clustermut_res <- c(all_clustermut_res_list1, all_clustermut_res_list2, all_clustermut_res_list3, all_clustermut_res_list4, all_clustermut_res_list5) %>% 
  reduce(bind_rows) %>% # time-consuming, a simpler method may needed to merge two datasets without considering column names and types 
  as_tibble() %>% 
  left_join(cancer_sam, by = c('Tumor_Sample_Barcode' = 'Sample')) %>% 
  select(ICGC_abbr_top, Tumor_Sample_Barcode, everything())
# save(all_clustermut_res, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/all_clustermut_res.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/all_clustermut_res.Rdata')

# kataegis (>=6 mutations in less than 1000bp)

all_kataegis_res_list1 <- list()
all_kataegis_res_list2 <- list()
all_kataegis_res_list3 <- list()
all_kataegis_res_list4 <- list()
all_kataegis_res_list5 <- list()

for(i in all_mutsig_sam_sub5){
  
  print(str_c('processing ', i))
  
  sin_mut <- comb_mutsig_mut_dat %>% 
    filter(Sample %in% i) %>% 
    select(Chromosome, Hugo_Symbol, Start_Position = Start_position, End_Position = Start_position, Tumor_Sample_Barcode = Sample, con.class = Attri_Signature) %>% 
    mutate(Chromosome = str_remove(Chromosome, 'chr'), 
           Chromosome = case_when(
             Chromosome %in% 'X' ~ '23', 
             Chromosome %in% 'Y' ~ '24', 
             TRUE ~ Chromosome
           ), 
           Chromosome = factor(Chromosome, levels = 1:24, labels = 1:24)) %>% 
    arrange(Chromosome, Start_Position) %>% 
    as.data.table() # it's a must because of many data.table operations following
  
  sin_res <- detect_mut_clusters(sin_mut)
  
  all_kataegis_res_list5[[i]] <- sin_res
  
}

save(all_kataegis_res_list5, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_kataegis_res_list5.Rdata')

# merge results

# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_kataegis_res_list1.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_kataegis_res_list2.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_kataegis_res_list3.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_kataegis_res_list4.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/all_kataegis_res_list5.Rdata')

all_kataegis_res <- c(all_kataegis_res_list1, all_kataegis_res_list2, all_kataegis_res_list3, all_kataegis_res_list4, all_kataegis_res_list5) %>% 
  reduce(bind_rows) %>% 
  as_tibble() %>% 
  left_join(cancer_sam, by = c('Tumor_Sample_Barcode' = 'Sample')) %>% 
  select(ICGC_abbr_top, Tumor_Sample_Barcode, everything())

# save(all_kataegis_res, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/all_kataegis_res.Rdata')


# CPA6 case for explanating the difference between cluster distance -------

# cluster distance = 2 or 6

chr8_dat <- comb_mutsig_clo_data %>% filter(Tumor_Sample_Barcode %in% 'SP101724', Hugo_Symbol %in% 'CPA6')

CPA6_dat <- chr8_dat %>% filter(Hugo_Symbol %in% 'CPA6') %>% mutate(mut_id = 1:16)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/SP101724_CPA6_mut_lolipop.pdf', width = 14, heigh = 3)

ggplot(CPA6_dat, aes(Start_Position, 1, color = con.class)) +
  geom_segment(aes(x = Start_Position, xend = Start_Position, y = 0, yend = 1, color = con.class), lwd = 0.5) +
  geom_point(size = 2) +
  annotate('rect', xmin = 68610000, xmax = 68650000, ymin = 0, ymax = 0.1) +
  scale_x_continuous(expand = c(0, NA)) +
  labs(x = 'genomic position', y = NULL, col = NULL) +
  geom_text_repel(aes(label = mut_id), show.legend = FALSE) +
  theme_classic() +
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

dev.off()

