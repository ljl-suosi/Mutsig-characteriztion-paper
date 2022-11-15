library(tidyverse)

load('/pub6/Temp/Liaojl/Data/SBS18_pancancer_research/Data/Processed/pcawg_wgs_cli.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_cnpp_data.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/pcawg_sub_cnpp_data.Rdata')

# PRAD metastasis donor sample type and purity

pcawg_sam_purity <- pcawg_cnpp_data %>% distinct(icgc_specimen_id, purity) %>% rename(cellularity = purity)

prost_met_donor <- pcawg_wgs_cli %>% 
  filter(histology_abbreviation %in% 'Prost-AdenoCA', !str_detect(dcc_specimen_type, 'Primary')) %>% 
  distinct(icgc_donor_id)

prost_met_cli <- pcawg_wgs_cli %>% semi_join(prost_met_donor, by = 'icgc_donor_id') # all PARD-UK

prost_met_samlab <- prost_met_cli %>% 
  mutate(subsample = map_chr(str_split(submitted_specimen_id, '_'), last)) %>% 
  select(icgc_donor_id, icgc_specimen_id, subsample) %>% 
  left_join(pcawg_sam_purity, by = 'icgc_specimen_id') %>% 
  filter(cellularity >= 0.4) # exclude samples with very low purity
# DO51955 T1: 0.263, M1: 0.32

save(prost_met_samlab, file = '/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/prost_met_samlab.Rdata')

# load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/prost_met_samlab.Rdata')


DPClust_donor <- prost_met_samlab %>% 
  filter(!((icgc_donor_id %in% 'DO51954' & str_detect(subsample, 'M')  & cellularity < 0.8) | (icgc_donor_id %in% 'DO51964' & str_detect(subsample, 'M') & cellularity < 0.9)))


# PRAD metastasis DPClust data preparation


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

prost_multi_cnpp <- pcawg_sub_cnpp_data %>% semi_join(DPClust_donor, by = 'icgc_specimen_id')

prost_multi_somatics <- prost_multi_cnpp %>% 
  mutate(sampleId = icgc_specimen_id, 
         totalReadCount = t_alt_count + t_ref_count, 
         normal_cn = 2, 
         adjustedVaf = t_alt_count / (purity * totalReadCount), 
         adjustedCopyNumber = normal_cn * (1 - purity) + total_cn * purity, 
         somaticPloidy = pmax(0, adjustedVaf * adjustedCopyNumber)) %>%  # mutation copy number equal to VAF/p * (Nnormal*(1-p) + Ntumor*p)
  rename(chromosome = Chromosome, alleleReadCount = t_alt_count) %>% 
  filter(!chromosome %in% c('X', 'Y'), !is.na(somaticPloidy), major_cn >= 1) %>%  # only retain chromosome 1:22 required by DPClust
  mutate(absolute_ccf = pmap_dbl(list(alleleReadCount, totalReadCount, purity, total_cn), absolute_ccf_estimation))


# no.chrs.bearing.mut

prost_met_somatics <- prost_multi_somatics %>% 
  mutate(no.chrs.bearing.mut = 1, 
         subclonal.fraction = absolute_ccf) %>% 
  select(sample = Donor_ID, icgc_specimen_id, chr = chromosome, end = Start_position, 
         WT.count = t_ref_count, mut.count = alleleReadCount, subclonal.CN = total_cn, 
         mutation.copy.number = absolute_ccf, subclonal.fraction, no.chrs.bearing.mut) # set mutation.copy.number = absolute_ccf * no.chrs.bearing.mut = absolute_ccf * 1

# retain prostate mutation data for combining cluster later in order to conduct signature analysis in each cluster
# mutations were combined based on donor but not sample

prost_met_unidonor_mut <- prost_multi_somatics %>% 
  distinct(Donor_ID, chromosome, Start_position, .keep_all = TRUE) %>% 
  select(Donor_ID, Hugo_Symbol:Start_position, Variant_Classification, Reference_Allele:Tumor_Seq_Allele2)

# vroom::vroom_write(prost_met_unidonor_mut, '/pub6/Temp/Liaojl/SBS18_pancancer_research/Data/Processed/prost_met_unidonor_mut.tsv')

vroom::vroom_write(prost_met_unidonor_mut, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/Run_DPClust_on_nchr1_ccf/prost_met_unidonor_mut.tsv')

# focus on union of mutations, mut.count of mutations were set to 0 if absent on some samples

prost_donor_default <- prost_met_somatics %>% 
  distinct(sample, chr, end, .keep_all = TRUE) %>% 
  select(-icgc_specimen_id) %>% 
  mutate(mut.count = 0, subclonal.CN = 2, mutation.copy.number = 0, subclonal.fraction = 0, no.chrs.bearing.mut = 1) %>%
  left_join(select(DPClust_donor, -cellularity), by = c('sample' = 'icgc_donor_id')) %>% 
  select(sample, icgc_specimen_id, subsample, everything())


prost_met_dpclust_data <- prost_donor_default %>% 
  left_join(select(prost_met_somatics, -sample), by = c('icgc_specimen_id', 'chr', 'end')) %>% 
  mutate(WT.count = coalesce(WT.count.y, WT.count.x), # Find first non-missing element
         mut.count = coalesce(mut.count.y, mut.count.x), 
         subclonal.CN = coalesce(subclonal.CN.y, subclonal.CN.x), 
         mutation.copy.number = coalesce(mutation.copy.number.y, mutation.copy.number.x), 
         subclonal.fraction = coalesce(subclonal.fraction.y, subclonal.fraction.x), 
         no.chrs.bearing.mut = coalesce(no.chrs.bearing.mut.y, no.chrs.bearing.mut.x)) %>% 
  select(sample, subsample:end, WT.count:no.chrs.bearing.mut)

# vroom::vroom_write(prost_met_dpclust_data, '/pub6/Temp/Liaojl/Data/SBS18_pancancer_research/Data/Processed/prost_met_dpclust_data.tsv')

vroom::vroom_write(prost_met_dpclust_data, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/Run_DPClust_on_nchr1_ccf/prost_met_dpclust_data.tsv')


prost_met_gdata <- prost_met_dpclust_data %>% 
  group_by(sample, subsample) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(datafile = str_c(sample, subsample, 'dp_input.txt', sep = '_')) %>% 
  arrange(sample, subsample)

prost_met_dpinput_summ <- prost_met_gdata %>% 
  select(sample, subsample, datafile) %>% 
  left_join(select(DPClust_donor, -icgc_specimen_id), by = c('sample' = 'icgc_donor_id', 'subsample'))

# write_tsv(prost_met_dpinput_summ, '/pub6/Temp/Liaojl/Run_DPClust/Input/prost_met_dpinput_summ.txt')
write_tsv(prost_met_dpinput_summ, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/Run_DPClust_on_nchr1_ccf/Input/prost_met_dpinput_summ.txt')

# write txt files

for(i in 1:nrow(prost_met_gdata)){
  
  write_tsv(prost_met_gdata$data[[i]], str_c('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/Run_DPClust_on_nchr1_ccf/Input/Data/', prost_met_gdata$datafile[[i]]))
  
}

