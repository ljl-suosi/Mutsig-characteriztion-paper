
library(tidyverse)
library(data.table)
library(patchwork)

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mut_sam_project.Rdata')
comb_mutsig_clo_data <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_alt_data.tsv')
comb_mutsig_we <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

comb_mutsig_we_filter <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% # signatures present in >= 5% samples
  ungroup()

comb_mutsig_mut_dat <- comb_mutsig_clo_data %>% 
  semi_join(comb_mutsig_we_filter, by = c('Sample', 'Attri_Signature' = 'Signature')) %>%  # Samples whose signature accuary >= 0.8
  inner_join(mut_sam_project, by = 'Sample')

# comb_mutsig_mut_dat <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_mutsig_mut_dat.tsv')

rf_cancer_colors <- tribble(
  ~cancer_type, ~color, 
  'Adrenal-neoplasm', '#0C487A', 
  'Biliary-AdenoCA', '#BAA133', 
  'Bladder-TCC', '#F8D0D6', 
  'Bone-Osteosarc', '#EE9028', 
  'Bone-Other', '#0C8D46', 
  'Breast', '#D92F86', 
  'Cervix', '#ECB065', 
  'CNS-GBM', '#A74C93', 
  'CNS-LGG', '#DD2929', 
  'CNS-Medullo', '#374F99', 
  'CNS-PiloAstro', '#CC98BF', 
  'ColoRect-AdenoCA', '#9DD6F0', 
  'Eso-AdenoCA', '#0E78AA', 
  'Eye-Melanoma', '#DC6B72', 
  'Head-SCC', '#93C9A3', 
  'Kidney-ChRCC', '#7A378B', #
  'Kidney-Papillary', '#B2212D', 
  'Kidney-RCC', '#EEAAAE', 
  'Liver-HCC', '#C7C9D7', 
  'Lung-AdenoCA', '#CEC0DC', 
  'Lung-SCC', '#997FB5', 
  'Lymph-BNHL', '#724C29', 
  'Lymph-CLL', '#502D80', 
  'Mesothelium-Mesothelioma', '#EDE33E', 
  'Myeloid-AML', '#D8EEFB', #
  'Myeloid-MDS/MPN', '#0FA094', 
  'Ovary-AdenoCA', '#CF7A29', 
  'Panc-AdenoCA', '#6B789B', 
  'Panc-Endocrine', '#DCBE29', 
  'Pheochromocytoma', '#CD6600', 
  'Prost-AdenoCA', '#7A1A1D', 
  'Sarcoma', '#CAA98D', 
  'Skin-Melanoma', '#B3CB47', 
  'Stomach-AdenoCA', '#2FA0D2', 
  'Testis-Ca', '#008B8B', 
  'Thy-AdenoCA', '#9370DB', 
  'Thymoma', '#E066FF', 
  'UCS', '#FF8C69', 
  'Uterus-AdenoCA', '#F7DEC5', 
  'other', 'grey60'
)


mutsig_sam_dat <- comb_mutsig_mut_dat %>% distinct(ICGC_abbr_top, Attri_Signature, Sample) %>% rename(Signature = Attri_Signature)
mutsig_sam_count <- mutsig_sam_dat %>% count(Signature, name = 'sam_count')
# 51 signatures

mutsig_sam_mut_count <- comb_mutsig_mut_dat %>% 
  count(Project, ICGC_abbr_top, Sample, Attri_Signature, name = 'mut_count') %>% 
  group_by(Sample) %>% 
  mutate(mut_frac = mut_count/sum(mut_count)) %>% 
  ungroup() %>% 
  filter(mut_count >= 30, mut_frac >= 0.1) %>%  # samples and signatures with sufficient statistical power (i.e. signatures with mut count >= 30 and mut frac >= 0.1)
  group_by(Attri_Signature) %>% 
  filter(n() >= 30) %>% # signature with sufficient statistical power (i.e. signataure with sample count >= 30)
  ungroup() %>% 
  select(Project:Sample, Signature = Attri_Signature, mut_count)

# mutation rate analysis -----------------------------------------

library(ggnewscale)

source('/boot3/bio_liaojl/Common_R_functions/custom_colors_v1_221017.R')

comb_mutsig_mut_dat <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/Temp/comb_mutsig_mut_dat.tsv')
comb_driver_gene <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')
comb_snv_strand <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_snv_strand.tsv')


mutsig_order <- comb_mutsig_mut_dat %>% distinct(Attri_Signature) %>% pull(Attri_Signature) %>% str_sort(numeric = TRUE)

var_region_dat <- tibble(Variant_Classification = c("Intron", 'IGR', 'RNA', 'Missense_Mutation', "3'UTR", 'lincRNA', "5'Flank", 'Silent', "5'UTR", 'Splice_Site', 'Nonsense_Mutation', 'De_novo_Start_OutOfFrame', 'De_novo_Start_InFrame', 'Start_Codon_SNP', 'Nonstop_Mutation', "3'Flank", 'Translation_Start_Site'), 
                         region = c('Intron', 'Nongenic', 'Nongenic', 'Exon', 'UTR', 'Nongenic', 'Nongenic', 'Exon', 'UTR', 'Intron', 'Exon', 'Exon', 'Exon', 'Exon', 'Exon', 'Nongenic', 'Nongenic'))

region_length <- tibble(region = c('Exon', 'Intron', 'UTR', 'Nongenic'), 
                        length = c(76.30, 1192.31, 49.56, 1780.66))

comb_mutsig_mut_region <- comb_mutsig_mut_dat %>% # 23,422,116
  anti_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top')) %>% # remove cancer gene for avoiding effects of positive selection
  left_join(var_region_dat, by = 'Variant_Classification') %>% 
  inner_join(comb_snv_strand, by = 'mut_id')
# 23,385,730

comb_mutsig_mut_region_tmb <- comb_mutsig_mut_region %>% 
  count(Project, ICGC_abbr_top, Sample, Attri_Signature, region) %>% 
  inner_join(region_length, by = 'region') %>% 
  mutate(tmb = n/length)

pcawg_mutsig_mut_region_tmb <- comb_mutsig_mut_region_tmb %>% filter(Project %in% 'PCAWG')

# data with sufficient statistical power

sam_mutsig_power <- comb_mutsig_mut_dat %>% 
  anti_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top')) %>% # remove cancer gene for avoiding effects of positive selection
  count(Sample, Attri_Signature) %>% 
  group_by(Sample) %>% 
  mutate(frac = n/sum(n)) %>% 
  ungroup() %>% 
  filter(n >= 30, frac >= 0.1) # samples and signatures with sufficient statistical power (i.e. signatures with mut count >= 30 and mut frac >= 0.1)

pcawg_mutsig_mut_region_tmb_pow <- pcawg_mutsig_mut_region_tmb %>% 
  semi_join(sam_mutsig_power, c('Sample', 'Attri_Signature')) %>% 
  group_by(Project, region, Attri_Signature) %>% 
  filter(n() >= 30) %>% # region signature with sufficient statistical power (i.e. signataure within a certain region with sample count >= 30)
  ungroup()

# significance test

mutsig_tmb_test <- ggpubr::compare_means(tmb ~ region, data = pcawg_mutsig_mut_region_tmb_pow, group.by = 'Attri_Signature', method = 'kruskal.test')
# all q_val < 0.0001

signif_label <- mutsig_tmb_test %>% mutate(y_val = 250)

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_mutation_analysis_revise/mutsig_region_tmb.pdf', width = 15, height = 3)

pcawg_mutsig_mut_region_tmb_pow %>% 
  mutate(Attri_Signature = factor(Attri_Signature, levels = mutsig_order), 
         region = factor(region, levels = c('Exon', 'UTR', 'Intron', 'Nongenic'))) %>% 
  ggplot() +
  geom_boxplot(aes(Attri_Signature, tmb, fill = region), outlier.color = 'grey60') +
  geom_text(aes(Attri_Signature, y_val, label = p.signif), size = 5, data = signif_label) +
  labs(x = NULL, fill = NULL, y = 'Tumor Mutation Burden (mut/Mb)') +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), label = c('0.01', '0.1', '1', '10', '100')) +
  scale_fill_brewer(palette = 'Set3') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# SBS1 exon vs. nongenenic, p_val: 6.304865e-06
# SBS1 exon vs. intron, p_val: 0.2227526
# SBS1 exon vs. UTR, p_val: 1.468078e-09
# SBS1 intron vs. UTR, p_val: 1.231836e-07

# DNA damage repair pathway proficient vs deficient -----------------------

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/DDR_pathway_sam_status.Rdata')

DDR_pathway_sam_status_filter <- DDR_pathway_sam_status %>% 
  semi_join(pcawg_mutsig_mut_region_tmb_pow, by = 'Sample')

DDR_pathway_sam_status_filter$DDR <- DDR_pathway_sam_status_filter %>% select(-Sample) %>% apply(1, function(x) any(x == TRUE))

DDR_pathway_sam_status_final <- DDR_pathway_sam_status_filter %>% 
  pivot_longer(-Sample, names_to = 'pathway', values_to = 'status') %>% 
  group_by(pathway) %>% 
  filter(sum(status) >= 15) %>% # DR was excluded because of small number of samples showing deficiency in this pathway
  ungroup() %>% 
  mutate(status = ifelse(status == TRUE, 'deficient', 'proficient'), 
         path_lab = str_c(pathway, status, sep = '-'))

pcawg_mutsig_mut_region_tmb_pow_DDR <- pcawg_mutsig_mut_region_tmb_pow %>% 
  inner_join(DDR_pathway_sam_status_final, by = 'Sample') %>% 
  select(-Project, -n, -length) %>% 
  group_by(Attri_Signature, region, pathway) %>% 
  filter(sum(status %in% 'deficient') >= 10) %>% 
  ungroup()

DDR_region_test <- function(region_1, region_2){
  
  # region_1 <- 'Exon'
  # region_2 <- 'Intron'
  # ratio <- median(region_1)/median(region_2)
  
  test_res <- pcawg_mutsig_mut_region_tmb_pow_DDR %>% 
    filter(region %in% c(region_1, region_2)) %>% 
    mutate(region = factor(region, levels = c(region_1, region_2))) %>% 
    group_by(path_lab, Attri_Signature) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(p_val = map_dbl(data, ~wilcox.test(tmb ~ region, data = .)$p.value), 
           q_val = p.adjust(p_val, method = 'BH'), 
           log2_ratio = map_dbl(data, ~{
             temp_dat <- group_by(., region) %>% summarise(med = median(tmb))
             log2(temp_dat$med[1]/temp_dat$med[2])
             })) %>% 
    select(-data) %>% 
    mutate(Attri_Signature = factor(Attri_Signature, levels = mutsig_order)) %>% 
    arrange(Attri_Signature, path_lab)
  
}

exon_intron_test <- DDR_region_test('Exon', 'Intron')
exon_UTR_test <- DDR_region_test('Exon', 'UTR')
exon_nongenic_test <- DDR_region_test('Exon', 'Nongenic')

stripe_dat <- tibble(xst = c(0.5, 1.5), 
                     xed = c(1.5, 2.5), 
                     reg_col = c('odd', 'even'))

signature_colors_12 <- tribble(
  ~signature, ~color, 
  'SBS1', '#8F6504', 
  'SBS5', '#927C64', 
  'SBS40', '#FFB50E', 
  'SBS2', '#FF6444', 
  'SBS13', '#FFFA00', 
  'SBS18', '#C41707', 
  'SBS3', '#0373C6', 
  'SBS12', '#9C8CC4', 
  'SBS17a', '#90AD1C', 
  'SBS17b', '#FF5FB3', 
  'SBS4', '#B10DA1', 
  'SBS29', '#CC6677'
)

exon_intron_plot_dat <- exon_intron_test %>% 
  separate(path_lab, into = c('pathway', 'status'), sep = '-') %>% 
  mutate(q_val_size = case_when(
    q_val >= 0.05 ~ 0.5, 
    q_val < 0.05 & q_val > 0.01 ~ 1, 
    q_val < 0.01 & q_val > 0.001 ~ 1.25, 
    q_val < 0.001 & q_val > 0.0001 ~ 1.5, 
    q_val < 0.0001 ~ 1.75
  ))

exon_UTR_plot_dat <- exon_UTR_test %>% 
  separate(path_lab, into = c('pathway', 'status'), sep = '-') %>% 
  mutate(q_val_size = case_when(
    q_val >= 0.05 ~ 0.5, 
    q_val < 0.05 & q_val > 0.01 ~ 1, 
    q_val < 0.01 & q_val > 0.001 ~ 1.25, 
    q_val < 0.001 & q_val > 0.0001 ~ 1.5, 
    q_val < 0.0001 ~ 1.75
  ))

exon_nongenic_plot_dat <- exon_nongenic_test %>% 
  separate(path_lab, into = c('pathway', 'status'), sep = '-') %>% 
  mutate(q_val_size = case_when(
    q_val >= 0.05 ~ 0.5, 
    q_val < 0.05 & q_val > 0.01 ~ 1, 
    q_val < 0.01 & q_val > 0.001 ~ 1.25, 
    q_val < 0.001 & q_val > 0.0001 ~ 1.5, 
    q_val < 0.0001 ~ 1.75
  ))


##### Exon-to-Intron plot

exon_intron_plot <- ggplot() +
  geom_point(data = exon_intron_plot_dat, aes(status, log2_ratio, col = Attri_Signature, size = q_val_size), position = position_jitterdodge(jitter.width = 0.5)) +
  geom_rect(data = stripe_dat, aes(xmin = xst, xmax = xed, ymin = -Inf, ymax = Inf, fill = reg_col), show.legend = FALSE) +
  geom_hline(yintercept = 0, col = 'grey60') +
  labs(x = NULL, y = 'Exon-to-Intron ratio\n(log2)', col = NULL) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3), label = c('-3', '-2', '-1', '0', '1', '2', '3')) +
  scale_fill_manual(breaks = c('odd', 'even'), values = c('#33333333', '#00000000')) + 
  scale_color_manual(breaks = signature_colors_12$signature, values = signature_colors_12$color) +
  coord_cartesian(ylim = c(-3, 3)) +
  facet_wrap(~ pathway, strip.position = 'bottom', nrow = 1) +
  scale_size(range = c(1, 4)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.spacing.x = unit(0,'pt'), 
        strip.placement = 'outside', 
        strip.background = element_blank(), 
        legend.background = element_blank())

##### Exon-to-UTR plot

exon_UTR_plot <- ggplot() +
  geom_point(data = exon_UTR_plot_dat, aes(status, log2_ratio, col = Attri_Signature, size = q_val_size), position = position_jitterdodge(jitter.width = 0.5)) +
  geom_rect(data = stripe_dat, aes(xmin = xst, xmax = xed, ymin = -Inf, ymax = Inf, fill = reg_col), show.legend = FALSE) +
  geom_hline(yintercept = 0, col = 'grey60') +
  labs(x = NULL, y = 'Exon-to-UTR ratio\n(log2)', col = NULL) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3), label = c('-3', '-2', '-1', '0', '1', '2', '3')) +
  scale_fill_manual(breaks = c('odd', 'even'), values = c('#33333333', '#00000000')) + 
  scale_color_manual(breaks = signature_colors_12$signature, values = signature_colors_12$color) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap(~ pathway, strip.position = 'bottom', nrow = 1) +
  scale_size(range = c(1, 4)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(), 
        panel.spacing.x = unit(0,'pt'), 
        strip.placement = 'outside', 
        strip.text.x = element_blank(), 
        strip.background = element_blank(), 
        legend.background = element_blank())

##### Exon-to-Nongenic plot

exon_nongenic_plot <- ggplot() +
  geom_point(data = exon_nongenic_plot_dat, aes(status, log2_ratio, col = Attri_Signature, size = q_val_size), position = position_jitterdodge(jitter.width = 0.5)) +
  geom_rect(data = stripe_dat, aes(xmin = xst, xmax = xed, ymin = -Inf, ymax = Inf, fill = reg_col), show.legend = FALSE) +
  geom_hline(yintercept = 0, col = 'grey60') +
  labs(x = NULL, y = 'Exon-to-Nongeneic ratio\n(log2)', col = NULL) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3), label = c('-3', '-2', '-1', '0', '1', '2', '3')) +
  scale_fill_manual(breaks = c('odd', 'even'), values = c('#33333333', '#00000000')) + 
  scale_color_manual(breaks = signature_colors_12$signature, values = signature_colors_12$color) +
  coord_cartesian(ylim = c(-4, 4)) +
  facet_wrap(~ pathway, strip.position = 'bottom', nrow = 1) +
  scale_size(range = c(1, 4)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(), 
        panel.spacing.x = unit(0,'pt'), 
        strip.placement = 'outside', 
        strip.text.x = element_blank(), 
        strip.background = element_blank(), 
        legend.background = element_blank())

##### Combine exon ratio plots

library(patchwork)

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_mutation_analysis_revise/exon_ratio_comb.pdf', width = 8, height = 7.5)

exon_nongenic_plot + exon_UTR_plot + exon_intron_plot + plot_layout(nrow = 3, guides = 'collect') & theme(legend.position = 'bottom')

dev.off()


# MSI vs MSS

pcawg_msi_sam <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/pcawg_MSI_sample.tsv')
# 34 samples

pcawg_mutsig_mut_region_tmb_pow_msi <- pcawg_mutsig_mut_region_tmb_pow %>% 
  left_join(pcawg_msi_sam, by = 'Sample') %>% 
  mutate(MSI_status = ifelse(MSI %in% 'positive', 'MSI', 'MSS'), 
         region_lump = ifelse(region %in% 'Exon', 'Exon', 'Non-Exon')) %>% 
  select(Sample, Attri_Signature, region_lump, n, length, MSI_status) %>% 
  group_by(Sample, MSI_status, Attri_Signature, region_lump) %>% 
  summarise(n = sum(n), length = sum(length)) %>% 
  ungroup() %>% 
  mutate(tmb = n/length) %>% 
  group_by(Attri_Signature, region_lump) %>% 
  filter(sum(MSI_status %in% 'MSI') >= 2) %>% # only SBS1 and SBS5 have more than 2 MSI samples
  ungroup()
# 1929
# only 9 common MSI samples

msi_test_res <- pcawg_mutsig_mut_region_tmb_pow_msi %>% 
  group_by(Attri_Signature, MSI_status) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~wilcox.test(tmb ~ region_lump, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH'), 
         q_label = case_when(
           q_val < 0.0001 ~ '****', 
           q_val >= 0.0001 & q_val < 0.001 ~ '***', 
           q_val >= 0.001 & q_val < 0.01 ~ '**', 
           q_val >= 0.01 & q_val < 0.05 ~ '*', 
           q_val >= 0.05 ~ 'n.s.'
         )) %>% 
  select(-data)
# A tibble: 4 × 4
# MSI_status Attri_Signature    p_val    q_val
# <chr>      <chr>              <dbl>    <dbl>
# MSS        SBS5            4.03e-93 1.61e-92
# MSS        SBS1            1.44e- 2 1.97e- 2
# MSI        SBS1            4.85e- 1 4.85e- 1
# MSI        SBS5            1.48e- 2 1.97e- 2

# SBS1, 6 MSI samples; SBS5, 8 MSI samples

library(ggnewscale)

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_mutation_analysis_revise/temp_msi.pdf', width = 6, height = 5)

ggplot() +
  geom_boxplot(data = pcawg_mutsig_mut_region_tmb_pow_msi, aes(MSI_status, tmb, col = region_lump), fill = NA, outlier.color = NA) +
  scale_color_manual(breaks = c('Exon', 'Non-Exon'), values = c('grey60', 'grey60')) +
  new_scale_colour() +
  geom_point(data = pcawg_mutsig_mut_region_tmb_pow_msi, aes(MSI_status, tmb, col = region_lump), position = position_jitterdodge(jitter.width = 0.25), alpha = 7/10) +
  scale_color_manual(breaks = c('Exon', 'Non-Exon'), values = c('#8DD3C7', '#FF7F00')) +
  geom_text(data = msi_test_res, aes(MSI_status, 20, label = q_label)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), label = c('0.01', '0.1', '1', '10', '100')) +
  labs(x = NULL, y = 'Tumor Mutation Burden (mut/Mb)') +
  facet_grid(. ~ Attri_Signature) +
  theme_bw() +
  theme(strip.text.x = element_text(margin = margin(0.4, 0, 0.4, 0, "cm")), 
        legend.position = 'bottom')


dev.off()

# similar results when confined to only retain cancer types in which MSI exist

# HRD score

comb_HRD_status <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_HRD_status.tsv')

pcawg_mutsig_mut_region_tmb_pow_HRD <- pcawg_mutsig_mut_region_tmb_pow %>% 
  left_join(select(comb_HRD_status, Sample, HRD_status), by = 'Sample') %>% 
  mutate(region_lump = ifelse(region %in% 'Exon', 'Exon', 'Non-Exon')) %>% 
  select(Sample, Attri_Signature, region_lump, n, length, HRD_status) %>% 
  group_by(Sample, HRD_status, Attri_Signature, region_lump) %>% 
  summarise(n = sum(n), length = sum(length)) %>% 
  ungroup() %>% 
  mutate(tmb = n/length) %>% 
  group_by(Attri_Signature, region_lump) %>% 
  filter(sum(HRD_status %in% 'HR-deficient') >= 10) %>% 
  ungroup()

HRD_status_test_res <- pcawg_mutsig_mut_region_tmb_pow_HRD %>% 
  group_by(Attri_Signature, HRD_status) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~wilcox.test(tmb ~ region_lump, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH'), 
         q_label = case_when(
           q_val < 0.0001 ~ '****', 
           q_val >= 0.0001 & q_val < 0.001 ~ '***', 
           q_val >= 0.001 & q_val < 0.01 ~ '**', 
           q_val >= 0.01 & q_val < 0.05 ~ '*', 
           q_val >= 0.05 ~ 'n.s.'
         )) %>% 
  select(-data)


stripe_dat_re <- tibble(xst = c(0.5, 1.5), 
                     xed = c(1.5, 2.5), 
                     reg_col = c('HRD-deficient', 'HRD-proficient'))


pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_mutation_analysis_revise/temp_HRD.pdf', width = 14, height = 5)

ggplot() +
  geom_boxplot(data = pcawg_mutsig_mut_region_tmb_pow_HRD, aes(HRD_status, tmb, col = region_lump), fill = NA, outlier.color = 'grey60', position = position_dodge(width = 1)) +
  geom_rect(data = stripe_dat_re, aes(xmin = xst, xmax = xed, ymin = 0, ymax = Inf, fill = reg_col)) +
  geom_text(data = HRD_status_test_res, aes(HRD_status, 40, label = q_label)) +
  labs(x = NULL, fill = NULL, col = NULL, y = 'Tumor Mutation Burden (mut/Mb)') +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), label = c('0.01', '0.1', '1', '10', '100')) +
  scale_color_manual(breaks = c('Exon', 'Non-Exon'), values = c('#8DD3C7', '#FF7F00')) +
  scale_fill_manual(breaks = c('HRD-deficient', 'HRD-proficient'), values = c('#33333333', '#00000000')) + 
  facet_wrap(~ factor(Attri_Signature, levels = mutsig_order), strip.position = 'bottom', nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.spacing.x = unit(0,'pt'), 
        strip.placement = 'outside', 
        strip.background = element_blank(), 
        legend.background = element_blank(), 
        legend.position = 'top')

dev.off()

# evolutionary selection analysis -----------------------------------------

library(dndscv)

pcawg_mutsig_exon_dndscv_prep <- comb_mutsig_mut_region %>% 
  filter(region %in% 'Exon') %>% 
  semi_join(pcawg_mutsig_mut_region_tmb_pow, by = c('Sample', 'Attri_Signature')) %>% 
  select(Attri_Signature, mut_id) %>% 
  separate(mut_id, into = c('sampleID', 'chr', 'pos', 'ref', 'mut'), sep = ':')
  
mutsig_exon_select_tres <- pcawg_mutsig_exon_dndscv_prep %>% 
  group_by(Attri_Signature) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(global_dnds = map(data, possibly(~dndscv(.)$globaldnds, NULL)))

mutsig_exon_select_tdnds_all <- mutsig_exon_select_tres %>% 
  select(-data) %>% 
  unnest(global_dnds) %>% 
  filter(name %in% 'wall') %>% 
  mutate(name = 'all')

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_mutation_analysis_revise/mutsig_exon_select.pdf', width = 6, height = 4)

mutsig_exon_select_tdnds_all %>% 
  mutate(Attri_Signature = factor(Attri_Signature, levels = mutsig_order), 
         signif_lab = ifelse(Attri_Signature %in% c('SBS1', 'SBS7a', 'SBS12', 'SBS17a'), '*', '')) %>% 
  ggplot(aes(Attri_Signature, mle, ymin = cilow, ymax = cihigh), col = 'grey60') +
  geom_pointrange(size = 1, col = 'grey60') +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(y = 1.3, label = signif_lab), size = 8) +
  labs(x = NULL, y = 'dN/dS') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'none')

dev.off()

# strand specific ---------------------------------------------------------

comb_mutsig_mut_region_ss <- comb_mutsig_mut_region %>% 
  filter(!strand %in% '-', !region %in% 'Nongenic') %>% 
  mutate(strand_abbr = ifelse(strand %in% 'transcribed', 'T', 'U'), 
         region_ss = str_c(region, strand_abbr, sep = ':'))

# comb_mutsig_mut_region_ss %>% filter(Project %in% 'PCAWG') %>% count(region_ss)

comb_mutsig_ss_tmb <- comb_mutsig_mut_region_ss %>% 
  count(Project, ICGC_abbr_top, Sample, Attri_Signature, region, region_ss) %>% 
  inner_join(region_length, by = 'region') %>% 
  mutate(tmb = n/length)

pcawg_mutsig_ss_tmb <- comb_mutsig_ss_tmb %>% filter(Project %in% 'PCAWG')

pcawg_mutsig_ss_tmb_pow <- pcawg_mutsig_ss_tmb %>% 
  semi_join(sam_mutsig_power, c('Sample', 'Attri_Signature')) %>% 
  group_by(Project, region_ss, Attri_Signature) %>% 
  filter(n() >= 30) %>% # region signature with sufficient statistical power (i.e. signataure within a certain region with sample count >= 30)
  ungroup()

# significance test

pcawg_mutsig_ss_tmb_test <- pcawg_mutsig_ss_tmb_pow %>% 
  select(-region, -(n:length)) %>% 
  separate(region_ss, into = c('region', 'type'), sep = ':') %>% 
  group_by(type, Attri_Signature) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(tmb ~ region, data = .)$p.value)) %>% 
  group_by(type) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH'), 
         q_label = case_when(
           q_val < 0.0001 ~ '****', 
           q_val >= 0.0001 & q_val < 0.001 ~ '***', 
           q_val >= 0.001 & q_val < 0.01 ~ '**', 
           q_val >= 0.01 & q_val < 0.05 ~ '*', 
           q_val >= 0.05 ~ 'n.s.'
         )) %>% 
  ungroup() %>% 
  select(-data)

# q_val
# SBS1,T: 0.00198, **; SBS1,U: 0.000129, ***; SBS13,T: 0.00758, **
# others < 0.0001 ****

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_mutation_analysis_revise/mutsig_region_tmb_ss.pdf', width = 18, height = 6)

pcawg_mutsig_ss_tmb_pow %>% 
  mutate(Attri_Signature = factor(Attri_Signature, levels = mutsig_order), 
         region_ss = factor(region_ss, levels = c('Exon:T', 'UTR:T', 'Intron:T', 'Exon:U', 'UTR:U', 'Intron:U'))) %>% 
  ggplot() +
  geom_boxplot(aes(Attri_Signature, tmb, fill = region_ss), outlier.color = 'grey60') +
  labs(x = NULL, fill = NULL, y = 'Tumor Mutation Burden (mut/Mb)') +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), label = c('0.01', '0.1', '1', '10', '100')) +
  scale_fill_manual(breaks = c('Exon:T', 'Exon:U', 'UTR:T', 'UTR:U', 'Intron:T', 'Intron:U'), values = c('#8DD3C7', '#009E80', '#FFFFB3', '#FFED6F', '#BEBADA', '#BC80BD')) +
  theme_bw() +
  theme(axis.text.x = element_text(vjust = -1))

dev.off()


# clustered mutation analysis -----------------------------------------------------

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/all_clustermut_res.Rdata')
# 378,337 clusters, 39 cancer types, 6295 samples

all_clustermut_res_pivot <- all_clustermut_res %>% 
  pivot_longer(starts_with('SBS'), names_to = 'Signature', values_to = 'count') %>% 
  filter(!is.na(count)) %>% 
  mutate(mut_cluster_id = str_c(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, sep = ':'), 
         frac = count/nMuts) %>% 
  rename(Sample = Tumor_Sample_Barcode)

sample_mut_clus_stat <- all_clustermut_res_pivot %>% 
  distinct(ICGC_abbr_top, Sample, mut_cluster_id) %>% 
  count(ICGC_abbr_top, Sample)

# remove Samples of which cluster count > 200

all_clustermut_res_filt <- all_clustermut_res_pivot %>% semi_join(filter(sample_mut_clus_stat, n <= 200), by = c('ICGC_abbr_top', 'Sample'))
# 78,818 clusters, 39 cancer types, 6122 samples

cluster_mutsig_sort <- all_clustermut_res_filt %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

clustermut_domi_mutsig <- all_clustermut_res_filt %>% filter(frac > 0.5) %>% select(Sample, mut_cluster_id, domi_sig = Signature, nMuts, count, frac)


mutsig_domi_cluster_sam_count <- clustermut_domi_mutsig %>% 
  distinct(domi_sig, Sample) %>% 
  count(domi_sig, name = 'domi_sam_count')
# 40 signatures

mutsig_domi_cluster_num <- clustermut_domi_mutsig %>% 
  count(domi_sig, Sample) %>% 
  group_by(domi_sig) %>% 
  summarise(mean_sc_num = mean(n)) %>% 
  ungroup()

mutsig_domi_cluster_stat <- mutsig_sam_count %>% 
  inner_join(mutsig_domi_cluster_sam_count, by = c('Signature' = 'domi_sig')) %>% 
  inner_join(mutsig_domi_cluster_num, by = c('Signature' = 'domi_sig')) %>% 
  mutate(domi_sam_frac = domi_sam_count/sam_count)

write_tsv(mutsig_domi_cluster_stat, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_domi_cluster_stat.tsv')

# only display sample count over 30

mutsig_domi_cluster_stat_rm <- mutsig_domi_cluster_stat %>% filter(sam_count >= 30)
# 32 signatures

# domi sam frac barplot

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure 2C/domi_clus_sam_frac_barplot.pdf', width = 13, height = 3)

mutsig_domi_cluster_stat_rm %>% 
  mutate(expression = map2(domi_sam_count, sam_count, ~bquote(over(.(.x), .(.y)))), 
         Signature = factor(Signature, levels = cluster_mutsig_sort)) %>% 
  ggplot(aes(x = Signature, domi_sam_frac, fill = mean_sc_num)) +
  geom_bar(stat = 'identity', color = 'grey60', lwd = 0.1) +
  geom_text(aes(label = expression), parse = TRUE, vjust = -0.15, size = 2) +
  labs(x = NULL, y = 'Proportion', fill = '# cluster/sample') +
  scale_fill_gradient(low = "white", high = "#006CB1", guide = guide_colorbar(label = TRUE,
                                                                              draw.ulim = TRUE, 
                                                                              draw.llim = TRUE,
                                                                              frame.colour = 'grey60',
                                                                              ticks = TRUE, 
                                                                              barwidth = 1,
                                                                              barheight = 4)) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank())

dev.off()


# domi cluster mutsig sam pie plot

domi_clus_donut_dat <- mutsig_sam_dat %>% 
  semi_join(clustermut_domi_mutsig, by = c('Sample', 'Signature' = 'domi_sig')) %>% 
  semi_join(mutsig_domi_cluster_stat_rm, by = 'Signature') %>% # 32 signatures
  count(Signature, ICGC_abbr_top, name = 'sam_num') %>% 
  mutate(Signature = factor(Signature, levels = cluster_mutsig_sort)) %>% 
  group_by(Signature) %>% 
  mutate(cancer_type_frac = sam_num/sum(sam_num)) %>% 
  ungroup()

# lump low frequency cancer types for better display

domi_clus_donut_dat_lump <- domi_clus_donut_dat %>% 
  group_by(Signature) %>% 
  mutate(ICGC_abbr_top_lump = ifelse((sam_num >= 15 & cancer_type_frac >= 0.05) | cancer_type_frac > 0.3, ICGC_abbr_top, 'other')) %>% 
  ungroup()

# merge other

domi_clus_donut_dat_lump_mg <- domi_clus_donut_dat_lump %>% 
  group_by(Signature, ICGC_abbr_top_lump) %>% 
  summarise(sam_num = sum(sam_num), cancer_type_frac = sum(cancer_type_frac)) %>% 
  ungroup()

pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure 2C/mutsig_domi_clus_sam_count_donut.pdf', width = 16, height = 10)

domi_clus_donut_dat_lump_mg %>% 
  ggplot(aes(x = 2, cancer_type_frac, fill = ICGC_abbr_top_lump)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(breaks = rf_cancer_colors$cancer_type, values = rf_cancer_colors$color) +
  coord_polar(theta = 'y') +
  labs(x = NULL, y = NULL, fill = NULL) +
  xlim(0.5, 2.5) + # the most important part for producing donut plot, otherwise a pie plot will be returned
  facet_grid(.~Signature) +
  theme_void()

dev.off()

lumped_mutsig_ttype <- domi_clus_donut_dat_lump %>% select(-(sam_num:cancer_type_frac))


# mutsig clustered mutation contribution

mutsig_sam_clusmut_contri_frac <- all_clustermut_res_filt %>% 
  group_by(Sample, Signature) %>% 
  summarise(clus_mut_count = sum(count)) %>% 
  ungroup() %>% 
  right_join(mutsig_sam_mut_count, by = c('Sample', 'Signature')) %>% 
  group_by(Sample) %>% 
  filter(any(!is.na(clus_mut_count))) %>% 
  ungroup() %>% 
  mutate(clus_mut_count = ifelse(is.na(clus_mut_count), 0, clus_mut_count), 
         clus_mut_frac = clus_mut_count/mut_count)
# 22 signatures, 5732 samples, 38 cancer types


mutsig_sam_clusmut_contri_frac_plot_dat <- mutsig_sam_clusmut_contri_frac %>% 
  mutate(Signature = factor(Signature, levels = cluster_mutsig_sort)) %>% 
  group_by(Signature) %>% 
  mutate(extreme_outlier_upper = quantile(clus_mut_frac)[4] + (quantile(clus_mut_frac)[4] - quantile(clus_mut_frac)[2]) * 3, 
         extreme_outlier_lower = quantile(clus_mut_frac)[2] - (quantile(clus_mut_frac)[4] - quantile(clus_mut_frac)[2]) * 3) %>% 
  ungroup() %>% 
  filter(clus_mut_frac <= extreme_outlier_upper, 
         clus_mut_frac >= extreme_outlier_lower) # Attention!!! remove extreme outlier for better displaying, some significant asterisks by ggpubr may need to modify 

# lump low frequency cancer types for better display

mutsig_sam_clusmut_contri_frac_plot_dat_lump <- mutsig_sam_clusmut_contri_frac_plot_dat %>% 
  left_join(lumped_mutsig_ttype, by = c('Signature', 'ICGC_abbr_top'))
# 93 NA 


pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure 2D/mutsig_clusmut_frac.pdf', width = 8, height = 5)

mutsig_sam_clusmut_contri_frac_plot_dat_lump %>% 
  ggplot(aes(Signature, clus_mut_frac)) +
  geom_jitter(aes(col = ICGC_abbr_top_lump), alpha = 0.5, width = 0.25, size = 1, shape = '.') +
  geom_boxplot(outlier.color = NA, fill = NA) +
  geom_hline(yintercept = median(mutsig_sam_clusmut_contri_frac$clus_mut_frac), linetype = 2) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = TRUE) +
  labs(x = NULL, y = 'Proportion of clustered mutation', col = NULL) +
  guides(color = guide_legend(ncol = 2)) +
  scale_color_manual(breaks = rf_cancer_colors$cancer_type, values = rf_cancer_colors$color) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust = 0.5))

dev.off()


# test for clustered mut fraction --------------------------------------------


# SBS7_frac_dat <- mutsig_sam_clusmut_contri_frac %>% filter(Signature %in% c('SBS7a', 'SBS7b'))
# SBS7_frac_dat_test <- ggpubr::compare_means(clus_mut_frac ~ Signature,  data = SBS7_frac_dat, method = 'wilcox.test')
# 0.031(SBS7a) vs 0.065(SBS7b), p = 3.61e-62


mutsig_clusmut_total_test <- ggpubr::compare_means(clus_mut_frac ~ Signature, data = mutsig_sam_clusmut_contri_frac, method = 'kruskal.test')
# p < 2e-16

mutsig_clusmut_base_test_ori <- ggpubr::compare_means(clus_mut_frac ~ Signature,  data = mutsig_sam_clusmut_contri_frac, ref.group = '.all.', method = 'wilcox.test')
mutsig_clusmut_base_test <- mutsig_clusmut_base_test_ori %>% mutate(p.signif = ifelse(p.adj < 0.05, p.signif, 'ns'))
# save(mutsig_clusmut_base_test, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_clusmut_base_test.Rdata')


# kataegis analysis -----------------------------------------------------

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/all_kataegis_res.Rdata')
# 41,743 kataegis, 30 cancer types, 1283 samples

all_kataegis_res_pivot <- all_kataegis_res %>% 
  pivot_longer(starts_with('SBS'), names_to = 'Signature', values_to = 'count') %>% 
  filter(!is.na(count)) %>% 
  mutate(mut_cluster_id = str_c(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, sep = ':'), 
         frac = count/nMuts) %>% 
  rename(Sample = Tumor_Sample_Barcode)

# remove Samples of which kataegis count > 200

sample_mut_kataegis_stat <- all_kataegis_res_pivot %>% 
  distinct(ICGC_abbr_top, Sample, mut_cluster_id) %>% 
  count(ICGC_abbr_top, Sample)

all_kataegis_res_filt <- all_kataegis_res_pivot %>% semi_join(filter(sample_mut_kataegis_stat, n <= 200), by = c('ICGC_abbr_top', 'Sample'))
# 8,222 kataegis, 30 cancer types, 1273 samples

kataegis_domi_mutsig <- all_kataegis_res_filt %>% filter(frac > 0.5) %>% select(Sample, mut_cluster_id, domi_sig = Signature, nMuts, count, frac)

mutsig_domi_kataegis_sam_count <- kataegis_domi_mutsig %>% 
  distinct(domi_sig, Sample) %>% 
  count(domi_sig, name = 'domi_sam_count')
# 25 signatures

mutsig_domi_kataegis_num <- kataegis_domi_mutsig %>% 
  count(domi_sig, Sample) %>% 
  group_by(domi_sig) %>% 
  summarise(mean_kataegis_num = mean(n)) %>% 
  ungroup()

mutsig_domi_kataegis_stat <- mutsig_sam_count %>% 
  inner_join(mutsig_domi_kataegis_sam_count, by = c('Signature' = 'domi_sig')) %>% 
  inner_join(mutsig_domi_kataegis_num, by = c('Signature' = 'domi_sig')) %>% 
  mutate(domi_sam_frac = domi_sam_count/sam_count)

# write_tsv(mutsig_domi_kataegis_stat, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_domi_kataegis_stat.tsv')


# domi kataegis sam frac barplot

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/domi_kataegis_sam_frac_barplot.pdf', width = 10, height = 3)

mutsig_domi_kataegis_stat %>% 
  mutate(expression = map2(domi_sam_count, sam_count, ~bquote(over(.(.x), .(.y)))), 
         Signature = factor(Signature, levels = cluster_mutsig_sort)) %>% 
  ggplot(aes(x = Signature, domi_sam_frac, fill = mean_kataegis_num)) +
  geom_bar(stat = 'identity', color = 'grey60', lwd = 0.1) +
  geom_text(aes(label = expression), parse = TRUE, vjust = -0.15, size = 2) +
  labs(x = NULL, y = 'Proportion', fill = '# kataegis/sample') +
  scale_fill_gradient(low = "white", high = "#006CB1", guide = guide_colorbar(label = TRUE,
                                                                              draw.ulim = TRUE, 
                                                                              draw.llim = TRUE,
                                                                              frame.colour = 'grey60',
                                                                              ticks = TRUE, 
                                                                              barwidth = 1,
                                                                              barheight = 4)) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank())

dev.off()


# domi kataegis mutsig sam donut plot

domi_kata_donut_dat <- mutsig_sam_dat %>% 
  semi_join(kataegis_domi_mutsig, by = c('Sample', 'Signature' = 'domi_sig')) %>% 
  count(Signature, ICGC_abbr_top, name = 'sam_num') %>% 
  mutate(Signature = factor(Signature, levels = cluster_mutsig_sort)) %>% 
  group_by(Signature) %>% 
  mutate(cancer_type_frac = sam_num/sum(sam_num)) %>% 
  ungroup()

# lump low frequency cancer types for better display

domi_kata_donut_dat_lump <- domi_kata_donut_dat %>% 
  group_by(Signature) %>% 
  mutate(ICGC_abbr_top_lump = ifelse((sam_num >= 10 & cancer_type_frac >= 0.01) | cancer_type_frac > 0.3, ICGC_abbr_top, 'other')) %>% 
  ungroup()


# 20 cancer types

# merge other

domi_kata_donut_dat_lump_mg <- domi_kata_donut_dat_lump %>% 
  group_by(Signature, ICGC_abbr_top_lump) %>% 
  summarise(sam_num = sum(sam_num), cancer_type_frac = sum(cancer_type_frac)) %>% 
  ungroup()

# pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_domi_kata_sam_count_donut.pdf', width = 20, height = 12)
pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure S3/mutsig_domi_kata_sam_count_donut.pdf', width = 20, height = 12)

domi_kata_donut_dat_lump_mg %>% 
  ggplot(aes(x = 2, cancer_type_frac, fill = ICGC_abbr_top_lump)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(breaks = rf_cancer_colors$cancer_type, values = rf_cancer_colors$color) +
  coord_polar(theta = 'y') +
  labs(x = NULL, y = NULL, fill = NULL) +
  xlim(0.5, 2.5) + # the most important part for producing donut plot, otherwise a pie plot will be returned
  facet_grid(.~Signature) +
  theme_void()

dev.off()

lumped_mutsig_ttype_kata <- domi_kata_donut_dat_lump %>% select(-(sam_num:cancer_type_frac))

# mutsig kataegis mutation contribution

mutsig_sam_kataegis_mut_contri_frac <- all_kataegis_res_filt %>% 
  group_by(Sample, Signature) %>% 
  summarise(kataegis_mut_count = sum(count)) %>% 
  ungroup() %>% 
  right_join(mutsig_sam_mut_count, by = c('Sample', 'Signature')) %>% 
  group_by(Sample) %>% 
  filter(any(!is.na(kataegis_mut_count))) %>% 
  ungroup() %>% 
  mutate(kataegis_mut_count = ifelse(is.na(kataegis_mut_count), 0, kataegis_mut_count), 
         clus_mut_frac = kataegis_mut_count/mut_count)
# 21 signatures, 1264 samples, 29 cancer types
# save(mutsig_sam_kataegis_mut_contri_frac, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_sam_kataegis_mut_contri_frac.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_sam_kataegis_mut_contri_frac.Rdata')


mutsig_sam_kataegis_mut_contri_frac_plot_dat <- mutsig_sam_kataegis_mut_contri_frac %>% 
  mutate(Signature = factor(Signature, levels = cluster_mutsig_sort)) %>% 
  group_by(Signature) %>% 
  mutate(extreme_outlier_upper = quantile(clus_mut_frac)[4] + (quantile(clus_mut_frac)[4] - quantile(clus_mut_frac)[2]) * 3, 
         extreme_outlier_lower = quantile(clus_mut_frac)[2] - (quantile(clus_mut_frac)[4] - quantile(clus_mut_frac)[2]) * 3) %>% 
  ungroup() %>% 
  filter(clus_mut_frac <= extreme_outlier_upper, 
         clus_mut_frac >= extreme_outlier_lower) # Attention!!! remove extreme outlier for better displaying, some significant asterisks by ggpubr may need to modify 

# lump low frequency cancer types for better display

mutsig_sam_kataegis_mut_contri_frac_plot_dat_lump <- mutsig_sam_kataegis_mut_contri_frac_plot_dat %>% 
  left_join(lumped_mutsig_ttype_kata, by = c('Signature', 'ICGC_abbr_top'))


# pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_kataegis_mut_frac.pdf', width = 12, height = 4)
pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/Figure S3/mutsig_kataegis_mut_frac.pdf', width = 12, height = 4)

mutsig_sam_kataegis_mut_contri_frac_plot_dat_lump %>% 
  mutate(Signature = factor(Signature, levels = cluster_mutsig_sort)) %>% 
  ggplot(aes(Signature, clus_mut_frac)) +
  geom_jitter(aes(col = ICGC_abbr_top_lump), alpha = 0.5, width = 0.25, size = 1, shape = '.') +
  geom_boxplot(outlier.color = NA, fill = NA) +
  geom_hline(yintercept = median(mutsig_sam_kataegis_mut_contri_frac$clus_mut_frac), linetype = 2) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = TRUE) +
  labs(x = NULL, y = 'Proportion of kataegis mutation', col = NULL) +
  guides(color = guide_legend(ncol = 2)) +
  scale_color_manual(breaks = rf_cancer_colors$cancer_type, values = rf_cancer_colors$color) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust = 0.5))

dev.off()


# test for kataegis mut fraction --------------------------------------------


mutsig_kataegis_mut_total_test <- ggpubr::compare_means(clus_mut_frac ~ Signature, data = mutsig_sam_kataegis_mut_contri_frac, method = 'kruskal.test')
# p < 2e-16

mutsig_kataegis_mut_base_test_ori <- ggpubr::compare_means(clus_mut_frac ~ Signature,  data = mutsig_sam_kataegis_mut_contri_frac, ref.group = '.all.', method = 'wilcox.test')
mutsig_kataegis_mut_base_test <- mutsig_kataegis_mut_base_test_ori %>% mutate(p.signif = ifelse(p.adj < 0.05, p.signif, 'ns'))
# save(mutsig_kataegis_mut_base_test, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/mutsig_kataegis_mut_base_test.Rdata')

# select gene for lollipop plot -------------------------------------------

library(fuzzyjoin)

kataegis_dat <- all_kataegis_res_pivot %>% 
  distinct(Sample, mut_cluster_id, Chromosome, Start_Position, End_Position) %>% 
  mutate(Chromosome = str_c('chr', Chromosome), 
         Chromosome = case_when(
           Chromosome %in% 'chr23' ~ 'chrX', 
           Chromosome %in% 'chr24' ~ 'chrY', 
           TRUE ~ Chromosome
         ))

exon_mut_dat <- comb_mutsig_mut_dat %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'De_novo_Start_OutOfFrame', 'De_novo_Start_InFrame', 'Start_Codon_SNP', 'Nonstop_Mutation'))

select_sample_gene_dat <- exon_mut_dat %>% 
  count(Sample, Hugo_Symbol) %>% 
  filter(n > 6, grepl('SP', Sample)) %>% # include genes with over 6 mutations in a PCAWG sample
  arrange(desc(n))

select_mut_dat <- exon_mut_dat %>% semi_join(select_sample_gene_dat, by = c('Sample', 'Hugo_Symbol'))
select_kata_dat <- kataegis_dat %>% semi_join(select_mut_dat, by = c('Sample', 'Chromosome'))

select_mut_kata_dat <- select_mut_dat %>% 
  fuzzy_inner_join(select_kata_dat, by = c('Sample', 'Chromosome', 'Start_position' = 'Start_Position', 'Start_position' = 'End_Position'), match_fun = list(`==`, `==`, `>=`, `<=`))
# time-consuming

# save(select_mut_kata_dat, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/select_mut_kata_dat.Rdata')

# Bone-Osteosarc, SP116478, MGAT4C, SP116478:12:86372801:86388330 cluster

