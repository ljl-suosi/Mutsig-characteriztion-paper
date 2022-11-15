
library(tidyverse)
library(survival)
library(survminer)
library(ggh4x)
library(patchwork)
library(forestmodel) # for cox forest plot

comb_mutsig_we <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')

comb_mutsig_we_sur <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>%
  mutate(group = ifelse(Weight > median(Weight), 'High', 'Low'), group = factor(group, levels = c('Low', 'High'))) %>% 
  ungroup() %>% 
  inner_join(comb_cli_data %>% select(Sample:OS.time, TMB), by = 'Sample') %>% 
  filter(!is.na(OS), !is.na(OS.time))
# 63,790

# Individual cancer type overall survival ---------------------------------

# univariate cox analysis

unicox_res_list <- list()

for(i in unique(comb_mutsig_we_sur$ICGC_abbr_top)){
  
  # i <- 'Skin-Melanoma'
  
  sur_data <- comb_mutsig_we_sur %>% 
    filter(ICGC_abbr_top %in% i) %>% 
    group_by(Signature) %>% 
    select(Signature, group, OS:OS.time) %>% 
    nest() %>% 
    ungroup()
  
  # univariate cox (for each mutation signature)
  
  unicox_res <- sur_data %>% 
    mutate(cox_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ group, data = .)), NA))) %>% 
    filter(!is.na(cox_res)) %>% 
    mutate(cox_mat = map(cox_res, ~tibble(coef = .$coefficients[, 1], 
                                          HR = .$coefficients[, 2], 
                                          HR_confint_lower = .$conf.int[, 3], 
                                          HR_confint_upper = .$conf.int[, 4], 
                                          p_val = .$coefficients[, 5]))) %>% 
    select(-((data:cox_res))) %>% 
    unnest(cox_mat) %>% 
    # mutate(ICGC_abbr_top = i, q_val = p.adjust(p_val, method = 'BH')) %>% 
    mutate(ICGC_abbr_top = i) %>% 
    select(ICGC_abbr_top, everything())
  
  unicox_res_list[[i]] <- unicox_res
  
  print(i)
}

mutsig_unicox_tib <- unicox_res_list %>% bind_rows() %>% filter_all(~!is.na(.))
mutsig_unicox_tib_signif <- mutsig_unicox_tib %>% filter(p_val < 0.05) # 35 correlations, 17 cancer types


# write_tsv(mutsig_unicox_tib, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/mutsig_unicox_tib.tsv')
# write_tsv(mutsig_unicox_tib_signif, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/mutsig_unicox_tib_signif.tsv')

# mutsig_unicox_tib_signif <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/mutsig_unicox_tib_signif.tsv')

# barplot

sur_sig_sort <- mutsig_unicox_tib_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

barplot_complex_dat <- mutsig_unicox_tib %>% 
  semi_join(mutsig_unicox_tib_signif, by = 'ICGC_abbr_top') %>% 
  semi_join(mutsig_unicox_tib_signif, by = 'Signature') %>% 
  mutate(sign = ifelse(coef > 0, 1, -1), log_p_val = sign * (-log10(p_val)), 
         Signature = factor(Signature, levels = sur_sig_sort))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/individual_mutsig_survival_association_barplot.pdf', width = 16, height = 5)

barplot_complex_dat %>% 
  mutate(log_p_color = case_when(
    log_p_val < log10(0.05) ~ log10(0.05), 
    log_p_val > -log10(0.05) ~ -log10(0.05), 
    TRUE ~ log_p_val
  )) %>% 
  # ggplot(aes(Signature, log_p_val, fill = as.character(log_p_val))) +
  ggplot(aes(Signature, log_p_val, fill = log_p_color)) +
  geom_bar(stat = "identity", color = 'black') +
  labs(x = NULL, y = 'survival association (log10 p-value)') +
  # scale_fill_manual(breaks = uni_logp_val, values = custom_colors) +
  scale_fill_gradient2(low = '#5082AF', mid = 'white', high = '#C9372E') +
  theme_bw() +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_grid2(.~ICGC_abbr_top, space = 'free', scales = 'free', strip = strip_vanilla(clip = 'off')) +
  theme(panel.grid = element_blank(), 
        strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 60, hjust = 1))

dev.off()

# km plot


sincan_kmplot_fun <- function(cancer_type, signature){
    
  sin_data <- comb_mutsig_we_sur %>% 
    filter(ICGC_abbr_top %in% cancer_type, Signature %in% signature) %>% 
    select(ICGC_abbr_top, Signature, group, OS:OS.time)
    
  sin_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group, data = sin_data),
                         pval = TRUE,
                         title = str_c(signature, cancer_type, sep = ' in '),
                         xlab = 'Survival time (Days)',
                         legend.labs = c('Low', 'High'),
                         palette = c('#5082AF', '#C9372E'),
                         risk.table = TRUE,
                         ggtheme = theme_bw(),
                         data = sin_data)$plot + 
    labs(col = NULL) + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(1, 1),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA), 
          panel.grid = element_blank())
  
  return(sin_plot)
}

# all km-plots

all_cox_signif_km_plots <- mutsig_unicox_tib_signif %>% 
  select(ICGC_abbr_top, Signature) %>% 
  mutate(plot_list = map2(ICGC_abbr_top, Signature, ~sincan_kmplot_fun(.x, .y)))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/individual_all_km_plots.pdf', width = 20, height = 10)

wrap_plots(all_cox_signif_km_plots$plot_list) + plot_layout(nrow = 4)

dev.off()

# selected km-plots

selected_km_plots <- tibble(ICGC_abbr_top = c('Breast', 'Panc-AdenoCA', 'CNS-Medullo', 'ColoRect-AdenoCA', 'Liver-HCC', 'Liver-HCC', 'Bladder-TCC', 'Ovary-AdenoCA', 'Lymph-CLL'), 
                            Signature = c('SBS18', 'SBS18', 'SBS18', 'SBS17b', 'SBS1', 'SBS5', 'SBS5', 'SBS3', 'SBS9')) %>% 
  mutate(km_plots_list = map2(ICGC_abbr_top, Signature, ~sincan_kmplot_fun(.x, .y)))


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/selected_km_plots.pdf', width = 10, height = 10)

wrap_plots(selected_km_plots$km_plots_list) + plot_layout(nrow = 3)

dev.off()


# multivariate cox analysis

# function: multivariate cox analysis for a specific signature in a certian cancer type

getCan_mutsig_multicox <- function(cli_sur_data, which_cancer_type, which_signature, factors = c('group', 'age', 'gender', 'stage', 'TMB'), return_forest_plot = FALSE){
    
  sin_cancer_data <- cli_sur_data %>% filter(ICGC_abbr_top %in% which_cancer_type, Signature %in% which_signature)
  
  sur_formu <- as.formula(str_c("Surv(OS.time, OS) ~ ", str_c(factors, collapse = ' + ')))
  
  multicox_res_ori <- coxph(sur_formu, data = sin_cancer_data)
  multicox_res <- summary(multicox_res_ori)
  multicox_tib <- tibble(variable = names(multicox_res$coefficients[, 1]),
                         coef = multicox_res$coefficients[, 1],
                         HR = multicox_res$coefficients[, 2],
                         HR_confint_lower = multicox_res$conf.int[, 3],
                         HR_confint_upper = multicox_res$conf.int[, 4],
                         p_val = multicox_res$coefficients[, 5])
  
  ######## Attention!!!!!!
  # aaa <- list(ggforest(multicox_res_ori, data = sin_cancer_data)) # store the plot in a list preahead to prevent the produce of blank page in PDF file
  # pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/temp.pdf')
  # aaa[[1]]
  # dev.off()
  
  if(return_forest_plot){
    
    cox_forest_plot <- forest_model(coxph(sur_formu, data = sin_cancer_data)) + 
      labs(title = str_c(which_signature, ', ', which_cancer_type)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
    
    return(cox_forest_plot)
  }else{
    return(multicox_tib)
  }
  
}

# adjusting age, gender, stage

mutsig_multicox_tib_ags <- mutsig_unicox_tib_signif %>% 
  select(ICGC_abbr_top, Signature) %>% 
  filter(!ICGC_abbr_top %in% c('Myeloid-AML', 'Sarcoma', 'Uterus-AdenoCA', 'Breast', 'Ovary-AdenoCA', 'Prost-AdenoCA')) %>% 
  mutate(multi_cox_res = map2(ICGC_abbr_top, Signature, ~getCan_mutsig_multicox(comb_mutsig_we_sur, .x, .y, factors = c('group', 'age', 'gender', 'stage', 'TMB')))) %>% 
  unnest(multi_cox_res)

# adjusting age, stage
# Actually, for SBS18 in Breast, it is inappropriate to adjust stage because there are only 43(43/678) SBS18-High samples among which only 
# 5 samples have stage information

mutsig_multicox_tib_as <- mutsig_unicox_tib_signif %>% 
  select(ICGC_abbr_top, Signature) %>% 
  filter(ICGC_abbr_top %in% c('Breast', 'Ovary-AdenoCA', 'Prost-AdenoCA')) %>% 
  mutate(multi_cox_res = map2(ICGC_abbr_top, Signature, ~getCan_mutsig_multicox(comb_mutsig_we_sur, .x, .y, factors = c('group', 'age', 'stage', 'TMB')))) %>% 
  unnest(multi_cox_res)

# adjusting age, gender

mutsig_multicox_tib_agen <- mutsig_unicox_tib_signif %>% 
  select(ICGC_abbr_top, Signature) %>% 
  filter(ICGC_abbr_top %in% c('Myeloid-AML', 'Sarcoma')) %>% 
  mutate(multi_cox_res = map2(ICGC_abbr_top, Signature, ~getCan_mutsig_multicox(comb_mutsig_we_sur, .x, .y, factors = c('group', 'age', 'gender', 'TMB')))) %>% 
  unnest(multi_cox_res)

# adjusting age, grade

mutsig_multicox_tib_agd <- mutsig_unicox_tib_signif %>% 
  select(ICGC_abbr_top, Signature) %>% 
  filter(ICGC_abbr_top %in% 'Uterus-AdenoCA') %>% 
  mutate(multi_cox_res = map2(ICGC_abbr_top, Signature, ~getCan_mutsig_multicox(comb_mutsig_we_sur, .x, .y, factors = c('group', 'age', 'grade', 'TMB')))) %>% 
  unnest(multi_cox_res)

mutsig_multicox_tib_all <- list(mutsig_multicox_tib_ags, mutsig_multicox_tib_as, mutsig_multicox_tib_agen, mutsig_multicox_tib_agd) %>% reduce(bind_rows)

mutsig_multicox_signif_signature <- mutsig_multicox_tib_all %>% 
  filter(str_detect(variable, 'group'), p_val < 0.05) %>% 
  select(-variable)
# 16 significant survival associations

mutsig_multicox_signif_signature_all <- mutsig_multicox_tib_all %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  semi_join(mutsig_multicox_signif_signature, by = c('ICGC_abbr_top', 'Signature')) %>% 
  ungroup()

# cox forest plot

selected_ggforest_plots <- tibble(ICGC_abbr_top = c('Panc-AdenoCA', 'Liver-HCC', 'CNS-Medullo', 'Lymph-CLL'), 
                                  Signature = c('SBS18', 'SBS1', 'SBS18', 'SBS9')) %>% 
  mutate(ggforest_plots_list = map2(ICGC_abbr_top, Signature, ~getCan_mutsig_multicox(comb_mutsig_we_sur, .x, .y, return_forest_plot = TRUE)))

# Attention!!! In linux environment, forest_model plot saved to a object will be changed!!!

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/selected_mutsig_cox_forest_plots.pdf', width = 12, height = 12)

wrap_plots(selected_ggforest_plots$ggforest_plots_list) + plot_layout(nrow = 2)

dev.off()


# Pan-carcinoma overall survival ------------------------------------------

pan_sur_data <- comb_mutsig_we %>% 
  inner_join(comb_cli_data %>% select(Sample:OS.time), by = 'Sample') %>% 
  group_by(Signature) %>% 
  # filter(mean(Weight > 0) >= 0.05) %>%
  mutate(group = ifelse(Weight > median(Weight), 'High', 'Low'), group = factor(group, levels = c('Low', 'High'))) %>% 
  filter(!is.na(OS), !is.na(OS.time)) %>% 
  select(Signature, ICGC_abbr_top, group, OS:OS.time) %>% 
  nest() %>% 
  ungroup()

pan_unicox_res <- pan_sur_data %>% 
  mutate(cox_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ group, data = .)), NA))) %>% 
  filter(!is.na(cox_res)) %>% 
  mutate(cox_mat = map(cox_res, ~tibble(coef = .$coefficients[, 1], 
                                        HR = .$coefficients[, 2], 
                                        HR_confint_lower = .$conf.int[, 3], 
                                        HR_confint_upper = .$conf.int[, 4], 
                                        p_val = .$coefficients[, 5]))) %>% 
  select(-((data:cox_res))) %>% 
  unnest(cox_mat) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH'))

pan_unicox_res_signif <- pan_unicox_res %>% filter(p_val < 0.05) # 24

# barplot

pan_sur_plot_data <- pan_unicox_res %>% 
  filter(!is.na(p_val)) %>% 
  mutate(sign = ifelse(coef > 0, 1, -1), 
         log_p_val = sign * (-log10(p_val)), 
         Signature = factor(Signature, levels = str_sort(Signature, numeric = TRUE)), 
         log_p_color = case_when(
           log_p_val < log10(0.05) ~ log10(0.05), 
           log_p_val > -log10(0.05) ~ -log10(0.05), 
           TRUE ~ log_p_val
         ))

# pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/pancan_mutsig_survival_association.pdf', height = 8, width = 12)

pan_sur_plot <- pan_sur_plot_data %>% 
  ggplot(aes(Signature, log_p_val, fill = log_p_color)) +
  geom_bar(stat = "identity", color = 'black') +
  labs(x = NULL, y = 'survival association (log10 p-value)', title = 'Overall survival association of Signature contribution') +
  # scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  scale_fill_gradient2(low = '#5082AF', mid = 'white', high = '#C9372E') +
  coord_cartesian(ylim = c(-10, 20)) +
  theme_bw() +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1))

# dev.off()


# Individual cancer type overall survival(driver) ------------------------------------------

comb_mutsig_clo_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_alt_data.tsv')
comb_driver_gene <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')

mutsig_clo_driver_data <- comb_mutsig_clo_data %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site')) %>% 
  filter_all(~!is.na(.)) %>% 
  semi_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top'))
# 20,872

sam_mutsig_driver <- mutsig_clo_driver_data %>% distinct(Sample, Attri_Signature) %>% mutate(driver_status = 'present')
# sam_mutsig_driver_clo <- mutsig_clo_driver_data %>% filter(clonality %in% 'clonal') %>% distinct(Sample, Attri_Signature) %>% mutate(driver_status = 'present')
# sam_mutsig_driver_sub <- mutsig_clo_driver_data %>% filter(clonality %in% 'subclonal') %>% distinct(Sample, Attri_Signature) %>% mutate(driver_status = 'present')


# cox analysis

individual_driver_sur_data <- comb_mutsig_we_sur %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% # the frequency of signature present require 5% above
  ungroup() %>% 
  distinct(ICGC_abbr_top, Sample, Signature, OS, OS.time) %>% 
  left_join(sam_mutsig_driver, by = c('Sample', 'Signature' = 'Attri_Signature')) %>% 
  mutate(driver_status = ifelse(is.na(driver_status), 'absent', driver_status)) %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(driver_status == 'present') >= 0.05) %>% # the frequency of signature driver present require 5% above
  select(ICGC_abbr_top, Signature, group = driver_status, OS:OS.time) %>% 
  nest() %>% 
  ungroup()


individual_driver_unicox_res <- individual_driver_sur_data %>% 
  mutate(cox_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ group, data = .)), NA))) %>% 
  filter(!is.na(cox_res)) %>% 
  mutate(cox_mat = map(cox_res, ~tibble(coef = .$coefficients[, 1], 
                                        HR = .$coefficients[, 2], 
                                        HR_confint_lower = .$conf.int[, 3], 
                                        HR_confint_upper = .$conf.int[, 4], 
                                        p_val = .$coefficients[, 5]))) %>% 
  select(-((data:cox_res))) %>% 
  unnest(cox_mat)

# write_tsv(individual_driver_unicox_res, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/individual_driver_unicox_res.tsv')

individual_driver_unicox_signif <- individual_driver_unicox_res %>% filter(p_val < 0.05) # 15
# write_tsv(individual_driver_unicox_signif, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/individual_driver_unicox_signif.tsv')


# barplot

sur_sig_sort_dri <- individual_driver_unicox_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

barplot_complex_dri_dat <- individual_driver_unicox_res %>% 
  semi_join(individual_driver_unicox_signif, by = 'ICGC_abbr_top') %>% 
  semi_join(individual_driver_unicox_signif, by = 'Signature') %>% 
  mutate(sign = ifelse(coef > 0, 1, -1), log_p_val = sign * (-log10(p_val)), Signature = factor(Signature, levels = sur_sig_sort_dri))

driver_sur_plot <- barplot_complex_dri_dat %>% 
  mutate(log_p_color = case_when(
    log_p_val < log10(0.05) ~ log10(0.05), 
    log_p_val > -log10(0.05) ~ -log10(0.05), 
    TRUE ~ log_p_val
  )) %>% 
  ggplot(aes(Signature, log_p_val, fill = log_p_color)) +
  geom_bar(stat = "identity", color = 'black') +
  labs(x = NULL, y = 'survival association (log10 p-value)', title = 'Overall survival association of Signature driver contribution') +
  scale_fill_gradient2(low = '#5082AF', mid = 'white', high = '#C9372E') +
  theme_bw() +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_grid2(.~ICGC_abbr_top, space = 'free', scales = 'free', strip = strip_vanilla(clip = 'off')) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank(), 
        # strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))

# Pan-carcinoma overall survival(contribution clonality) ------------------------------------------

comb_sam_clo_stat <- comb_mutsig_clo_data %>% 
  filter_at(c('Attri_Signature', 'clonality'), ~!is.na(.)) %>% 
  count(ICGC_abbr_top, Sample, Attri_Signature, clonality) %>% 
  pivot_wider(names_from = clonality, values_from = n, values_fill = 0) %>% 
  group_by(Sample) %>% 
  mutate(clonal_weight = clonal/sum(clonal), subclonal_weight = subclonal/sum(subclonal)) %>% 
  ungroup() %>% 
  group_by(ICGC_abbr_top, Attri_Signature) %>% 
  mutate(clonal_group = ifelse(clonal_weight > median(clonal_weight), 'High', 'Low'), 
         clonal_group = factor(clonal_group, levels = c('Low', 'High')), 
         subclonal_group = ifelse(subclonal_weight > median(subclonal_weight), 'High', 'Low'), 
         subclonal_group = factor(subclonal_group, levels = c('Low', 'High'))) %>% 
  ungroup()

comb_mutsig_clo_sur <- comb_sam_clo_stat %>% 
  select(Sample, Attri_Signature, clonal_group, subclonal_group) %>% 
  inner_join(comb_mutsig_we_sur, by = c('Sample', 'Attri_Signature' = 'Signature')) %>% 
  select(-(Exposure:group))
  
# clonal mutations contributed by signature

# cox analysis

clonal_unicox_res <- comb_mutsig_clo_sur %>% 
  select(ICGC_abbr_top, Signature = Attri_Signature, group = clonal_group, OS:OS.time) %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(cox_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ group, data = .)), NA))) %>% 
  filter(!is.na(cox_res)) %>% 
  mutate(cox_mat = map(cox_res, ~tibble(coef = .$coefficients[, 1], 
                                        HR = .$coefficients[, 2], 
                                        HR_confint_lower = .$conf.int[, 3], 
                                        HR_confint_upper = .$conf.int[, 4], 
                                        p_val = .$coefficients[, 5]))) %>% 
  select(-((data:cox_res))) %>% 
  unnest(cox_mat) %>% 
  filter_all(~!is.na(.))

clonal_unicox_signif <- clonal_unicox_res %>% filter(p_val < 0.05) # 13

# boxplot

sur_sig_sort_clo <- clonal_unicox_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

barplot_complex_clo_dat <- clonal_unicox_res %>% 
  semi_join(clonal_unicox_signif, by = 'ICGC_abbr_top') %>% 
  semi_join(clonal_unicox_signif, by = 'Signature') %>% 
  mutate(sign = ifelse(coef > 0, 1, -1), log_p_val = sign * (-log10(p_val)), Signature = factor(Signature, levels = sur_sig_sort_clo))

# pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/individual_mutsig_clo_survival_association_barplot.pdf', width = 14, height = 6)

clo_sur_plot <- barplot_complex_clo_dat %>% 
  mutate(log_p_color = case_when(
    log_p_val < log10(0.05) ~ log10(0.05), 
    log_p_val > -log10(0.05) ~ -log10(0.05), 
    TRUE ~ log_p_val
  )) %>% 
  ggplot(aes(Signature, log_p_val, fill = log_p_color)) +
  geom_bar(stat = "identity", color = 'black') +
  labs(x = NULL, y = 'survival association (log10 p-value)', title = 'Overall survival association of Signature clonal contribution') +
  scale_fill_gradient2(low = '#5082AF', mid = 'white', high = '#C9372E') +
  theme_bw() +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_grid2(.~ICGC_abbr_top, space = 'free', scales = 'free', strip = strip_vanilla(clip = 'off')) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank(), 
        # strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))

# dev.off()

# subclonal mutations contributed by signature

# cox analysis

subclonal_unicox_res <- comb_mutsig_clo_sur %>% 
  select(ICGC_abbr_top, Signature = Attri_Signature, group = subclonal_group, OS:OS.time) %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(cox_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ group, data = .)), NA))) %>% 
  filter(!is.na(cox_res)) %>% 
  mutate(cox_mat = map(cox_res, ~tibble(coef = .$coefficients[, 1], 
                                        HR = .$coefficients[, 2], 
                                        HR_confint_lower = .$conf.int[, 3], 
                                        HR_confint_upper = .$conf.int[, 4], 
                                        p_val = .$coefficients[, 5]))) %>% 
  select(-((data:cox_res))) %>% 
  unnest(cox_mat) %>% 
  filter_all(~!is.na(.))

subclonal_unicox_signif <- subclonal_unicox_res %>% filter(p_val < 0.05) # 8

# barplot

sur_sig_sort_sub <- subclonal_unicox_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

barplot_complex_sub_dat <- subclonal_unicox_res %>% 
  semi_join(subclonal_unicox_signif, by = 'ICGC_abbr_top') %>% 
  semi_join(subclonal_unicox_signif, by = 'Signature') %>% 
  mutate(sign = ifelse(coef > 0, 1, -1), log_p_val = sign * (-log10(p_val)), Signature = factor(Signature, levels = sur_sig_sort_sub))

sub_sur_plot <- barplot_complex_sub_dat %>% 
  mutate(log_p_color = case_when(
    log_p_val < log10(0.05) ~ log10(0.05), 
    log_p_val > -log10(0.05) ~ -log10(0.05), 
    TRUE ~ log_p_val
  )) %>% 
  ggplot(aes(Signature, log_p_val, fill = log_p_color)) +
  geom_bar(stat = "identity", color = 'black') +
  labs(x = NULL, y = 'survival association (log10 p-value)', title = 'Overall survival association of Signature subclonal contribution') +
  scale_fill_gradient2(low = '#5082AF', mid = 'white', high = '#C9372E') +
  theme_bw() +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  facet_grid2(.~ICGC_abbr_top, space = 'free', scales = 'free', strip = strip_vanilla(clip = 'off')) +
  coord_cartesian(ylim = c(-3.5, 3.5)) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank(), 
        # strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))

# supplementary figures 

layout <- "
AAAAAA
AAAAAA
BBBCCC
BBBDDD
"

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/prognostic_association/combined_supp_plots.pdf', width = 18, height = 10)

pan_sur_plot + driver_sur_plot + clo_sur_plot + sub_sur_plot + plot_layout(design = layout)

dev.off()

cli_signif_cat_stat <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/cli_signif_cat_stat.tsv')


TMB_stat <- cli_signif_cat_stat %>% 
  filter(Clinical_Variable %in% 'TMB', Signature %in% c('SBS1', 'SBS5', 'SBS2', 'SBS13', 'SBS40')) %>% 
  select(ICGC_abbr_top, Signature, variable_value, median_Weight) %>% 
  pivot_wider(names_from = variable_value, values_from = median_Weight) %>% 
  mutate(`High > Low` = High_TMB > Low_TMB)

# TMB_stat %>% count(Signature, `High > Low`)
# # A tibble: 10 × 3
# Signature `High > Low`     n
# <chr>     <lgl>        <int>
# SBS1      FALSE           16
# SBS1      TRUE             2
# SBS13     FALSE            4
# SBS13     TRUE             5
# SBS2      FALSE            4
# SBS2      TRUE             5
# SBS40     FALSE            6
# SBS40     TRUE             1
# SBS5      FALSE           14
# SBS5      TRUE             4
