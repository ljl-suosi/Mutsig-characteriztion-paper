
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggh4x)

comb_mutsig_clo_data <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_data.tsv')
comb_mutsig_we <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')
cancer_mutsig <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/cancer_mutsig.tsv') # mutsig freq >= 0.05

all_mutsig_level <- comb_mutsig_we %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

# signature fitting

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/clonality_sig_fit_strict.Rdata')

# test

cancer_sam_count_prev <- clonality_sig_fit_strict %>% distinct(ICGC_abbr_top, Sample) %>% count(ICGC_abbr_top, name = 'total_n') %>% filter(total_n > 10)

prev_cancer_mutsig <- clonality_sig_fit_strict %>% 
  filter(!Signature %in% 'Unknown', Exposure >= 1) %>% 
  distinct(ICGC_abbr_top, Sample, Signature) %>% 
  count(ICGC_abbr_top, Signature) %>% 
  inner_join(cancer_sam_count_prev, by = 'ICGC_abbr_top') %>% 
  filter(n/total_n >= 0.05) %>% 
  select(ICGC_abbr_top, Signature)

mutsig_clo_fit_strict_test <- clonality_sig_fit_strict %>% 
  semi_join(prev_cancer_mutsig, by = c('ICGC_abbr_top', 'Signature')) %>% 
  select(-Exposure) %>% 
  pivot_wider(names_from = clonality, values_from = Weights, values_fill = 0) %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~wilcox.test(.$subclonal, .$clonal, paired = TRUE)$p.value, NA))) %>% # paired test
  select(-data) %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
  ungroup()
# 41 signatures, 27 cancer types

mutsig_clo_fit_strict_test_signif <- mutsig_clo_fit_strict_test %>% filter(q_val < 0.05)
# 30 signatures, 26 cancer types

# write_tsv(mutsig_clo_fit_strict_test_signif, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_clo_fit_strict_test_signif.tsv')

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')
mutsig_cli_relate_signif <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clinical_genomic_association/mutsig_cli_relate_signif.tsv')
comb_mutsig_we <- vroom::vroom('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

# Figure S4A ---------------------------------------------------------------
# lineplot
# Sigature contribution change across stages or grades

stage_grade_signif <- mutsig_cli_relate_signif %>% filter(Clinical_Variable %in% c('stage', 'grade'))

# 10 cancer types, 13 cancer stage or grade associations

# line plots between stages and grades

comb_mutsig_we_filter <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% 
  ungroup()

stage_sort <- comb_cli_data %>% distinct(stage) %>% pull(stage) %>% str_sort(numeric = TRUE)
grade_sort <- comb_cli_data %>% distinct(grade) %>% pull(grade) %>% str_sort(numeric = TRUE)

comb_mutsig_sg <- comb_mutsig_we_filter %>% 
  left_join(select(comb_cli_data, Sample, stage, grade), by = 'Sample') %>% 
  semi_join(stage_grade_signif, by = 'ICGC_abbr_top') %>% 
  mutate(stage = factor(stage, levels = stage_sort), 
         grade = factor(grade, levels = grade_sort))

forggline_plot <- function(cancer_type, stage_var){
  
  # cancer_type <- 'Panc-AdenoCA'
  # stage_var <- 'stage'
  
  # cancer_type <- 'Panc-AdenoCA'
  # cancer_type <- 'Breast'
  # stage_var <- 'stage'
  
  sin_data <- comb_mutsig_sg %>% 
    filter(ICGC_abbr_top %in% cancer_type) %>% 
    select(all_of(stage_var), Signature, Weight) %>% 
    filter_all(~!is.na(.)) %>% 
    set_names(c('stage_var', 'Signature', 'Weight'))
  
  sin_signif_dat <- stage_grade_signif %>% 
    filter(ICGC_abbr_top %in% cancer_type, Clinical_Variable %in% stage_var) %>% 
    mutate(label = case_when(
      p_val < 0.0001 ~ '****', 
      p_val >= 0.0001 & p_val < 0.001 ~ '***', 
      p_val >= 0.001 & p_val < 0.01 ~ '**', 
      p_val >= 0.01 & p_val < 0.05 ~ '*', 
    ))

  sin_data_plot_dat <- sin_data %>% 
    group_by(stage_var, Signature) %>% 
    summarise(mean_weight = mean(Weight), se_weight = sd(Weight)/sqrt(length(Weight))) %>% # Standard Error
    ungroup() %>% 
    mutate(ymin_weight = mean_weight - se_weight, 
           ymin_weight = ifelse(ymin_weight < 0, 0, ymin_weight), 
           ymax_weight = mean_weight + se_weight)
  
  sin_signif_label <- sin_data_plot_dat %>% 
    arrange(stage_var) %>% 
    distinct(Signature, .keep_all = TRUE) %>% 
    left_join(sin_signif_dat, by = 'Signature')
  
  sin_signature_label <- sin_data_plot_dat %>% 
    arrange(desc(stage_var)) %>% 
    distinct(Signature, .keep_all = TRUE)
  
  stage_number <- n_distinct(sin_data_plot_dat$stage_var)
  
  sin_data_plot_dat %>% 
    ggplot(aes(stage_var, mean_weight, group = Signature, col = Signature)) +
    geom_line(show.legend = FALSE) +
    geom_point() +
    labs(x = stage_var, col = NULL, y = 'Signature contribution', title = cancer_type) +
    geom_pointrange(aes(ymin = ymin_weight, ymax = ymax_weight)) +
    geom_text_repel(aes(label = label), data = sin_signif_label, xlim = c(NA, 0.8), direction = 'y', segment.alpha = 0.5, segment.linetype = "dotted") +
    geom_text_repel(aes(label = Signature), data = sin_signature_label, size = 2.5, xlim = c(stage_number+0.1, NA), direction = 'y', segment.alpha = 0.5, segment.linetype = "dotted") +
    ggsci::scale_color_futurama() +
    theme_classic() +
    theme(legend.position = 'none', 
          plot.title = element_text(hjust = 0.5))
    
}

# if one cancer has both grade and stage, select only one to display

selected_sg_plots <- tibble(ICGC_abbr_top = c('Liver-HCC', 'Breast', 'ColoRect-AdenoCA', 'Eso-AdenoCA', 'Lymph-CLL', 'Panc-AdenoCA', 'Panc-Endocrine', 'Prost-AdenoCA', 'Stomach-AdenoCA', 'Uterus-AdenoCA'), 
                            Clinical_Variable = c('grade', 'grade', 'stage', 'stage', 'stage', 'stage', 'stage', 'stage', 'stage', 'grade')) %>% 
  mutate(plot_list = map2(ICGC_abbr_top, Clinical_Variable, ~forggline_plot(.x, .y)))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_stage_grade_errorbar_selected.pdf', width = 22, height = 8)

wrap_plots(selected_sg_plots$plot_list) + plot_layout(nrow = 2)

dev.off()

# Figure 4A ---------------------------------------------------------------
# point plot
# Fold change of signature contribution between subclonal and clonal

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/clonality_sig_fit_strict.Rdata')
mutsig_clo_fit_strict_test_signif <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_clo_fit_strict_test_signif.tsv')

cancer_sam_count_prev <- clonality_sig_fit_strict %>% distinct(ICGC_abbr_top, Sample) %>% count(ICGC_abbr_top, name = 'total_n') %>% filter(total_n > 10)

prev_cancer_mutsig <- clonality_sig_fit_strict %>% 
  filter(!Signature %in% 'Unknown', Exposure >= 1) %>% 
  distinct(ICGC_abbr_top, Sample, Signature) %>% 
  count(ICGC_abbr_top, Signature) %>% 
  inner_join(cancer_sam_count_prev, by = 'ICGC_abbr_top') %>% 
  filter(n/total_n >= 0.05) %>% 
  select(ICGC_abbr_top, Signature)


mutsig_clo_fc_strict <- clonality_sig_fit_strict %>% 
  semi_join(prev_cancer_mutsig, by = c('ICGC_abbr_top', 'Signature')) %>% 
  group_by(Signature, ICGC_abbr_top, clonality) %>% 
  summarise(mean_weight = mean(Weights)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = clonality, values_from = mean_weight, values_fill = 0) %>% 
  mutate(mutsig_FC = subclonal/clonal)

# mutsig_clo_fit_strict_test_signif_fc <- mutsig_clo_fit_strict_test_signif %>% left_join(mutsig_clo_fc_strict, by = c('Signature', 'ICGC_abbr_top'))
# write_tsv(mutsig_clo_fit_strict_test_signif_fc, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_clo_fit_strict_test_signif_fc.tsv')


fc_sig_sort <- mutsig_clo_fc_strict %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

fc_cancer_colors <- tribble(
  ~cancer_type, ~color, 
  'Biliary-AdenoCA', '#BAA133', 
  'Bladder-TCC', '#F8D0D6', 
  'Breast', '#D92F86', 
  'Cervix', '#ECB065', 
  'ColoRect-AdenoCA', '#9DD6F0', 
  'CNS-Medullo', '#374F99', 
  'Eso-AdenoCA', '#0E78AA', 
  'CNS-GBM', '#A74C93', 
  'Head-SCC', '#93C9A3', 
  'Kidney-RCC', '#EEAAAE', 
  'Lymph-BNHL', '#724C29', 
  'CNS-PiloAstro', '#CC98BF', 
  'Liver-HCC', '#C7C9D7', 
  'Lung-AdenoCA', '#CEC0DC', 
  'Lung-SCC', '#997FB5', 
  'Lymph-CLL', '#502D80', 
  'Ovary-AdenoCA', '#CF7A29', 
  'Panc-AdenoCA', '#6B789B', 
  'Panc-Endocrine', '#DCBE29', 
  'Prost-AdenoCA', '#7A1A1D', 
  'Myeloid-AML', '#D8EEFB', 
  'Myeloid-MDS/MPN', '#0FA094', 
  'Skin-Melanoma', '#B3CB47', 
  'Stomach-AdenoCA', '#2FA0D2', 
  'Uterus-AdenoCA', '#F7DEC5', 
  'Bone-Osteosarc', '#EE9028', 
  'Bone-Other', '#0C8D46'
)

mutsig_fc_plot_dat <- mutsig_clo_fc_strict %>% 
  left_join(mutsig_clo_fit_strict_test_signif, by = c('Signature', 'ICGC_abbr_top')) %>% 
  group_by(Signature) %>% 
  filter(any(!is.na(p_val))) %>% # only display signatures with more than one significant clonality contribution change result
  ungroup() %>% 
  mutate(Signature = factor(Signature, levels = fc_sig_sort), 
         p_size_level = case_when(
           is.na(p_val) ~ 0.5, 
           p_val < 0.05 & p_val > 0.01 ~ 1, 
           p_val < 0.01 & p_val > 0.001 ~ 1.25, 
           p_val < 0.001 & p_val > 0.0001 ~ 1.5, 
           p_val < 0.0001 ~ 1.75
  ))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_clo_fc_pointplot.pdf', width = 18, height = 6)

mutsig_fc_plot_dat %>% 
  ggplot(aes(Signature, mutsig_FC, col = ICGC_abbr_top, size = p_size_level)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.3)) + 
  scale_color_manual(breaks = fc_cancer_colors$cancer_type, values = fc_cancer_colors$color) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_classic() +
  labs(x = NULL, y = 'Fold Change', col = NULL, size = NULL) +
  guides(col = guide_legend(nrow = 4), size = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = c(0.3, 0.8), 
        legend.background = element_blank())

dev.off()


# Figure 4B ---------------------------------------------------------------
# lineplot
# Signature contribution change between clonal and subclonal

signature_colors <- tribble(
  ~Signature, ~color, 
  'SBS1', '#8F6504', 
  'SBS5', '#927C64', 
  'SBS2', '#F17E27', 
  'SBS13', '#EB9B31', 
  'SBS40', '#FFB50E', 
  'SBS15', '#CC6677', 
  'SBS18', '#C41707', 
  'SBS3', '#0373C6', 
  'SBS45', '#B10DA1', 
  'SBS4', '#C1A72F', 
  'SBS10b', '#FF5FB3', 
  'SBS17a', '#70CBEF', 
  'SBS17b', '#91E4A6', 
  'SBS30', '#525975', 
  'SBS7a', '#9C8CC4', 
  'SBS49', '#ED2891', 
  'SBS6', '#B2509E', 
  'SBS7b', '#D49DC7', 
  'SBS7c', '#DEA0FD', 
  'SBS8', '#90AD1C', 
  'SBS10a', '#E8C51D', 
  'SBS19', '#F9ED32', 
  'SBS36', '#104A7F', #
  'SBS7d', 'grey60', 
  'SBS9', 'grey60', 
  'SBS12', '#16FF32', 
  'SBS16', 'grey60', 
  'SBS29', 'grey60', 
  'SBS32', 'grey60', 
  'SBS37', 'grey60', 
  'SBS38', 'grey60', 
  'SBS39', 'grey60', 
  'SBS44', 'grey60', 
  'SBS58', 'grey60', 
  'SBS59', 'grey60', #
  'SBS11', 'grey60', 
  'SBS14', 'grey60', 
  'SBS20', 'grey60', 
  'SBS22', 'grey60', 
  'SBS23', 'grey60', 
  'SBS27', 'grey60', 
  'SBS28', 'grey60', 
  'SBS33', 'grey60', 
  'SBS34', 'grey60', 
  'SBS41', 'grey60', 
  'SBS42', 'grey60', 
  'SBS48', 'grey60', 
  'SBS52', 'grey60', 
  'SBS56', 'grey60', 
  'SBS57', 'grey60', 
  'SBS60', 'grey60'
)

# 'SBS2', '#FF6444', 'darkred'
# 'SBS2', '#7D3526'
# 'SBS13', '#614527'


mutsig_clo_fit_plot_dat <- clonality_sig_fit_strict %>% filter(!Signature %in% 'Unknown')

forggline_clo_fit_plot <- function(cancer_type){
  
  # cancer_type <- 'Bladder-TCC'
  
  sin_data <- mutsig_clo_fit_plot_dat %>% filter(ICGC_abbr_top %in% cancer_type)
  
  sin_signif_dat <- mutsig_clo_fit_strict_test_signif %>% 
    filter(ICGC_abbr_top %in% cancer_type) %>% 
    mutate(label = case_when(
      p_val < 0.0001 ~ '****', 
      p_val >= 0.0001 & p_val < 0.001 ~ '***', 
      p_val >= 0.001 & p_val < 0.01 ~ '**', 
      p_val >= 0.01 & p_val < 0.05 ~ '*', 
    ))
  
  sin_data_plot_dat <- sin_data %>% 
    group_by(clonality, Signature) %>% 
    summarise(mean_weight = mean(Weights), se_weight = sd(Weights)/sqrt(length(Weights))) %>% # Standard Error
    ungroup() %>% 
    mutate(ymin_weight = mean_weight - se_weight, 
           ymin_weight = ifelse(ymin_weight < 0, 0, ymin_weight), 
           ymax_weight = mean_weight + se_weight)
  
  sin_signif_label <- sin_data_plot_dat %>% filter(clonality %in% 'clonal') %>% left_join(sin_signif_dat, by = 'Signature')
  
  sin_signature_label <- sin_data_plot_dat %>% filter(clonality %in% 'subclonal')
  
  # pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_clonality_errorbar.pdf')
  
  sin_data_plot_dat %>% 
    ggplot(aes(clonality, mean_weight, group = Signature, col = Signature)) +
    geom_line(show.legend = FALSE) +
    geom_point() +
    labs(x = NULL, col = NULL, y = 'Signature contribution', title = cancer_type) +
    geom_pointrange(aes(ymin = ymin_weight, ymax = ymax_weight)) +
    geom_text_repel(aes(label = label), data = sin_signif_label, xlim = c(NA, 0.8), direction = 'y', segment.alpha = 0.5, segment.linetype = "dotted") +
    geom_text_repel(aes(label = Signature), data = sin_signature_label, size = 2.5, xlim = c(2.2, NA), direction = 'y', segment.alpha = 0.5, segment.linetype = "dotted") +
    # ggsci::scale_color_futurama() +
    scale_color_manual(breaks = signature_colors$Signature, values = signature_colors$color) +
    theme_classic() +
    theme(legend.position = 'none', 
          plot.title = element_text(hjust = 0.5))
  
  # dev.off()
  
}

mutsig_clo_fit_errorbar_plots_selected <- c('Prost-AdenoCA', 'Skin-Melanoma') %>% map(~forggline_clo_fit_plot(.))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_stm_clonality_fit_errorbar_selected.pdf', width = 6, height = 4)

wrap_plots(mutsig_clo_fit_errorbar_plots_selected) + plot_layout(ncol = 2)

dev.off()


mutsig_clo_fit_errorbar_plots_other <- mutsig_clo_fit_strict_test_signif %>% 
  distinct(ICGC_abbr_top) %>% 
  filter(!ICGC_abbr_top %in% c('Prost-AdenoCA', 'Skin-Melanoma')) %>% 
  arrange(ICGC_abbr_top) %>% 
  mutate(plot_list = map(ICGC_abbr_top, ~forggline_clo_fit_plot(.)))
# 24

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/mutsig_stm_clonality_fit_errorbar_other.pdf', width = 21, height = 9)

wrap_plots(mutsig_clo_fit_errorbar_plots_other$plot_list) + plot_layout(nrow = 3)

dev.off()

# Figure 4C ---------------------------------------------------------------
# barplot
# comparison of mutational signature contribution between trunck and branch

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/PRAD_met_sig_we.Rdata') # PRAD_cluster_sig_we, PRAD_tb_sig_we

signature_colors_tree <- tribble(
  ~Signature, ~color, 
  'SBS1', '#8F6504', 
  'SBS5', '#927C64', 
  'SBS40', '#FFB50E', 
  'SBS2', '#FF6444', 
  'SBS13', '#FFFA00', 
  'SBS18', '#C41707', 
  'SBS3', '#0373C6', 
  'SBS41', '#9C8CC4', 
  'SBS8', '#90AD1C', 
  'SBS37', '#FF5FB3', 
  'SBS28', '#B10DA1', 
  'SBS33', '#CC6677', 
  'SBS58', '#70CBEF',
  'SBS45', '#DEA0FD', 
  'SBS52', '#9C8CC4', 
  'Unknown', 'grey'
)

p_anno_dat <- PRAD_tb_sig_we %>% distinct(icgc_donor_id) %>% mutate(start = 'trunk', end = 'branch', label = '****', y = 1.05)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/trunk_bran_clonality_barplot.pdf', width = 8, height = 4)

PRAD_tb_sig_we %>% 
  mutate(timing = factor(timing, levels = c('trunk', 'branch'))) %>% 
  ggplot() +
  geom_bar(aes(timing, Weights, fill = Signature), stat = 'identity') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(breaks = signature_colors_tree$Signature, values = signature_colors_tree$color) +
  ggsignif::geom_signif(data = p_anno_dat, aes(xmin = start, xmax = end, annotations = label, y_position = y), manual = TRUE) +
  theme_classic() +
  facet_grid2(. ~ icgc_donor_id, strip = strip_vanilla(clip = 'off')) +
  labs(x = NULL, y = 'Fraction of signatures', fill = NULL) +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text.x = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
        strip.background = element_blank(), 
        panel.grid = element_blank())

dev.off()

