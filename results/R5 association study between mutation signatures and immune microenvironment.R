
library(tidyverse)
library(patchwork)
library(ggh4x)
library(GSVA)
library(scatterpie)

# library(enrichplot)
# Attention!!! some packages exist in a removed environment

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_tpm_expr_final.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mutsig_expr_group.Rdata') # mean(Weight > 0) >= 0.05

# caculation of CYT score

comb_CYT_score_data <- comb_tpm_expr_final %>% 
  filter(Hugo_Symbol %in% c('GZMA', 'PRF1')) %>% 
  pivot_longer(-Hugo_Symbol, names_to = 'Sample_ID', values_to = 'tpm_val') %>% 
  group_by(Sample_ID) %>% 
  mutate(CYT_score = exp(mean(log(tpm_val)))) %>% # geometric mean of GZMA and PRF1
  ungroup() %>% 
  distinct(Sample_ID, CYT_score) %>% 
  filter(!is.na(CYT_score))

save(comb_CYT_score_data, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_CYT_score_data.Rdata')

# caculation of T-cell-inflamed GEP score

source('/pub6/Temp/Liaojl/Code/common_R_functions/T.cell.inflamed.GEP.R')

comb_gene_expr_m <- comb_tpm_expr_final %>% column_to_rownames('Hugo_Symbol') %>% mutate_all(~log2(. + 1))

comb_GEP_score_data_list <- predictorScore(exp.matrix = comb_gene_expr_m, geneType = "Name", extraFileDir = "/pub6/temp/xlw/RCodes/Jobs/Job2_MutationReshapeTME/0.Function/Immunophenotype")
comb_GEP_score_data <- comb_GEP_score_data_list$score %>% tibble(patient_id = names(.), GEP_score = .)

save(comb_GEP_score_data, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_GEP_score_data.Rdata')


# immune check points -----------------------------------------------------

# PMID: 34643725
# https://ngdc.cncb.ac.cn/cancerscem/search
# Immune Checkpoint
# BTLA, CD244, CD274, CD276, CD40, CD40LG, CD48, CD80, CD86, CEACAM1, CLEC4G, CTLA4, FGL1, HAVCR2, HLA-A, HLA-B, HLA-C, HMGB1, KIR3DL1, LAG3, LGALS3, LGALS9, NECTIN2, PDCD1, PDCD1LG2, PVR, SNCA, TIGIT, TNFRSF14, TNFRSF18, TNFRSF4, TNFRSF9, TNFSF18, TNFSF4, TNFSF9, VSIR
# 36


immune_check_expr <- comb_tpm_expr_final %>% 
  filter(Hugo_Symbol %in% c('BTLA', 'CD244', 'CD274', 'CD276', 'CD40', 'CD40LG', 'CD48', 'CD80', 'CD86', 'CEACAM1', 'CLEC4G', 'CTLA4', 'FGL1', 'HAVCR2', 'HLA-A', 'HLA-B', 'HLA-C', 'HMGB1', 'KIR3DL1', 'LAG3', 'LGALS3', 'LGALS9', 'NECTIN2', 'PDCD1', 'PDCD1LG2', 'PVR', 'SNCA', 'TIGIT', 'TNFRSF14', 'TNFRSF18', 'TNFRSF4', 'TNFRSF9', 'TNFSF18', 'TNFSF4', 'TNFSF9', 'VSIR')) %>% 
  column_to_rownames('Hugo_Symbol') %>% 
  mutate_all(~log2(. + 1)) %>% 
  rownames_to_column('Hugo_Symbol') %>% 
  as_tibble() %>% 
  pivot_longer(-Hugo_Symbol, names_to = 'Sample_ID', values_to = 'log2(TPM+1)') %>% 
  inner_join(mutsig_expr_group, by = c('Sample_ID' = 'Patient'))

immune_check_test_res <- immune_check_expr %>% 
  group_by(ICGC_abbr_top, Signature, Hugo_Symbol) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~wilcox.test(`log2(TPM+1)` ~ group, data = .)$p.value, NA)), 
         `High > Low` = map_lgl(data, ~group_by(., group) %>% summarise(med_score_value = median(`log2(TPM+1)`)) %>% mutate(High_higher_score = (med_score_value[2] - med_score_value[1]) > 0) %>% pull(High_higher_score) %>% .[1]), 
         `High > Low` = ifelse(p_val < 0.05, `High > Low`, NA)) %>% 
  group_by(ICGC_abbr_top, Hugo_Symbol) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
  ungroup() %>% 
  select(-data)

immune_check_test_res_signif <- immune_check_test_res %>% filter(q_val < 0.05)

# heatmap

immune_check_gene_color <- tribble(
  ~Hugo_Symbol, ~color,
  'BTLA', '#0C487A',
  'CD40', '#BAA133',
  'CD40LG', '#F8D0D6',
  'CD48', '#EE9028',
  'CD80', '#0C8D46',
  'CD86', '#D92F86',
  'CD244', '#ECB065',
  'CD274', '#A74C93', 
  'CD276', '#DD2929', 
  'CEACAM1', '#374F99',
  'CLEC4G', '#90AD1C',
  'CTLA4', '#CC98BF',
  'FGL1', '#9DD6F0',
  'HAVCR2', '#0E78AA', 
  'HLA-A', '#DC6B72',
  'HLA-B', '#93C9A3',
  'HLA-C', '#7A378B',
  'HMGB1', '#B2212D',
  'KIR3DL1', '#EEAAAE',
  'LAG3', '#C7C9D7',
  'LGALS3', '#CEC0DC',
  'LGALS9', '#997FB5', 
  'PDCD1', '#724C29', 
  'PDCD1LG2', '#502D80',
  'PVR', '#EDE33E',
  'SNCA', '#D8EEFB',
  'TIGIT', '#0FA094',
  'TNFRSF4', '#CF7A29',
  'TNFRSF9', '#6B789B',
  'TNFRSF14', '#DCBE29',
  'TNFRSF18', '#CD6600',
  'TNFSF4', '#7A1A1D',
  'TNFSF9', '#CAA98D', 
  'TNFSF18', '#B3CB47'
)

immune_check_test_res_signif %>% distinct(Hugo_Symbol) %>% transmute(Hugo_Symbol = str_sort(Hugo_Symbol, numeric = TRUE)) %>% print(n=34)


immue_signature_sort <- immune_check_test_res_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)
immue_cancer_sort <- immune_check_test_res_signif %>% distinct(ICGC_abbr_top) %>% pull(ICGC_abbr_top) %>% str_sort(numeric = TRUE)

# only display cancer types with 4 or more significant genes

immune_check_retained_cancer <- immune_check_test_res_signif %>% distinct(ICGC_abbr_top, Hugo_Symbol) %>% count(ICGC_abbr_top) %>% filter(n >= 6)

immune_check_heatmap_dat <- immune_check_test_res_signif %>% 
  semi_join(immune_check_retained_cancer, by = 'ICGC_abbr_top') %>% 
  mutate(sign = ifelse(`High > Low`, 1, -1), 
         log_p_val = sign * (-log10(p_val))) %>% 
  select(ICGC_abbr_top, Hugo_Symbol, Signature, log_p_val)

# only display signatures with 2 or more significant associations

immune_check_retained_signature <- immune_check_heatmap_dat %>% count(Signature) %>% filter(n >= 2)

immune_check_heatmap_plot_dat <- immune_check_heatmap_dat %>% 
  semi_join(immune_check_retained_signature, by = 'Signature') %>% 
  mutate(Signature = factor(Signature, levels = immue_signature_sort), 
         ICGC_abbr_top = factor(ICGC_abbr_top, levels = immue_cancer_sort))

immune_check_tile_dat <- immune_check_heatmap_plot_dat %>% mutate(ttype_gene_id = str_c(ICGC_abbr_top, Hugo_Symbol, sep = ':')) %>% distinct(ICGC_abbr_top, Hugo_Symbol, ttype_gene_id)


immune_check_heatmap_plot <- immune_check_heatmap_plot_dat %>% 
  ggplot(aes(Hugo_Symbol, Signature)) +
  geom_tile(aes(fill = log_p_val)) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = c('midnightblue', 'white', 'darkred')) +
  ggthemes::theme_few() +
  facet_grid2(. ~ ICGC_abbr_top, scales = 'free', space = 'free', strip = strip_vanilla(clip = 'off')) +
  theme(axis.ticks = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(), 
        panel.background = element_rect(fill = "#C7C4C1"), 
        panel.border = element_rect(size = 0.5),
        # strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3),
        strip.background = element_blank(), 
        panel.spacing = unit(0, 'cm'))


immune_check_anno_plot <- immune_check_tile_dat %>% 
  ggplot(aes(ttype_gene_id, 1, fill = Hugo_Symbol)) +
  geom_bar(stat = 'identity') +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(breaks = immune_check_gene_color$Hugo_Symbol, values = immune_check_gene_color$color) +
  ggthemes::theme_few() +
  facet_grid( ~ ICGC_abbr_top, scales = 'free', space = 'free') +
  guides(fill = guide_legend(nrow = 4)) +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing = unit(0, 'cm'), 
        strip.background = element_blank(), 
        panel.border = element_rect(size = 0.5), 
        strip.text = element_blank())


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/mutsig_immune_check_heatmap.pdf', width = 16, height = 8)

immune_check_heatmap_plot / immune_check_anno_plot + plot_layout(heights = c(20, 1), guides = 'collect') & theme(legend.position = 'bottom')
# immune_check_anno_plot / immune_check_heatmap_plot + plot_layout(heights = c(1, 20), guides = 'collect') & theme(legend.position = 'top')

dev.off()


# violin plots

forImmuneBoxplot <- function(which_signature, which_cancer_type, which_gene){
  
  # which_signature <- 'SBS1'
  # which_cancer_type <- 'Panc-AdenoCA'
  # which_gene <- 'LAG3'

  immune_check_expr %>% 
    filter(Signature %in% which_signature, ICGC_abbr_top %in% which_cancer_type, Hugo_Symbol %in% which_gene) %>% 
    ggplot(aes(group, `log2(TPM+1)`, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = 'white', outlier.color = NA) +
    labs(x = NULL, title = str_c(which_signature, which_cancer_type, sep = ', '), fill = NULL) +
    scale_fill_manual(values = c('#5873A8', '#C98035'), labels = c('Low contribution', 'High contribution')) +
    annotate('text', x = 1.5, y = -Inf, label = which_gene, vjust = -0.5, hjust = 0.5) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.grid = element_blank())
  
}

Bladder_SBS1_LAG3_plot <- forImmuneBoxplot('SBS1', 'Bladder-TCC', 'LAG3') + annotate('text', x = 1.5, y = Inf, label = 'p=4.83e-9', vjust = 1.5)
Breast_SBS13_CTLA4_plot <- forImmuneBoxplot('SBS13', 'Breast', 'CTLA4') + annotate('text', x = 1.5, y = Inf, label = 'p=1.25e-14', vjust = 1.5)
ColoRect_SBS18_CD274_plot <- forImmuneBoxplot('SBS18', 'ColoRect-AdenoCA', 'CD274') + annotate('text', x = 1.5, y = Inf, label = 'p=6.25e-3', vjust = 1.5)
Stomach_SBS20_PDCD1_plot <- forImmuneBoxplot('SBS20', 'Stomach-AdenoCA', 'PDCD1') + annotate('text', x = 1.5, y = Inf, label = 'p=9.93e-5', vjust = 1.5)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/immune_check_expr_violin_plot.pdf', width = 6, height = 6)

wrap_plots(list(Bladder_SBS1_LAG3_plot, Breast_SBS13_CTLA4_plot, ColoRect_SBS18_CD274_plot, Stomach_SBS20_PDCD1_plot)) + plot_layout(nrow = 2, guides = 'collect') & theme(legend.position = 'bottom')

dev.off()

# CYT, T-cell-inflamed GEP score and immune score ---------------------------------------

cyt_inflamed_test_dat <- comb_CYT_score_data %>% 
  left_join(comb_GEP_score_data, by = c('Sample_ID' = 'patient_id')) %>% 
  pivot_longer(-Sample_ID, names_to = 'metrics', values_to = 'score') %>% 
  inner_join(mutsig_group, by = c('Sample_ID' = 'Patient'))

cyt_inflamed_test_res <- cyt_inflamed_test_dat %>% 
  group_by(ICGC_abbr_top, Signature, metrics) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~wilcox.test(score ~ group, data = .)$p.value, NA)), 
         `High > Low` = map_lgl(data, ~group_by(., group) %>% summarise(med_score_value = median(score)) %>% mutate(High_higher_score = (med_score_value[2] - med_score_value[1]) > 0) %>% pull(High_higher_score) %>% .[1]), 
         `High > Low` = ifelse(p_val < 0.05, `High > Low`, NA)) %>% 
  group_by(ICGC_abbr_top, metrics) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
  ungroup() %>% 
  select(-data)

cyt_inflamed_test_res_signif <- cyt_inflamed_test_res %>% filter(q_val < 0.05)


# CYT, GEP, immune score dot plot

# immune score

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/xcell_test_res_signif.Rdata')

immune_score_test_res_signif <- xcell_test_res_signif %>% filter(cell_type %in% 'immune score') %>% rename(metrics = cell_type)

immune_metrics_res_signif <- cyt_inflamed_test_res_signif %>% bind_rows(immune_score_test_res_signif)

immune_metrics_sig_idx <- immune_metrics_res_signif %>% mutate(Signature = str_sort(Signature, numeric = TRUE)) %>% distinct(Signature) %>% rowid_to_column('sig_idx')
immune_metrics_ttype_idx <- immune_metrics_res_signif %>% distinct(ICGC_abbr_top) %>% arrange(ICGC_abbr_top) %>% rowid_to_column('ttype_idx')


immune_metrics_plot_dat <- immune_metrics_res_signif %>% 
  left_join(immune_metrics_sig_idx, by = 'Signature') %>% 
  left_join(immune_metrics_ttype_idx, by = 'ICGC_abbr_top') %>% 
  mutate(sign = ifelse(`High > Low`, 1, -1), 
         sign_logP = sign * (-log10(p_val)))


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/immune_metrics_dotplot.pdf', width = 6, height = 8)

ggplot() +
  geom_point(data = filter(immune_metrics_plot_dat, metrics %in% 'CYT_score'), aes(sig_idx, ttype_idx, col = sign_logP), pch = 1, size = 5) +
  geom_point(data = filter(immune_metrics_plot_dat, metrics %in% 'GEP_score'), aes(sig_idx, ttype_idx, col = sign_logP), pch = 2, size = 5) +
  geom_point(data = filter(immune_metrics_plot_dat, metrics %in% 'immune score'), aes(sig_idx, ttype_idx, col = sign_logP), pch = 6, size = 5) +
  scale_color_gradientn(colors = c('midnightblue', 'white', 'darkred')) +
  scale_x_continuous(breaks = immune_metrics_sig_idx$sig_idx, labels = immune_metrics_sig_idx$Signature) +
  scale_y_continuous(breaks = immune_metrics_ttype_idx$ttype_idx, labels = immune_metrics_ttype_idx$ICGC_abbr_top) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# ssGSEA on immune pathways -------------------------------------------------

load('/pub6/Temp/Liaojl/Data/MutSig_Driver/Processed_Data/Transcriptomic_Data/GRCh37_protein_coding_gene.Rdata')

immune_pathway_gene <- qusage::read.gmt('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/immune_pathway_gene.gmt')

immune_pathway_symbol <- list()

for(i in names(immune_pathway_gene)){
  
  path_gene <- immune_pathway_gene[[i]]
  immune_pathway_symbol[[i]] <- GRCh37_protein_coding_gene %>% filter(Gene_ID %in% path_gene) %>% pull(Gene_Symbol)
  
}


comb_tpm_expr_final_mat <- comb_tpm_expr_final %>% column_to_rownames('Hugo_Symbol') %>% as.matrix()

mutsig_immune_ssgsea_score <- gsva(comb_tpm_expr_final_mat, immune_pathway_symbol, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)

immune_ssgsea_group_data <- mutsig_immune_ssgsea_score %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Patient') %>% 
  as_tibble() %>% 
  pivot_longer(-Patient, names_to = 'immune_pathway', values_to = 'ssGSEA_score') %>% 
  left_join(mutsig_expr_group, by = 'Patient') %>% 
  dplyr::select(-Project)


immune_ssgsea_score_test <- immune_ssgsea_group_data %>% 
  group_by(ICGC_abbr_top, Signature, immune_pathway) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~wilcox.test(ssGSEA_score ~ group, data = .)$p.value, NA)), 
         `High > Low` = map_lgl(data, ~group_by(., group) %>% summarise(med_score_value = median(ssGSEA_score)) %>% mutate(High_higher_score = (med_score_value[2] - med_score_value[1]) > 0) %>% pull(High_higher_score) %>% .[1]), 
         `High > Low` = ifelse(p_val < 0.05, `High > Low`, NA)) %>% 
  dplyr::select(-data) %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
  ungroup()

immune_ssgsea_score_test_signif <- immune_ssgsea_score_test %>% filter(q_val < 0.05)
# 371

# heatmap

immue_ssgsea_signature_sort <- immune_ssgsea_score_test_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)
immue_ssgsea_cancer_sort <- immune_ssgsea_score_test_signif %>% distinct(ICGC_abbr_top) %>% pull(ICGC_abbr_top) %>% str_sort(numeric = TRUE)

# only display cancer types with 3 or more significant genes

immune_ssgsea_retained_cancer <- immune_ssgsea_score_test_signif %>% distinct(ICGC_abbr_top, immune_pathway) %>% count(ICGC_abbr_top) %>% filter(n >= 3)

immune_ssgsea_score_plot_im <- immune_ssgsea_score_test_signif %>% 
  semi_join(immune_ssgsea_retained_cancer, by = 'ICGC_abbr_top') %>% 
  mutate(sign = ifelse(`High > Low` == TRUE, 1, -1), 
         sign_log_p = sign * (-log10(p_val))) %>% 
  dplyr::select(-sign)

# only display signatures with 2 or more significant associations

immune_ssgsea_retained_signature <- immune_ssgsea_score_plot_im %>% count(Signature) %>% filter(n >= 2)

immune_ssgsea_score_plot_dat <- immune_ssgsea_score_plot_im %>% 
  semi_join(immune_ssgsea_retained_signature, by = 'Signature') %>% 
  mutate(Signature = factor(Signature, levels = immue_ssgsea_signature_sort), 
         ICGC_abbr_top = factor(ICGC_abbr_top, levels = immue_ssgsea_cancer_sort))

ssgsea_pathway_tile_dat <- immune_ssgsea_score_plot_dat %>% mutate(ttype_path_id = str_c(ICGC_abbr_top, immune_pathway, sep = ':')) %>% distinct(ICGC_abbr_top, immune_pathway, ttype_path_id)

Immune_pathways_colors <- tribble(
  ~pathway, ~color,
  'Antigen_Processing_and_Presentation', '#FEDB05',
  'Antimicrobials', '#FFF385',
  'BCRSignalingPathway', '#7F7F7D',
  'Chemokine_Receptors', '#644469',
  'Chemokines', '#EF3536',
  'Cytokine_Receptors', '#C10504',
  'Cytokines', '#B6B6B6',
  'Interferon_Receptor', '#E56905', 
  'Interferons', '#FBA629', 
  'Interleukins', '#CFD5C9',
  'Interleukins_Receptor', '#C692C0',
  'NaturalKiller_Cell_Cytotoxicity', '#2E3336',
  'TCRsignalingPathway', '#0C54C3',
  'TGFb_Family_Member', '#7DAFE2',
  'TGFb_Family_Member_Receptor', '#A98D0F',
  'TNF_Family_Members', '#C17D12',
  'TNF_Family_Members_Receptors', '#E9B96D'
)

mutsig_immune_gsea_hp <- immune_ssgsea_score_plot_dat %>% 
  ggplot(aes(immune_pathway, Signature)) +
  geom_tile(aes(fill = sign_log_p)) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = c('midnightblue', 'white', 'darkred')) +
  ggthemes::theme_few() +
  facet_grid2(. ~ ICGC_abbr_top, scales = 'free', space = 'free', strip = strip_vanilla(clip = 'off')) +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        panel.background = element_rect(fill = "#C7C4C1"), 
        panel.border = element_rect(color = 'grey60', size = 0.5), 
        strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
        strip.background = element_blank(), 
        panel.spacing = unit(0, 'cm'))

immune_ssgsea_anno_plot <- ssgsea_pathway_tile_dat %>% 
  ggplot(aes(ttype_path_id, 1, fill = immune_pathway)) +
  geom_bar(stat = 'identity') +
  labs(x = NULL, y = NULL, fill = 'Immune pathways') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(breaks = Immune_pathways_colors$pathway, values = Immune_pathways_colors$color) +
  ggthemes::theme_few() +
  facet_grid( ~ ICGC_abbr_top, scales = 'free', space = 'free') +
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing = unit(0, 'cm'), 
        strip.background = element_blank(), 
        panel.border = element_rect(color = 'grey60', size = 0.5), 
        strip.text = element_blank())

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/mutsig_immune_ssgsea_heatmap.pdf', width = 14.5, height = 7)

mutsig_immune_gsea_hp / immune_ssgsea_anno_plot + plot_layout(heights = c(20, 1), guides = 'collect')

dev.off()

# violin plots

forImmunessgSEAplot <- function(which_signature, which_cancer_type, which_immune_pathway){
  
  immune_ssgsea_group_data %>% 
    filter(Signature %in% which_signature, ICGC_abbr_top %in% which_cancer_type, immune_pathway %in% which_immune_pathway) %>% 
    ggplot(aes(group, ssGSEA_score, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = 'white', outlier.color = NA) +
    labs(x = NULL, y = 'ssGSEA score', title = str_c(which_signature, which_cancer_type, sep = ', '), fill = NULL) +
    scale_fill_manual(values = c('#5873A8', '#C98035'), labels = c('Low contribution', 'High contribution')) +
    annotate('text', x = 1.5, y = -Inf, label = which_immune_pathway, vjust = -0.5, hjust = 0.5) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.grid = element_blank())
  
}

Breast_SBS1_Chemokines_plot <- forImmunessgSEAplot('SBS1', 'Breast', 'Chemokines') + annotate('text', x = 1.5, y = Inf, label = 'p=4.10e-11', vjust = 1.5)
Cervix_SBS13_Antimi_plot <- forImmunessgSEAplot('SBS13', 'Cervix', 'Antimicrobials') + annotate('text', x = 1.5, y = Inf, label = 'p=3.38e-3', vjust = 1.5)
ColoRect_SBS15_Interferons_plot <- forImmunessgSEAplot('SBS15', 'ColoRect-AdenoCA', 'Interferons') + annotate('text', x = 1.5, y = Inf, label = 'p=6.71e-6', vjust = 1.5)
Stomach_SBS20_TNF_plot <- forImmunessgSEAplot('SBS20', 'Stomach-AdenoCA', 'TNF_Family_Members') + annotate('text', x = 1.5, y = Inf, label = 'p=2.86e-6', vjust = 1.5)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/immune_pathway_ssGSEA_violin_plot.pdf', width = 6, height = 6)

wrap_plots(list(Breast_SBS1_Chemokines_plot, Cervix_SBS13_Antimi_plot, ColoRect_SBS15_Interferons_plot, Stomach_SBS20_TNF_plot)) + plot_layout(nrow = 2, guides = 'collect') & theme(legend.position = 'bottom')

dev.off()


# intersection of immune genes --------------------------------------------

# immune-related pathway gene

sub_len <- immune_pathway_symbol %>% map_dbl(~length(.))
path_vec <- rep(names(sub_len), times = sub_len)
immune_pathway_gene_dat <- tibble(Hugo_Symbol = unlist(immune_pathway_symbol), pathway = path_vec)

immune_pathway_idx <- immune_pathway_gene_dat %>% distinct(pathway) %>% arrange(pathway) %>% rowid_to_column('x_idx')

immune_gene_common_count <- immune_pathway_gene_dat %>% count(Hugo_Symbol, name = 'type') %>% mutate(type = as.character(type))

immune_pathway_gene_stat <- immune_pathway_gene_dat %>% 
  left_join(immune_gene_common_count, by = 'Hugo_Symbol') %>% 
  count(pathway, type) %>% 
  group_by(pathway) %>% 
  mutate(gene_num = sum(n), 
         log_gene_num = log2(gene_num), 
         frac = n/gene_num) %>% 
  ungroup() %>% 
  left_join(immune_pathway_idx, by = 'pathway') %>% 
  select(-n) %>% 
  pivot_wider(names_from = type, values_from = frac) %>% 
  mutate_all(~ifelse(is.na(.), 0, .))

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/immune_pathway_gene_count_scatterpie.pdf', width = 8, height = 4)

immune_pathway_gene_stat %>% 
  ggplot() +
  geom_scatterpie(aes(x_idx, log_gene_num, r = log_gene_num * 0.1), cols = unique(immune_gene_common_count$type), data = immune_pathway_gene_stat) +
  coord_equal() +
  scale_x_continuous(expand = c(0.01, 0.01), breaks = immune_pathway_idx$x_idx, labels = immune_pathway_idx$pathway) +
  labs(x = NULL, y = 'log2(gene count)', fill = NULL) +
  ggsci::scale_fill_jco() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        legend.position = 'top')

dev.off()

