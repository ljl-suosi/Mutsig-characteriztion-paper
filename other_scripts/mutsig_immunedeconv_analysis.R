
library(immunedeconv)
library(tidyverse)
library(patchwork)
library(ggh4x)

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_tpm_expr_final.Rdata')
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mutsig_expr_group.Rdata') # mean(Weight > 0) >= 0.05

comb_tpm_expr_final_mat <- comb_tpm_expr_final %>% column_to_rownames('Hugo_Symbol') %>% as.matrix()

res_mcp_counter <- deconvolute(comb_tpm_expr_final_mat, "mcp_counter")
# save(res_mcp_counter, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/immune_deconv/res_mcp_counter.Rdata')

res_xcell <- deconvolute(comb_tpm_expr_final_mat, "xcell")
# save(res_xcell, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/immune_deconv/res_xcell.Rdata')

# test function

mutsig_cell_test_fun <- function(cell_score_dat, mutsig_expr_group){
  
  cell_test_dat <- cell_score_dat %>% 
    column_to_rownames('cell_type') %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('Patient') %>% 
    as_tibble() %>% 
    pivot_longer(-Patient, names_to = 'cell_type', values_to = 'frac') %>% 
    inner_join(mutsig_expr_group, by = 'Patient')
  
  cell_test_res <- cell_test_dat %>% 
    group_by(ICGC_abbr_top, Signature, cell_type) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(p_val = map_dbl(data, possibly(~wilcox.test(frac ~ group, data = .)$p.value, NA)), 
           `High > Low` = map_lgl(data, ~group_by(., group) %>% summarise(med_score_value = median(frac)) %>% mutate(High_higher_score = (med_score_value[2] - med_score_value[1]) > 0) %>% pull(High_higher_score) %>% .[1]), 
           `High > Low` = ifelse(p_val < 0.05, `High > Low`, NA)) %>% 
    group_by(ICGC_abbr_top, cell_type) %>% 
    mutate(q_val = p.adjust(p_val, method = 'BH')) %>% 
    ungroup() %>% 
    select(-data)
  
  cell_test_res_signif <- cell_test_res %>% filter(q_val < 0.05)

}


# plot function


cell_heatmap_fun <- function(cell_test_res_signif){
  
  signature_sort <- cell_test_res_signif %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)
  cancer_sort <- cell_test_res_signif %>% distinct(ICGC_abbr_top) %>% pull(ICGC_abbr_top) %>% str_sort(numeric = TRUE)
  
  heatmap_dat <- cell_test_res_signif %>% 
    mutate(sign = ifelse(`High > Low` == TRUE, 1, -1), 
           sign_log_p = sign * (-log10(p_val))) %>% 
    dplyr::select(-sign)
  
  # tile_dat <- mcp_counter_test_res_signif %>% mutate(ttype_ctype_id = str_c(ICGC_abbr_top, cell_type, sep = ':')) %>% distinct(ICGC_abbr_top, cell_type, ttype_ctype_id)
  tile_dat <- cell_test_res_signif %>% mutate(ttype_ctype_id = str_c(cell_type, ICGC_abbr_top, sep = ':')) %>% distinct(ICGC_abbr_top, cell_type, ttype_ctype_id)
  
  
  heatmap_plot <- heatmap_dat %>% 
    mutate(Signature = factor(Signature, levels = signature_sort), 
           ICGC_abbr_top = factor(ICGC_abbr_top, levels = cancer_sort)) %>% 
    ggplot(aes(ICGC_abbr_top, Signature)) +
    geom_tile(aes(fill = sign_log_p)) +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(low = 'midnightblue', mid = 'white', high = 'darkred') +
    ggthemes::theme_few() +
    facet_grid2(. ~ cell_type, scales = 'free', space = 'free', strip = strip_vanilla(clip = 'off')) +
    theme(axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_blank(), 
          panel.background = element_rect(fill = "#C7C4C1"), 
          panel.border = element_rect(colour = "grey", size = 0.5), 
          strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
          strip.background = element_blank(), 
          panel.spacing = unit(0, 'cm')) +
    my_legend_theme
  
  anno_plot <- tile_dat %>% 
    ggplot(aes(ttype_ctype_id, 1, fill = ICGC_abbr_top)) +
    geom_bar(stat = 'identity') +
    labs(x = NULL, y = NULL, fill = 'Cancer types') +
    scale_x_discrete(expand = c(0, 0)) +
    # scale_fill_manual(breaks = LM22_colors$cell_type, values = LM22_colors$color) +
    ggthemes::theme_few() +
    facet_grid( ~ cell_type, scales = 'free', space = 'free') +
    guides(fill = guide_legend(ncol = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing = unit(0, 'cm'), 
          strip.background = element_blank(), 
          panel.border = element_rect(colour = "grey", size = 0.5), 
          strip.text = element_blank()) +
    my_legend_theme  
  
  comb_plots <- heatmap_plot / anno_plot + plot_layout(heights = c(20, 1), guides = 'collect')
  
}

# mcp_counter -------------------------------------------------------------

load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/immune_deconv/res_mcp_counter.Rdata') # 190

mcp_counter_test_res_signif <- mutsig_cell_test_fun(res_mcp_counter, mutsig_expr_group)
# save(mcp_counter_test_res_signif, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/mcp_counter_test_res_signif.Rdata')

# heatmap

mcp_counter_cell_heatmap <- cell_heatmap_fun(epic_test_res_signif)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/immune_cell_infiltration/mutsig_mcp_counter_heatmap.pdf', width = 22, height = 6)

dev.off()

# radical plots

retain_cancer_type <- mcp_counter_test_res_signif %>% count(ICGC_abbr_top) %>% filter(n > 2) # cancer types with more than 2 significant relations
retain_signature <- mcp_counter_test_res_signif %>% count(Signature) %>% filter(n > 2)

mcp_counter_test_res_signif_cut <- mcp_counter_test_res_signif %>% 
  semi_join(retain_cancer_type, by = 'ICGC_abbr_top') %>% 
  semi_join(retain_signature, by = 'Signature') %>% 
  filter(cell_type != 'cytotoxicity score')
# 149
# 10 cell types, 16 cancer types, 16 signatures

cancer_colors <- tribble(
  ~ICGC_abbr_top, ~color, 
  'Bladder-TCC', '#658FB1', 
  'Breast', '#E90500', 
  'Cervix', '#3E5262', 
  'ColoRect-AdenoCA', '#4CAFD2',  
  'Eso-AdenoCA', '#BAE5EE',  
  'Head-SCC', '#E66304',  
  'Kidney-RCC', '#F49E54',  
  'Liver-HCC', '#FEDBBA',  
  'Lung-AdenoCA', '#DCE4B0',  
  'Lung-SCC', '#799525',  
  'Lymph-BNHL', '#F0DBB4',  
  'Lymph-CLL', '#2EC8C3',  
  'Panc-AdenoCA', '#8978A1',  
  'Skin-Melanoma', '#53437D',  
  'Stomach-AdenoCA', '#54B672',  
  'Thy-AdenoCA', '#F3E97B'
)

# radical plots

mcp_ctype_sort_hh <- mcp_counter_test_res_signif_cut %>% 
  distinct(cell_type) %>% 
  pull(cell_type) %>% 
  str_c('hh', ., sep = '_') %>% 
  str_sort(numeric = TRUE)

mcp_ctype_sort_hl <- mcp_counter_test_res_signif_cut %>% 
  distinct(cell_type) %>% 
  pull(cell_type) %>% 
  str_c('hl', ., sep = '_') %>% 
  str_sort(numeric = TRUE, decreasing = TRUE)

mcp_signature_sort <- mcp_counter_test_res_signif_cut %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

mcp_radical_dat <- mcp_counter_test_res_signif_cut %>% 
  mutate(cor = ifelse(`High > Low`, 'hh', 'hl'), 
         log_p_val = -log10(p_val)) %>% 
  select(ICGC_abbr_top, cell_type, Signature, cor, log_p_val) %>% 
  pivot_wider(names_from = cor, values_from = log_p_val) %>% 
  pivot_longer(hh:hl, names_to = 'cor', values_to = 'log_p_val') %>% 
  mutate(ctype_id = factor(str_c(cor, cell_type, sep = '_'), levels = c(mcp_ctype_sort_hh, mcp_ctype_sort_hl)))

########## for selected Signatures

mcp_radical_dat_selected <- mcp_radical_dat %>% filter(Signature %in% c('SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS40'))
cancer_colors_selected <- cancer_colors %>% semi_join(mcp_radical_dat_selected, by = 'ICGC_abbr_top')
mcp_signature_sort_selected <- mcp_radical_dat_selected %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

selected_radical_plots <- list()

for(i in mcp_signature_sort_selected){
  
  print(str_c('processing ', i, ' ...'))
    
  sincan_dat <- mcp_radical_dat_selected %>% filter(Signature %in% i) %>% arrange(ctype_id)
  sincan_uni_ctypeid <- unique(sincan_dat$ctype_id)
  
  sincan_y_val <- sincan_dat %>% 
    mutate(log_p_val = ifelse(is.na(log_p_val), 0, log_p_val)) %>% 
    group_by(ctype_id) %>% 
    summarise(y_sum = sum(log_p_val)) %>% 
    ungroup()
  
  sincan_ymax <- sincan_y_val %>% pull(y_sum) %>% max()
  
  background_coord_red <- tibble(xmin = 0.5, xmax = 0.5 + (length(sincan_uni_ctypeid)/2), ymin = 0, ymax = sincan_ymax)
  background_coord_blue <- tibble(xmin = 0.5 + (length(sincan_uni_ctypeid)/2), xmax = 0.5 + (length(sincan_uni_ctypeid)), ymin = 0, ymax = sincan_ymax)
  
  sincan_ctype_labels <- sincan_y_val %>% 
    left_join(distinct(sincan_dat, ctype_id, cell_type), by = 'ctype_id') %>% 
    mutate(label = ifelse(y_sum > 0, cell_type, '')) %>% 
    pull(label)
  
  sincan_plot <- ggplot() +
    geom_bar(data = sincan_dat, aes(ctype_id, log_p_val, fill = ICGC_abbr_top), stat = 'identity', col = 'grey60', size = 0.1) +
    geom_blank(aes(fill = ICGC_abbr_top), data = cancer_colors_selected) + # geom_blank for extral legends
    geom_vline(xintercept = c(0.5, 0.5 + length(sincan_uni_ctypeid)/2), colour = "black") +
    geom_rect(data = background_coord_red, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "red", alpha = 0.08) +
    geom_rect(data = background_coord_blue, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.08) +
    coord_polar(theta = "x") +
    scale_fill_manual(breaks = cancer_colors_selected$ICGC_abbr_top, values = cancer_colors_selected$color) +
    scale_x_discrete(labels = sincan_ctype_labels) +
    labs(x = NULL, y = NULL, fill = NULL, title = i) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          panel.grid = element_line(linetype = 2),
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          legend.background = element_blank())
  
  selected_radical_plots[[i]] <- sincan_plot
  
}

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/mutsig_mcp_counter_radical_plot_selected.pdf', width = 16, height = 8)

wrap_plots(selected_radical_plots) + guide_area() + plot_layout(nrow = 2, guides = 'collect')

dev.off()


########## for other Signatures

mcp_radical_dat_other <- mcp_radical_dat %>% filter(!Signature %in% mcp_signature_sort_selected)
cancer_colors_other <- cancer_colors %>% semi_join(mcp_radical_dat_other, by = 'ICGC_abbr_top')
mcp_signature_sort_other <- mcp_radical_dat_other %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

other_radical_plots <- list()

for(i in mcp_signature_sort_other){
  
  print(str_c('processing ', i, ' ...'))
    
  sincan_dat <- mcp_radical_dat_other %>% filter(Signature %in% i) %>% arrange(ctype_id)
  sincan_uni_ctypeid <- unique(sincan_dat$ctype_id)
  
  sincan_y_val <- sincan_dat %>% 
    mutate(log_p_val = ifelse(is.na(log_p_val), 0, log_p_val)) %>% 
    group_by(ctype_id) %>% 
    summarise(y_sum = sum(log_p_val)) %>% 
    ungroup()
  
  sincan_ymax <- sincan_y_val %>% pull(y_sum) %>% max()
  
  background_coord_red <- tibble(xmin = 0.5, xmax = 0.5 + (length(sincan_uni_ctypeid)/2), ymin = 0, ymax = sincan_ymax)
  background_coord_blue <- tibble(xmin = 0.5 + (length(sincan_uni_ctypeid)/2), xmax = 0.5 + (length(sincan_uni_ctypeid)), ymin = 0, ymax = sincan_ymax)
  
  sincan_ctype_labels <- sincan_y_val %>% 
    left_join(distinct(sincan_dat, ctype_id, cell_type), by = 'ctype_id') %>% 
    mutate(label = ifelse(y_sum > 0, cell_type, '')) %>% 
    pull(label)
  
  sincan_plot <- ggplot() +
    geom_bar(data = sincan_dat, aes(ctype_id, log_p_val, fill = ICGC_abbr_top), stat = 'identity', col = 'grey60', size = 0.1) +
    geom_blank(aes(fill = ICGC_abbr_top), data = cancer_colors_other) + # geom_blank for extral legends
    geom_vline(xintercept = c(0.5, 0.5 + length(sincan_uni_ctypeid)/2), colour = "black") +
    geom_rect(data = background_coord_red, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "red", alpha = 0.08) +
    geom_rect(data = background_coord_blue, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.08) +
    coord_polar(theta = "x") +
    scale_fill_manual(breaks = cancer_colors_other$ICGC_abbr_top, values = cancer_colors_other$color) +
    scale_x_discrete(labels = sincan_ctype_labels) +
    labs(x = NULL, y = NULL, fill = NULL, title = i) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          panel.grid = element_line(linetype = 2),
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          legend.background = element_blank())
  
  other_radical_plots[[i]] <- sincan_plot
  
}

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_immune_analysis/mutsig_mcp_counter_radical_plot_other.pdf', width = 24, height = 8)

wrap_plots(other_radical_plots) + guide_area() + plot_layout(nrow = 2, guides = 'collect')

dev.off()

