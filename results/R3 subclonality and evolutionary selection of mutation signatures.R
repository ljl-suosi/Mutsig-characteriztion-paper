
library(tidyverse)
library(patchwork)
library(ggh4x)
library(gggibbous) # moon charts
library(scatterpie)

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mut_sam_project.Rdata')

comb_driver_gene <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_driver_gene.tsv')
comb_mutsig_clo_data <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_clo_alt_data.tsv')
comb_mutsig_we <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_mutsig_we.tsv')

cancer_mutsig <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/cancer_mutsig.tsv') # frequency >= 0.05


comb_mutsig_we_filter <- comb_mutsig_we %>% 
  group_by(ICGC_abbr_top, Signature) %>% 
  filter(mean(Weight > 0) >= 0.05) %>% 
  ungroup()


# Attention!!! comb_mutsig_clo_data must be limited to Samples of which signature accuary >= 0.8

comb_mutsig_clo_data_lim <- comb_mutsig_clo_data %>% 
  semi_join(comb_mutsig_we_filter, by = c('Sample', 'Attri_Signature' = 'Signature'))
# 7993 samples, 23,422,116 muts

comb_mutsig_clo_ns_data <- comb_mutsig_clo_data_lim %>% filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site'))
comb_mutsig_clo_ns_dri_data <- comb_mutsig_clo_ns_data %>% semi_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top'))

all_mutsig_sort <- comb_mutsig_clo_data_lim %>% filter(!is.na(Attri_Signature)) %>% distinct(Attri_Signature) %>% pull(Attri_Signature) %>% str_sort(numeric = TRUE)


# Fraction of driver mutation (limited to PCAWG samples) ---------------------------------------------

pcawg_sam <- mut_sam_project %>% filter(Project %in% 'PCAWG')

pcawg_mutsig_mut_count <- comb_mutsig_clo_data_lim %>% 
  semi_join(pcawg_sam, by = 'Sample') %>% 
  filter(!is.na(Attri_Signature)) %>% 
  count(ICGC_abbr_top, Attri_Signature, name = 'Mut_Count')

pcawg_mutsig_ns_dri_frac <- comb_mutsig_clo_data_lim %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site')) %>% 
  semi_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top')) %>% 
  semi_join(pcawg_sam, by = 'Sample') %>% 
  filter(!is.na(Attri_Signature)) %>% 
  count(ICGC_abbr_top, Attri_Signature, name = 'NS_driCount') %>% 
  right_join(pcawg_mutsig_mut_count, by = c('ICGC_abbr_top', 'Attri_Signature')) %>% 
  mutate(NS_driCount = ifelse(is.na(NS_driCount), 0, NS_driCount), 
         Dri_frac = NS_driCount/Mut_Count)

# write_tsv(mutsig_ns_dri_frac, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_ns_dri_frac.tsv')
# mutsig_ns_dri_frac <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_ns_dri_frac.tsv')


# only prevalent signatures

mutsig_ns_dri_frac_prev <- pcawg_mutsig_ns_dri_frac %>% semi_join(cancer_mutsig, by = c('ICGC_abbr_top', 'Attri_Signature' = 'Signature'))

prev_sig_sort <- mutsig_ns_dri_frac_prev %>% filter(NS_driCount > 1) %>% distinct(Attri_Signature) %>% pull(Attri_Signature) %>% str_sort(numeric = TRUE)
cancer_type_sort <- mutsig_ns_dri_frac_prev %>% filter(NS_driCount > 1) %>% count(ICGC_abbr_top) %>% arrange(n) %>% pull(ICGC_abbr_top)

prev_plot_dat <- mutsig_ns_dri_frac_prev %>% 
  filter(NS_driCount > 1) %>% # each signature must contribute > 1 driver mutations
  mutate(ICGC_abbr_top = factor(ICGC_abbr_top, levels = cancer_type_sort), 
         Attri_Signature = factor(Attri_Signature, levels = prev_sig_sort)) %>% 
  arrange(ICGC_abbr_top, Attri_Signature)

# dot circos plot

## data preparation

selected_prev_plot_dat <- prev_plot_dat %>% 
  group_by(ICGC_abbr_top) %>% 
  mutate(NS_driFrac = NS_driCount/sum(NS_driCount)) %>% 
  ungroup()

# mutsig driver contribution

# write_tsv(selected_prev_plot_dat, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_driver_contribution.tsv')
# selected_prev_plot_dat <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_driver_contribution.tsv')

cancer_colors <- tribble(
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
  'Uterus-AdenoCA', '#F7DEC5'
)


mutsig_num <- n_distinct(selected_prev_plot_dat$Attri_Signature)

stripe_dat <- tibble(xst = c(seq(0.5, mutsig_num + 0.5, by = 1)), 
                     xed = c(xst[-1], Inf), 
                     reg_col = rep_len(c('odd', 'even'), length.out = length(xst))) %>% 
  filter(xed != Inf)


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_dri_dot_circos_plots.pdf')

ggplot() +
  geom_jitter(data = selected_prev_plot_dat, aes(Attri_Signature, NS_driFrac, col = ICGC_abbr_top, size = Dri_frac), width = 0.3) +
  geom_rect(data = stripe_dat, aes(xmin = xst, xmax = xed, ymin = -Inf, ymax = Inf, fill = reg_col), show.legend = FALSE) +
  scale_color_manual(breaks = cancer_colors$cancer_type, values = cancer_colors$color) +
  scale_fill_manual(breaks = c('odd', 'even'), values = c('#33333333', '#00000000')) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_reverse(limit = c(1, 0)) +
  coord_polar(theta = 'x') +
  labs(x = NULL, y = NULL, col = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())

dev.off()

# Total subclonal contribution of signatures ------------------------------

sub_ns_dri_data <- comb_mutsig_clo_ns_dri_data %>% filter(clonality %in% 'subclonal')
# subclonal: 13.4%(2556/19146)
# n_distinct(comb_mutsig_clo_data_lim$Attri_Signature): 51 signatures
# n_distinct(sub_ns_dri_data$Attri_Signature): 36 signatures
# 36/51 signatures contribute at least one subclonal driver mutation

sig_sub_frac <- sub_ns_dri_data %>% count(Attri_Signature) %>% mutate(frac = n/sum(n)) %>% arrange(desc(frac))

# barplot

count_label <- sub_ns_dri_data %>% count(Attri_Signature, name = 'count')

library(gg.gap)

sub_ns_dri_plot <- sub_ns_dri_data %>% 
  mutate(Attri_Signature = factor(Attri_Signature, levels = sig_sub_frac$Attri_Signature), 
         Variant_Classification = str_remove(Variant_Classification, '_Mutation')) %>% 
  ggplot() +
  geom_bar(aes(Attri_Signature, fill = Variant_Classification)) +
  geom_text(data = count_label, aes(Attri_Signature, count, label = count), vjust = -0.5, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, .03))) +
  scale_fill_manual(breaks = c('Missense', 'Nonsense', 'Splice_Site', 'Nonstop'), values = c('#699F40', '#5E398D', '#D26820', 'black')) + 
  labs(x = NULL, y = 'subclonal driver number', fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank())


pdf('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Revise_1/mutsig_driver_contribution_analysis/mutsig_subclonl_barplot.pdf', width = 8, height = 6)

gg.gap(plot = sub_ns_dri_plot,
       segments = list(c(180, 240), c(250, 400), c(440, 800)),
       ylim = c(0, 850),
       rel_heights = c(1, 0, 0.1, 0, 0.1, 0, 0.1))

dev.off()

# Fraction of clonal/subclonal driver mutation --------------------------------------

mutsig_dri_clo_frac <- comb_mutsig_clo_ns_dri_data %>% 
  semi_join(cancer_mutsig, by = c('ICGC_abbr_top', 'Attri_Signature' = 'Signature')) %>% 
  filter_at(c('Attri_Signature', 'clonality'), ~!is.na(.)) %>% 
  count(ICGC_abbr_top, Attri_Signature, clonality) %>% 
  pivot_wider(names_from = clonality, values_from = n, values_fill = 0) %>% 
  mutate(clo_frac = clonal/(clonal + subclonal), sub_frac = 1 - clo_frac)
# Ovary: SBS2 (3/3) clonal, SBS13 (2/2) subclonal

write_tsv(mutsig_dri_clo_frac, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_dri_clo_frac.tsv')

# moon chart

clo_frac_plot_sig_sort <- mutsig_dri_clo_frac %>% 
  distinct(Attri_Signature) %>% 
  mutate(Attri_Signature = str_sort(Attri_Signature, numeric = TRUE), 
         x_id = 1:length(Attri_Signature))

clo_frac_plot_cancer_sort <- mutsig_dri_clo_frac %>% 
  distinct(ICGC_abbr_top) %>% 
  mutate(ICGC_abbr_top = str_sort(ICGC_abbr_top, numeric = TRUE), 
         y_id = 1:length(ICGC_abbr_top))

mutsig_dri_clo_frac_plot_dat <- mutsig_dri_clo_frac %>% 
  mutate(total_count = clonal + subclonal) %>% 
  select(-(clonal:subclonal)) %>% 
  pivot_longer(clo_frac:sub_frac, names_to = 'Clonality', values_to = 'fraction') %>% 
  mutate(right = ifelse(Clonality %in% 'sub_frac', TRUE, FALSE)) %>% 
  left_join(clo_frac_plot_sig_sort, by = 'Attri_Signature') %>% 
  left_join(clo_frac_plot_cancer_sort, by = 'ICGC_abbr_top')


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_dri_clofrac_moon.pdf', width = 10, height = 7)

mutsig_dri_clo_frac_plot_dat %>%
  ggplot() +
  geom_moon(aes(x_id, y_id, ratio = fraction, right = right, fill = right, size = log2(total_count)), col = NA) +
  scale_fill_manual(breaks = c(FALSE, TRUE), values = c('#C9372E', '#5082AF')) +
  scale_x_continuous(expand = c(0.02, 0.02), breaks = clo_frac_plot_sig_sort$x_id, labels = clo_frac_plot_sig_sort$Attri_Signature) +
  scale_y_continuous(expand = c(0.02, 0.02), breaks = clo_frac_plot_cancer_sort$y_id, labels = clo_frac_plot_cancer_sort$ICGC_abbr_top) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank())

dev.off()


# Fraction of clonal/subclonal driver mutation (view of genes) --------------------------------------

signature_colors <- tribble(
  ~signature, ~color, 
  'SBS1', '#8F6504', 
  'SBS2', '#F17E27', 
  'SBS3', '#0373C6', 
  'SBS4', '#C1A72F', 
  'SBS5', '#927C64', 
  'SBS6', '#B2509E', 
  'SBS7a', '#9C8CC4', 
  'SBS7b', '#D49DC7', 
  'SBS7c', '#DEA0FD', 
  'SBS7d', '#0C487A', 
  'SBS8', '#90AD1C', 
  'SBS9', '#BAA133', 
  'SBS10a', '#E8C51D', 
  'SBS10b', '#FF5FB3', 
  'SBS11', '#F9ED32', 
  'SBS12', '#16FF32', 
  'SBS13', '#EB9B31', 
  'SBS14', '#ED2891', 
  'SBS15', '#CC6677', 
  'SBS16', '#F8D0D6', 
  'SBS17a', '#70CBEF', 
  'SBS17b', '#91E4A6', 
  'SBS18', '#C41707', 
  'SBS20', '#EE9028', 
  'SBS22', '#F7DEC5', 
  'SBS27', '#FF8C69', 
  'SBS28', '#E066FF', 
  'SBS29', '#008B8B', 
  'SBS30', '#525975', 
  'SBS36', '#104A7F', 
  'SBS37', '#2FA0D2', 
  'SBS38', '#B3CB47', 
  'SBS39', '#CAA98D', 
  'SBS40', '#FFB50E', 
  'SBS44', '#7A1A1D', 
  'SBS45', '#B10DA1'
)


# only retain cancer types which carry at least 10 subclonal mutations

cancer_type_retain <- comb_mutsig_clo_ns_dri_data %>% 
  filter(clonality %in% 'subclonal') %>% 
  count(ICGC_abbr_top) %>% 
  filter(n >= 10) %>% 
  arrange(desc(n))
# 23 cancer types

# deleted cancer types
# UCS                          7
# Myeloid-MDS/MPN              5
# CNS-PiloAstro                4
# Thy-AdenoCA                  4
# Eye-Melanoma                 3
# Mesothelium-Mesothelioma     3
# Thymoma                      3
# Kidney-ChRCC                 2
# Kidney-Papillary             2
# Adrenal-neoplasm             1
# Biliary-AdenoCA              1

# only retain genes which carry at least 20 subclonal mutations

gene_sub_retain <- comb_mutsig_clo_ns_dri_data %>% 
  filter(clonality %in% 'subclonal') %>% 
  semi_join(cancer_type_retain, by = 'ICGC_abbr_top') %>% # this condition removed one gene, 35 -> 34
  count(Hugo_Symbol) %>% 
  filter(n >= 20) %>% 
  arrange(desc(n))
# 34 genes

mutsig_gene_dri_clo_frac <- comb_mutsig_clo_ns_dri_data %>% 
  semi_join(gene_sub_retain, by = 'Hugo_Symbol') %>% 
  semi_join(cancer_type_retain, by = 'ICGC_abbr_top') %>% 
  semi_join(cancer_mutsig, by = c('ICGC_abbr_top', 'Attri_Signature' = 'Signature')) %>% 
  filter_at(c('Attri_Signature', 'clonality'), ~!is.na(.)) %>% 
  count(Hugo_Symbol, ICGC_abbr_top, clonality, Attri_Signature) %>% 
  group_by(Hugo_Symbol, ICGC_abbr_top, clonality) %>% 
  mutate(radius = sum(n), 
         frac = n/radius, 
         radius_scale = log2(radius)/6) %>% 
  ungroup()

# write_tsv(mutsig_gene_dri_clo_frac, '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_gene_dri_clo_frac.tsv')
# mutsig_gene_dri_clo_frac <- read_tsv('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_gene_dri_clo_frac.tsv')

sp_sig_sort <- mutsig_gene_dri_clo_frac %>% distinct(Attri_Signature) %>% pull(Attri_Signature) %>% str_sort(numeric = TRUE)

sp_gc_sort <- mutsig_gene_dri_clo_frac %>% 
  mutate(gc_comb = str_c(Hugo_Symbol, clonality, sep = '_'), text_col = ifelse(clonality %in% 'clonal', '#C9372E', '#5082AF')) %>% 
  distinct(Hugo_Symbol, clonality, gc_comb, text_col) %>% 
  arrange(gc_comb) %>% 
  rowid_to_column('x_idx') %>% 
  mutate(x_idx = x_idx + lag(x_idx, default = 0)) # index: 1, 3, 5, 7, ...

sp_cancer_sort <- mutsig_gene_dri_clo_frac %>% 
  distinct(ICGC_abbr_top) %>% 
  arrange(ICGC_abbr_top) %>% 
  rowid_to_column('y_idx') %>% 
  mutate(y_idx = y_idx + lag(y_idx, default = 0))

mutsig_gene_scatterpie_dat <- mutsig_gene_dri_clo_frac %>% 
  mutate(gc_comb = str_c(Hugo_Symbol, clonality, sep = '_')) %>% 
  left_join(sp_gc_sort, by = 'gc_comb') %>% 
  left_join(sp_cancer_sort, by = 'ICGC_abbr_top') %>% 
  select(-n) %>% 
  pivot_wider(names_from = Attri_Signature, values_from = frac, values_fill = 0)

radius_legend <- seq(min(mutsig_gene_dri_clo_frac$radius_scale), max(mutsig_gene_dri_clo_frac$radius_scale), length.out = 5)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_gene_driver_frac.pdf', width = 15, height = 8)

ggplot() +
  geom_scatterpie(aes(x_idx, y_idx, r = radius_scale), cols = sp_sig_sort, data = mutsig_gene_scatterpie_dat) +
  geom_scatterpie_legend(radius_legend, x = 29, y = 31, labeller = function(x) round(2^(6*x))) + # labeller convert back to log2(radius)/6
  coord_equal() +
  scale_x_continuous(expand = c(0.01, 0.01), breaks = sp_gc_sort$x_idx, labels = sp_gc_sort$Hugo_Symbol) +
  scale_y_continuous(expand = c(0.01, 0.01), breaks = sp_cancer_sort$y_idx, labels = sp_cancer_sort$ICGC_abbr_top) +
  labs(x = NULL, y = NULL, fill = NULL) +
  # scale_fill_discrete(guide = guide_legend(nrow = 3)) +
  scale_fill_manual(breaks = signature_colors$signature, values = signature_colors$color, guide = guide_legend(nrow = 3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = sp_gc_sort$text_col),
        panel.grid.major = element_blank(),
        legend.position = 'top')

dev.off()

# Subclonal selection among mutational signatures -------------------------

library(dndscv)

# total

dndscv_prep_data <- comb_mutsig_clo_data_lim %>% 
  semi_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top')) %>% 
  filter(!is.na(Attri_Signature), !is.na(clonality)) %>% 
  select(Attri_Signature, clonality, mut_id) %>% 
  separate(mut_id, into = c('sampleID', 'chr', 'pos', 'ref', 'mut'), sep = ':')

mutsig_select_tres <- dndscv_prep_data %>% 
  group_by(Attri_Signature, clonality) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(global_dnds = map(data, possibly(~dndscv(.)$globaldnds, NULL)))

mutsig_select_tdnds <- mutsig_select_tres %>% 
  select(-data) %>% 
  unnest(global_dnds) %>% 
  # filter(cihigh != Inf, name != 'wtru', cihigh < 2000) %>% 
  filter(name != 'wtru') %>% 
  mutate(name = case_when(
    name %in% 'wmis' ~ 'Missense', 
    name %in% 'wnon' ~ 'Nonsense', 
    name %in% 'wspl' ~ 'Splice site', 
    name %in% 'wall' ~ 'All'
  ), status = ifelse(cilow > 1, name, 'n.s.'))

save(mutsig_select_tdnds, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mutsig_select_tdnds.Rdata')
# load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mutsig_select_tdnds.Rdata')

# add NA_val to distingush not significant value () and not available value

mutsig_select_tdnds_trans <- mutsig_select_tdnds %>%
  select(-(mle:cihigh)) %>%
  pivot_wider(names_from = Attri_Signature, values_from = status, values_fill = 'NA_val') %>% 
  pivot_longer(-(clonality:name), names_to = 'Attri_Signature', values_to = 'status')

# heatmap

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_selection_heatmap.pdf', width = 9, height = 2)

mutsig_select_tdnds_trans %>% 
  mutate(Attri_Signature = factor(Attri_Signature, levels = all_mutsig_sort), 
         name = factor(name, levels = c('Missense', 'Nonsense', 'Splice site', 'All'))) %>% 
  ggplot(aes(Attri_Signature, name)) +
  geom_tile(aes(fill = status), color = 'white') +
  labs(x = NULL, y = NULL, fill = 'dN/dS>1') +
  scale_fill_manual(breaks = c('Missense', 'Nonsense', 'Splice site', 'All', 'n.s.', 'NA_val'), values = c('#699F40', '#5E398D', '#D26820', '#5787A5', 'white', '#858585')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggthemes::theme_few() +
  facet_grid(clonality ~ ., switch = 'y') +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text.y = element_blank(), 
        axis.text = element_text(angle = 45, hjust = 1),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12))

dev.off()

# errorbar plots

sin_errobar_plot <- function(which_signature){
  
  sin_plot_dat <- mutsig_select_tdnds %>% 
    filter(Attri_Signature %in% which_signature) %>% 
    mutate(name = factor(name, levels = c('Missense', 'Nonsense', 'Splice site', 'All')))
  
  sin_plot_dat %>% 
    ggplot(aes(name, mle, ymin = cilow, ymax = cihigh, col = name)) +
    geom_point(size = 2) +
    geom_errorbar(size = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = NULL, y = 'dN/dS', col = NULL, title = which_signature) +
    scale_color_manual(breaks = c('Missense', 'Nonsense', 'Splice site', 'All'), values = c('#699F40', '#5E398D', '#D26820', '#5787A5')) +
    facet_wrap(~clonality, scales='free_x', strip.position = 'bottom') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.spacing.x = unit(10,'pt'), 
          axis.text.x = element_blank(), 
          strip.placement = 'outside', 
          strip.background = element_blank(), 
          axis.ticks.x = element_blank())
  
}

SBS3_errorbar_plot <- sin_errobar_plot('SBS3')
SBS5_errorbar_plot <- sin_errobar_plot('SBS5')
SBS40_errorbar_plot <- sin_errobar_plot('SBS40')
SBS44_errorbar_plot <- sin_errobar_plot('SBS44')

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_selection_errorbar.pdf', width = 10, height = 3)

wrap_plots(list(SBS3_errorbar_plot, SBS5_errorbar_plot, SBS40_errorbar_plot, SBS44_errorbar_plot)) + plot_layout(guides = 'collect', nrow = 1)

dev.off()


# individual cancer types

dndscv_prep_data2 <- comb_mutsig_clo_data_lim %>% 
  semi_join(comb_driver_gene, by = c('Hugo_Symbol', 'ICGC_abbr_top')) %>% 
  filter(!is.na(Attri_Signature), !is.na(clonality)) %>% 
  select(ICGC_abbr_top, Attri_Signature, clonality, mut_id) %>% 
  separate(mut_id, into = c('sampleID', 'chr', 'pos', 'ref', 'mut'), sep = ':')

mutsig_select_ires <- dndscv_prep_data2 %>% 
  group_by(ICGC_abbr_top, Attri_Signature, clonality) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(global_dnds = map(data, possibly(~dndscv(.)$globaldnds, NULL)))

mutsig_select_idnds <- mutsig_select_ires %>% 
  select(-data) %>% 
  unnest(global_dnds) %>% 
  # filter(cihigh != Inf, name != 'wtru', cihigh < 2000) %>% 
  filter(name != 'wtru') %>% 
  mutate(name = case_when(
    name %in% 'wmis' ~ 'Missense', 
    name %in% 'wnon' ~ 'Nonsense', 
    name %in% 'wspl' ~ 'Splice site', 
    name %in% 'wall' ~ 'All'
  ), status = ifelse(cilow > 1, name, 'n.s.'))

# save(mutsig_select_idnds, file = '/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/mutsig_select_idnds.Rdata')
# 27 cancer types, 32 signatures

mutsig_select_sig_sort <- mutsig_select_idnds %>% distinct(Attri_Signature) %>% pull(Attri_Signature) %>% str_sort(numeric = TRUE)

# heatmap


# if a cancer type lack clonal or subclonal rows, add info to prevent the apperance of NA values in geom_tile ploting
# case 1: for example, Liver-HCC only has clonal result, we must add NA_val to subclonal result

mutsig_sup_data <- mutsig_select_idnds %>% 
  select(-(mle:status)) %>%
  mutate(clo_name = str_c(clonality, name, sep = ':'), clo_name_status = 1) %>% 
  select(-(clonality:name)) %>% 
  pivot_wider(names_from = clo_name, values_from = clo_name_status, values_fill = 0) %>% 
  pivot_longer(-(ICGC_abbr_top:Attri_Signature), names_to = 'clo_name', values_to = 'clo_name_status') %>% 
  select(-clo_name_status) %>% 
  separate(clo_name, into = c('clonality', 'name'), sep = ':')


# add NA_val to distingush not significant value () and not available value
# case 2: for example, Liver-HCC clonal only has Missense and Nonsense result, we must add NA-val to Splice and All result

mutsig_select_idnds_trans <- mutsig_sup_data %>%
  left_join(mutsig_select_idnds, by = c('ICGC_abbr_top', 'Attri_Signature', 'clonality', 'name')) %>% 
  mutate(status = ifelse(is.na(status), 'NA_val', status)) %>% 
  select(-(mle:cihigh)) %>%
  pivot_wider(names_from = Attri_Signature, values_from = status, values_fill = 'NA_val') %>% 
  pivot_longer(-(ICGC_abbr_top:name), names_to = 'Attri_Signature', values_to = 'status')


icancer_heatmap <- function(mutation_type){
  
  # mutation_type <- 'All'
  
  plot_dat <- mutsig_select_idnds_trans %>% 
    filter(name %in% mutation_type) %>% 
    mutate(Attri_Signature = factor(Attri_Signature, levels = mutsig_select_sig_sort), 
           clonality = factor(clonality, levels = c('subclonal', 'clonal')), 
           status_code = case_when(
             status %in% 'NA_val' ~ 0,
             status %in% 'n.s.' ~ 1,
             TRUE ~ 2))
  
  # for clearer display, delete signatures/cancer types which are non-significant in either clonal or subclonal and meanwhile the other is not available across all cancer types/signatures
  
  retained_dat <- plot_dat %>%
    group_by(Attri_Signature, ICGC_abbr_top) %>%
    filter(sum(status_code) > 1) %>%
    ungroup() %>%
    distinct(ICGC_abbr_top, Attri_Signature)
  
  plot_dat %>% 
    semi_join(retained_dat, by = 'Attri_Signature') %>%
    semi_join(retained_dat, by = 'ICGC_abbr_top') %>%
    ggplot(aes(Attri_Signature, clonality)) +
    geom_tile(aes(fill = status), color = 'white') +
    labs(x = NULL, y = NULL, fill = 'dN/dS>1') +
    scale_fill_manual(breaks = c('Missense', 'Nonsense', 'Splice site', 'All', 'n.s.', 'NA_val'), values = c('#699F40', '#5E398D', '#D26820', '#5787A5', 'white', '#858585')) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), position = 'right') +
    ggthemes::theme_few() +
    facet_grid(ICGC_abbr_top ~ ., switch = 'y') +
    theme(axis.ticks = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12))

}

All_icancer_heatmap <- icancer_heatmap('All')

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_selection_heatmap_icancer_all.pdf', width = 7, height = 8)

All_icancer_heatmap

dev.off()


# supplemetal figures

icancer_heatmap_sup <- function(mutation_type){
  
  # mutation_type <- 'Missense'
  
  plot_dat <- mutsig_select_idnds_trans %>% 
    filter(name %in% mutation_type) %>% 
    mutate(Attri_Signature = factor(Attri_Signature, levels = mutsig_select_sig_sort), 
           clonality = ifelse(clonality %in% 'clonal', 'C', 'S'), 
           status_code = case_when(
             status %in% 'NA_val' ~ 0,
             status %in% 'n.s.' ~ 1,
             TRUE ~ 2))
  
  # for clearer display, delete signatures/cancer types which are non-significant in either clonal or subclonal and meanwhile the other is not available across all cancer types/signatures
  
  retained_dat <- plot_dat %>%
    group_by(Attri_Signature, ICGC_abbr_top) %>%
    filter(sum(status_code) > 1) %>%
    ungroup() %>%
    distinct(ICGC_abbr_top, Attri_Signature)
  
  # pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_selection_heatmap_icancer_missense.pdf', width = 12, height = 5)
  
  plot_dat %>% 
    semi_join(retained_dat, by = 'Attri_Signature') %>%
    semi_join(retained_dat, by = 'ICGC_abbr_top') %>%
    ggplot(aes(clonality, Attri_Signature)) +
    geom_tile(aes(fill = status), color = 'white') +
    labs(x = NULL, y = NULL, fill = 'dN/dS>1') +
    scale_fill_manual(breaks = c('Missense', 'Nonsense', 'Splice site', 'All', 'n.s.', 'NA_val'), values = c('#699F40', '#5E398D', '#D26820', '#5787A5', 'white', '#858585')) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    ggthemes::theme_few() +
    facet_grid2(. ~ ICGC_abbr_top, strip = strip_vanilla(clip = 'off')) +
    theme(axis.ticks = element_blank(),
          strip.text = element_text(angle = 45, hjust = 0.3, vjust = 0.3), 
          strip.background = element_blank())
  
  # dev.off()
  
}

Missense_icancer_heatmap <- icancer_heatmap_sup('Missense')
Nonsense_icancer_heatmap <- icancer_heatmap_sup('Nonsense')
Splice_icancer_heatmap <- icancer_heatmap_sup('Splice site')


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_driver_clonality_analysis/mutsig_selection_heatmap_icancer_other_type.pdf', width = 12, height = 16)

Missense_icancer_heatmap + Nonsense_icancer_heatmap + Splice_icancer_heatmap + plot_layout(ncol = 1)

dev.off()


