library(tidyverse)
library(trackViewer)
library(maftools)

# Bone-Osteosarc, SP116478, MGAT4C
load('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/MGAT4C_HGSV_anno.Rdata')
pro_struc_dat_ori <- readRDS(system.file('extdata', 'protein_domains.RDs', package = 'maftools'))

rainfall_sig_colors <- tribble(
  ~Signature, ~color, 
  'SBS1', '#D8EEFB', 
  'SBS2', '#0E78AA', 
  'SBS3', '#DCBE29', 
  'SBS5', '#C7C9D7', 
  'SBS13', '#374F99'
)

MGAT4C_prot_struc <- pro_struc_dat_ori %>% filter(HGNC %in% 'MGAT4C')
# only one transcript: NM_013244

MGAT4C_AA_len <- unique(MGAT4C_prot_struc$aa.length)
MGAT4C_mut_hgvs <- MGAT4C_HGSV_anno %>% 
  filter(Refseq_ID %in% 'NM_013244') %>% 
  mutate(aa_pos = as.numeric(str_extract(HGVS_short, '[0-9]+'))) %>% 
  left_join(rainfall_sig_colors, by = 'Signature')

# AA change setting

aa_change_pos <- GRanges('chr12', IRanges(MGAT4C_mut_hgvs$aa_pos, width = 1, names = MGAT4C_mut_hgvs$HGVS_short))
aa_change_pos$color <- MGAT4C_mut_hgvs$color
aa_change_pos$border <- MGAT4C_mut_hgvs$color
aa_change_pos$label.parameter.rot <- 60
aa_change_pos$label.parameter.gp <- gpar(col = MGAT4C_mut_hgvs$color)
### Attention!!! find it hard to change label color ideally, need to manually change in AI

# domain setting

domain_pos <- GRanges('chr12', IRanges(MGAT4C_prot_struc$Start, MGAT4C_prot_struc$End, names = MGAT4C_prot_struc$Label))
domain_pos$height <- 0.08
domain_pos$fill <- '#FF8833'


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/SP116478_MGAT4C_lollipop_plot.pdf', width = 10, height = 3)

lolliplot(aa_change_pos, domain_pos, ranges = GRanges('chr12', IRanges(0, MGAT4C_AA_len)))

dev.off()





