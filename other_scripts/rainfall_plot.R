
library(data.table)

# transform data from segment level to whole genome chromosome level
#--- Change segment sizes into linear scale
transformSegments = function(segmentedData, build = 'hg19'){
  
  build.opts = c('hg19', 'hg18', 'hg38')
  
  if(!build %in% build.opts){
    stop('Available reference builds: hg18, hg19, hg38')
  }
  
  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){ #hg38
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  } else{
    stop('Available reference builds: hg18, hg19, hg38')
  }
  
  seg.spl = split(segmentedData, segmentedData$Chromosome)
  
  seg.spl.transformed = seg.spl[[1]]
  if(nrow(seg.spl.transformed) > 0){
    seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
  }
  
  chr.lens.sumsum = cumsum(chr.lens)
  
  for(i in 2:length(seg.spl)){
    
    x.seg = seg.spl[[i]]
    if(nrow(x.seg) > 0){
      x.seg$Start_Position_updated = x.seg$Start_Position + chr.lens.sumsum[i-1]
      x.seg$End_Position_updated = x.seg$End_Position + chr.lens.sumsum[i-1]
    }
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }
  
  return(seg.spl.transformed)
}


hg19_chr_lens <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                   146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                   102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                   155270560, 59373566)
hg19_lens_cumsum <- cumsum(hg19_chr_lens)


sam_rainfall_plot_dat <- function(sample){
  
  # sample <- 'SP116478'
  
  sin_sam_dat <- comb_mutsig_mut_dat %>% filter(Sample %in% sample)
  
  sin_mut <- sin_sam_dat %>% 
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
  
  sin_trans_mut <- transformSegments(sin_mut)
  sin_trans_mut$diff <- suppressWarnings( log10(c(0, diff(sin_trans_mut$Start_Position_updated))+1) )
  sin_trans_mut <- sin_trans_mut[complete.cases(diff)]
  
}


sam_rainfall_plot <- function(plot_dat){
  
  # base_breaks_x <- function(x){
  #   b <- pretty(x)
  #   d <- data.frame(y=-Inf, yend=-Inf, x=min(b), xend=max(b))
  #   list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE),
  #        scale_x_continuous(breaks=b))
  # }
  
  # plot_dat <- SP116478_rainplot_dat
  
  base_breaks_y <- function(x){
    b <- pretty(x)
    d <- data.frame(x=-Inf, xend=-Inf, y=min(b), yend=max(b))
    list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE),
         scale_y_continuous(breaks=b, expand = c(0, 0)))
  }
  
  # pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/rainfall_test.pdf', width = 12, height = 3)
  
  plot_dat %>% 
    ggplot(aes(Start_Position_updated, diff, col = con.class)) +
    geom_vline(xintercept = hg19_lens_cumsum, size = 0.25, color = 'grey60') +
    geom_point(shape = '.') +
    labs(x = NULL, y = 'log10(inter event distance)', col = NULL, title = plot_dat$Tumor_Sample_Barcode[1]) +
    scale_x_continuous(breaks = hg19_lens_cumsum, labels = c(1:22, 'X', 'Y')) +
    scale_color_manual(breaks = rainfall_sig_colors$Signature, values = rainfall_sig_colors$color) +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    base_breaks_y(plot_dat$diff)
  
  # dev.off()
  
}

# rainfall_sig_colors <- tribble(
#   ~Signature, ~color, 
#   'SBS1', '#BAA133', 
#   'SBS2', '#0E78AA', 
#   'SBS5', '#DCBE29', 
#   'SBS7a', '#2FA0D2', 
#   'SBS7b', '#9DD6F0', 
#   'SBS7d', '#D8EEFB',
#   'SBS13', '#374F99', 
#   'SBS18', '#C7C9D7', 
#   'SBS38', '#CAA98D', 
#   'SBS40', '#D92F86', 
# )


rainfall_sig_colors <- tribble(
  ~Signature, ~color, 
  'SBS1', '#D8EEFB', 
  'SBS2', '#0E78AA', 
  'SBS3', '#DCBE29', 
  'SBS5', '#C7C9D7', 
  'SBS13', '#374F99'
)

SP116478_rainplot_dat <- sam_rainfall_plot_dat('SP116478')
selected_region <- SP116478_rainplot_dat %>% as_tibble() %>% filter(Chromosome == 12, Start_Position >= 86372801, Start_Position <= 86388330)

pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/SP116478_rainfall.pdf', width = 13, height = 4)

sam_rainfall_plot(SP116478_rainplot_dat) + annotate('rect', xmin = 1950914406 + 86372801, xmax = 1950914406 + 86388330, ymin = min(selected_region$diff), ymax = max(selected_region$diff), fill = '#EC7014')

dev.off()


pdf('/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_mutation_analysis/SP116478_rainfall_amplified.pdf', width = 13, height = 4)

sam_rainfall_plot(SP116478_rainplot_dat) + coord_cartesian(xlim = c(1950914406, 2084766301)) + annotate('rect', xmin = 1950914406 + 86372801, xmax = 1950914406 + 86388330, ymin = min(selected_region$diff), ymax = max(selected_region$diff), fill = '#EC7014')

dev.off()
