
library(intervals)

genecode_v19_anno <- read_tsv('/boot3/bio_liaojl/Common_data/processed_data/genecode_v19_anno.tsv')

# exon length -------------------------------------------------------------

genecode_v19_exon <- genecode_v19_anno %>% 
  filter(gene_type %in% 'protein_coding', transcript_type %in% 'protein_coding', type %in% 'exon', seqnames %in% str_c('chr', c(1:22, 'X', 'Y'))) %>% 
  select(chrom = gene_id, exonStart = start, exonEnd = end, strand) %>% 
  group_by(chrom) %>% 
  filter(length(unique(strand)) == 1) %>% # exclude genes of which exons located at both strands within or among transcripts; no genes were excluded
  nest() %>% 
  ungroup()
# some genes are protein_coding genes (e.g., ENSG00000255275.3), but they don't have protein coding transcripts


# merge all transcripts' exons

genecode_v19_exon_comb <- genecode_v19_exon %>% 
  mutate(comb_res = map(data, ~as_tibble(as.matrix(interval_union(Intervals(.)))))) %>% 
  select(-data) %>% 
  unnest(comb_res) %>% 
  rename(gene_id = chrom, start = V1, end = V2) %>% 
  mutate(width = end - start)
# sum(genecode_v19_exon_comb$width), >76M; 76299328 bp

# other methods
# retain transcript with maximum width, > 67M
# retain transcript with minimum width, > 22M


# intron length -----------------------------------------------------------

# for each gene, [(end of last exon) - (start of first exon)] - exon width


genecode_v19_exon_stend <- genecode_v19_exon_comb %>% 
  group_by(gene_id) %>% 
  summarise(start_min = min(start), end_max = max(end)) %>% 
  mutate(stend_width = end_max - start_min)
# 20,156 genes

gene_exon_width <- genecode_v19_exon_comb %>% 
  group_by(gene_id) %>% 
  summarise(exon_full_len = sum(width))

intron_len <- genecode_v19_exon_stend %>% 
  inner_join(gene_exon_width, by = 'gene_id') %>% 
  mutate(intron_width = stend_width - exon_full_len)

# sum(intron_len$intron_width), > 1192M; 1192307199 bp


# UTR length --------------------------------------------------------------

# for each gene, gene length - [(end of last exon) - (start of first exon)]


# exclude genes of which exons located at both strands within or among transcripts; though no genes were excluded

genecode_v19_protein_coding_gene <- genecode_v19_anno %>% 
  #  separate(gene_id, into = c('gene_id', NA), sep = '\\.') %>% 
  filter(type %in% 'gene', gene_type %in% 'protein_coding') %>% 
  select(chrom = seqnames, start:strand, Gene_ID = gene_id, Hugo_Symbol = gene_name)

# write_tsv(genecode_v19_protein_coding_gene, '/boot3/bio_liaojl/Common_data/processed_data/genecode_v19_protein_coding_gene.tsv')

utr_len <- genecode_v19_protein_coding_gene %>% 
  select(gene_len = width, Gene_ID) %>% 
  inner_join(intron_len, by = c('Gene_ID' = 'gene_id')) %>% 
  mutate(UTR_width = gene_len - stend_width) %>% 
  select(-(start_min:end_max))

# sum(utr_len$UTR_width), > 49M; 49558686 bp


# non-genic region length -------------------------------------------------

# chromosome length

# https://en.wikipedia.org/wiki/Human_genome
# https://archive.ph/6kkzz

# chrom_len <- tibble(chrom = str_c('chr', c(1:22, 'X', 'Y')), 
#                     size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 
#                              141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 
#                              81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566))
# total: 3095677412


# https://grch37.ensembl.org/Homo_sapiens/Location/Genome
# total: 3098825702

# gene length
# sum(utr_len$gene_len), > 1318M; 1318165213 bp

# non-genic length
# 3098825702 - 1318165213 = 1780660489; > 1780M


