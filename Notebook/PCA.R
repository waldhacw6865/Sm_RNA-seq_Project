
counts_df <- read_tsv("/data/users/waldhacw6865/Sm_RNA-seq_Project/Notebook/counts/star_counts.tsv",
                      comment = "#") |>
             mutate(across(where(is.numeric), as.integer))

counts_summary <- counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "/data/users/waldhacw6865/Sm_RNA-seq_Project/Notebook/dedup/star.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)

sample_summary <- counts_df |>
    select(Geneid, contains('/data/users/waldhacw6865/Sm_RNA-seq_Project/Notebook/dedup/star.bam:')) |>
    rename_with(~str_remove(., "/data/users/waldhacw6865/Sm_RNA-seq_Project/Notebook/dedup/star.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)

genes_to_remove = sample_summary$Geneid

counts_filt <- counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)
counts_filt

counts_m <- counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(counts_m) <- counts_filt$Geneid
counts_m

dists <- dist(t(counts_m))
dists

dists_df <- as.matrix(dists) |>
    as_tibble(rownames = 'sample')

dist_plot <- dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
dist_plot

pca_fit <- t(log10(counts_m + 1)) |> 
  prcomp(scale = TRUE)
pca_fit

library(broom)

pca_fit |>
  augment(t(counts_m)) |>
  dplyr::rename(sample = .rownames) |>
  separate(sample, into = c('tissue', 'age'), sep = '_') |>
  mutate(age = str_remove(age, '[0-9]')) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = tissue, shape = age)) + 
  geom_point(size = 4)

  metadata <- data.frame(sample_id = colnames(counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           age = str_sub(sample_id, 5, 6),
           rep = str_sub(sample_id, 7))
rownames(metadata) <- metadata$sample_id
metadata <- select(metadata, -sample_id)
metadata

all(rownames(metadata) == colnames(counts_m))

dds <- DESeqDataSetFromMatrix(countData = counts_m,
                              colData = metadata,
                              design = ~ tissue)
dds <- DESeq(dds)
dds

res <- results(dds)
res

as_tibble(res, rownames = 'gene_id') |>
    filter(gene_id == 'Smp_329140')

filter(counts_filt, Geneid == 'Smp_329140')

as_tibble(res, rownames = 'gene_id') |>
    filter(gene_id == 'Smp_333930')

metadata2 <- metadata |>
    mutate(tissue_age = str_c(tissue, "_", age)) |>
    select(tissue_age, rep)
metadata2

dds2 <- DESeqDataSetFromMatrix(countData = counts_m,
                               colData = metadata2,
                               design = ~ tissue_age)
dds2 <- DESeq(dds2)
dds2

LIV_ma_vs_INT_ma <- results(dds2, contrast = c('tissue_age', 'LIV_ma', 'INT_ma'))
LIV_ma_vs_INT_ma

as_tibble(LIV_ma_vs_INT_ma, rownames = 'gene_id') |>
    filter(gene_id == 'Smp_333930')

volcano_data <- as_tibble(LIV_ma_vs_INT_ma, rownames = 'gene_id')

volcano_plot <- volcano_data |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2) +
    theme_minimal()
volcano_plot

vsd <- vst(dds)
vsd

plotPCA(vsd, intgroup = c("tissue", "age"))

vsd_dists <- dist(t(assay(vsd)))

vsd_dists_df <- as.matrix(vsd_dists) |>
    as_tibble(rownames = 'sample')

vsd_dist_plot <- vsd_dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
vsd_dist_plot


