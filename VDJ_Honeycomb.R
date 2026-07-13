# From scRepertoir contig.list
# Honeycomb from AIRR and rds ---------------------
# Get contigs
library(scRepertoire)
library(dplyr)
library(ggplot2)
library(ggforce)
library(Seurat)

# contig.list already loaded via loadContigs(format="BD")
combined <- combineTCR(
  contig.list,
  samples = samples,
  ID = samples,
  cells = "T-AB"
)

# Extract clonotype table
clones <- combinedTCR[[1]] %>%
  group_by(clone_id = CTaa) %>%
  summarise(
    size = n(),
    v = first(v_call),
    j = first(j_call),
    cdr3 = first(CDR3.aa),
    sample = first(sample)
  ) %>%
  filter(!is.na(clone_id))

# Assign honeycomb coordinates
clones <- clones %>%
  mutate(
    row = floor((row_number() - 1) / 20),
    col = (row_number() - 1) %% 20,
    x = col + ifelse(row %% 2 == 0, 0, 0.5),
    y = -row
  )

# Plot
ggplot(clones, aes(x, y)) +
  geom_regon(
    aes(fill = sample, r = sqrt(size) * 0.1, sides = 6, angle = 0),
    alpha = 0.9,
    color = "black",
    size = 0.2
  ) +
  scale_fill_brewer(palette = "Set2") +
  coord_equal() +
  theme_void() +
  ggtitle("Honeycomb-style clonotype plot (BD AIRR)")


# UMAP Honeycomb ---------------------------
library(dplyr)

umap <- Embeddings(seu, "umap") %>% as.data.frame()
umap$cell_id <- rownames(umap)

meta <- seu@meta.data %>%
  select(cell_id = barcode, clonotype = CTaa, sample)

df <- umap %>%
  left_join(meta, by = "cell_id") %>%
  filter(!is.na(clonotype))

centroids <- df %>%
  group_by(clonotype) %>%
  summarise(
    x = mean(UMAP_1),
    y = mean(UMAP_2),
    size = n(),
    sample = first(sample)
  )

library(ggforce)
library(ggplot2)

hex_plot <- ggplot(centroids, aes(x, y)) +
  geom_regon(
    aes(
      r = sqrt(size) * 0.05,   # scale hex size by clone size
      fill = sample
    ),
    sides = 6,
    angle = 0,
    color = "black",
    alpha = 0.9,
    size = 0.2
  ) +
  scale_fill_brewer(palette = "Set2") +
  coord_equal() +
  theme_void() +
  ggtitle("UMAP‑anchored honeycomb clonotype plot")


umap_points <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey80", size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_void()

umap_points + 
  geom_regon(
    data = centroids,
    aes(x, y, r = sqrt(size) * 0.05, fill = sample),
    sides = 6,
    angle = 0,
    color = "black",
    alpha = 0.9,
    size = 0.2
  ) +
  scale_fill_brewer(palette = "Set2") +
  ggtitle("UMAP + Honeycomb Clonotypes")

