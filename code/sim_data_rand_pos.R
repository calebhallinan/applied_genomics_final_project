library(Seurat)
library(tidyverse)
library(sf)
library(anndata)
set.seed(42)



# Read data ---------------------------------------------------------------

## read sc counts from MS1 (there are three)
sc_mtx <- Read10X_h5('original_data/single_cell/Kidney_MS_1.h5')

## read srt counts
srt_mtx <- Matrix::readMM('original_data/male_cartana/male_ctx.mtx')
rownames(srt_mtx) <- read.csv('original_data/male_cartana/male_genes.csv', 
                              header = F) %>% 
  pull(V1)
colnames(srt_mtx) <- read.csv('original_data/male_cartana/male_cells.csv', 
                              header = F) %>% 
  pull(V1)
## read srt metadata
srt_meta <- read.csv('original_data/male_cartana/male_meta.csv') %>% 
  select(!X)



# Explore data ------------------------------------------------------------

## check n of cells
dim(sc_mtx)
Matrix::colSums(sc_mtx) %>% hist()

## check correlation of area and transcripts
cor_val <- cor.test(srt_meta$area, 
                    srt_meta$n_transcripts, use = 'complete.obs')
srt_meta %>% 
  ggplot() +
  geom_point(aes(x = area, y = n_transcripts), size = .01) +
  geom_smooth(aes(x = area, y = n_transcripts), method = "lm") +
  labs(title = paste('cor', round(cor_val$estimate, 2), 
                     'p-val', cor_val$p.value))


# Assign locations --------------------------------------------------------

## Select subset ---------------------------------------------------------

## select smaller sample to be compatible with the number of scs
sub_srt_meta <- srt_meta %>% 
  filter(x < 15000, y < 20000,
         x > 5000, y > 13000)

## check area nas, 134 cells
sum(is.na(sub_srt_meta$area))
## check if it is evenly spread
sub_srt_meta[is.na(sub_srt_meta$area), ] %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .01) +
  scale_color_manual(values = rainbow(14)) +
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  labs(title = 'Cartana Male Sham') + 
  theme_minimal()
## remove area nas
sub_srt_meta <- sub_srt_meta %>% 
  drop_na(area)

## check area sizes, some cells have low values
sub_srt_meta$area %>% hist(breaks = 100)
sub_srt_meta %>% 
  filter(area >= 10) %>% 
  pull(area) %>% 
  hist(breaks = 100)
## 18 cells
sum(sub_srt_meta$area < 10)
sub_srt_meta %>% 
  filter(area < 10) %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .01) +
  scale_color_manual(values = rainbow(14)) +
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  labs(title = 'Cartana Male Sham') + 
  theme_minimal()
## remove small areas
sub_srt_meta <- sub_srt_meta %>% 
  filter(area >= 10)

## plot data
sub_srt_meta %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .01) +
  scale_color_manual(values = rainbow(14)) +
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  labs(title = 'Cartana Male Sham') + 
  theme_minimal()


## Assign positions ------------------------------------------------------

## number of cells in srt subset
sub_srt_ncell <- dim(sub_srt_meta)[1]

## select random cells from sc data
sc_cells <- colnames(sc_mtx)
cells <- sample(sc_cells, sub_srt_ncell)
mtx <- sc_mtx[, cells]
## save mtx
# adata <- AnnData(X = mtx, 
#                  obs = data.frame(cell_id = colnames(mtx)), 
#                  var = data.frame(gene_id = rownames(mtx)))
# write_h5ad(ad, "sim_data/ad.h5ad")
saveRDS(mtx, 'sim_data/mtx.RDS')

## add cell names from the sc data to the srt metadata in order by area
meta <- sub_srt_meta %>% 
  arrange(area) %>% 
  mutate(cell = names(sort(colSums(mtx))))
write.csv(meta, 'sim_data/sc_srt_meta.csv')



## Filter by srt genes ---------------------------------------------------

srt_genes <- rownames(srt_mtx)
sc_genes <- rownames(sc_mtx)
common_genes <- intersect(srt_genes, sc_genes)
saveRDS(filtered_srt_mtx, 'sim_data/filtered_srt_mtx.RDS')


## Create cell boundaries ------------------------------------------------

## create circles representing the cell boundaries
circles <- lapply(1:dim(meta)[1], 
                  FUN = function(i) {
                    cell <- meta[i, ]
                    st_buffer(st_point(c(cell$x, cell$y)), 
                              dist = sqrt(cell$area / pi))
                  })
names(circles) <- meta$cell
meta <- meta %>% 
  mutate(boundary = circles)

circles_sf <- st_sfc(circles) %>% st_sf()
ggplot(data = circles_sf) +
  geom_sf(fill = "lightblue", color = "blue") +
  theme_minimal() +
  labs(title = "Cell areas")
## some cells overlap

## Simulate diffusion 10 ----------------------------------------------------

### Simulate transcript positions ----------------------------------------

## percentage out of the cell
pct_out <- .1
## distance from the cell centroid as a function of area
min_dist <- 1.01
max_dist <- 2

## empty matrix to get the transcripts later
ld10_mtx <- sparseMatrix(dim(filtered_srt_mtx)[1], 
                         dim(filtered_srt_mtx)[2], 
                         x = 0)
rownames(ld10_mtx) <- rownames(filtered_srt_mtx)
colnames(ld10_mtx) <- colnames(filtered_srt_mtx)
## create df with transcript coordinates represented by the cell center
transcripts <- data.frame()
## for every cell
for (i in 1:dim(meta)[1]) {
  print(paste(i, dim(meta)[1]))
  cell <- meta[i, ]
  cell_id <- cell$cell
  ## get expressed transcripts 
  ts <- filtered_srt_mtx[which(filtered_srt_mtx[, cell_id] > 0), cell_id]
  ## calculate transcripts outside and inside the cell
  ts_out <- round(ts * pct_out)
  ts_in <- ts - ts_out
  ## populate exp matrix with the transcripts inside the cell
  for (ts_name in names(ts_in)) {
    ld10_mtx[ts_name, cell_id] <- ts_in[ts_name]
  }
  ## calculate positions of the transcripts outside
  r <- sqrt(cell$area / pi)
  for (ts_name in names(ts_in)) {
    ## generate random position for transcript given cell centroid
    ## calculate angle
    theta <- runif(1, 0, 2*pi)
    ## calculate distance
    d <- runif(1, r * min_dist, r * max_dist)
    x <- cell$x + d * cos(theta)
    y <- cell$y + d * sin(theta)
    ## add to dataframe
    transcripts <- rbind(transcripts, c(ts_name, x, y, cell_id))
  }
}
colnames(transcripts) <- c('gene', 'x', 'y', 'orig_cell')
transcripts <- transcripts %>% 
  mutate(x = as.numeric(x), y = as.numeric(y))

### Reassign transcripts to cells ----------------------------------------

## convert transcripts to st
sf_transcripts <- sf::st_as_sf(transcripts, coords = c('x','y'))

tmp <- ld10_mtx
## for each cell, calculate intersection
for (cell_id in meta$cell) {
  ## get cell geometry
  cell_geo <- circles[[cell_id]]
  ## intersect genes with cell
  points_within <- st_within(sf_transcripts, cell_geo, sparse = FALSE)
  genes_in <- sf_transcripts[points_within, ]$gene
  if (length(genes_in) > 0) {
    print(paste(cell_id, length(genes_in)))
    count_genes_in <- table(genes_in)
    for (ts_name in names(count_genes_in)) {
      ld10_mtx[ts_name, cell_id] <- ld10_mtx[ts_name, cell_id] + 
        count_genes_in[ts_name]
    }
  }
}

## save
saveRDS(ld10_mtx, file = 'filtered_ld10_mtx.RDS')

## check differences
added_genes <- colSums(ld10_mtx - tmp)
mean(added_genes) ## 2.220984
median(added_genes) ## 1
max(added_genes) ## 22

## check original data
mean(colSums(filtered_srt_mtx)) ## 67.43361
median(colSums(filtered_srt_mtx)) ## 23
max(colSums(filtered_srt_mtx)) ## 1118

## check diffusion data
mean(colSums(ld10_mtx)) ## 63.98208
median(colSums(ld10_mtx)) ## 24
max(colSums(ld10_mtx)) ## 1024

## check data without adding genes (10 pcs loss)
mean(colSums(tmp)) ## 61.7611
median(colSums(tmp)) ## 22
max(colSums(tmp)) ## 1007


## Simulate diffusion 20 ----------------------------------------------------

### Simulate transcript positions ----------------------------------------

## percentage out of the cell
pct_out <- .2
## distance from the cell centroid as a function of area
min_dist <- 1.01
max_dist <- 2

## empty matrix to get the transcripts later
ld_mtx <- sparseMatrix(dim(filtered_srt_mtx)[1], 
                         dim(filtered_srt_mtx)[2], 
                         x = 0)
rownames(ld_mtx) <- rownames(filtered_srt_mtx)
colnames(ld_mtx) <- colnames(filtered_srt_mtx)
## create df with transcript coordinates represented by the cell center
transcripts <- data.frame()
## for every cell
for (i in 1:dim(meta)[1]) {
  print(paste(i, dim(meta)[1]))
  cell <- meta[i, ]
  cell_id <- cell$cell
  ## get expressed transcripts 
  ts <- filtered_srt_mtx[which(filtered_srt_mtx[, cell_id] > 0), cell_id]
  ## calculate transcripts outside and inside the cell
  ts_out <- round(ts * pct_out)
  ts_in <- ts - ts_out
  ## populate exp matrix with the transcripts inside the cell
  for (ts_name in names(ts_in)) {
    ld_mtx[ts_name, cell_id] <- ts_in[ts_name]
  }
  ## calculate positions of the transcripts outside
  r <- sqrt(cell$area / pi)
  for (ts_name in names(ts_in)) {
    ## generate random position for transcript given cell centroid
    ## calculate angle
    theta <- runif(1, 0, 2*pi)
    ## calculate distance
    d <- runif(1, r * min_dist, r * max_dist)
    x <- cell$x + d * cos(theta)
    y <- cell$y + d * sin(theta)
    ## add to dataframe
    transcripts <- rbind(transcripts, c(ts_name, x, y, cell_id))
  }
}
colnames(transcripts) <- c('gene', 'x', 'y', 'orig_cell')
transcripts <- transcripts %>% 
  mutate(x = as.numeric(x), y = as.numeric(y))

### Reassign transcripts to cells ----------------------------------------

## convert transcripts to st
sf_transcripts <- sf::st_as_sf(transcripts, coords = c('x','y'))

## for each cell, calculate intersection
for (cell_id in meta$cell) {
  ## get cell geometry
  cell_geo <- circles[[cell_id]]
  ## intersect genes with cell
  points_within <- st_within(sf_transcripts, cell_geo, sparse = FALSE)
  genes_in <- sf_transcripts[points_within, ]$gene
  if (length(genes_in) > 0) {
    print(paste(cell_id, length(genes_in)))
    count_genes_in <- table(genes_in)
    for (ts_name in names(count_genes_in)) {
      ld_mtx[ts_name, cell_id] <- ld_mtx[ts_name, cell_id] + 
        count_genes_in[ts_name]
    }
  }
}

## save
saveRDS(ld_mtx, file = 'filtered_ld20_mtx.RDS')
