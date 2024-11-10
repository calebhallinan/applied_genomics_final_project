library(Seurat)
library(tidyverse)
library(sf)
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
scrattch.io::write_dgCMatrix_h5(mtx, cols_are = 'cell_names', 
                                'sim_data/mtx.h5', gene_ids = NULL)
# write.csv(mtx, 'sim_data/mtx.csv')

## add cell names from the sc data to the srt metadata in order by area
meta <- sub_srt_meta %>% 
  arrange(area) %>% 
  mutate(cell = names(sort(colSums(mtx))))
write.csv(meta, 'sim_data/sc_srt_meta.csv')



## Filter by srt genes ---------------------------------------------------

srt_genes <- rownames(srt_mtx)
sc_genes <- rownames(sc_mtx)
common_genes <- intersect(srt_genes, sc_genes)

filtered_srt_mtx <- mtx[common_genes, ]
scrattch.io::write_dgCMatrix_h5(filtered_srt_mtx, cols_are = 'cell_names', 
                                'sim_data/filtered_srt_mtx.h5')


## Create cell boundaries ------------------------------------------------

## create circles representing the cell boundaries
circles <- lapply(1:dim(meta)[1], 
                  FUN = function(i) {
                    cell <- meta[i, ]
                    st_buffer(st_point(c(cell$x, cell$y)), 
                              dist = sqrt(cell$area / pi))
                  })
meta <- meta %>% 
  mutate(boundary = circles)


## Simulate diffusion ----------------------------------------------------

## create df with transcript coordinates represented by the cell center
transcript_list <- list()
## for every cell
for (c in 1:dim(meta)[1]) {
  cell <- meta[i, 'cell']
  ## get expressed transcripts
  ts <- which(mtx[, cell] > 0)
  for (t in 1:dim(mtx)[1]) {
    
    append(transcript_list, c)
  }
}

  