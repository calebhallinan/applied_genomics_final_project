library(tidyverse)
library(Seurat)
library(sf)
library(anndata)
library(biomaRt)
library(Matrix)
library(BiocParallel)

set.seed(42)



### Read data
sc_meta <- read.csv('sim_data/prox_assign/sc_meta.csv', row.names = 1)
meta <- sc_meta %>%
    mutate(radius = sqrt(area / pi))
head(meta)
mtx <- readRDS("sim_data/prox_assign/mtx_500_genes.RDS")



### Sim transcript pos
ncells <- dim(meta)[1]

## percentage out of the cell
pct_out <- .05
## distance from the cell centroid as a function of area
min_dist <- 1
max_dist <- 2

## populate matrix with transcripts that will not suffer diffusion
ld_mtx <- mtx - round(mtx * pct_out)

## for every chunk
simulate_transcripts <- function(chunks) {
  mtx_chunk <- chunks[[1]]
  meta_chunk <- chunks[[2]]
  ## create df with transcript coordinates represented by the cell center
  transcripts <- data.frame()
  ## for every cell
  for (i in 1:dim(meta_chunk)[1]) {
    cell <- meta_chunk[i, ]
    cell_id <- cell$cell
    ## get expressed transcripts 
    ts <- mtx_chunk[which(mtx_chunk[, cell_id] > 0), cell_id]
    ## calculate transcripts outside and inside the cell
    ts_out <- round(ts * pct_out)
    ## calculate positions of the transcripts outside
    r <- sqrt(cell$area / pi)
    for (ts_name in names(ts_out)) {
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
  return(transcripts)
}

## separate mtx and meta into chunks
n_chunks <- 50
cells_per_chunck <- ceiling(ncol(mtx) / n_chunks)
## mtx
submatrices <- list()
for (i in 1:(n_chunks-1)) {
  start_col <- (i - 1) * cells_per_chunck + 1
  end_col <- i * cells_per_chunck
  submatrix <- mtx[, start_col:end_col]
  submatrices[[i]] <- submatrix
}
submatrices[[n_chunks]] <- mtx[, ((n_chunks-1)*cells_per_chunck) : ncol(mtx)]
## meta
subdfs <- list()
for (i in 1:(n_chunks-1)) {
  start_col <- (i - 1) * cells_per_chunck + 1
  end_col <- i * cells_per_chunck
  subdf <- meta[start_col:end_col, ]
  subdfs[[i]] <- subdf
}
subdfs[[n_chunks]] <- meta[((n_chunks-1)*cells_per_chunck) : nrow(meta), ]
## join into list
arguments <- lapply(1:n_chunks, FUN = function(i) {list(submatrices[[i]], subdfs[[i]])})

## run parallelization
bpparam <- MulticoreParam(workers = n_chunks)
transcripts_list <- bplapply(arguments, simulate_transcripts, BPPARAM = bpparam)
## kachow

## join dataframes
transcripts <- bind_rows(transcripts_list)


## intersect transcripts with cells
## for each cell, calculate intersection
interesect_transcripts <- function(meta_chunk){
  ## create empty mtx
  sim_mtx <- Matrix::sparseMatrix(dim(mtx)[1],
    dim(mtx)[2],
    x = 0
  )
  rownames(sim_mtx) <- rownames(mtx)
  colnames(sim_mtx) <- colnames(mtx)
  ## for each cell
  for (i in 1:dim(meta_chunk)[1]) {
    cell <- meta_chunk[i, ]
    ## intersect genes with cell
    ## filter genes near to speed the process
    near_transcripts <- transcripts %>% 
      mutate(distance = sqrt((x - cell$x)^2 + (y - cell$y)^2)) %>% 
      filter(distance <= cell$radius)
    ## if no genes are near, just skip
    if (dim(near_transcripts)[1] > 0) {
      genes_in <- near_transcripts$gene
      ## add the number of copies to the sim matrix
      count_genes_in <- table(genes_in)
      for (ts_name in names(count_genes_in)) {
        sim_mtx[ts_name, cell$cell] <- count_genes_in[ts_name]
      }
    }
  }
  return(sim_mtx)
}
## run parallelization
bpparam <- MulticoreParam(workers = n_chunks)
mtxs_list <- bplapply(subdfs, interesect_transcripts, BPPARAM = bpparam)
## kachow
sum_mtx <- Reduce('+', mtxs_list)

## sum ld genes with genes that did not ld
ld_mtx <- ld_mtx + sum_mtx
## check if there was any ld
all(ld_mtx == mtx)

## save sim data
saveRDS(ld_mtx, 'sim_data/prox_assign/ld_mtx_5pct_parallelized.RDS')
