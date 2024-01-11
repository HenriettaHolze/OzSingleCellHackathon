MixingScoreSingleClone <-
  function(cells,
           clone, 
           k,
           cells_to_clones_map,
           clone_size,
           dist_matrix_clones) {
    col_dist <- colnames(dist_matrix_clones)
    k <- floor(k)
    sum_ki <- sapply(cells, function(cell) {
      # get nearest neighbours (not cell itself)
      kneigh <-
        col_dist[order(dist_matrix_clones[cell,])[seq(2, k + 1)]]
      # count how many neighbours are of same clone
      sum(cells_to_clones_map[kneigh] == clone) / (2 * k)
    })
    sum_i <- sum(sum_ki) / clone_size
    
    return(sum_i)
  }

#' @description
#' Calculates a mixing score between 2 clones
#'
#' @returns the mixing score between 2 clones as float
#'
#' @param clones_to_cells_map list with clones as names and vector of cell IDs as values
#' @param clones vector with 2 clone names
#' @param clone_sizes table with clone sizes
#' @param k_frac fraction of cells of a clone to use as nearest neighbours
#' @param dist_matrix distance matrix between all cells of 2 clones
#' @param cells_to_clones_map list with cell IDs as names and clones as values
#'
MixingScoreTwoClones <-
  function(clones_to_cells_map,
           clones,
           clone_sizes,
           k_frac,
           dist_matrix,
           cells_to_clones_map) {
    k_x <- k_frac * clone_sizes[clones[1]]
    k_y <- k_frac * clone_sizes[clones[2]]
    cells_both_clones <-
      as.vector(unlist(clones_to_cells_map[c(clones[1], clones[2])]))

    dist_matrix_clones <-
      dist_matrix[cells_both_clones, cells_both_clones]
    
    mixing_score_x <-
      MixingScoreSingleClone(
        cells = clones_to_cells_map[[clones[1]]],
        clone = clones[1],
        k = k_x,
        cells_to_clones_map = cells_to_clones_map,
        clone_sizes[clones[1]],
        dist_matrix_clones
      )
    mixing_score_y <-
      MixingScoreSingleClone(
        cells = clones_to_cells_map[[clones[2]]],
        clone = clones[2],
        k = k_y,
        cells_to_clones_map = cells_to_clones_map,
        clone_size = clone_sizes[clones[2]],
        dist_matrix_clones = dist_matrix_clones
      )
    
    mixing_score <- 2 - 2 * (mixing_score_x + mixing_score_y)
    
    return(mixing_score)
  }

#' @description
#' Calculates a symmetrical mixing score between a set of clones given an embedding.
#'
#' @returns a dataframe with clone names as row and column names and mixing scores as values.
#'
#' @param embedding dataframe or matrix with embedding. rows are cells and columns are dimensions. 
#' @param dims dimensions to use from embedding
#' @param k_frac fraction of cells of a clone to use as nearest neighbours
#' @param clones_df dataframe with cells as rownames and a column with clone IDs
#' @param clone_id name of column with clone IDs
#' @param cores compute cores to use
#'
MixingScoreClones <- function(embedding,
                              k_frac = 0.3,
                              clone_df,
                              cores = 1,
                              clone_id) {
  clone_sizes <- table(clone_df[clone_id])
  
  clones <- names(clone_sizes)
  clones_to_cells_map <-
    lapply(clones, function(c)
      rownames(clone_df)[clone_df[clone_id] == c])
  names(clones_to_cells_map) <- clones
  
  cells_to_clones_map <- clone_df[[clone_id]]
  names(cells_to_clones_map) <- rownames(clone_df)
  
  dist_matrix <- stats::dist(embedding)
  dist_matrix <- as.matrix(dist_matrix)
  
  # prepare empty results dataframe
  df_scores <-
    data.frame(matrix(nrow = length(clones), ncol = length(clones)))
  
  clone_index_pairs <-
    as.list(as.data.frame(t(
      gtools::permutations(
        n = length(clones),
        r = 2,
        repeats.allowed = F
      )
    )))
  
  clone_index_pairs <-
    clone_index_pairs[sapply(clone_index_pairs, function(x)
      x[1] < x[2])]
  
  # iterate over all clone pairs
  mixing_scores <-
    parallel::mclapply(clone_index_pairs, mc.cores = cores, function(clone_pair) {
      MixingScoreTwoClones(
        clones_to_cells_map,
        clones[clone_pair],
        clone_sizes,
        k_frac,
        dist_matrix,
        cells_to_clones_map
      )
    })
  
  for (i in seq(1, length(clone_index_pairs))) {
    df_scores[clone_index_pairs[[i]][1], clone_index_pairs[[i]][2]] <-
      mixing_scores[[i]]
    df_scores[clone_index_pairs[[i]][2], clone_index_pairs[[i]][1]] <-
      mixing_scores[[i]]
  }
  
  rownames(df_scores) <- clones
  colnames(df_scores) <- clones
  
  return(df_scores)
}

#' @param object Seurat object
#' @param clones vector of clone names to calculate mixing score for
#' @param clone_id column in metadata containing clone information
MixingScoreSeurat <- function(object,
                              clones,
                              clone_id,
                              reduction = "pca",
                              dims = seq(1, 20),
                              k = 0.3,
                              cores = 1) {
  clones <- as.character(clones)
  
  # dataframe with clone IDs and cell IDs as rownames
  clone_df <- object@meta.data %>%
    filter(!!as.name(clone_id) %in% clones) %>%
    select(!!as.name(clone_id))
  
  # get embedding with specified dimensions from Seurat object
  # subset embedding to cells from selected clones
  embedding <-
    Embeddings(object = object[[reduction]])[rownames(clone_df), dims]
  
  scores <- MixingScoreClones(
    embedding = embedding,
    k_frac = k,
    clone_df = clone_df,
    cores = cores,
    clone_id = clone_id
  )
  return(scores)
}


#' @description
#' Old mixing score function, requires to downsample the larger clone.
#' https://rdrr.io/github/satijalab/seurat/man/MixingMetric.html
#' https://rdrr.io/github/satijalab/seurat/src/R/integration.R
#' 
mixingCoefficient <- function(object, 
                              reduction,
                              dims,
                              k,
                              clones,
                              cores = 1,
                              max_cells = 50,
                              rep = 1,
                              clone_id = "clone_id") {
  clones <- as.character(clones)
  
  # get embedding with specified dimensions from Seurat object
  embeddings <- Embeddings(object = object[[reduction]])[, dims]
  
  # dataframe of cell IDs and clone IDs
  cells_df <- object@meta.data %>%
    rownames_to_column(var = "cell_id_rowname") %>%
    filter(!!as.name(clone_id) %in% clones) %>%
    select(!!as.name(clone_id), cell_id_rowname)
  
  # subset embedding to cells from selected clones
  embeddings <- embeddings[cells_df$cell_id_rowname, ]
  
  # calculate distance matrix between all cells of all selected clones
  dist_matrix <- dist(embeddings)
  dist_matrix <- as.matrix(dist_matrix)
  
  # create clone to cell_ids dict
  clone_to_cell_dics <- split(cells_df$cell_id_rowname, cells_df[[clone_id]])
  
  # prepare empty results dataframe
  df_scores <-
    data.frame(matrix(nrow = length(clones), ncol = length(clones)))
  
  clone_index_pairs <-
    as.list(as.data.frame(t(
      gtools::permutations(
        n = length(clones),
        r = 2,
        repeats.allowed = F
      )
    )))
  
  clone_index_pairs <-
    clone_index_pairs[sapply(clone_index_pairs, function(x)
      x[1] < x[2])]
  
  # cells_vector
  cells_vector <- cells_df$clone_id
  names(cells_vector) <- cells_df$cell_id_rowname
  
  # iterate over all clone pairs
  mixing_scores <-
    mclapply(clone_index_pairs, mc.cores = cores, function(clone_pair) {
      c = clones[clone_pair[1]]
      b = clones[clone_pair[2]]
      
      cells1 <- as.vector(unlist(clone_to_cell_dics[c]))
      cells2 <- as.vector(unlist(clone_to_cell_dics[b]))
      
      # iterate over sampling repetitions
      mixing_scores <- sapply(rep, function(r) {
        # sample to smaller clone size
        if (length(cells1) < length(cells2)) {
          cells2 <- sample(cells2, size = length(cells1), replace = F)
        } else if (length(cells1) > length(cells2)) {
          cells1 <- sample(cells1, size = length(cells2), replace = F)
        }
        
        cell_list <- c(cells1, cells2)
        
        dist_matrix_clones <-
          dist_matrix[cell_list, cell_list]
        
        # k is half the number of cells of one clone
        k_clone_pair <- floor(length(cell_list) * 0.5 * k)
        col_dist <- colnames(dist_matrix_clones)
        # calculate average identity
        s <- sapply(cell_list, function(cell) {
          # get nearest neighbours (2-21 to not take same cell)
          kneigh <-
            col_dist[order(dist_matrix_clones[cell, ])[seq(2, k_clone_pair + 1)]]
          sum(cells_vector[kneigh] == cells_vector[cell])
        })
        
        # calculate mixing score as in
        mixing_score <-
          1 - (sum(s) / length(cell_list) - k_clone_pair / 2) / (k_clone_pair - k_clone_pair / 2)
        return(mixing_score)
      })
      
      return(mean(mixing_scores))
    })
  
  for (i in seq(1, length(clone_index_pairs))) {
    df_scores[clone_index_pairs[[i]][1], clone_index_pairs[[i]][2]] <-
      mixing_scores[[i]]
    df_scores[clone_index_pairs[[i]][2], clone_index_pairs[[i]][1]] <-
      mixing_scores[[i]]
  }
  
  rownames(df_scores) <- clones
  colnames(df_scores) <- clones
  return(df_scores)
}

#' @description
#' Calculates intra-clonal heterogeneity metrics: mean, median and standard deviation of cell-cell distances.
#'
intra_clonal_distance <-
  function(object, reduction, dims, clones, clone_id = "clone_id") {
    intra_clone_distances <- list()
    
    embeddings <- Seurat::Embeddings(object = object[[reduction]])[, dims]
    
    for (c_i in seq(length(clones))) {
      c = clones[c_i]
      
      # identify which cells belong to cone
      cells <-
        object@meta.data %>%
        filter(!!as.name(clone_id) %in% c) %>%
        rownames_to_column("cell") %>%
        pull(cell)
      
      # calculate distance matrix between cells in a clone
      dist_matrix <- dist(embeddings[cells, ])
      # calculate mean, median and sd of distances
      intra_clone_distances[[as.character(c)]] <-
        c(mean(dist_matrix), median(dist_matrix), sd(dist_matrix))
    }
    
    # returns a names list
    return(intra_clone_distances)
  }
