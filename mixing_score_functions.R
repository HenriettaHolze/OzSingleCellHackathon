MixingScoreSingleClone <- function(cells, k, cells_to_clones_map, clone_size, dist_matrix_clones) {
  col_dist <- colnames(dist_matrix_clones)
  k <- floor(k)
  sum_ki <- sapply(cells, function(cell) {
    # get nearest neighbours (not cell itself)
    kneigh <-
      col_dist[order(dist_matrix_clones[cell,])[seq(2, k + 1)]]
    # count how many neighbours are of same clone
    sum(cells_to_clones_map[kneigh] == cells_to_clones_map[cell]) / (2 * k)
  })
  sum_i <- sum(sum_ki) / clone_size
  
  return(sum_i)
}

MixingScoreTwoClones <-
  function(clones_to_cells_map, clones, clone_sizes, k_frac, dist_matrix, cells_to_clones_map) {
    k_x <- k_frac * clone_sizes[clones[1]]
    k_y <- k_frac * clone_sizes[clones[2]]
    
    cells_both_clones <-
      as.vector(unlist(clones_to_cells_map[c(clones[1], clones[2])]))
    dist_matrix_clones <-
      dist_matrix[cells_both_clones, cells_both_clones]
    
    mixing_score_x <-
      MixingScoreSingleClone(cells = clones_to_cells_map[[clones[1]]],
                             k = k_x,
                             cells_to_clones_map = cells_to_clones_map,
                             clone_sizes[clones[1]],
                             dist_matrix_clones)
    mixing_score_y <-
      MixingScoreSingleClone(
        cells = clones_to_cells_map[[clones[2]]],
        k = k_y,
        cells_to_clones_map = cells_to_clones_map,
        clone_size = clone_sizes[clones[2]],
        dist_matrix_clones = dist_matrix_clones
      )
    
    mixing_score <- 2 - 2 * (mixing_score_x + mixing_score_y)
    
    return(mixing_score)
  }


#' @param embedding dataframe with embedding. Rows match rows in clone_df
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
                              clone_id = "clone_id") {
  clone_sizes <- table(clone_df[clone_id])
  
  clones <- names(clone_sizes)
  clones_to_cells_map <-
    lapply(clones, function(c)
      rownames(clone_df)[clone_df[clone_id] == c])
  names(clones_to_cells_map) <- clones
  
  cells_to_clones_map <- clone_df[[clone_id]]
  names(cells_to_clones_map) <- rownames(clone_df)
  
  dist_matrix <- stats::dist(embeddings)
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
      MixingScoreTwoClones(clones_to_cells_map,
                           clones[clone_pair],
                           clone_sizes,
                           k_frac,
                           dist_matrix,
                           cells_to_clones_map)
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
  embeddings <-
    Embeddings(object = object[[reduction]])[rownames(clone_df), dims]
  
  scores <- MixingScoreClones(
    embedding = embeddings,
    k_frac = k,
    clone_df = clone_df,
    cores = cores,
    clone_id = clone_id
  )
  return(scores)
}