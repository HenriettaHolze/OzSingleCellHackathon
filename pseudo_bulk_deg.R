#' @description
#' Aggregate expression in seurat object by clones
#' 
#' @param seurat_object Seurat object with clone metadata.
#' @param clone_id Metadata column with clone information
#' @param assays Assay to aggregate and perform DEG analysis on. Always using counts slot.

AggregateClones <-
  function(seurat_object,
           clone_id = "clone_id",
           assays = "originalexp") {
    seurat_object_aggr <-
      AggregateExpression(
        seurat_object,
        assays = assays,
        return.seurat = T,
        group.by = c(clone_id),
        slot = "counts",
      )
    
    colnames(seurat_object_aggr) <- gsub("-", "_", colnames(seurat_object_aggr))
    return(seurat_object_aggr)
  }

#' @description
#' Pseudo bulk DEG analysis of clone groups.
#' 
#' @param seurat_object_aggr Seurat object aggregated by clones.
#' @param cluster_anno Vector with clones as names and clusters as values.
#' @param assays Assay to aggregate and perform DEG analysis on. Always using counts slot.
#' 
CloneGroupDeg <- function(seurat_object_aggr,
                          cluster_anno) {
  
  # add clone metadata 
  seurat_object_aggr <-
    AddMetaData(seurat_object_aggr, metadata = cluster_anno, col.name = "clone_groups")
  
  # get markers for all groups
  Idents(seurat_object_aggr) <- "clone_groups"
  cluster_markers <-
    FindAllMarkers(seurat_object_aggr,
                   group.by = "clone_groups",
                   test.use = "DESeq2",
    )
  
  return(cluster_markers)
}


