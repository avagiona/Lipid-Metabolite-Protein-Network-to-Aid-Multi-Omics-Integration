library("dplyr")
library("igraph")
library("ggplot2")
library("NetHypGeom")


#set pathway
setwd("~/R-4.2.1/multi_omics_embedding")
options(scipen = 999)

#import file 
coor_multi_omics<-readRDS("coordinates_protein_lipids_new.RData")

nodes_polar<-coor_multi_omics$polar

find_clusters_by_top_gaps <- function(theta, top_n_gaps = 10) {
  # Sort the angular coordinates
  sorted_theta <- sort(theta)
  
  # Compute the differences between successive theta values
  theta_diff <- diff(sorted_theta)
  
  # Identify the top_n largest gaps
  if (top_n_gaps > length(theta_diff)) {
    stop("Number of top gaps exceeds the number of available gaps.")
  }
  
  # Get the indices of the top_n largest gaps
  top_gaps_indices <- order(theta_diff, decreasing = TRUE)[1:top_n_gaps]
  top_gaps <- theta_diff[top_gaps_indices]
  
  # Use these top gaps to define cluster boundaries
  gap_indices <- sort(top_gaps_indices)  # Ensure gaps are processed in order
  
  # Initialize list to hold clusters
  clusters <- list()
  
  # Define the start and end of each cluster
  start_idx <- 1
  for (gap_idx in gap_indices) {
    end_idx <- gap_idx + 1
    cluster_size <- end_idx - start_idx
    clusters[[length(clusters) + 1]] <- list(
      start = start_idx,
      end = end_idx,
      size = cluster_size
    )
    start_idx <- end_idx + 1
  }
  
  # Add the final cluster if there are remaining points
  if (start_idx <= length(sorted_theta)) {
    clusters[[length(clusters) + 1]] <- list(
      start = start_idx,
      end = length(sorted_theta),
      size = length(sorted_theta) - start_idx + 1
    )
  }
  
  # Return clusters and top gaps used
  return(list(clusters = clusters, top_gaps = top_gaps))
}

# Find clusters using the top 5 largest gaps
result <- find_clusters_by_top_gaps(nodes_polar$theta, top_n_gaps = 15)

get_cluster_sizes <- function(clustering_result) {
  # Ensure that the input is a list containing the clusters
  if (!"clusters" %in% names(clustering_result)) {
    stop("The input must be the result of find_clusters_by_top_gaps function.")
  }
  
  # Extract cluster sizes into a vector
  cluster_sizes <- sapply(clustering_result$clusters, function(cluster) cluster$size)
  
  return(cluster_sizes)
}

cluster_sizes<-get_cluster_sizes(result)
cluster_sizes

get_cluster_ends <- function(clustering_result) {
  # Ensure that the input is a list containing the clusters
  if (!"clusters" %in% names(clustering_result)) {
    stop("The input must be the result of find_clusters_by_top_gaps function.")
  }
  
  # Extract the end indices into a vector
  cluster_ends <- sapply(clustering_result$clusters, function(cluster) cluster$end)
  
  return(cluster_ends)
}

cluster_ends <- get_cluster_ends(result)



























