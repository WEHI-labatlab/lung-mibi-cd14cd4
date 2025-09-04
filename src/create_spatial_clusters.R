#' @author Big K
#' 
#' 
#'
#' @param 
#'
#' @return
#' @examples
#' 

library(dplyr)
library(tidyr)



create_spatial_clusters <- function(cell_df, 
                                    image_col,
                                    phenotype_col,
                                    x_coordinates, 
                                    y_coordinates, 
                                    k_closest=10,
                                    max_dist=50, 
                                    number_clusters_kmeans=7,
                                    keep_proportions=TRUE,
                                    random.state = 123){
  
  proportions_df <- find_nearest_neighbour_proportions(cell_df, 
                                                      image_col,
                                                      phenotype_col,
                                                      x_coordinates, 
                                                      y_coordinates, 
                                                      k_closest,
                                                      max_dist)
  
  set.seed(random.state)
  cluster_objs <- kmeans(proportions_df %>% subset(select=-c(Image)), number_clusters_kmeans)
  cluster <- cluster_objs$cluster
  if (keep_proportions){
    cluster <- cbind(proportions_df, cluster)
  }
  return(cluster)
}

find_nearest_neighbour_proportions <- function(cell_df, 
                                               image_col,
                                               phenotype_col,
                                               x_coordinates, 
                                               y_coordinates, 
                                               k_closest=10,
                                               max_dist=50){
    cell_df$temp_cell_id <- rownames(cell_df)
    
    proportions_df <- cell_df %>% 
        group_by(!!as.name(image_col)) %>%
        group_modify(~ {
            
            
            find_nearest_neighbour_proportions_(
            cell_df=.x,
            phenotype_col=phenotype_col,
            x_coordinates=x_coordinates, 
            y_coordinates=y_coordinates, 
            k_closest=k_closest,
            max_dist=max_dist
        )}, .keep = TRUE) %>%
        ungroup()
    
    
    proportions_df <- proportions_df[match(cell_df$temp_cell_id, proportions_df$temp_cell_id),]
    proportions_df <- proportions_df %>% replace(is.na(.), 0)
    proportions_df <- proportions_df %>% subset(select=-c(temp_cell_id))
    
    return(proportions_df)
}


find_nearest_neighbour_proportions_ <- function(cell_df, 
                                    phenotype_col,
                                    x_coordinates, 
                                    y_coordinates, 
                                    k_closest=10,
                                    max_dist=50
                                    ){
  if (nrow(cell_df) == 1){
      return(cell_df[, c("temp_cell_id")])
  }
    
    # add one since if you want 10 nearest neighbours, you want 11 nearest
  # including the cell itself

    
  k_closest <- k_closest + 1
  
  if (k_closest > nrow(cell_df)){
    warning(paste0(cell_df$Image[1], " has less rows than the set nearest neighbours. NN will be calculated with nrow()"))
    k_closest <- nrow(cell_df)
  }
  
  
  cell_df <- data.frame(cell_df, check.names = FALSE)
  
  # Run a nearest neighbor analysis to determine the closest cells to the cells of interest
  closest <- RANN::nn2(
    data = cell_df[, c(x_coordinates, y_coordinates)], 
    query = cell_df[, c(x_coordinates, y_coordinates)], 
    k = k_closest
  )
  
  # Data frame of: rows = query cells and columns = top k closest cell indexes
  dists <- data.frame(closest$nn.dists)
  colnames(dists) <- paste("Dist", 1:k_closest-1, sep = ":")
  
  # Data frame of index of the closest cells
  nn.idx <- data.frame(closest$nn.idx)
  colnames(nn.idx) <- paste("NN:idx", 1:k_closest-1, sep = ":")
  
  # remove the first column since it is itself
  dists <- dists %>% select(-1)
  nn.idx <- nn.idx %>% select(-1)
  
  
  # Initialise empty dataframe to store the string nearest neighbours
  nn <- nn.idx
  colnames(nn) <- paste("NN", 1:(k_closest-1), sep=":")
  nn[] <- cell_df[match(as.integer(unlist(nn.idx)), as.integer(rownames(cell_df))), phenotype_col]
  
  # rename cell type if too far away
  null_cell_type <- "TOO_FAR"
  
  # remove neighbours which are too far away
  nn[dists > max_dist] <- null_cell_type
  
  # calculate the proportions. For each cell, count the number of a certain cell
  # type being its neighbour. 
  nn$temp_id <- as.numeric(rownames(nn))
  proportions <- nn %>% 
    pivot_longer(cols = -temp_id, names_to = "Hierarchy:Level", 
                 values_to = "Phenotype") %>%
    mutate(Phenotype = as.character(Phenotype)) %>%
    group_by(temp_id) %>%
    count(Phenotype) %>%
    pivot_wider(names_from = Phenotype, values_from = n) %>%
    replace(is.na(.), 0)
  proportions <- proportions %>% subset(select=-c(temp_id))
  
  proportions$temp_cell_id <- cell_df$temp_cell_id
  
  
  f <- function(x) x / (k_closest - 1)
  
  if ("TOO_FAR" %in% colnames(proportions)){
    proportions <- proportions %>% subset(select=-c(TOO_FAR))
  }
  proportions <- proportions %>% 
      mutate_if(is.numeric, f)
  
  return(proportions)
}


