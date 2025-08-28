# file with functdions used in HistoMapST
# author: Roman Laddach
# date: 14-Aug-2025
# licence: MIT

library("rdist")

read_output_from_QuPath <- function(path_to_file, file_name){
  # function to read csv output from QuPath based on the path and file name
  
  QuPath_full <- read.csv(paste(path_to_file, "/", file_name, sep=""))
  # adding rownames
  rownames(QuPath_full) <- paste("hist", 1:dim(QuPath_full)[1], sep="_")
  
  return(QuPath_full)
}


select_coordinates_from_QuPath <- function(object_name, x_coord, y_coord){
  # function to subset the QuPath output to only include coordinates
  
  QuPath_coordinates = object_name[, c(x_coord, y_coord)]
  # adding rownames from original file
  rownames(QuPath_coordinates) =  rownames(QuPath_full)
  # selecting and renaming relevant columns
  colnames(QuPath_coordinates) = c("x_coord", "y_coord")
  # flipping y coordinates to match Cartesian coordinate system
  QuPath_coordinates$y_coord = -QuPath_coordinates$y_coord
  
  return(QuPath_coordinates)
}


read_coordinates_from_SR <- function(path_to_file){
  # function to retrieve coordinates from Seurat object
  
  SR_coordinates <- read.csv(paste(path_to_file, "tissue_positions_list.csv", sep=""), 
                             col.names = c("tissue", "row", "col", "imagerow", "imagecol"),
                             row.names = 1)
  # filtering to only include spots under the tissue
  SR_coordinates <- SR_coordinates[SR_coordinates$tissue == 1, ]
  # flipping the y coordinatesto match Cartesian coordinate system and representation in Seurat outputs
  SR_coordinates$imagerow = -SR_coordinates$imagerow
  # selecting and renaming relevant columns
  SR_coordinates_subset = SR_coordinates[, c("imagecol", "imagerow")]
  colnames(SR_coordinates_subset) = c("x_coord", "y_coord")
  # adding rownames 
  rownames(SR_coordinates_subset) = rownames(SR_coordinates)
  
  return(SR_coordinates_subset)
}


identify_cells_within_spots <- function(SR_coordinates_df, QuPath_coordinates_df){
  # function to convert spot distances from um to pixels and subsequently calculate the distance between cell centroids from QuPath
  # and spaceranger output and subsequently return a dataframe containing cells within a spot 
  
  # conversion of spot distances in um to pixels
  spot_distance_dist = dist(SR_coordinates[,1:2])
  spot_2_spot_dist_px = round(min(as.vector(spot_distance_dist)))
  radius_um = 55/2
  spot_2_spot_dist_um = 100
  radius_px = radius_um * spot_2_spot_dist_px / spot_2_spot_dist_um
  
  # calculating the distance between centroids of spots and cells
  distances_combined_df = as.data.frame(cdist(SR_coordinates[,1:2], QuPath_coordinates[,1:2]))
  rownames(distances_combined_df) = rownames(SR_coordinates)
  colnames(distances_combined_df) = rownames(QuPath_coordinates)
  distances_combined_df$spot_ID = rownames(distances_combined_df)
  
  # filtering for cells within the radius of a spot
  cells_within_spots = melt(distances_combined_df)
  cells_within_spots = cells_within_spots[cells_within_spots$value <= radius_px, ]
  
  colnames(cells_within_spots)[2:3] = c("cell_ID", "distance")
  
  return(cells_within_spots)
}