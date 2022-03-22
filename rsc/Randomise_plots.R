packages <- c("raster", "rgeos", "rgdal", "sp")
for(i in 1:NROW(packages)){
  if(!require(packages[i], character.only = TRUE)){
    install.packages(packages[i])
    library(packages[i], character.only = TRUE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
block <- args[[1]]
plot <- args[[2]]
n <- args[[3]]

shp_dir <- "/home/manuel/Nextcloud/LELE_2021/QGIS/Shapefile"
shapefile_dir <- file.path(shp_dir, "Grid_10m.shp")

Grid_10m <- rgdal::readOGR(shapefile_dir)

cat("\n\nRandomised sample for block", block, "plot", plot, "contains the following subplots:\n")
for(i in 1:n){
  row <- sample(seq(2, 9), 1)
  column <- sample(LETTERS[seq(2, 9)], 1)
  
  cell <- Grid_10m[Grid_10m@data$Quadrat == paste0(column, as.character(row)) &
                     Grid_10m@data$Block == block &
                     Grid_10m@data$Plot == plot,]
  cell <- cell@polygons[[1]]
  cell <- cell@Polygons[[1]]
  
  x_coords <- cell@coords[, 1]
  y_coords <- cell@coords[, 2]
  
  corners <- data.frame(x = x_coords, y = y_coords)[-5, ]
  corner_points_UTMz35S <- sp::SpatialPoints(corners, proj4string = CRS("+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs"))
  
  corner_points_WGS84 <- sp::spTransform(corner_points_UTMz35S, CRS("+proj=longlat +datum=WGS84 +no_defs"))
  cat("\nSubplot ID:", paste0(column, as.character(row), "\n"))
  cat("Corner coordinates:\n")
  print(corner_points_WGS84@coords)
}
