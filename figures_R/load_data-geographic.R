library(leaflet) # For interactive maps
library(sf) # General map package
library(spatstat) # For risk maps
library(sparr) # For risk maps

if(!exists("samples")){
  stop("Need to run load_data-sequencing.R first!")
}

# Load geographic data ####

## |- Load and name datasets ####

# Manually designed with MapHub
gdf <- st_read("../market_map/market.geojson")
# `contour` is needed to identify elements corresponding to contour in the market map
contour <- st_read("../market_map/market-contour.geojson")

# Project geographic data (3857 needed with leaflet)
epsg_WGS84 <- 4326 # Initial code for this projection
epsg <- 3857 # Code needed for leaflet

gdf <- st_transform(gdf, epsg) 
contour <- st_transform(contour, epsg)

# From the MapHub map, retrieve names using number of datapoints (hard-coded)
# table(gdf$group) # Show number of points in each group and thereby identify their codes
dicoMap <- c("market_map", "cases", "cases_koopmans", "boundary", "positive_sites", "negative_sites", "unc_negative_sites", "negative_drains", "positive_drains", "stalls", "wildlife_stalls", "landmarks", "boundary_withDrains", "street_labels", "stall_labels", "street_labels_Wside")
names(dicoMap) <- c("3450322459", "3932803488", "3183518500", "3324846982", "1511299674", "4250682362", "524707348", "4097735875", "3869397518", "3957334658", "1134260820", "1847080433", "2949713191", "1616031262", "2778039285", "3869706544")#"570381184")
# indices of groups to rename
irename <- which(is.element(gdf$group, names(dicoMap)))
# Rename the groups for which we have a name
gdf$group[irename] <- dicoMap[as.character(gdf$group[irename])]
rm(irename)
table(gdf$group)

# Get the market boundary
market_boundary <- gdf[which(gdf$group=="boundary"), ]
market_boundary_withDrains <- gdf[which(gdf$group=="boundary_withDrains"), ]

# Landmarks 
landmarks <- gdf[which(gdf$group == "landmarks"), ]
landmarkNames <- unname(unlist(st_drop_geometry(gdf[gdf$group == "landmarks", "title"])))
landmarks$address <- landmarks$title # Add column with name matching what is used elsewhere

# Wildlife stalls
wildlifeStalls <- gdf[which(gdf$group == "wildlife_stalls"), ]

# Street labels
streetLabels <- gdf[which(gdf$group == "street_labels"), ]
streetLabelsWside <- gdf[which(gdf$group == "street_labels_Wside"), ]


# Stall labels
stallLabels <- gdf[which(gdf$group == "stall_labels"), ]

# Separate stalls and pillars
tmp <- gdf[gdf$group == "market_map",]

# Identify polygons - pillars
ipg <- which(substr(tmp$title, 1, 5) == "Polyg")

# Remove pillars
tmp2 <- tmp[-c(ipg, ipg - 1), ]

# Identify elements that correspond to the contour
idContour <- contour$title
icontour <- which(is.element(tmp2$title, idContour))

map_stalls <- tmp2[-icontour, ]
map_contours <- tmp2[icontour, ]

# Remove empty data
map_stalls <- map_stalls[!is.na(map_stalls$group), ]

rm(tmp, tmp2, icontour, idContour, ipg)

# Make sure we have locations
# Check that we could assign a position to all samples
stopifnot(all(is.element(unique(samples$address), landmarkNames)))
  
# Check that we have all positions for the cases
stopifnot(all(is.element(cases$location, landmarkNames)))





## |- Rotate geographic data ####

# We want a horizontal version of the geographic data for plotting

# Load rotation functions and parameters
source("map_computeRotation.R") 

# Compute common centroid
# Combine the shapes
x_combined <- sf::st_combine(gdf)
# Get centroid
center_coords <- sf::st_centroid(x_combined)

for(shapeName in c("market_boundary", "landmarks", "wildlifeStalls", "map_stalls", "market_boundary_withDrains", "streetLabels", "stallLabels", "streetLabelsWside")){
  cat("Rotating", shapeName, "...")
  # Initialize new shape
  assign(paste0(shapeName, "_r"), get(shapeName))
  # Rotate coordinates and translate East wing
  newCoords <- st_translate(st_rotate(get(shapeName), radians = angle2, center_coords = center_coords))
  # Put the new coordinates in the new shape
  assign(paste0(shapeName, "_r"), newCoords)
  cat("done.\n")
}
rm(newCoords, x_combined, center_coords)


## |- Set windows for spatial analyses ####

# For spatial analyses: Set the window: contour of the market, within which spatial analysis done
window <- as.owin(market_boundary)
window_withDrains <- as.owin(market_boundary_withDrains)

window_r <- as.owin(market_boundary_r)
window_withDrains_r <- as.owin(market_boundary_withDrains_r)



## Add geographic information to the other datasets ####

# Samples
samples_geo <- merge(landmarks[, c("address", "geometry")], samples)
samples_geo_r <- merge(landmarks_r[, c("address", "geometry")], samples)

# Cases
cases$address <- cases$location
cases_geo <- merge(landmarks[, c("address", "geometry")], cases)
cases_geo_r <- merge(landmarks_r[, c("address", "geometry")], cases)

# By stall (one line per stall address)
allStalls_geo <- samples_geo[!duplicated(samples_geo$address), ]
allStalls_geo_r <- samples_geo_r[!duplicated(samples_geo_r$address), ]




