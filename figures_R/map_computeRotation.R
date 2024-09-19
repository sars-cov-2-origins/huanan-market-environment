#
# We want the market map to be horizontal. We therefore rotate the coordinates. 
# To compute the rotation angle, we take two points in the same alley, 
# and compute the angle between the straight line linking these two points, and the horizontal axis. 
#
# We also want the two wings to be closer so that the maps take up less space 
# (this does not affect the calculations, done on each wing separately)
#

# Coordinates of two points in the same alley
A <- c(114.256381, 30.619407)
B <- c(114.257078, 30.619572)

# Compute the angle between the points
theta <- -acos((B[1]-A[1])/(sqrt((B[1] - A[1])^2 + (B[2] - A[2])^2)))

# Rotation angle
angle2 <- -theta + 0.03

# Define the rotation matrix
rotationMatrix <- function(angle = theta) matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), byrow = TRUE, ncol = 2)

# Vectors of coordinates
makeVec <- function(x){
  matrix(x, ncol = 1)
}

# Test
rotationMatrix() %*% makeVec(A) 
rotationMatrix() %*% makeVec(B) 

# Function to compute a point's coordinates after the rotation
rotatePoint <- function(x, angle = theta){
  c(rotationMatrix(angle) %*% makeVec(x))
}

# Repeat this on a matrix
rotateCoordinates <- function(m, angle = theta){
  matrix(apply(m, 1, rotatePoint, angle = theta), byrow = TRUE, ncol = 2)
}

# About a point
# Modified from 
# https://stackoverflow.com/questions/31873151/how-rotate-map-in-r
# Rotate an sf geom around a center point. If no center is
# specified then it rotates around the center of the geom.
st_rotate <- function(x, radians = NULL, degrees = NULL, center_coords = NULL){
  # Get angle in rad
  if(!is.null(degrees)){ 
    if(degrees < -360 | degrees > 360) stop('Degrees must be in the range -360 to 360')
    radians <- degrees * pi/180
  }else{
    if(is.null(radians)){
    stop("need to enter radians or degrees!")
    }
  }
  # Compute midpoint if not provided
  if(is.null(center_coords)){
    # Combine all in a single shape
    x_combined <- sf::st_combine(x)
    # Get centroid
    center_coords <- sf::st_centroid(x_combined)
  }
  # Compute transformation matrix
  transform_matrix <- matrix(c(cos(radians), sin(radians), -sin(radians), cos(radians)), 2, 2)
  # Compute new coordinates
  newdata <- x
  for (i in seq_len(nrow(x))){
    st_geometry(newdata[i, ]) <- (x[i, ]$geometry - center_coords) * transform_matrix + center_coords
  }
  newdata
}

# Clean
rm(A, B)

# Function to translate the right wing a bit closer to the left wing
st_translate <- function(geodata, separationPoint = 12719050, deltaX = 35){
  # geodata: the dataset of geographic data
  # separationPoint: we only move points to the right of this value
  # deltaX: distance by which we move the points
  
  newdata <- geodata
  
  for(iline in seq_len(nrow(geodata))){
    # Get x coordinates
    xx <- st_coordinates(geodata[iline, ])[, 1]
    # Identify the ones we want to move
    iToMove <- which(xx > separationPoint)
    # Move matrix (yes, horizontally!)
    mTranslate <- matrix(0, ncol = length(xx), nrow = 2)
    mTranslate[1, iToMove] <- - deltaX
    # Move the points
    st_geometry(newdata[iline, ]) <- (geodata[iline, ]$geometry + mTranslate)
  }
  newdata
}
