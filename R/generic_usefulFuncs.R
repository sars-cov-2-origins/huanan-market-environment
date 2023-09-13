# Function to make a color transparent
# Source https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
makeTransparent<-function(someColor, alpha=100, noRescaling = FALSE)
{
  if(alpha<=1 & !noRescaling){
    cat("I am assuming that you entered alpha on a [0, 1] scale, and converting it to [0, 255]. Add `noRescaling = TRUE` to prevent rescaling.\n")
    alpha <- 255 * alpha
  }
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#..........................................................
# Function to plot common and lating names
# via https://stackoverflow.com/questions/29943251/displaying-values-from-a-character-vector-as-italic-labels-in-boxplot-in-r

commonAndLatin <- function(vecCommon, vecLatin){
  # atop to stack
  x <- as.expression(lapply(seq_along(which(ii)), function(i) bquote(.(vecCommon[i]) ~ "(" * italic(.(gsub("\\.", " ", vecLatin[i]))) * ")")))
  x
}

#...........................................................
# Function to plot a palette
plotpal <- function(pal){
  n <- length(pal)
  plot(1:n, rep(0, n), pch = 15, col = pal, cex = 4)
}

#...........................................................
# Add an image to a plot and preserve aspect ratio
# via https://stackoverflow.com/questions/27800307/adding-a-picture-to-plot-in-r
addImg <- function(
    obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
    x = NULL, # mid x coordinate for image
    y = NULL, # mid y coordinate for image
    width = NULL, # width of image (in x coordinate units)
    interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}

#...........................................................
# Sum ignoring NAs
sumIgnoreNA <- function(v){
  sum(v, na.rm = TRUE)
}

#...........................................................
# Length ignoring NAs
lengthIgnoreNA <- function(v){
  length(v[!is.na(v)])
}

#...........................................................
# Number of non zero and non NA
sumPos <- function(v){
  sum(v > 0, na.rm = TRUE)
}

#...........................................................
# Unique ignoring NAs
uniqueIgnoreNA <- function(v){
  unique(v[!is.na(v)])
}

#...........................................................
# Define function to convert geometry coordinates to matrix
geometry2matrix <- function(data){
  matrix(unlist(data$geometry), ncol = 2, byrow = TRUE)
}

#...........................................................
# Clear viewer
# https://stackoverflow.com/questions/49217915/how-to-clear-all-charts-from-viewer-pane-in-rstudio
# install.packages("rstudioapi")
clear_viewer_pane <- function() {
  dir <- tempfile()
  dir.create(dir)
  TextFile <- file.path(dir, "blank.html")
  writeLines("", con = TextFile)
  rstudioapi::viewer(TextFile) 
}
