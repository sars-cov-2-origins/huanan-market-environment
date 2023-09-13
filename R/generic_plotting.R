library(MetBrewer) # For colors as well

if(!exists("samples")){
  stop("Need to run load_data-sequencing.R first!")
}

# Define colors ####
# Species colors

# Base palette
thepal <- met.brewer("Cross", 9, "discrete"); thepal

# Define palette
# (hard-coded numbers based on the chosen palette to maximize contrast)
# Human
colHuman <- thepal[6]
# SARS-CoV-2
colSC2 <- gray(0.3)
# Key species
colKey <- colorRampPalette(thepal[1:3])(length(keySpecies)-1)[c(1, 5, 2, 6, 3, 7, 4, 8, 9)] # Shuffle for contrast
plotpal(colKey)
# Other Mammals
colOMam <- colorRampPalette(thepal[4:5])(length(spMammals) - (length(keySpecies)))
plotpal(colOMam)
# Other species
colOth <- colorRampPalette(thepal[7:9])(length(spMZ) - (length(spMammals)))
plotpal(colOth)
# Join them all in a single palette
colSp <- c(colSC2, colHuman, colKey, colOMam, colOth)
names(colSp) <- c("SC2", keySpecies, spOMam, spOth)
#  Show the resulting palette
plotpal(colSp)

#.........................................................

# Colors for SC2 positivity results (qualitative)
tmppal <- MetBrewer::met.brewer("Hiroshige", 11)
col2 <- tmppal[1]
col1 <- tmppal[3]
col0 <- tmppal[7]
colPCRSeq <- c(col0, col1, col2)
names(colPCRSeq) <- as.character(0:2)
colPCRSeq
rm(tmppal)

#..........................................................
# Market map details
colMarketBackground <- gray(0.95)
colWildlifeStall <- gray(0.65)
colStallBorder <- gray(0.975)

# Color of the border of bubbles
borderBubble <- gray(0)
opacityBubble <- 0.8

marBubblePlot <- c(2, 2, 3, 5) # Margins of bubble plot (leaving space for legends)

insetxval <- -0.13 # x position of the legend
insetyval <- 0.2 # y position of the second legend
txtwd <- 20 # Width of the legend text (to calibrate box widths)
ttladj <- 0.5 # Adj ot the legend title
cexLegend <- 0.9 # Relative size of the legend text

lwdMarketBoundary <- 1.5

opWildlifeStall <- 1 # Opacity wildlife stalls
weightWildlifeStall <- 1 # Line weight wildlife stalls

# Palette of the p value heatmap in the market
col0 <- colMarketBackground
baseColorHeatmap = "orange"
baseColorBubbles <- "red"


# Bubble plot pdf size
wpdf <- 6.5 # Width
hpdf <- 6 # Height

wpdf.maps <- 3.31
hpdf.maps <- 3.1
pointsize.maps <- 8
dLeg <- 0.42 # mai added if legend is plotted
mai.maps <- c(0.01, 0.01, 0.17, 0.01) # mai for market map


## |- Parameters for figures over time

#  Colors
colsType <- met.brewer("Java", 3)#[c(1, 2, 4)]
names(colsType) <- unique(samples$Sample.type)[c(1, 3, 2)]
op.overtime <- 100 # Opacity
#  Pch
pchType <- c(21, 22, 23)
names(pchType) <- names(colsType)

# Plot size and options
wpdf.overtime <- 7
hpdf.overtime <- 6
mai.overtime <- 0.25*c(6, 4.5, 2, 3.25)
mgp.overtime <- c(3, 0.5, 0)
tck.overtime <- -0.01
cexa.overtime <- 0.9


# Functions ####

## |- Common plotting functions

# Max number of samples, for the combinations that we consider (Figure 2)
#  Compute max numbers of samples in the different combinations
#  Jan 01-12 samples
tmp <- samples[which(samples$Sampling.date <= "2020-01-12"), ]
agg <- aggregate(tmp$Sample.ID, by = list(tmp$address), FUN = length)
m1 <- max(agg$x)
#  Drains ans sewage
tmp <- samples[which(samples$Sample.type != "Environmental swab"), ]
agg <- aggregate(tmp$Sample.ID, by = list(tmp$address), FUN = length)
m2 <- max(agg$x)
# Max of max, to have the same scale in all figures
maxSamples <- max(c(m1, m2))


# Prepare to plot bubbles
# Function to get their color, when plotting proportion positives
colProp <- function(p, baseCol, ncol, opacity = opacityBubble){
  # p: value of proportion (single value)
  # baseCol: base color for the bubbles
  # ncol: number of different colors
  # opacity: global opacity of the bubbles
  
  # Define palette vector from white to base color, with ncol colors, and opacity as chosen
  thepal <- makeTransparent(colorRampPalette(c("white", baseCol))(ncol+2)[-c(2, 3)], opacity * 255)
  # Make the white color more transparent (and the white color only)
  thepal[1] <- makeTransparent(thepal[1], 0.5*255)
  # Vector of boundary values
  seqCol <- seq(0, 1, length.out = ncol) 
  # Place p on this scale
  ip <- sum(p > seqCol) + 1
  # Return the corresponding color
  list(col = thepal[ip], pal = thepal)
}

# Set min and max values (-log10(p)) as global variables
minpp <- 3.8
maxpp <- 7.75

# Function to get their color, when plotting proportion reads or total reads
colFreq <- function(p, lightestGrayLevel = 1, pal = NA, opacity = 1, plotTotalReads = FALSE){
  # p: proportion reads or absolute number of reads
  # lightestGrayLevel: value of the level of the lightest gray color
  # pal: type of palette
  # minpp: min value on the plot
  # maxpp: max value on the plot
  # opacity: as the name indicates
  # plotTotalReads: whether we plot total reads (TRUE) or proportion reads (FALSE)
  
  # If no palette is provided, define one
  if(all(is.na(pal))) pal <- met.brewer("VanGogh3", 101, direction = -1)
  # Other: gray(seq(0, 1, length.out = 101))
  
  # Initialize output: white when p = 0
  out <- rep("white", length(p))
  # Compute log-transformed values, for p>0
  if(!plotTotalReads){ # Proportions -> -log10
    pp <- -log(p[p > 0], 10)
  }else{ # Absolute numbers -> log10
    pp <- log(p[p > 0], 10) 
  }
  stopifnot(all(pp >= minpp) & all(pp <= maxpp))
  
  # Compute color for p>0, based on (-)log10 values
  out[p > 0] <- pal[1 + round(100 * (pp - minpp)/(maxpp - minpp) * lightestGrayLevel)]
  # Apply opacity change
  out <- makeTransparent(out, 255 * opacity)
  out
}


# Vectorialize function
getSampleCol <- function(v, baseCol = "red", ncol = 21){
  # v: vector of proportions
  vapply(seq_along(v), function(i) colProp(v[i], baseCol, ncol)$col, "c")
}

# Function to get bubble size
# Parameters for point sizes
# NB: cex changes the *radius* of a point. We work with sqrt for it to control surface instead.
#minCexSample <- 0.75^2 # Min sample point size (1 sample)
#maxCexSample <- 3^2 # Max sample point size (maxSamples) - big to fit the large number of samples
# Compute cex
getSampleCex <- function(v, maxsamp = maxSamples, 
                         minCexSample = 0.75^2, # Min sample point size (1 sample)
                         maxCexSample = 2.75^2 # Max sample point size (maxSamples) - big to fit the large number of samples
){
  out <- sqrt(minCexSample + (maxCexSample - minCexSample) * ((v - 1) / (maxsamp - 1)))
  # Make sure to remove those for which we do not have samples
  out[v == 0] <- NA
  c(out)
}

getStallData <- function(spl){
  spl[!duplicated(st_drop_geometry(spl$address)), ]
}

#........................................................
# Define color of the positivity heatmap
getPalHeatmap <- function(col0, baseColorHeatmap, opacityFirstColor = 0.2, doubleFirst = TRUE, ncol = 7){
  #  First gradient
  tmp <- colorRampPalette(c(col0, baseColorHeatmap), alpha = TRUE)(ncol)
  if(doubleFirst){
    #  Palette doubling the background color, made transparent
    colorRampPalette(c(makeTransparent(tmp[1], 255 * opacityFirstColor), tmp[1:ncol]), alpha = TRUE)
  }else{
    colorRampPalette(c(tmp[1:ncol]), alpha = TRUE)
  }
}

palHeatmap <- getPalHeatmap(makeTransparent(colMarketBackground, 1*255), baseColorHeatmap, opacityFirstColor = 1)


# Load animal and SC2 pictures ####

# Load SC2 picture (Font Awesome)
SC2pic <- png::readPNG("../figs/icons/virus.png")

# Load saved phylopics
load("../figs/icons/phylopics.RData")

