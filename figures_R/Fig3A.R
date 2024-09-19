# Plot composition heatmap for samples collected on Jan 12, 
# focusing on animals described in Xiao et al. (2021)

# Load data if necessary
if(!exists("samples")){
  source("load_data-sequencing.R")
}

if(!exists("thepal")){
  source("generic_plotting.R")
}

addPhylopics <- TRUE # Whether to add images
addHomo <- TRUE # Whether to plot Homo sapiens
homoSeparate <- FALSE # Whether Homo sapiens is kept separate
orderByReads <- TRUE # Whether to order by total reads (vs. by nb of stalls in which species is present)
colorReadNb <- TRUE # Whether the color scale shows absolute read numbers (instead of proportion)
addGrayBG <- FALSE # Whether to add gray backgroun to better see lines and columns
commonNames <- TRUE # Whether to put common names (TRUE) or latin names (FALSE)
removeBR <- TRUE # Whether to remove the other bamboo rats (Xiao may have misindentified the species)
usePseudo <- TRUE # Whether to pseudonymize the stalls


thresholdAbundance <- 200 # Threshold number of total reads to be considered for inclusion

if(addPhylopics){
  library(rphylopic) # For animal silhouettes
}

# Subset of the samples data for 12-Jan
s12 <- samples[which(samples$Sampling.date == "2020-01-12" & !is.na(samples$SARS2.reads)), ]

# Sort samples by stall and then by SC2 positivity (horizontal order)
iorder <- order(s12$address, -s12$PCRSeqPos)
ss12 <- s12[iorder, ]
n <- nrow(ss12)

# Load Xiao species list
xiaoSp <- read.csv("../metadata/Xiao_etal_animalList.csv")
orderedSpecies <- gsub(" ", "\\.", xiaoSp$Latin_name)
# Check which species are absent
missingSp <- orderedSpecies[!is.element(orderedSpecies, names(ss12))]
missingSp
# Get the genera
generaX <- vapply(orderedSpecies, function(x) strsplit(x, "\\.")[[1]][1], "x")
generaMZ <- vapply(spMZ, function(x) strsplit(x, "\\.")[[1]][1], "x")
missingG <- generaX[!is.element(generaX, generaMZ)]
missingG
sort(spMZ[is.element(generaMZ, generaX)])

# Xiaos that we have
xiaoSp2 <- xiaoSp$species[is.element(xiaoSp, spMZ)]
# Including species for which we have genus
tmp <- generaX[is.element(generaX, generaMZ)]
tmp2 <- (spMZ[is.element(generaMZ, tmp)])
genera2 <- vapply(tmp2, function(x) strsplit(x, "\\.")[[1]][1], "x")
xiaoSp3 <- names(genera2[order(match(genera2, generaX))])

# xiaoSp2 has Xiao's species, stricly
# xiaoSp3 additionally includes species for which the genus is in Xiao
orderedSpecies <- xiaoSp3

if(removeBR){
  orderedSpecies <- orderedSpecies[!is.element(orderedSpecies, c("Rhizomys.sinensis", "Rhizomys.sumatrensis"))]
}

# Identify species with a minimum number of reads over the samples we are considering
vv <- apply(ss12[, orderedSpecies], 2, sumIgnoreNA) # Compute total reads in the sample
sort(vv)
length(vv)
keepSpecies <- names(vv[which(vv > thresholdAbundance)]) # Names of the species that we are keeping
length(keepSpecies)

if(addHomo & !homoSeparate){
  keepSpecies <- c("Homo.sapiens", keepSpecies)
}

# Reorder by abundance
if(orderByReads){
  # Compute total number of reads by species
  vv <- apply(ss12[, keepSpecies], 2, sumIgnoreNA) 
}else{
  # Compute number of samples in which the species is detected
  vv <- apply(ss12[, keepSpecies], 2, sumPos) 
}

iorder <- order(-vv)
# Reorder species
keepSpecies <- keepSpecies[iorder]

if(addHomo & homoSeparate){
  keepSpecies <- c("Homo.sapiens", keepSpecies)
}

# Min and max values, to compute color scale
# Max proportion (1)
maxp <- max(ss12[, keepSpecies] / ss12$Read_pairs_after_trimming)
if(!colorReadNb){ # Work with proportions
  # All proportions
  tmp <- ss12[, keepSpecies] / ss12$Read_pairs_after_trimming
  # Log-transformed positive values
  tmpl <- -log(tmp[tmp>0], 10)
  rm(tmp)
}else{ # Work with absolute values
  # All values 
  tmp <- ss12[, keepSpecies] 
  # Log-transformed positive values
  tmpl <- log(tmp[tmp>0], 10)
  rm(tmp)
}
# Define color scale using minpp and maxpp
colFreq <- function(p, lightestGrayLevel = 1, pal = NA, minpp = min(tmpl), maxpp = max(tmpl), opacity = 1, plotTotalReads = colorReadNb){
  # p: proportion reads or absolute number of reads
  # lightestGrayLevel: value of the level of the lightest gray color
  # pal: type of palette
  
  if(all(is.na(pal))) pal <- met.brewer("VanGogh3", 101, direction = ifelse(colorReadNb, 1, -1))
  # Other: gray(seq(0, 1, length.out = 101))
  
  # Initialize output: white when p = 0
  out <- rep("white", length(p))
  # Compute log-transformed values, for p>0
  if(!plotTotalReads){ # Proportions -> -log10
    pp <- -log(p[p > 0], 10)
  }else{ # Absolute numbers -> log10
    pp <- log(p[p > 0], 10) 
  }
  # Compute color for p>0, based on (-)log10 values
  out[p > 0] <- pal[1 + round(100 * (pp - minpp)/(maxpp - minpp) * lightestGrayLevel)]
  # Apply opacity change
  out <- makeTransparent(out, 255 * opacity)
  out
}

# x positions of the rectangles
dx <- 0.4 # Width of the rectangles
dstall <- 0.8 # Space between rectangles of different stalls
dy <- 1 # Height of the rectangles
dSC <- 0.4 # Space between rectangles of different types of species
dSp <- 0.15 # Space between rectangles of different species
lwdRec <- 0.25 # Line width of the rectangles
spacing <- "       " # Space before species names
cexLeg <- 0.675

# x values: 
xvals <- dx * (1:n) # x values
#   Add space for stalls (hard-coded for Jan 12)
xvals <- xvals + rep(dstall * ((1:7) - 1), each = 10)

# y values: 
yvals <- -seq_along(keepSpecies) * dy
if(addHomo & homoSeparate){
  #   Add space for species of different groups
  ddy <- c(0, rep(- dSC, length(keepSpecies) - 1))
}else{
  ddy <- cumsum(rep(- dSp, length(keepSpecies)))
}
yvals <- yvals + ddy


fname <- "../figs/Fig3A.pdf"
pdf(fname, width = 5.5, height = 3.1)
par(mar = c(0.5, 0.1, 0.5, 7.5))
# Initialize plot
plot(0, xlim = c(xvals[1], ceiling(max(xvals))), ylim = c(-dy * length(keepSpecies), 2), axes = FALSE, xlab = "", ylab = "", type = "n", asp = 1)

if(addGrayBG){
  # Add gray background to better see lines and columns
  colGray <- gray(0.9)
  col2 <- gray(1)
  colsX <- rep(c(colGray, col2), ceiling(length(xvals)/2))[1:length(xvals)]
  colsY <- rep(c(colGray, col2), ceiling(length(yvals)/2))[1:length(yvals)]
  #for(i in seq_along(xvals)){
  rect(xleft = xvals, xright = xvals + dx, 
       ybottom = min(yvals), ytop = dy + dSC, 
       col = colsX, lwd = -1)
  #}
}

# Plot rectangles for SC2
rect(xleft = xvals, xright = xvals + dx, 
     ybottom = rep(0, n) + dSC, ytop = rep(dy, n) + dSC, col = colPCRSeq[as.character(ss12$PCRSeqPos)], lwd = lwdRec)
# Add read numbers of SC2
iR <- which(ss12$SARS2.reads > 0) # Indices of samples with non zero SC2 reads
cexReads <- 0.15 # size of reads labels
text(x = xvals[iR] + dx/2, y = dSC + dy/2, labels = ss12$SARS2.reads[iR], 
     col = "white", cex = cexReads, adj = c(0.5, 0.5), 
     srt = 90)

# Add text legend for SC2
par(xpd = TRUE)
text(x = max(xvals), y = dy/2 + dSC, labels = paste0(spacing, "SARS-CoV-2"), cex = cexLeg, adj = c(0, 0.5))

# Stall labels
#   x positions
ixlbl <- ((1:7)-1) * 10 + 5 
xlbl <- (xvals[ixlbl] + xvals[ixlbl + 1])/2
if(usePseudo){
  pseudos <- read.csv("../metadata/stall_pseudonyms.csv")
  tmp <- pseudos$pseudonym
  names(tmp) <- pseudos$address
  lbl <- paste("Stall", tmp[unique(ss12$address)]) 
}else{
  lbl <- gsub("West\\|", "", unique(ss12$address))
}
#   text
text(x = xlbl, y = 3.5*dy + dSC, adj = c(0.5, 0), labels = lbl, cex = cexLeg)

# Sample labels
tt <- ss12$Lab.code
ii <- which(nchar(tt) > 4)
tt[ii] <- s12[ii, "Sample.ID"]
text(x = xvals + dx/2, y = dy + dSC, adj = c(0, 0.5), labels = paste0("   ", tt), cex = 0.4*cexLeg, srt = 90)

par(xpd = FALSE)

# Plot rectangles for species
for(i in seq_along(keepSpecies)){
  # Get the vector of what we are going to plot as colors
  if(!colorReadNb){ # Proportion reads
    pcol <- ss12[, keepSpecies[i]] / ss12$Read_pairs_after_trimming
  }else{ # Absolute number of reads
    pcol <- ss12[, keepSpecies[i]]
  }
  # Plot the rectangles
  rect(xleft = xvals, xright = xvals + dx, 
       ybottom = rep(yvals[i], n), ytop = rep(yvals[i] + dy, n), 
       #col = c("white", "black")[1 + (ss12[, keepSpecies[i]] > 0)], 
       col = colFreq(pcol),
       lwd = lwdRec)
  # Add read numbers
  if(colorReadNb){
    iR <- which(ss12[, keepSpecies[i]] > 0)
    cols <- rep(gray(0.7), length(iR))
    cols[ss12[, keepSpecies[i]][iR] < 1000] <- gray(0.2)
    text(x = xvals[iR] + dx/2, y = rep(yvals[i] + dy/2, length(iR)), 
         labels = (ss12[, keepSpecies[i]])[iR],
         col = cols, cex = cexReads, adj = c(0.5, 0.5), srt = 90)
  }
  
  # Add legend for species
  par(xpd = TRUE)
  if(commonNames == TRUE){
    tmp <- dicoSp[keepSpecies[i]] # Species name, which we will capitalize
    splbl <- paste0(spacing, toupper(substr(tmp, 1, 1)), substr(tmp, 2, nchar(tmp)))
  }else{
    splbl <- paste0(spacing, gsub("\\.", " ", keepSpecies[i]))
  }
  text(x = max(xvals), y = yvals[i] + 0.5*dy, labels = splbl, adj = c(0, 0.5), cex = cexLeg, font = ifelse(commonNames, 1, 3)) # italized if latin names
  par(xpd = FALSE)
}

if(addPhylopics){
  # Add images
  phypics <- phylopics[keepSpecies] 
  # Add phylopics
  sizepic <- 2.5*dx
  maxheight <- 1*dy
  xpic <- max(xvals) + dstall + sizepic/2
  for(i in seq_along(keepSpecies)){
    cat(i, "")
    sp <- keepSpecies[i]
    # Compute size of the picture, based on its aspect ratio, 
    xratio <- phylopics[[sp]]@summary@xscale[2]/phylopics[[sp]]@summary@yscale[1]
    #  and of the maximal height
    ysz <- ifelse(sizepic/xratio > maxheight, maxheight, sizepic/xratio)
    # Add the image
    add_phylopic_base(img = phylopics[[sp]], x = xpic, y = yvals[i] + dy/2, fill = "black", color = NA, ysize = ysz)
  }
  # Add SC2 image
  par(xpd = TRUE)
  addImg(SC2pic, x = xpic, y = dSC + dy/2, width = sizepic)
  par(xpd = FALSE)
}

dev.off()
system(paste0("open ", fname))

