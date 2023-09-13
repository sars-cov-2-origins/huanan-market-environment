# Plot: Barplots of species composition

# Load data if necessary
if(!exists("samples")){
  source("load_data-sequencing.R")
}

if(!exists("thepal")){
  source("generic_plotting.R")
}

library("rphylopic", quietly = TRUE)

# Function to extract the values of the key species and of ntop top non-key species
getTops <- function(sampleID, data = samples, ntop = 3, refProp = "MZ", keySp = keySpecies){
  # sampleID: sampleID of the sample to plot
  # data: dataset to use
  # ntop: number of "other" species to plot
  # refProp: reference relative to which proportions are calculated (MZ or Mammals); also affects the other species plotted
  # keySp: vector of names of key species, i.e. species which are always plotted
  
  # Subset of the data with this lab code
  subdat <- data[which(data$Sample.ID == sampleID), ]
  
  # Get numbers of reads as vector
  # Define reference set
  if(refProp == "Mammals"){
    spp <- spMammals
  }
  if(refProp == "MZ"){
    spp <- spMZ
  }
  
  abundances <- unlist(subdat[spp]) #for metazoa use spMZ
  # Corresponding proportions
  proportions <- abundances / sum(abundances)
  
  if(ntop > 0){
    # Find the top ntop among other MZ
    #propFilt <- proportions[proportions>0]
    #propFilt <- propFilt[!sapply(propFilt, is.na)]
    others <- spp[!is.element(spp, keySp)] # Species among spp which are not the keyspecies
    topOthers <- names(sort(unlist(subdat[others]), decreasing = TRUE))
    
    topOthers <- topOthers[!duplicated(topOthers)][1:ntop]
    topOthers <- topOthers[!sapply(topOthers, is.na)]
    
    #for plotting:
    
    # From this we keep: the key species and the top others
    skeep <- c(keySp, topOthers) 
  }else{
    skeep <- keySp
  }
  
  skeep <- skeep[!duplicated(skeep)]
  # Output
  list(proportions = proportions[skeep],
       readcounts = abundances[skeep],
       totreads = sum(abundances), 
       SC2 = subdat$SARS2.reads, 
       ntop = ntop)  
}

# Get the max numbers of reads
maxread <- max(samples[, spMZ], na.rm = TRUE) # Among MZ
maxtotread <- max(c(maxread, max(samples$SARS2.reads, na.rm = TRUE))) # Including SC2


# Function to plot the output as barplot, read counts
plotBarReads <- function(ll){
  # ll: output of the getTops function
  
  # Combine the results in a single vector, and rev to put SC2 on top
  ov <- v <- rev(c("SC2" = ll$SC2, ll$readcounts))
  ov[is.na(v)] <- "0"
  # Parameters
  xmin <- 0.7 # Min x value, affects how the results look on the log scale
  v[v == 0] <- xmin # Change 0 into xmin because we plot on log scale
  names(v)[is.na(v)] <- "" # Rename NAs
  v[is.na(v)] <- xmin
  
  cexleg <- 0.9 # Cex of the text of the legends
  
  par(mgp = c(0.3, 0.3, 0), 
      tck = - 0.0125, cex.axis = 1, cex.lab = 1, cex = 1, 
      #mar = c(1.5, 11.5, 2.5, 2.5), 
      las = 1, xpd = FALSE)
  
  colpal <- colSp[names(v)] # Color palette
  
  # Plot results as horizontal barplot, with same scale for all figures via maxtotread
  bp <- barplot(v, names.arg = NA, col = colSp[names(v)], 
                horiz = TRUE, xlim = c(xmin, maxtotread), log = "x", border = 1, axes = FALSE, 
                space = 1.25)
  # Add axes and text legends
  axis(1, pos = 0, lwd = 0, lwd.ticks = 1, at = c(10^(0:6)), labels = c("1", "10", "100", "1000", "10000", "100000", "10^6"), cex.axis = cexleg)
  # Add 0 separately so that we have both 0 and 1 (otherwise R removes one)
  axis(1, pos = 0, lwd = 0, lwd.ticks = 1, at = xmin, labels = "0", cex.axis = cexleg)
  par(xpd = TRUE)
  spaceleg <- "  " # Space between name and bar
  # Species names, in the left margin
  text(x = xmin, y = bp[-length(bp)], labels = paste0(gsub("\\.", " ", names(v)[-length(bp)]), spaceleg), adj = 1, font = 3, cex = cexleg, col = colSp[names(v)[-length(bp)]])
  text(x = xmin, y = bp[length(bp)], labels = paste0("SARS-CoV-2", spaceleg), adj = 1, font = 1, col = colSp["SC2"], cex = cexleg)
  # Add read counts, in the right of the bars
  text(x = v, y = bp, labels = paste0(" ", ov), adj = 0, font = 1, cex = 0.75, col = colSp[names(v)])
  par(xpd = FALSE)
  
  out <- list(speciesName = names(v), ypos = c(bp), readCount = v, colpal = colpal)
}  

# Get max read proportions
maxreadProp <- 0. # Among mammals
maxtotreadProp <- 0. # Including SC2

# Function to plot the output as barplot
plotBarReadsProp <- function(ll){
  # ll: output of the getTops function
  
  # Combine the results in a single vector, and rev to put SC2 on top
  ov <- v <- rev(c("SC2" = ll$SC2, ll$readcounts))
  ovProp <- vProp <- log10(rev(c("SC2" = ll$SC2/ll$totreads, ll$proportions)))
  
  offset <- 4.5
  # Parameters
  xmin <- -offset #0.7 # Min x value, affects how the results look on the log scale
  vProp[v == 0] <- xmin
  v[v == 0] <- 0.7 # Change 0 into xmin because we plot on log scale
  cexleg <- 1 # Cex of the text of the legends
  
  #skinny
  par(mgp = c(0.3, 0.3, 0), 
      tck = - 0.0125, cex.axis = 1, cex.lab = 1, cex = 1, 
      mar = c(1.5, 11.5, 2.2, 2.5), 
      las = 1, xpd = FALSE)
  
  # Define color palette
  if(ll$ntop > 0){
    colpal <- colSp[names(v)]
  }else{
    # Palette with alternating colors
    #    colpal <- c(rep(thepal[1:2], length(v))[1:(length(v)-2)], colSp[c("Homo.sapiens", "SC2")]) # SC2, Human, and then the mix of species
    #    names(colpal) <- names(v)
    colpal <- colSp[names(v)]
  }
  
  # Plot results as horizontal barplot, with same scale for all figures via maxtotread
  bp <- barplot(vProp + offset, names.arg = NA, col = colpal, 
                horiz = TRUE, xlim = c(xmin + offset, maxtotreadProp + offset), border = 1, axes = FALSE, 
                space = .6)
  # Add axes and text legends
  xx <- c((xmin + offset + 0.5):(0 + offset))
  axis(1, pos = 0, lwd = 0, lwd.ticks = 1, at = xx, labels = c(as.expression(sapply(head(xx, -1), function(x) bquote(italic(10^.(x - offset))))), 1), cex.axis = cexleg)
  # Add 0 separately so that we have both 0 and 1 (otherwise R removes one)
  axis(1, pos = 0, lwd = 0, lwd.ticks = 1, at = xmin + offset, labels = "0", cex.axis = cexleg)
  par(xpd = TRUE)
  spaceleg <- "  " # Space between name and bar
  # Species names, in the left margin
  tmpnames <- dicoSp[names(v)[-length(bp)]] # Species names (will be capitalized)
  tmpcol <- colpal[-length(v)] # Define colors
  tmpcol[vProp[-length(vProp)] == xmin] <- gray(0, 0) # Remove names when no reads
  text(x = xmin+offset, y = bp[-length(bp)], labels = paste0(vapply(tmpnames, function(v){paste0(toupper(substr(v, 1, 1)), substr(v, 2, nchar(v)))}, "x"), spaceleg), adj = 1, font = 1, cex = cexleg, col = 'black')
  text(x = xmin+offset, y = bp[length(bp)], labels = paste0("SARS-CoV-2", spaceleg), adj = 1, font = 1, col = 'black', cex = cexleg)
  # Add read counts, in the right of the bars
  text(x = vProp+offset, y = bp, labels = paste0(" ", ov), adj = 0, font = 1, cex = 0.9*cexleg, col = colpal)
  par(xpd = FALSE)
  
  out <- list(speciesName = names(v), ypos = c(bp), readCount = v, proportions = vProp, colpal = colpal)
}  

plotSpeciesComposition <- function(sID, png = FALSE, pdf = TRUE, format = 'readCount', ntop = 3, refProp = "MZ", keySp = keySpecies){
  # Metadata for this sample ID
  md <- samples[which(samples$Sample.ID == sID), ]
  # Plot
  fname <- paste0("../figs/Fig2_", gsub("#", "h", sID), ".png")
  if(format=='prop'){
    if(png) png(fname, width = 5, height = 4, units = "in", res = 200)
    if(pdf) pdf(gsub("png", "pdf", fname), width = 3.35, height = 2.7, pointsize = 8)
    pbr <- plotBarReadsProp(getTops(sID, ntop = ntop, refProp = refProp, keySp = keySp))
  }else{
    if(png) png(fname, width = 8, height = 4, units = "in", res = 200)
    if(pdf) pdf(gsub("png", "pdf", fname), width = 8, height = 4)
    pbr <- plotBarReads(getTops(sID, ntop = ntop, refProp = refProp, keySp = keySp))
  }
  # Add title 
  # Text for PCR information
  if(is.na(md$PCRPos)){
    pcrInfo <- "PCR information unavailable"
  }else{
    if(md$PCRPos == 1){
      ct <- md$CT
      ctInfo <- ifelse(ct == "+", "not available", paste0("= ", ct))
      pcrInfo <- paste0("PCR-positive (Ct ", ctInfo, ").")
    }else{
      pcrInfo <- "PCR-negative."
    }
  }
  if(format == "readCount"){
    ln <- 0.5 # Line at which title is plotted
    mtext(side = 3, paste0(md$Stall_corrected, "\n", md$Sample.ID, " (", md$Lab.code, "), ", md$Sample.information, ""), line = ln)
    mtext(side = 3, paste0("Sampled ", md$Sampling.date, ", ", pcrInfo), line = ln - 0.9, cex = 0.8)
  }else{
    ln <- 0.8 # Line at which title is plotted
    mtext(side = 3, paste0(md$Lab.code, " (", md$Sample.information, ")"), line = ln)
    mtext(side = 3, paste0(ifelse(md$PCRPos == 1, "PCR+", "PCR-"), ", ", "NGS+"), line = ln - 1, cex = 1, col = ifelse(md$PCRPos == 1, colPCRSeq["2"], colPCRSeq["1"]))
  }
  # Add annotations
  par(xpd = TRUE)
  xpic <- -1.55*10^5
  xsize <- 2
  
  # Get vector of uuids for the species we are considering
  # (-1: removing SC2 as it is not plotting with phylopic)
  vuuid <- uuids[head(pbr$speciesName, -1)]
  # Keep only those for which read count is strictly positive
  ikeep <- which(floor(pbr$readCount[names(vuuid)]) > 0)
  # If there are images to plot
  if(length(ikeep) > 0){
    vuuid <- vuuid[ikeep]
    
    # Create new plot for annotations
    ypos <- pbr$ypos[ikeep]
    
    ylm <- par("usr")[3:4] # Extract ylim of the current plot
    par(new = TRUE)
    # Empty new plot
    plot(seq_along(ypos), ypos, xlim = c(0, 1), axes = FALSE, xlab = "", ylab = "", ylim = ylm, xaxs = "i", yaxs = "i", type = "n")
    pusr <- par("usr")
    xyscale <- 1.65*(pusr[4] - pusr[3])/(pusr[2] - pusr[1])
    
    if(format=='prop'){
      xpic <- -.77 
      xsize <- 0.045
      
    }else{
      xpic <- -0.41
      xsize <- 0.035
    }
    ysize <- xsize * xyscale
    if(pbr$readCount["SC2"] >= 1){
      # Add SC2 picture
      addImg(SC2pic, x = xpic, y = tail(pbr$ypos, 1), width = xsize*1.5)
    }
    
    # Add phylopics
    # We can only set its height, but we need max width as well
    sizepic <- 2.5
    maxheight <- pbr$ypos[2] - pbr$ypos[1]
    cat("Plotting phylopics... ")
    for(i in seq_along(vuuid)){
      cat(i, "")
      sp <- names(vuuid[i])
      # xy ratio of the phylopic
      xratio <- phylopics[[sp]]@summary@xscale[2]/phylopics[[sp]]@summary@yscale[1] 
      # Define y size based on the width, xyratio, and set maximal height (could also have used min!)
      ysz <- ifelse(sizepic/xratio > maxheight, maxheight, sizepic/xratio)
      # Plot phylopic
      add_phylopic_base(img = phylopics[[sp]], x = xpic, y = ypos[i], color = NA, fill = pbr$colpal[sp], ysize = ysz)
    }
    
  }
  
  par(xpd = FALSE)
  par(new = FALSE)
  if(png | pdf) dev.off()
}

# Special version for Q61-70 samples
cutoffReadCount <- 300
cutoffDetects <- 2
tmp <- samples[which(samples$Lab.code == 'Q61' | samples$Lab.code == 'Q64'| samples$Lab.code == 'Q68' | samples$Lab.code == 'Q69' | samples$Lab.code == 'Q70'), ]

# at least 2 detections
tmpDetect <- tmp[spMammals]>0
speciesOfInterestDetects <- apply(tmpDetect, 2, sum, na.rm=TRUE)

# or total 300 reads across the Q61-Q70 samples
speciesOfInterest <- sapply(tmp[spMammals], sum, na.rm=TRUE)
speciesOfInterest  <- speciesOfInterest [speciesOfInterest >=cutoffReadCount | speciesOfInterestDetects >= cutoffDetects]
speciesOfInterest  <- speciesOfInterest [order(speciesOfInterest, decreasing=TRUE)]
speciesOfInterest  <- names(speciesOfInterest )



# Include all key species, plus all other mammals identified in Q61-70 samples
# speciesOfInterest = speciesOfInterest[!(speciesOfInterest %in% keySpecies)]
#store a temp version while running this
keySpecies0 <- keySpecies
keySpecies2 <- c("Homo.sapiens", speciesOfInterest[!(speciesOfInterest %in% c("Homo.sapiens"))])


tmp <- samples$Sample.ID[which(is.element(samples$Lab.code, c('Q61', 'Q64', 'Q68', 'Q69', 'Q70')))] # Sample ID for the lab codes we want
for(sID in tmp){
  if(!is.na(sID)){
    cat(sID, "")
    plotSpeciesComposition(sID, ntop = 0, format = 'prop', refProp = "Mammals", keySp = keySpecies2)
    cat("\n")
  }
}
keySpecies <- keySpecies0
rm(tmp)

