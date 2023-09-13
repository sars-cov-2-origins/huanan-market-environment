# Barplot of sample positivity (PCR and NGS) on Jan 01 and Jan 12

# Load data if necessary
if(!exists("samples")){
  source("load_data-sequencing.R")
}

if(!exists("thepal")){
  source("generic_plotting.R")
}

# Define function to plot the result
plotSamplesPos <- function(date.range, location = "address", merge4 = TRUE){
  # date.range: either "01-12Jan", or "alldates"
  # location: whether we use "address" or "Stall_corrected_merged"
  # merge4: whether to merge 4|26 and 4|28
  
  newsamples <- samples
  
  if(merge4 == TRUE){
    newsamples[is.element(newsamples$address, c("West|4|26", "West|4|28")), "address"] <- "West|4|26-28"
  }
  
  # Subset the data according to the date range we want
  if(date.range == "01-12Jan"){
    stl <- newsamples[newsamples$Sampling.date <= "2020-01-12", ]
    # Aggregate the results by stall, combining the two dates
    m <- agg1 <- aggregate(cbind(nbPos = stl$isPos, nbPCRPos = stl$PCRPos, nbTot = !is.na(stl$isPos)), by = list(address = stl[, location]), FUN = sum)
    # Define the order in which the samples are displayed
    orderPos <- order(-m$nbPCRPos, -m$nbPos, m$nbTot, m$address)
    theorder <- orderPos
    # Define labels
    lbls <- paste0(m$address[theorder], " ")
  }
  
  if(date.range == "alldates"){
    stl <- newsamples
    # Aggregate the results by stall and date
    m <- agg1 <- aggregate(cbind(nbPos = stl$isPos, nbPCRPos = stl$PCRPos, nbTot = !is.na(stl$isPos)), by = list(address = stl[, location], date = stl$Sampling.date, Sample.type = stl$Sample.type), FUN = sum)
    # Define the order in which the samples are displayed
    orderPos <- order(-m$nbPCRPos, -m$nbPos, m$nbTot, m$address, m$date)
    theorder <- orderPos
    # Define labels
    lbls <- paste0(m$address[theorder], ", ", format.Date(as.Date(m$date[theorder]), "%d%b"), " (", substr(m$Sample.type[theorder], 1, 1), ")", "", " ")
    # Change name of Water drain into drain
    m$Sample.type[m$Sample.type == "Water drain"] <- "Drain"
  }
  
  # Plotting
  fname <- paste0("../figs/positivities_perstall_", date.range, ifelse(location == "Stall_corrected_merged", "_merged", ""), ifelse(merge4, "_stall4merged", ""),".pdf")
  
  if(date.range == "01-12Jan" & location == "address" & merge4 == TRUE){
    fname <- "../figs/FigS1.pdf"
  }
  pdf(file = fname, width = 7.5 * ifelse(date.range == "01-12Jan", 17/15, 1), height = 30 * ifelse(date.range == "01-12Jan", 0.5, 1))
  par(las = 1, mar = c(1, 8, 1.5, 0.5), mgp = c(2, 0.2, 0), tck = -0.01)
  lwdH <- 6
  # Initialize plot
  plot(m$nbTot[theorder], seq_along(theorder), axes = FALSE, type = "n", 
       xlim = c(0, max(m$nbTot)), ylim = range(theorder),
       xlab = "", ylab = "")
  
  xx <- rev(seq_along(theorder)) # Positions of the stalls
  # Add bars for total number of reads, using previous color code
  segments(y0 = xx, y1 = xx, 
           x0 = 0, x1 = m$nbTot[theorder], 
           lwd = lwdH, col = colPCRSeq["0"], 
           lend = "butt")
  # Add bars for number of positive samples (PCR or sequencing)
  segments(y0 = xx, y1 = xx, 
           x0 = 0, x1 = m$nbPos[theorder], 
           lwd = lwdH, col = colPCRSeq["1"], 
           lend = "butt")
  # Add bars for number of positive samples (PCR only)
  segments(y0 = xx, y1 = xx, 
           x0 = 0, x1 = m$nbPCRPos[theorder], 
           lwd = lwdH, col = colPCRSeq["2"], 
           lend = "butt")
  # Graduations
  for(i in 1:20) abline(v = i, col = "white")
  
  par(xpd = TRUE)
  # Stall labels
  text(y = xx, x = 0, adj = 1, srt = 0, cex = 0.7, labels = lbls)
  # Legend
  lgd <- c("PCR+", "PCR-, NGS+", "PCR-, NGS- or no NGS")
  legend(x = 0, y = length(theorder) + 3, fill = rev(colPCRSeq), legend = lgd, 
         horiz = TRUE, pt.cex = 1.5, cex = 1, bty = "n", xjust = 0, yjust = 0, text.width = 0.35*nchar(lgd))
  par(xpd = FALSE)
  # Add axes
  axis(1, pos = 0.5, lwd = 0, lwd.ticks = 1, at = 1:20)
  axis(3, pos = 0.5 + max(theorder), lwd = 0, lwd.ticks = 1, at = 1:20)
  dev.off()
}

# Plot
plotSamplesPos("01-12Jan", "address", merge4 = TRUE)
