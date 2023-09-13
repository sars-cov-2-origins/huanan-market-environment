# Load data if necessary
if(!exists("samples")){
  source("load_data-sequencing.R")
}

if(!exists("thepal")){
  source("generic_plotting.R")
}

if(!exists("risk_PCR_Jan112_r")){
  source("map_spatialAnalyses.R")
}

# Color scales ####



# Plot the base market map ####

plotBaseMap <- function(plotStreetNumbers = FALSE, plotStallNumbers = FALSE, colNumbers = NA, factorCex = 1){
  # plotStreetNumbers: whether to plot street numbers
  # plotStallNumbers: whether to plot stall numbers
  # colNumbers: color of the number
  # factorCex: multiplying factor for cex, because not the same aspect with pdf
  
  par(new = FALSE, xpd = FALSE, 
      mar = marBubblePlot)
  # Plot background market
  plot(market_boundary_r$geometry, col = colMarketBackground, lwd = lwdMarketBoundary, border = NA)
  # Plot stalls
  plot(map_stalls_r$geometry, col = colStallBorder, add = TRUE)
  # Add wildlife stalls
  plot(wildlifeStalls_r$geometry, col = colWildlifeStall, border = colStallBorder, 
       weight = weightWildlifeStall, opacity = opWildlifeStall, add = TRUE)
  
  if(is.na(colNumbers)){
    colL <- colStallBorder
  }else{
    colL <- colNumbers
  }
  
  if(plotStreetNumbers){
    # Add street numbers
    text(geometry2matrix(streetLabels_r)[, 1], geometry2matrix(streetLabels_r)[, 2], labels = streetLabels_r$title, col = colL, cex = factorCex*0.85, font = 2)
  }
  if(plotStallNumbers){
    # Add stall numbers
    text(geometry2matrix(stallLabels_r)[, 1], geometry2matrix(stallLabels_r)[, 2], labels = stallLabels_r$title, col = colL, cex = factorCex*0.7)
  }
}

# Function to plot bubbles ####


plotBubbleMap <- function(stallNumerators, stallDenominators, stallPositions, plotRR = FALSE, rr = NA, logProp = FALSE, nlg = 11, colorBubbles = baseColorBubbles, titlemap, plotLegend = TRUE, sampleSizes = NA){
  
  # stallNumerators: vector of values at the numerator of the proportions
  # stallDenominator: vector of values at the denominator of the proportions (sample size)
  # stallPositions: vector of geometries
  # plotRR: whether to plot risk map
  # rr: risk output
  # logProp: whether to plot the proportions on log scale (if TRUE, continuous scale)
  # nlg: number of discrete colors, if discrete scale
  # colorBubbles: base color for the bubbles
  # titlemap: plot title
  # plotLegend: whether to plot legend
  # sampleSizes: if plotting reads, need to have information on sample sizes
  
  risk <- rr # Rename it to be able to change it
  
  # Point sizes
  if(all(is.na(sampleSizes))){
    sampleSizes <- stallDenominators
  }
  cexs <- getSampleCex(sampleSizes)
  
  # Compute proportions
  fracs <- stallNumerators / stallDenominators
  
  # Order 
  iorder <- order(fracs, decreasing = FALSE)
  
  # Compute background colors depending on the type of data
  if(!logProp){
    # Linear scale (proportion detections)
    nlg <- 11 # Number of discrete values used
    # (NB: legend hard-coded for 11 values)
    bgs <- getSampleCol(fracs, baseCol = colorBubbles, ncol = nlg)
    
  }else{
    # log scale (proportion reads)
    # Define colors
    bgs <- colFreq(fracs, opacity = opacityBubble)
  }
  
  # Start plotting
  par(mar = marBubblePlot)
  
  # Plot base: risk surface, or just market background
  if(plotRR){
    par(mfrow = c(1, 1))
    # Define ticks for the legend
    axtcks <- c(1, 0.5, 0.1, 0.05, 0.01, 0.005)
    # Initialize plot
    plot(market_boundary_r$geometry, col = colMarketBackground, lwd = lwdMarketBoundary, border = colMarketBackground)
    
    # same p-value plotting approach as for SC2
    zmax <- ceiling(maxl10P * 2) / 2 # Round maximal value
    risk$rrs$rr <- -log(risk$rrs$P, 10) # Trick to be able to plot -log10
    plot(risk$rrs, add = TRUE, 
         zlim = c(0, zmax), col = palHeatmap, 
         tol.type = c("upper"), tol.args = list(levels = c(0.05,0.01), lty = c(2, 1), lwd = 1, drawlabels = TRUE, col = gray(0.5, 0)),  # transparency 0 does not plot
         auto.axes = FALSE, 
         addcontour = FALSE, 
         riblab = "", ribbon = FALSE, 
         ribargs = list(at = -log(axtcks, 10), labels = axtcks, lwd = 0, lwd.ticks = 1), 
         ribsep = 0.2, ribwid = 0.025, 
         asp = 1)
    
    # Add again the market boundary, now just for the contour 
    # (risk plot does not put nice contour)
    #   Remove drain contour if present
    plot(market_boundary_withDrains_r$geometry, col = gray(0, 0), lwd = 2, border = "white", add = TRUE)
    #   Market contour if present
    plot(market_boundary_r$geometry, col = gray(0, 0), lwd = lwdMarketBoundary, border = col0, add = TRUE)
    
    # Plot stalls
    plot(map_stalls_r$geometry, col = colStallBorder, add = TRUE)
    
    # Add wildlife stalls
    plot(wildlifeStalls_r$geometry, col = colWildlifeStall, border = colStallBorder, 
         weight = weightWildlifeStall, opacity = opWildlifeStall, add = TRUE)
    
  }else{
    plotBaseMap()
  }
  
  # Add the bubbles
  plot(stallPositions[iorder], cex = cexs[iorder], pch = 21, bg = bgs[iorder], add = TRUE, col = borderBubble, lwd = 0.3)
  
  if(plotLegend){
    
    # Add Legend
    par(xpd = TRUE)
    # Legends
    #  Relative position of the bubble legends
    #  Sample size legend
    if(max(sampleSizes) > 15){
      svals <- c(1, 5, 10, 15) # Values of the points we are labeling
    }else{
      svals <- c(1, 2, 5, 10) # Values of the points we are labeling
    }
    
    legend("bottomleft", pt.cex = getSampleCex(svals), pch = 21, col = 1, legend = svals, inset = c(-0.03, 0.015), xjust = 0, yjust = 0, 
           cex = 0.9, bty = "n", title = "Number of\nsamples", 
           text.width = txtwd, title.adj = ttladj)
    
    if(!logProp){ # Discrete scale
      #  Positivity legend
      lgt <- rep("", nlg)
      lgt[c(1, 2, 6, 11)] <- c("0", "0-0.1", "0.4-0.5", "0.9-1")
      legend(x = "topright",
             legend = rev(lgt),
             fill = rev(getSampleCol(seq(0, 1, length.out = nlg), baseCol = colorBubbles)),
             border = borderBubble,
             y.intersp = 0.7,
             x.intersp = 0.3,
             cex = cexLegend, text.font = 1, box.lwd = 0, horiz = FALSE, bty = 'n', 
             inset = c(insetxval, 0), xjust = 1, yjust = 0, 
             pt.cex = 1, 
             title = "Proportion \ndetection", 
             text.width = txtwd, title.adj = ttladj)
      
    }else{
      # Vector of values of proportions
      tmp <- seq(minpp, maxpp, length.out = 101)
      # Get round values from it
      tmpp <- pretty(tmp)
      # Limit these values to the ones included in our vector
      tmpp <- tmpp[tmpp >= minpp & tmpp <= maxpp]
      # Identify the indices of the round values
      itmp <- sapply(tmpp, function(x) which(x <= tmp)[1])
      # Prepare the legend text
      lgs <- rep("", length(tmp))
      lgs[itmp] <- as.expression(sapply(tmpp, function(x) bquote(italic(10^.(-x)))))
      # Plot the legend
      legend("topright", fill = colFreq(p = 10^-tmp, opacity = opacityBubble), y.intersp = 0.1, border = NA,
             x.intersp = 0.5,
             legend = lgs, bg = "white", box.lwd = 0, cex = cexLegend, pt.cex = 1, 
             title = "Proportion\nreads\n", 
             inset = c(insetxval, insetyval), xjust = 1, yjust = 0, bty = "n", text.width = txtwd, title.adj = ttladj
      )
      rm(tmp, tmpp)
    }
    if(plotRR){
      # Heatmap legend
      # Vector of values of proportions
      nn <- 101
      tmp <- seq(0, zmax, length.out = nn)
      # Identify the indices of the round values
      itmp <- sapply(axtcks, function(x) which(x >= 10^(-tmp))[1])
      # Prepare the legend text
      lgs <- rep("", length(tmp))
      lgs[itmp] <- axtcks
      # Plot the legend
      legend("topright", fill = palHeatmap(nn), y.intersp = 0.115, border = NA,
             x.intersp = 0.5,
             legend = lgs, bg = "white", box.lwd = 0, cex = cexLegend, pt.cex = 1,
             bty = "n",
             title = "Relative risk, \np-value\n",
             inset = c(insetxval - 0.015, 0.45), xjust = 1, yjust = 0,
             text.width = txtwd, title.adj = ttladj
      )
    }
    par(xpd = FALSE)
    
  }
  title(main = titlemap)
  
}









