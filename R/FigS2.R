# Plot Ct values over time 

# Load data if necessary
if(!exists("samples")){
  source("load_data-sequencing.R")
}

if(!exists("thepal")){
  source("generic_plotting.R")
}


# Subset of the data with CT values
cts <- samples[!is.na(samples$CT), c("address", "Sample.ID", "Lab.code", "Sampling.date", "CT", "PCRPos", "Sample.type")]
cts
# New column to select CT values that are OK
cts$CT_OK <- FALSE
cts[which(is.element(nchar(cts$CT), c(2, 4, 5))), "CT_OK"] <- TRUE


# Number of samples on Jan12
nJ12 <- nrow(samples[samples$Sampling.date == "2020-01-12" & samples$PCRPos == 1 & samples$CT == "+", ])

# Add column for date and convert it to numerical value
cts$date <- as.Date(cts$Sampling.date)
cts$numDate <- cts$date - as.Date("2020-01-01")

# Extract CT values for stair1-2
CTstair12 <- cts[which(cts$address == "West|5|stair1-2"), ]
CTstair12 <- rbind(CTstair12, CTstair12)
CTstair12$CT_OK <- TRUE
CTstair12[1, "CT"] <- strsplit(CTstair12[1, "CT"], "/")[[1]][1]
CTstair12[2, "CT"] <- substr(strsplit(CTstair12[2, "CT"], "/")[[1]][3], 1, 5)
CTstair12

plotCT <- function(pdf = TRUE){
  
  par(las = 1)
  if(pdf){
    fname <- "../figs/FigS2.pdf"
    pdf(fname, width = wpdf.overtime, heigh = hpdf.overtime)
  }
  par(mgp = mgp.overtime, tck = tck.overtime, mai = mai.overtime)
  plot(jitter(as.numeric(cts$numDate[cts$CT_OK]), factor = 1), - as.numeric(cts$CT[cts$CT_OK]), 
       pch = pchType[cts$Sample.type[cts$CT_OK]], col = colsType[cts$Sample.type[cts$CT_OK]], bg = makeTransparent(colsType[cts$Sample.type[cts$CT_OK]], op.overtime), 
       axes = FALSE, 
       xlab = "", ylab = "", 
       xlim = range(cts$numDate))
  # Add stair1-2
  points(as.numeric(CTstair12$numDate), -as.numeric(CTstair12$CT), 
         pch = pchType[cts$Sample.type], col = colsType[CTstair12$Sample.type], bg = makeTransparent(colsType[CTstair12$Sample.type], op.overtime), 
         type = "o"
  )
  
  # Axes
  #  Y
  yvals <- -38:-23
  for(side in c(2, 4)) axis(side, at = yvals, labels = -yvals, las = 1)
  mtext(side = 2, line = 3, text = "CT", las = 1)
  #  X
  days <- table(cts$date[cts$CT_OK | cts$address == "West|5|stair1-2"])
  theline <- 0
  axis(1, at = c(as.Date(names(days)) - as.Date("2020-01-01")), labels = paste0(format(as.Date(names(days)), "%d %b"), " (n=", days, ")"), las = 3, line = theline, lwd = 0, lwd.ticks = 1, cex = cexa.overtime)
  axis(1, at = c(as.Date("2020-01-12") - as.Date("2020-01-01")), labels = paste0(format(as.Date("2020-01-12"), "%d %b"), "* (n=", nJ12, ")"), las = 3, line = theline, lwd = 0, lwd.ticks = 1, cex = cexa.overtime)
  mtext(side = 1, line = 5.5, text = "Sampling date")
  legend("topright",
         pch = c(16, 15, 18), col = makeTransparent(colsType, 100), 
         legend = names(colsType), title = "Type of sample", bty = "n")
  
  if(pdf){
    dev.off()
  }
  
}

plotCT(pdf = TRUE)
