# Plot species detection in Stall West 6-29

# Load data if necessary
if(!exists("samples")){
  source("load_data-sequencing.R")
}

if(!exists("thepal")){
  source("generic_plotting.R")
}


# Subset of data in this stall
dat629 <- samples[which(substr(samples$address, 1, 9) == "West|6|29"), ]

# Subselect samples with sequencing data
dat629 <- dat629[!is.na(dat629$Homo.sapiens), ]
nrow(dat629)

# Numbers of samples per date
nsamp629 <- aggregate(dat629$address, by = list(Sampling.date = dat629$Sampling.date), FUN = length)
names(nsamp629)[2] <- "nbSamples"

dat629 <- merge(dat629, nsamp629, all = TRUE)
rm(nsamp629)

# Check address and sample types
table(dat629$Sampling.date, dat629$address, useNA = "ifany")

# Species detected in the stall
abs629 <- apply(dat629[, spMZ], 2, sum, na.rm = TRUE) # Abundances
sp629 <- spMZ[which(abs629>0)]

# Proportions of total reads
props629 <- as.data.frame(dat629[, c(spMZ, "SARS2.reads")] / matrix(rep(dat629$Read_pairs_after_trimming, length(spMZ)), ncol = length(spMZ)))
names(props629) <- c(spMZ, "SARS2.reads")

# For plotting: get range of values
xmin <- 0.5 * 10^-9
xmax <- max(props629, na.rm = TRUE)

# Convert date into numerical version
dat629$numDate <- as.numeric(as.Date(dat629$Sampling.date) - as.Date("2020-01-01")) + 1
# Add horizontal jitter to the points
xx <- jitter(dat629$numDate, factor = 0.5)


# Define function for plotting
plot629 <- function(sp, pdf = TRUE){
  # sp: species name (column name)
  # pdf: whether to save as pdf
  
  # Get values, and change 0 to xmin to plot on log scale
  yy <- props629[, sp]
  yy[which(yy == 0)] <- xmin
  
  if(pdf){
    fname <- paste0("../figs/FigS6_", sp, ".pdf")
    pdf(fname, width = wpdf.overtime, height = hpdf.overtime)
  }
  par(mai = mai.overtime)
  # Plot the points
  plot(xx, yy, log = "y", axes = FALSE, 
       xlab = "", ylab = "Proportion of total reads", 
       pch = pchType[dat629$Sample.type], bg = makeTransparent(colsType[dat629$Sample.type], op.overtime), col = colsType[dat629$Sample.type], 
       ylim = c(xmin, xmax))
  # Title: common species name, first letter capitalized
  spname <- ifelse(sp == "SARS2.reads", "SARS-CoV-2", dicoSp[sp])
  title(main = paste0(toupper(substr(spname, 1, 1)), substr(spname, 2, nchar(spname))), font.main = 1)
  # X axis
  tmp <- dat629[!duplicated(dat629$Sampling.date), ]
  axis(1, at = tmp$numDate, labels = paste0(format(as.Date(tmp$Sampling.date), "%d %b"), " (n=", tmp$nbSamples, ")"), las = 3, cex.axis = cexa.overtime, lwd = 0, lwd.ticks = 0, line = -0.5)
  mtext(side = 1, line = 5.5, text = "Sampling date")
  
  # Y axis
  yvals <- 4:9
  ylabs <- as.expression(sapply(yvals, function(x) bquote(10^-.(x))))
  for(side in c(2, 4)){
    axis(side, las = 1, at = 10^-yvals, labels = ylabs)
    axis(side, at = xmin, labels = "0", las = 1)
  }
  
  # Legend  
  legend("topright", yjust = 0, 
         pch = c(16, 15), col = makeTransparent(colsType[1:2], 100), 
         legend = names(colsType)[1:2], #title = "Type of sample", 
         bty = "n")
  
  #  x1 <- yy[dat629$Sampling.date == "2020-01-12"]
  #  x2 <- yy[dat629$Sampling.date >= "2020-02-20"]
  #  print(paste(sp, t.test(log(x1), log(x2))$p.value))
  if(pdf){
    dev.off()
  }
}

# Plot for Xiao species and SC2
for(sp in c("Nyctereutes.procyonoides", "Paguma.larvata")){
  plot629(sp, pdf = TRUE)
}
