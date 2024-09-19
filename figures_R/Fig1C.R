#
# Plot tMRCA distributions
#

# Note: check samples sizes

# Initializations

# Load data
tMRCA <- read.csv("../sars2_phylogenetics/combined_tMRCA.csv")
head(tMRCA)


# NB: each analysis has a multiple of 9001 samples
#     (1–4 chains were needed per analysis)
cat("Chains needed for each analysis: \n")
print(table(tMRCA$analysis)/9001)

# Colors
# (2 and 3 colors are obtained from Paletton)
colMarket <- "#D40000"
colMarket2 <- "#AA0060"
colWuhan <- "#D4AA00"
colWuhan2 <- "#D4CD00"
colAll <- "#5F8DD3"
colAll2 <- colAll  #"#7766D7" 
colAll3 <- "#53D0B8"

# Markers for mean and median
pchMean <- 5
pchMedian <- 3

# Function to make color transparent
# from <https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color>
makeTransparent = function(..., alpha = 0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  if(alpha == 0){
    return(gray(0, 0))
  }else{
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    .makeTransparent = function(col, alpha) {
      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    return(newColor)
  }
}

# Check the market data
sub <- tMRCA[tMRCA$analysis == "market_unconstrained", ]
mean(sub$tMRCA < 2019.00)
# This is a very minor proportion of the data
# Change to min date to simplify plotting
mindate <- 2019.00
tMRCA[tMRCA$tMRCA < mindate, "tMRCA"] <- mindate
# Check that this does not alter the results
sub2 <- tMRCA[tMRCA$analysis == "market_unconstrained", ]
cat("Change in the median tMRCA after masking the few very very early dates")
print(c(median(sub$tMRCA), median(sub2$tMRCA)))


# Convert times into dates
num2date <- function(numdate){
  # numdate can be a vector
  yr <- floor(numdate)
  out <- as.Date(paste0(as.character(yr), "-01-01")) + round(365 * (numdate - yr))
  out
}

tMRCA$tMRCA.date <- num2date(tMRCA$tMRCA)

# Inverse function: convert dates into numerical
date2num <- function(dates){
  # dates can be a vector
  yr <- as.numeric(substr(dates, 1, 4))
  deltat <- dates - as.Date(paste0(as.character(yr), "-01-01"))
  out <- yr + deltat / 365
}

tMRCA$tMRCA2 <- date2num(tMRCA$tMRCA.date)

# Define dictionary to typeset names
datasets <- sort(unique(tMRCA$analysis))
datasets

if(length(datasets) == 7){
  dico <- c("All, recCA", "All, unconstrained", 
            "All, original", 
            "Market, unconstrained", "Market, recCA", 
            "Wuhan, unconstrained", "Wuhan, recCA")
}
if(length(datasets) == 4){
  dico <- c("All, recCA", "All, unconstrained", 
            "Market, unconstrained", "Wuhan, unconstrained")
}

if(length(datasets) != length(dico)){
  stop("Check `dico`; its length needs to match the number of analyses in the file.")
}
names(dico) <- datasets
# Check consistency
dico

# Vector of all days for histogram
days <- seq(as.Date(min(tMRCA$tMRCA.date)), as.Date(max(tMRCA$tMRCA.date)), by = "day")
# Save with offset to ensure proper centering
days.hist <- c(days - 0.5, days[length(days)] + 0.5)
days.hist.num <- date2num(days.hist) # Numerical version

# Compute distributions per day
h <- list()
for(aa in datasets){
  sub <- tMRCA[tMRCA$analysis == aa, ]
  print(aa)
  print(range(sub$tMRCA.date))
}


for(aa in datasets){
  sub <- tMRCA[tMRCA$analysis == aa, ]
  
  h[[aa]] <- hist(sub$tMRCA.date, breaks = days.hist[days.hist > mindate], main = aa)
#  lo <- smooth.spline(h[[aa]]$density ~ h[[aa]]$mids)
#  # <https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r>
#  lines(predict(lo, h[[aa]]$mids)$x, predict(lo, h[[aa]]$mids)$y, col = "red", lwd = 5)
}

# Define colors
subdatasets <- c("all_Lv_recCA", "all_Lv_unconstrained", "market_unconstrained", "wuhan_unconstrained",
                 "wuhan", "market", 
                 "all_old", 
                 "wuhan_recCA", "market_recCA")
colsDS <- c(colAll2, colAll, colWuhan, colMarket, 
            colWuhan, colMarket, 
            colAll3, 
            colWuhan2, colMarket2)
names(colsDS) <- subdatasets


# Function to plot a violin plot
# is used in the plotVioMRCA function defined below
# There is no smoothing here, basically just plots a polygon
plotvio <- function(x, y, pos, width = 0.4, ...){
  # x: x values of the distribution
  # y: y values of the distribution, densities
  # pos: position of the violin plot
  # width: max width of the violin
  # ...: parameters of the polygon
  
  # Just to make sure we are not plotting value with 0 densities
  i0 <- which(y > 0)
  # Note: The violins are horizontal
  
  # Width of the violins
  scale <- width / max(y)
  # Prepare polygon coordinates
  xx <- c(x[i0], rev(x[i0]), x[i0][1])
  yy <- c(pos + scale * y[i0], pos - scale * rev(y[i0]), pos + scale * y[i0][1])

  polygon(xx, yy, ...)
}

# Function to plot tMRCA violin plots
plotVioMRCA <- function(data, datasets = NULL, hval = 0.01, viowidth = 0.4, positions = NULL, cols = NULL, add = FALSE, xmin = NULL, alpha95 = 0.6, alphaOut = 0, save = TRUE, filename = "plot", ext = "png"){
  # data: data containing tMRCAs
  # datasets: analyses to plot
  # hval: smoothing factor
  # viowidth: width of the violin plots
  # positions: positions of the plots (automatic if NULL)
  # cols: colors of the violin plots
  # add: whether to add to an existing plot
  # xmin: minimum x value
  # alpha95: alpha value in the 95% interval
  # alphaOut: alpha value outside
  # ext: file type if figure is saved
  
  # Checks and assignments
  if(is.null(datasets)) datasets <- sort(unique(data$analysis))
  
  if(is.null(xmin)){
    # Define xmin as the min of the lowest 95% quantile
    xmin <- min(unlist(lapply(datasets, function(x){
      ss <- data[data$analysis == x, ]
      quantile(ss$tMRCA, 0.025)
    })))
  }
  xmax <- max(data$tMRCA)
  
  if(length(viowidth) == 1) viowidth <- rep(viowidth, length(datasets))
  
  if(all(is.null(positions))) positions <- seq_along(datasets)
  
  # Colors
  if(all(is.null(cols))){
    cols <- seq_along(datasets)
    names(cols) <- datasets
  }
  
  # Initialize plot
  deltax <- 0.05 * (xmax - xmin) # Slack around the min max positions
  if(!add){
    ww <- 6
    hh <- 0.5 + 0.5 * length(datasets)
    if(save){
      if(ext == "png"){
        png(paste0(filename, ".png"), width = ww, height = hh, units = "in", res = 250)
      }
      if(ext == "pdf"){
        pdf(paste0(filename, ".pdf"), width = ww, height = hh)
      }
      par(mai = c(0.35, 0.25, 0.15, 1.5), mgp = c(1.5, 0.5, 0))
    }else{
      par(mar = c(3, 1, 1, 6))
    }
    ymin <- min(positions) - 1.2*viowidth[1]
    ymax <- max(positions) + 1.2*viowidth[1]
    plot(0, xlim = c(xmin - deltax, xmax + deltax), ylim = c(ymin, ymax), 
         axes = FALSE, xlab = "", ylab = "")
    # Locations of the x ticks
    tmp <- seq(as.Date("2019-01-01"), as.Date("2020-01-01"), by = "month")
    tckstimes <- sort(c(tmp, tmp+14))
    # Time axis
    axis(1, at = date2num(tckstimes), labels = format(tckstimes, "%b %d"))
  }
  
  # For each dataset / violin plot
  for(i in seq_along(datasets)){
    # Subset of the data for this analysis
    sub <- data[data$analysis == datasets[i], ]
    # Compute smoothed distribution, on the points at which there are data
    yy <- sm::sm.density(sub$tMRCA, xlim = range(sub$tMRCA), h = hval, 
                         eval.points = days.hist.num[days.hist.num >= min(sub$tMRCA) & days.hist.num <= max(sub$tMRCA)], 
                         display = "none")
    
    # Identify the 95% interval in the original data
    qq <- quantile(sub$tMRCA, c(0.025, 0.975))
    # Select the days in the interval
    ii <- which(yy$eval.points >= qq[1] & yy$eval.points <= qq[2])
    
    thepos <- positions[i] # Position of the violin plot
    thecol <- unname(c(cols[datasets[i]])) # Main color for the violin

    # Plot violin with all data
    plotvio(yy$eval.points, yy$estimate, pos = thepos, col = makeTransparent(thecol, alpha = alphaOut), border = thecol, width = viowidth[i], lwd = 0.25)
    # Plot violin with data in the 95% interval
    plotvio(yy$eval.points[ii], yy$estimate[ii], pos = thepos, col = makeTransparent(thecol, alpha = alpha95), border = thecol, width = viowidth[i])
    # Add median and mean
    med <- median(sub$tMRCA)
    mn <- mean(sub$tMRCA)
    lwdpts <- 2
    points(med, y = thepos, pch = pchMedian, col = thecol, lwd = lwdpts)
    points(mn, y = thepos, pch = pchMean, col = thecol, lwd = lwdpts)
    
    # Add text label
    par(xpd = TRUE)
    text(x = xmax, y = thepos, labels = paste0(" ", dico[datasets[i]]), adj = c(0, 0.5))
    par(xpd = FALSE)
  }
  
  # Add legend for the markers
  colLeg <- ifelse(length(unique(colsDS[datasets])) == 1, unique(colsDS[datasets]), "black") 
  legend("topleft", pch = c(pchMean, pchMedian), legend = c("Mean", "Median"), lwd = 2, col = colLeg, lty = 0, 
         horiz = TRUE, bty = "n", 
         cex = 0.8, yjust = 0, y.intersp = 0.5, x.intersp = 0)
  if(save){
    dev.off()
  }
}

# Compute summary stats and put them in a table
# (there is probably a cleaner way to do it)

# Median
med <- aggregate(tMRCA$tMRCA, by = list(analysis = tMRCA$analysis), FUN = median)
med$x <- num2date(med$x)
names(med)[names(med) == "x"] <- "median"
# Mean 
mn <- aggregate(tMRCA$tMRCA.date, by = list(analysis = tMRCA$analysis), FUN = mean)
mn$x <- num2date(mn$x)
names(mn)[names(mn) == "x"] <- "mean"
tmp <- merge(med, mn)
# Lower value of the 95% interval
q1 <- aggregate(tMRCA$tMRCA, by = list(analysis = tMRCA$analysis), FUN = function(z) quantile(z, probs = c(0.025)))
q1$x <- num2date(q1$x)
names(q1)[names(q1) == "x"] <- "2.5%"
tmp <- merge(tmp, q1)
# Higher value of the 95% interval
q2 <- aggregate(tMRCA$tMRCA, by = list(analysis = tMRCA$analysis), FUN = function(z) quantile(z, probs = c(0.975)))
q2$x <- num2date(q2$x)
names(q2)[names(q2) == "x"] <- "97.5%"
summarystats <- merge(tmp, q2)

# Print values
cat("Summary stats:\n")
print(summarystats)

# Plot! ##### Plot! ####median()

# Check
plotVioMRCA(data = tMRCA, datasets = subdatasets[rev(1:4)], cols = colsDS, xmin = date2num(as.Date("2019-09-01")), save = FALSE)

# Plot and save figures

plotVioMRCA(data = tMRCA, datasets = datasets[rev(c(1, 2))], cols = colsDS, xmin = date2num(as.Date("2019-09-01")), filename = "../figs/Fig1C_alldata", ext = "pdf")

plotVioMRCA(data = tMRCA, datasets = datasets[rev(c(1, 2, 4, 3))], cols = colsDS, xmin = date2num(as.Date("2019-08-01")), filename = "../figs/Fig1C_Supp_alldata-allanalyses")



