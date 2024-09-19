# Plot SARS-CoV-2 positivity maps

# Load data if necessary
if(!exists("samples")){
  source("load_data-sequencing.R")
}

if(!exists("thepal")){
  source("generic_plotting.R")
}

if(!exists("landmarks")){
  source("load_data-geographic.R")
}

if(!exists("plotBubbleMap")){
  source("map_plottingfunctions.R")
}

# Fig 2A - PCR positive samples from Jan 01 and Jan 12 ####

tmp <- samples112_r
tmp$present <- 1
# Aggregate by stall
agg <- aggregate(tmp[, c("PCRPos", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
rm(agg, tmp)

# Plot
pdf(paste0("../figs/Fig2A.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]) + dLeg, height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, 
              plotRR = TRUE, rr = risk_PCR_Jan112_r, 
              logProp = FALSE, titlemap = "01-12 Jan, PCR positivity")

dev.off()
rm(num, denom, geom)


# Fig 2B - NGS positive samples from Jan 12 ####

tmp <- samples12_r
# Subselect the ones for which we have sequencing data
tmp <- tmp[!is.na(tmp$Homo.sapiens), ]
tmp$present <- 1
# Aggregate by stall
agg <- aggregate(tmp[, c("seqPos", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
rm(agg, tmp)

# Plot 
# Note the different margins, because we are not plotting the legend again
pdf(paste0("../figs/Fig2B.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, 
              plotRR = TRUE, rr = risk_seq_Jan12_r, 
              logProp = FALSE, titlemap = "12 Jan, sequencing", plotLegend = FALSE)

dev.off()
rm(num, denom, geom)

# Fig 2B - NGS positive samples from Jan 12 ####

tmp <- samplesDrains_r
tmp$present <- 1
# Aggregate by stall
agg <- aggregate(tmp[, c("PCRPos", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
rm(agg, tmp)

# Plot 
# Note the different margins, because we are not plotting the legend again
pdf(paste0("../figs/Fig2C.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, 
              plotRR = TRUE, rr = risk_PCR_drains_r, 
              logProp = FALSE, titlemap = "27 Jan - 15 Feb, Drains PCR", plotLegend = FALSE)

dev.off()
rm(num, denom, geom)


