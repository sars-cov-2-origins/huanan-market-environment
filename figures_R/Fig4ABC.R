# Plot other viruses detection maps

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

# Fig 4A - Raccoon dog amdovirus in samples collected on Jan 01 and 12 ####

# Base number of samples
tmp <- samples112_r
# Subset of data for which we have sequencing data
tmp <- tmp[!is.na(tmp$Homo.sapiens), ]
tmp$present <- 1
tmp <- tmp[, c("address", "Sample.ID", "present")]
dim(tmp)

tmpvir <- viral[which(viral$Virus_shorter == "Raccoon dog amdovirus" & viral$Sampling.date <= "2020-01-12"), c("address", "Sample.ID", "read_count")]
dim(tmpvir)

tmpvir <- merge(tmp, tmpvir, all = TRUE)
tmpvir$detected <- !is.na(tmpvir$read_count)

dim(tmpvir)

# Aggregate by stall
agg <- aggregate(tmpvir[, c("detected", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
rm(agg, tmpvir)

# Plot
pdf(paste0("../figs/Fig4A.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, 
              plotRR = TRUE, rr = risk_RDAmdo_r, 
              logProp = FALSE, titlemap = "Raccoon dog amdovirus", plotLegend = TRUE,  use_large_points=TRUE)

dev.off()

rm(num, denom, geom)


# Fig 4B - Civet kobuvirus in samples collected on Jan 01 and 12 ####

tmpvir <- viral[which(viral$Virus_shorter == "Civet kobuvirus" & viral$Sampling.date <= "2020-01-12"), c("address", "Sample.ID", "read_count")]
dim(tmpvir)

tmpvir <- merge(tmp, tmpvir, all = TRUE)
tmpvir$detected <- !is.na(tmpvir$read_count)

dim(tmpvir)

# Aggregate by stall
agg <- aggregate(tmpvir[, c("detected", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
rm(agg, tmpvir)

# Plot
pdf(paste0("../figs/Fig4B.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, 
              plotRR = TRUE, rr = risk_CvKo_r, 
              logProp = FALSE, titlemap = "Civet kobuvirus", plotLegend = TRUE, use_large_points=TRUE)

dev.off()

rm(num, denom, geom)

# Fig 4C - Bamboo rat coronavirus in samples collected on Jan 01 and 12 ####

tmpvir <- viral[which(viral$Virus_shorter == "Bamboo rat coronavirus" & viral$Sampling.date <= "2020-01-12"), c("address", "Sample.ID", "read_count")]
dim(tmpvir)

tmpvir <- merge(tmp, tmpvir, all = TRUE)
tmpvir$detected <- !is.na(tmpvir$read_count)

dim(tmpvir)

# Aggregate by stall
agg <- aggregate(tmpvir[, c("detected", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
rm(agg, tmpvir)

# Plot
pdf(paste0("../figs/Fig4C.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, 
              plotRR = TRUE, rr = risk_BRCoV_r, 
              logProp = FALSE, titlemap = "Bamboo rat coronavirus", plotLegend = TRUE, use_large_points=TRUE)

dev.off()

rm(num, denom, geom)

