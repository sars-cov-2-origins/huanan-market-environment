# Plot animal detection maps
source("load_data-sequencing.R")

source("generic_plotting.R")
# Load data if necessary
if(!exists("landmarks")){
  source("load_data-geographic.R")
}

if(!exists("plotBubbleMap")){
  source("map_plottingfunctions.R")
}

# Fig 3B - Human ####

tmp <- samples112_r
# Subselect the ones for which we have sequencing data
tmp <- tmp[!is.na(tmp$Homo.sapiens), ]

# Here we are going to plot the proportions of reads
tmp$tot <- tmp$Read_pairs_after_trimming
tmp$detect <- (tmp$Homo.sapiens)
tmp$present <- 1

# Aggregate by stall
agg <- aggregate(tmp[, c("detect", "tot", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
ssize <- as.numeric(unlist(st_drop_geometry(agg[, 4])))
rm(agg)

# Plot
pdf(paste0("../figs/Fig3B.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, sampleSizes = ssize,
              plotRR = FALSE,  
              logProp = TRUE, titlemap = "Human", plotLegend = TRUE, use_large_points=TRUE)

dev.off()
rm(num, denom, geom, ssize)


# Fig 3C - Raccoon dog ####

tmp$detect <- (tmp$Nyctereutes.procyonoides)

# Aggregate by stall
agg <- aggregate(tmp[, c("detect", "tot", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
ssize <- as.numeric(unlist(st_drop_geometry(agg[, 4])))
rm(agg)

# Plot
pdf(paste0("../figs/Fig3C.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, sampleSizes = ssize,
              plotRR = FALSE,  
              logProp = TRUE, titlemap = "Raccoon dog", plotLegend = TRUE, use_large_points=TRUE)

dev.off()
rm(num, denom, geom, ssize)


# Fig 3D - Bamboo rat ####

tmp$detect <- (tmp$Rhizomys.pruinosus)

# Aggregate by stall
agg <- aggregate(tmp[, c("detect", "tot", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
ssize <- as.numeric(unlist(st_drop_geometry(agg[, 4])))
rm(agg)

# Plot
pdf(paste0("../figs/Fig3D.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, sampleSizes = ssize,
              plotRR = FALSE,  
              logProp = TRUE, titlemap = "Hoary bamboo rat", plotLegend = TRUE, use_large_points=TRUE)

dev.off()
rm(num, denom, geom, ssize)

# Fig 3E - Masket palm civet ####

tmp$detect <- (tmp$Paguma.larvata)

# Aggregate by stall
agg <- aggregate(tmp[, c("detect", "tot", "present")], by = list(address = tmp$address), FUN = sumIgnoreNA)

num <- as.numeric(unlist(st_drop_geometry(agg[, 2])))
denom <- as.numeric(unlist(st_drop_geometry(agg[, 3])))
geom <- agg$geometry
ssize <- as.numeric(unlist(st_drop_geometry(agg[, 4])))
rm(agg)

# Plot
pdf(paste0("../figs/Fig3E.pdf"), width = wpdf.maps + sum(mai.maps[c(2, 4)]), height = hpdf.maps + sum(mai.maps[c(1, 3)]), pointsize = pointsize.maps)

plotBubbleMap(stallNumerators = num, stallDenominators = denom, stallPositions = geom, sampleSizes = ssize,
              plotRR = FALSE,  
              logProp = TRUE, titlemap = "Masqued palm civet", plotLegend = TRUE, use_large_points=TRUE)

dev.off()
rm(num, denom, geom, ssize)