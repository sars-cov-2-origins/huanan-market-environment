source("generic_usefulFuncs.R")

# Load sample metadata ####

samples <- read.csv("../metadata/Liu_etal_2023_with_sequencing.csv")
head(samples)
dim(samples)

# Column to identify positive samples according to Liu et al.
samples$isPositive <- 1*(samples$SARS.CoV.2.qPCR.result == "Positive")
# Note: this is not just qPCR, this includes positivity via sequencing

# Create dictionaries for the different names

#  Sample ID to Lab code
dicoID2Code <- samples$Lab.code
names(dicoID2Code) <- samples$Sample.ID

#  Lab code to Sample ID
dicoCode2ID <- samples$Sample.ID
names(dicoCode2ID) <- samples$Lab.code

# Remove columns that we will not use
ii <- which(names(samples) == "Sequencing_Run")
samples <- samples[, -ii]

# Load data on PCR results ####

# Add information about PCR
# From Liu et al. 2023 Supp Table 2
l22 <- read.csv("../metadata/Liu_PCR_results.csv")

# Write PCR result in binary format
l22$PCRPos <- 1 * (l22$PCR == "+")

# Add the information to the sample metadata dataset
samples <- merge(samples, l22[, c("Sample.ID", "PCRPos", "CT")], by = "Sample.ID", all.x = TRUE)

# NAs are actually 0, because only positive samples were listed
samples[is.na(samples$PCRPos), "PCRPos"] <- 0

# Compare the results
sum(samples$PCRPos)
sum(samples$isPositive) # isPositive includes the 3 samples that were PCR- 
                        # but in which SARS-CoV-2 was detected by sequencing
rm(l22)

# Load sequencing results ####

## |- mtDNA mappings ####

# Read counts
dataMZ <- read.csv("../mitochondrial_mappings/mitochondrial_metazoa_counts_93.tsv", sep = "\t")
head(dataMZ)
dim(dataMZ)
stopifnot(!any(duplicated(names(dataMZ)))) # Make sure that there are no duplicates

# Covered bases
covMZ <- read.csv("../mitochondrial_mappings/mitochondrial_metazoa_coveredbases_93.tsv", sep = "\t")
# Information about contig lengths
mt.len <- read.csv("../mitochondrial_mappings/mitochondrial_lengths.tsv", sep = "\t", header = FALSE)
names(mt.len) <- c("Contig", "Contig_length")

# Load species information
dataSpecies <- read.csv("../mitochondrial_mappings/species_descriptions_with_common_name.csv", sep = ",")
dim(dataSpecies)
head(dataSpecies)
dataSpecies <- merge(dataSpecies, mt.len)

# There are duplicates, remove the ones with least coverage
dataSpecies[duplicated(dataSpecies$Species) | duplicated(dataSpecies$Species, fromLast = TRUE), ]
dataSpecies <- dataSpecies[!is.element(dataSpecies$Contig, c("NC_025503.1", "OM736804.1")), ]
dataMZ <- dataMZ[, !is.element(names(dataMZ), c("NC_025503.1.Liposcelis.entomophila", "OM736804.1.Monopterus.albus"))]
covMZ <- covMZ[, !is.element(names(covMZ), c("NC_025503.1.Liposcelis.entomophila", "OM736804.1.Monopterus.albus"))]

# Extract species names 
# /!\ Position of the first one is hard-coded
firstNamePos <- 6
# Check position; middle one should be the first species name
names(dataMZ)[c(-1, 0, 1) + firstNamePos]
names(covMZ)[c(-1, 0, 1) + firstNamePos]
speciesNames <- tail(names(dataMZ), -(firstNamePos - 1))
stopifnot(!any(duplicated(speciesNames))) # Check that no duplicated name
# Function to remove accession from name
remove.accession <- function(nm){
  # Get position of the second . in the name
  i2 <- (gregexpr("\\.", nm)[[1]])[2]
  # Return part of the column name after that second point
  substr(nm, i2+1, nchar(nm))
}
tmp <- sapply(speciesNames, remove.accession)
stopifnot(!all(duplicated(tmp)))

# Assign these names to the dataset
stopifnot(length(tmp) == length(firstNamePos:length(names(dataMZ))))
stopifnot(length(tmp) == length(firstNamePos:length(names(covMZ))))
names(dataMZ)[firstNamePos:length(names(dataMZ))] <- tmp
names(covMZ)[firstNamePos:length(names(covMZ))] <- tmp
speciesNames <- unname(tmp)
rm(tmp)
stopifnot(!any(duplicated(speciesNames)))

# Check kinds of animals we have
unique(dataSpecies$Group)
unique(dataSpecies$Class)

## |- Process sequencing results ####

# Identify column at which names start
firstColSp <- firstNamePos 
# For visual inspection

# Please check that the first printed value is not a species name, and that the second one is
names(dataMZ)[c(-1, 0, 1) + firstColSp]
# Other check
stopifnot(!grepl("\\.", names(dataMZ)[firstColSp-1]) & grepl("\\.", names(dataMZ)[firstColSp]))

# Extract species names (column names)
spMZ <- names(dataMZ)[firstColSp:ncol(dataMZ)]
stopifnot(all(spMZ == speciesNames))

# Add indicator variable for presence of sequencing results
dataMZ$mtResults <- 1 

# Make sure that all the lab codes are correctly written
all(is.element(dataMZ$Lab.code, samples$Lab.code))
dataMZ$Lab.code[!(is.element(dataMZ$Lab.code, samples$Lab.code))] # Missing codes, if any

# Simplify the dataset before merging: keep only sample ID and sequencing results
dataMZ <- dataMZ[, c("Lab.code", "Run", "Sample_category", "mtResults", spMZ)]

## |- Metadata for the animal species ####

### |- Taxonomic information ####
# Define key species which will be plotted separately
keySpecies <- c("Homo.sapiens", "Nyctereutes.procyonoides", "Erinaceus.amurensis", "Hystrix.brachyura", "Rhizomys.pruinosus", "Marmota.himalayana", "Paguma.larvata", "Mustela.sibirica", "Arctonyx.collaris", "Muntiacus.reevesi")

#  Names of other mammals
spMammals <- sub(" ", "\\.", dataSpecies[dataSpecies$Group == "Mammal", ]$Species)
spOMam <- spMammals[!is.element(spMammals, keySpecies)]
#  Names of other species (not Mammals)
spOth <- spMZ[!is.element(spMZ, spMammals)]


# More refined groups

# This file is created via `species_taxonomy.R`
#  using NCBI taxonomic information. 

taxo <- read.csv("../mitochondrial_mappings/speciesTaxonomy.csv")
taxo <- taxo[!is.na(taxo$species) & !duplicated(taxo$species), ]
row.names(taxo) <- taxo$species

# Define typeSpecies recursively
typeSpecies <- taxo[, "phylum"]
table(typeSpecies)

typeSpecies[which(typeSpecies == "Chordata" & is.element(taxo$species, spOth))] <- "Other Chordates"
typeSpecies[which(taxo[, "class"] == "Aves")] <- "Birds"
typeSpecies[which(is.element(taxo$species, spOMam))] <- "Other mammals"
typeSpecies[which(is.element(taxo$species, keySpecies))] <- "Key mammals"
typeSpecies[which(is.element(taxo$species, "Homo.sapiens"))] <- "Human"
tmp <- c("Human", "Key mammals", "Other mammals", "Birds", "Other Chordates", "Arthropoda")
tmp2 <- unique(typeSpecies)
# Define code to order the different types of species
dicoSpecies <- 1:length(tmp2)
names(dicoSpecies) <- c(tmp, setdiff(tmp2, tmp))
dicoSpecies
rm(tmp, tmp2)

metaSpecies <- data.frame(species = taxo$species, type = typeSpecies, code = dicoSpecies[typeSpecies])
row.names(metaSpecies) <- metaSpecies$species

### |- Other species metadata: UUIDs for plotting (phylopic), common names ####

# Load UUIDs
# file created via `speciesImages.R`
tmp <- read.csv("../metadata/species_uuid.csv")
# Format them as dictionary
uuids <- tmp$uuid
names(uuids) <- tmp$species

# Add the information to the taxo table
taxo <- merge(taxo, tmp, all = TRUE)
rm(tmp)

# Dictionnary for Common names
tmp <- dataSpecies[, c("Species", "Common_name", "Category")]
names(tmp) <- c("displayLatinName", "commonName", "Category")
tmp$latinName <- gsub(" ", "\\.", tmp$displayLatinName)
dicoSp <- tmp$commonName
names(dicoSp) <- tmp$latinName
tmp$species <- tmp$latinName

# Add the information to the taxo table
taxo <- merge(taxo, tmp[, c("species", "commonName", "displayLatinName")], all.x = TRUE)
rm(tmp)
# Load SARS-CoV-2 results and other sequencing information ####

# Load SC2 specific results
sc2 <- read.csv("../sars2_mappings/sars2_reads_post_trimming.tsv", sep = "\t")
head(sc2)

# Check if there are duplicates in Run names
any(duplicated(sc2$Sample))

# Check that all data are there
stopifnot(all(is.element(sc2$Sample, dataMZ$Run)))
stopifnot(all(is.element(dataMZ$Run, sc2$Sample)))

# Rename columns
names(sc2)[names(sc2) == "Sample"] <- "Run"
names(sc2)[names(sc2) == "Read_count"] <- "SARS2.reads"
names(sc2)[names(sc2) == "Covered_bases"] <- "SARS2.coverage"

# Add the information to the MZ results
tmp <- merge(dataMZ, sc2, by = "Run", all = TRUE)
stopifnot(nrow(tmp) == nrow(dataMZ))
dataMZ <- tmp

tmp <- merge(covMZ, sc2, by = "Run", all = TRUE)
stopifnot(nrow(tmp) == nrow(covMZ))
covMZ <- tmp

rm(sc2)

# Missing a column for sample ID, create it using the lab/ID dictionary
dataMZ$Sample.ID <- dicoCode2ID[dataMZ$Lab.code]
covMZ$Sample.ID <- dicoCode2ID[covMZ$Lab.code]

# Load sequencing metadata to link sequencing results and samples ####
sdat <- read.csv("../metadata/Sequencing_run_info.tsv", sep = "\t")
dim(sdat)
head(sdat)
dim(dataMZ)

table(sdat$Sample_category)
table(sdat$Sample_Type)

# Check composition
addmargins(table(sdat$LibrarySource, sdat$Sample_Type))
# -> we also have non-metagenomic results, to be separated from the rest.

# Check match between Run and Sample.ID
tmp1 <- dataMZ[order(dataMZ$Run), c("Run", "Sample.ID")]
tmp2 <- sdat[order(sdat$Run), c("Run", "LibraryNameFixed")]
stopifnot(all(tmp1 == tmp2))
rm(tmp1, tmp2)

# Merge with other sequencing results
i1 <- which(names(sdat) == "Sample_Type")
dataSeq <- merge(dataMZ, sdat[, i1:ncol(sdat)], all = TRUE)
dim(dataSeq)
head(dataSeq)
rm(dataMZ, i1)

# Separate out non-metagenomic results
iSC2 <- which(dataSeq$LibrarySource == "VIRAL RNA" | dataSeq$Sample_Type == "sars2_amplicon")
dataViral <- dataSeq[iSC2, ]
dim(dataViral)
dataViral$Run

# Keep only metagenomic results
iMetagenomic <- which(dataSeq$LibrarySource != "VIRAL RNA" & dataSeq$Sample_Type != "sars2_amplicon")
dataMetagenomic <- dataSeq[iMetagenomic, ]
dim(dataMetagenomic)

rm(dataSeq, iSC2)

# Check the types of sequencing data
table(dataMetagenomic$LibraryLayout, avgLength = dataMetagenomic$avgLength, useNA = "ifany")
table(dataMetagenomic$Platform, useNA = "ifany")

# Check duplicates
any(duplicated(dataMetagenomic$Run))
any(duplicated(dataMetagenomic$LibraryName))
any(duplicated(dataMetagenomic$LibraryNameFixed))
any(duplicated(dataMetagenomic$Sample.ID))

# Duplicates in LibraryNameFixed: are duplicated because Single & Paired sequencing
idupl <- which(duplicated(dataMetagenomic$LibraryNameFixed) | duplicated(dataMetagenomic$LibraryNameFixed, fromLast = TRUE))
#print(paste("Duplicated: ", sort(dataMetagenomic$LibraryNameFixed[idupl])))
tmp <- dataMetagenomic[idupl, ]
sort(names(tmp))
duplSamples <- tmp[order(tmp$LibraryName), c("LibraryName", "LibraryNameFixed", "LibraryLayout", "Run")]
# Duplicated:
duplSamples
rm(idupl, tmp)

# Duplicates in Sample.ID: Single&Paired sequencing
idupl <- which(duplicated(dataMetagenomic$Sample.ID) | duplicated(dataMetagenomic$Sample.ID, fromLast = TRUE))
tmp <- dataMetagenomic[idupl, ]
tmp[order(tmp$Sample.ID), c("Sample.ID", "LibraryName", "LibraryNameFixed", "LibraryLayout")]
rm(idupl, tmp)

# Samples for which we are missing mtDNA results, if any
missingMTDNA <- dataMetagenomic[is.na(dataMetagenomic$Sample.ID), c("Run", "LibraryName", "Lab.code")]
missingMTDNA

stopifnot(!any(duplicated(names(dataMetagenomic))))

# Combine results by sample
# -> This is summing results when we have multiple runs for a given sample
dataMetagenomicBySample <- aggregate(dataMetagenomic[, c("SARS2.reads", spMZ, "Read_pairs_after_trimming", "spots", "mtResults")], by = list(Sample.ID = dataMetagenomic$LibraryNameFixed), FUN = sumIgnoreNA)
head(dataMetagenomicBySample)
dim(dataMetagenomicBySample)

# Add these results to the samples dataset
samples <- merge(samples, dataMetagenomicBySample, all = TRUE)
rm(dataMetagenomicBySample)
dim(samples)

# Are we missing samples? No
samples[is.na(samples$Sampling.date), ]$Sample.ID

# `dataMetagenomic` is by runs, 
# `samples` is by sample

## Process sample data ####

# New column for whether the sample has sequencing results
samples$isSeq <- 1 * !is.na(samples$SARS2.reads)

# New column to add sequencing results in binary format
samples$seqPos <- 1 * (samples$SARS2.reads > 0)

# Combine PCR result and sequencing result
samples$PCRSeqPos <- samples$PCRPos + samples$seqPos
# Note: has NAs when no sequence

# Compare Liu et al.'s isPositive with what we have
isPos <- samples$PCRPos > 0
isPos[which(samples$seqPos > 0)] <- TRUE # Cannot to it at once because of NAs

all(isPos == samples$isPositive)
# Sample for which Liu and our conclusions differ
ii <- which(samples$isPositive != isPos)
samples[ii, c("Sample.ID", "Lab.code", "Sampling.date", "Sample_location", "SARS2.reads")]

# Add our information
samples$isPos <- isPos


# Add information about available sequences
#   Subset of samples for which we have sequencing data and an address
samplesSeq <- samples[which(!is.na(samples$SARS2.reads) & !is.na(samples$address)), ]
dim(samplesSeq)

# Add total number of samples
tmp <- aggregate(samplesSeq$address, by = list(address = samplesSeq$address), FUN = length)
names(tmp)[2] <- "nbSeq"
samplesSeq <- merge(samplesSeq, tmp, all.x = TRUE)
rm(tmp)

# Add back this information into the main dataset
samples <- merge(samples, samplesSeq[, c("Sample.ID", "address", "nbSeq")], all.x = TRUE)
dim(samples)
rm(samplesSeq)



# Define column combining results
samples$positivity <- paste0("PCR", samples$PCRPos, "_seq", samples$seqPos)
table(samples$positivity)



# Load case data ####
cases <- read.csv("../metadata/cases.csv")

head(cases)


