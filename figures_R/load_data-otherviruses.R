
if(!exists("samples")){
  stop("Need to run load_data-sequencing.R first!")
}


## Load results on other viruses

# Load the data
viral_metadata <- read.csv("../viruses/viral_metadata.tsv", sep = "\t")
head(viral_metadata)

viral_rc <- read.csv("../viruses/filtered_viral_counts_97_95_20_200.tsv", sep = "\t")
head(viral_rc)
dim(viral_rc)

viral_names <- read.csv("../viruses/viral_metadata_names.csv")
head(viral_names)

# Merge metadata files
viral_metadata <- merge(viral_metadata, viral_names, all.x = TRUE)

# Create Reference column to link the two datasets
refs <- unname(vapply(viral_metadata$Virus_full_name, function(x) strsplit(x, " ")[[1]][1], FUN.VALUE = "x"))
viral_metadata$Reference <- refs
viral_metadata$contig = viral_metadata$Reference
# Add metadata to viral results
viral_rc <- merge(viral_rc, viral_metadata, all.x = TRUE)
rm(viral_metadata)

# Rename columns to match
names(viral_rc)[which(names(viral_rc) == "file")] <- "Run"

if(!exists("landmarks")){
  # Are the addresses OK? No, so fix them
  all(is.element(unique(viral_rc$Stall), landmarks$title))
  setdiff(viral_rc$Stall, landmarks$title)
}

all(is.element(unique(viral_rc$Lab_code), samples$Lab.code))
setdiff(viral_rc$Lab_code, samples$Lab.code)

# Retrieve corresponding lab codes (excluding amplicon seq)
viral_rc <- merge(viral_rc, dataMetagenomic[, c('Run', 'Lab.code')])
dim(viral_rc)

viral_rc <- merge(viral_rc, samples[, c("Lab.code", "Sample.ID", "address", "Sampling.date")], all.x = TRUE)
dim(viral_rc)

# Are the addresses OK now?
all(is.element(unique(viral_rc$address), landmarks$title))
setdiff(viral_rc$address, landmarks$title)
unique(viral_rc[is.na(viral_rc$address), "Lab.code"])
# Sewers outside of the market -> we are fine!

# New column for detection
viral_rc$nDetect <- 1
viral_rc$keep <- 1 * (viral_rc$Interesting == "Yes")
# Keep mottle viruses
viral_rc[which(grepl("mottle", viral_rc$Virus_name)), "keep"] <- 1
# Keep SC2
viral_rc[which(viral_rc$Host_general == "SARS2"), "keep"] <- 1

# Keep Human viruses
viral_rc[which(viral_rc$Host_general == "Human"), "keep"] <- 1

# Simplify names
viral_rc$Virus_simpler_name <- viral_rc$Virus_shorter

# # Columns to keep
cls <- c("address", "keep", "Virus_simpler_name", "Virus_shorter", "Host_general", "read_count", "Reference", "Sample.ID", "Lab.code", "nDetect", "viral_length", "covered_bases", "Sampling.date")
viral <- viral_rc[, cls]
# Subset the data to keep only the ones we want
viral <- viral[viral$keep == 1, ]
head(viral)
unique(viral$Host_general)
sort(unique(viral$Virus_simpler_name))
dim(viral)

# Order by host then by name
viral <- viral[order(viral$Host_general, viral$Virus_simpler_name), ]

# Add geographic information
viral_geo_r <- merge(landmarks_r[, c("address", "geometry")], viral)

