library(sf) # General map package
library(spatstat) # For risk maps
library(sparr) # For risk maps

# Load data ####

## |- SARS-CoV-2 ####

# Load the data if they are not loaded already
if(!exists("samples")){
  source("load_data-sequencing.R")
}
if(!exists("landmarks")){
  source("load_data-geographic.R")
}
if(!exists("viral")){
  source("load_data-otherviruses.R")
}

# Check that the rotated and unrotated data are ordered the same way
stopifnot(all((samples_geo_r$Precise.location == 1) == (samples_geo$Precise.location == 1)))

conditionLocation <- (samples_geo_r$Precise.location == 1)
cat("We could precisely locate the sample: ")
addmargins(table(precise = samples$Precise.location, hsm=samples$FromHSM & samples$Sample_location != "Warehouse", useNA = "ifany"))

# Rotated data
# Samples collected on or before Jan 12
samples112_r <- samples_geo_r[samples_geo_r$Sampling.date <= "2020-01-12" & conditionLocation, ]
# Samples collected on Jan 1
samples1_r <- samples_geo_r[samples_geo_r$Sampling.date == "2020-01-01" & conditionLocation, ]
# Samples collected on Jan 12
samples12_r <- samples_geo_r[samples_geo_r$Sampling.date == "2020-01-12" & conditionLocation, ]
# Drains only, all dates
samplesDrains_r <- samples_geo_r[samples_geo_r$Sample.type == "Water drain" & conditionLocation, ]

# Identify positive and negative samples
# Jan 1-12, PCR+
pos_PCR_Jan112_r <- samples112_r[samples112_r$PCRPos == 1, ]
neg_PCR_Jan112_r <- samples112_r[samples112_r$PCRPos == 0, ]

# Jan 12, Sequencing +
pos_seq_Jan12_r <- samples12_r[samples12_r$seqPos == 1, ]
neg_seq_Jan12_r <- samples12_r[samples12_r$seqPos == 0, ]

# Drains, PCR+
pos_PCR_drains_r <- samplesDrains_r[samplesDrains_r$PCRPos == 1, ]
neg_PCR_drains_r <- samplesDrains_r[samplesDrains_r$PCRPos == 0, ]

## |- Other viruses ####

viruses112 <- viral[which(viral$Sampling.date <= "2020-01-12"), ]
unique(viruses112$Virus_shorter)
dim(viruses112)

# Subselect data for which we have sequences
seq112_r <- samples112_r[!is.na(samples112_r$Homo.sapiens), c("Sample.ID", "Lab.code", "address", "Sampling.date", "geometry")]
# Compute total number of sequences on these dates
tmp <- aggregate(seq112_r$Sample.ID, by = list(address = seq112_r$address), FUN = length)
names(tmp)[2] <- "nbSeq112"
seq112_r <- merge(seq112_r, tmp, all = TRUE)

BRCoV <- viruses112[which(viruses112$Virus_shorter == "Bamboo rat coronavirus"), ]
CvKo <- viruses112[which(viruses112$Virus_shorter == "Civet kobuvirus"), ]
RDAmdo <- viruses112[which(viruses112$Virus_shorter == "Raccoon dog amdovirus"), ]

any(duplicated(BRCoV$Sample.ID))
any(duplicated(CvKo$Sample.ID))
any(duplicated(RDAmdo$Sample.ID))

# Add the samples in which the virus was not found

BRCoV <- merge(BRCoV, seq112_r, all = TRUE)
CvKo <- merge(CvKo, seq112_r, all = TRUE)
RDAmdo <- merge(RDAmdo, seq112_r, all = TRUE)

# Identify positive and negative samples
# BRCoV
pos_BRCoV_r <- BRCoV[which(BRCoV$nDetect == 1), ]
neg_BRCoV_r <- BRCoV[is.na(BRCoV$nDetect), ]

# CvKo
pos_CvKo_r <- CvKo[which(CvKo$nDetect == 1), ]
neg_CvKo_r <- CvKo[is.na(CvKo$nDetect), ]

# RDAmdo
pos_RDAmdo_r <- RDAmdo[which(RDAmdo$nDetect == 1), ]
neg_RDAmdo_r <- RDAmdo[is.na(RDAmdo$nDetect), ]


# Define functions for spatial analysis ####

# Function to create point pattern
# Define function for ppp, including jitter 
#  (set eps = 0 to remove jitter)
#  Jitter needed because issues arise when the points are the absolute exact same location.
myppp <- function(dataset, win, epsilon, convert2matrix = TRUE){
  # epsilon: jitter amplitude, set it to 0 to remove it. 
  
  if(convert2matrix){
    # Convert the dataset's geometry to a matrix
    dat <- geometry2matrix(dataset)
  }else{ # Already as matrix
    dat <- dataset 
  }
  
  # Add or not jitter
  if(epsilon > 0){
    # Compute ppp with jitter
    ppp(x = jitter(dat[,1], amount = epsilon),
        y = jitter(dat[,2], amount = epsilon), 
        window = win, 
        check = TRUE)
  }else{
    # Compute ppp without jitter
    ppp(x = dat[,1],
        y = dat[,2], 
        window = win, 
        check = TRUE)
  }
}

# Define LSCV function to standardize what we do
computeLSCVrisk <- function(pos, control){
  LSCV.risk(pos, control, type = "fixed", seqres = 100, 
            method = "hazelton", auto.optim = TRUE, hlim = c(0, 20))
}

# Define risk function to have the same parameters for all
runRisk <- function(pos, control, hzero){
  sparr::risk(pos, control, 
              h0 = hzero,
              resolution = 1024, tolerate = TRUE, 
              log = FALSE, adapt = FALSE)
}

# Function to compute risk, combining the functions previously defined
computeRisk <- function(pos, neg, window, epsilon = 1E-6){
  # Transform as ppp
  pos_ppp <- myppp(pos, window, epsilon)
  #neg_ppp <- myppp(neg, window, epsilon) # not needed
  control_ppp <- myppp(rbind(pos, neg), window, epsilon)
  
  # Compute h0
  h0 <- computeLSCVrisk(pos_ppp, control_ppp)
  
  # Compute risk
  rsk <- runRisk(pos_ppp, control_ppp, hzero = h0)
  
  out <- list(rrs = rsk, h0 = h0)
}


# Computations! (Takes time) ####

## |- SARS-CoV-2 ####

cat("Spatial risk analyses for SARS-CoV-2 running...\n")
# Jan 1-12, PCR+
cat("Rotated, 1-12 Jan\n")
risk_PCR_Jan112_r <- computeRisk(pos_PCR_Jan112_r, neg_PCR_Jan112_r, window_r)

# Jan 12, Seq+
cat("Rotated, 12 Jan\n")
risk_seq_Jan12_r <- computeRisk(pos_seq_Jan12_r, neg_seq_Jan12_r, window_r)

# Drains, PCR+ (use larger window)
cat("Rotated, drains\n")
risk_PCR_drains_r <- computeRisk(pos_PCR_drains_r, neg_PCR_drains_r, window_withDrains_r)

cat("... done!\n")


## |- Other viruses ####
cat("Spatial risk analyses for other viruses running...\n")

cat("Bamboo rat CoV\n")
risk_BRCoV_r <- computeRisk(pos_BRCoV_r, neg_BRCoV_r, window_r)

cat("Civet kobuvirus\n")
risk_CvKo_r <- computeRisk(pos_CvKo_r, neg_CvKo_r, window_r)

cat("Raccoon dog amdovirus\n")
risk_RDAmdo_r <- computeRisk(pos_RDAmdo_r, neg_RDAmdo_r, window_r)

cat("... done!\n")


# Get the maximum -log10(P) in our plots to homogenize the scale
# (which is why we have to run everything even though this takes time!)
maxl10P <- max(c(max(-log(risk_PCR_Jan112_r$rrs$P, 10)), max(-log(risk_seq_Jan12_r$rrs$P, 10)), max(-log(risk_PCR_drains_r$rrs$P, 10)), max(-log(risk_BRCoV_r$rrs$P, 10)), max(-log(risk_CvKo_r$rrs$P, 10)), max(-log(risk_RDAmdo_r$rrs$P, 10))))

