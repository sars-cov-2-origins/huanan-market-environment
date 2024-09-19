#
# Format tables showing BEAST results
# with recCA constraint, and unconstrained
#

# Initializations ####

# Change directory if needed
if(!file.exists("BEAST/global.Lv2024.masked.GLM.HKY.asr.states.combined.results.tsv")){
  setwd("../sars2_phylogenetics/") 
}
system("ls BEAST")

# Load data
tb_unc <- read.csv("BEAST/global.Lv2024.masked.GLM.HKY.asr.states.combined.results.tsv", sep = "\t")
head(tb_unc)
tb_rec <- read.csv("BEAST/global.Lv2024.masked.GLM.HKY.recCA_GARD.states.combined.results.tsv", sep = "\t")
head(tb_rec)

tb_unc$type <- "Unconstrained"
tb_rec$type <- "recCA"
tb_combined <- rbind(tb_unc, tb_rec)

# Check sample sizes
tmp <- aggregate(tb_combined$count, by = list(tb_combined$type), FUN = sum); tmp
nS <- unique(tmp$x)
stopifnot(length(nS) == 1) # Stop if we have multiple lengths

# 1) Format table, focus on roots shown in the main text ####

# Define the roots that we want to select
roots1B <- c("None", "C8782T_T28144C", "T28144C", "C8782T", "C8782T_C18060T_T28144C", "C8782T_T28144C_C29095T")
names(roots1B) <- c("B", "A", "C/C", "T/T", "A+18060T", "A+29095T")
tmp <- tb_combined[is.element(tb_combined$root, roots1B), ]
dim(tmp)
# Dictionary to switch between names
invdico <- names(roots1B)
names(invdico) <- roots1B

# Format output table
outtb <- as.data.frame(matrix(NA, ncol = length(roots1B)+1), nrow = 2)
names(outtb) <- c("Type", roots1B)
outtb
# Fill in the output table
i <- 1
for(tp in unique(tmp$type)){
  outtb[i, "Type"] <- tp
  for(j in roots1B){
    z <- tb_combined[tb_combined$type == tp & tb_combined$root == j, ]
    if(nrow(z) == 0){
      outtb[i, j] <- 0
    }else{
      outtb[i, j] <- z$count
    }
  }
  i <- i+1
}
names(outtb)[2:ncol(outtb)] <- invdico[names(outtb)[2:ncol(outtb)]]
outtb

# Output table with proportions instead of counts
outtb_prop <- outtb
outtb_prop[2:ncol(outtb_prop)] <- outtb_prop[2:ncol(outtb_prop)] / nS
outtb_prop

# Output table with percentages (rounded)
outtb_percent <- outtb_prop
outtb_percent[2:ncol(outtb_prop)] <- round(outtb_percent[2:ncol(outtb_prop)] * 100, 1)
outtb_percent

# Export the table
write.csv(outtb_percent, row.names = FALSE, file = "roots_percents.csv")

# 2) Full table ####
allroots <- unique(tb_combined$root)

# Reformat
tmp1 <- tb_rec[, 1:3]
names(tmp1)[names(tmp1) == "count"] <- "recCA"
tmp2 <- tb_unc[, 1:3]
names(tmp2)[names(tmp2) == "count"] <- "unconstrained"
tb_c2 <- merge(tmp1, tmp2, all = TRUE)
# Remove NAs, replace them by 0s
tb_c2$recCA[is.na(tb_c2$recCA)] <- 0
tb_c2$unconstrained[is.na(tb_c2$unconstrained)] <- 0
# Reorder by decreasing recCA counts
tb_c2 <- tb_c2[order(tb_c2$recCA, decreasing = TRUE), ]
# Reformat columns to show percentages as well
tb_c2$recCA <- paste0(tb_c2$recCA, " (", round(100* tb_c2$recCA / nS, 1), "%)")
tb_c2$unconstrained <- paste0(tb_c2$unconstrained, " (", round(100* tb_c2$unconstrained / nS, 1), "%)")
# Reformat root, remove "_" and replace them by ", "
tb_c2$root <- gsub("_", ", ", tb_c2$root)
tb_c2

write.csv(tb_c2, row.names = FALSE, file = "roots_percents_all.csv")

