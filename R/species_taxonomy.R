# Obtain taxonomic information for our species

getInfo <- function(sp){
  taxinfo <- myTAI::taxonomy(organism = sp, db = "ncbi", output = "classification")
  out <- taxinfo[is.element(taxinfo$rank, c("kingdom", "phylum", "class")), "name"]
  out <- c(out, rep(NA, 3))[1:3] # Fill in with NA if necessary
  Sys.sleep(1)
  out
}

res <- matrix(NA, ncol = 4, nrow = length(spMZ))

cat("out of", nrow(res), ":")
for(i in 1:nrow(res)){
  print(i)
  res[i, ] <- c(spMZ[i], getInfo(spMZ[i]))
}
#rres <- t(res)
tmp <- as.data.frame(res)
names(tmp) <- c("species", "kingdom", "phylum", "class")
write.csv(tmp, file = "../../data/speciesTaxonomy.csv", row.names = FALSE)
