#' Carri J. LeRoy, PhD 3-12-23
#' 
#' 2016 Elwha study. Calculate species richness and Simpson's diversity for fungi on leaf litter (ITS OTUs)
#' Best to use transposed data: "phy_Elwha_Integer.rds"
#' 
#' # Update R and RStudio: 
install.packages("installr")
library(installr)
updateR()
# To update RStudio, go to "Help" and click "Check for Updates"
# To update packages, go to "Tools" and "Check for Package Updates"

if(!require(vegan)){install.packages("vegan")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(magrittr)){install.packages("magrittr")}
if(!require(here)){install.packages("here")}
if(!require(lme4)){install.packages("lme4")}
if(!require(ggthemes)){install.packages("ggthemes")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(iNEXT)){install.packages("iNEXT")}
if(!require(gt)){install.packages("gt")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require("BiocManager")){install.packages("BiocManager")}
BiocManager::install("phyloseq")
a
devtools::install_github("vmikk/metagMisc")

library(doMC)
library(vegan)
library(tidyverse)
library(phyloseq)
library(magrittr)
library(here)
library(lme4)
library(ggthemes)
library(ggpubr)
library(iNEXT)
library(gt)
library(ggplot2)
library(metagMisc)
#
### Read in phyloseq object
phyF <- readRDS("phy_Elwha_Integer.rds") # skip this step if data are already loaded as phyF

# remove empty OTUs from phyloseq object
sweepOTUs <- function(phyF){
  prune_taxa(taxa_sums(phyF) > 0,phyF)
}

phyF %<>% sweepOTUs

# Output phyB phyloseq .rds as two csv files: 
main = as(otu_table(phyF), "matrix")
samp1 = as(sample_data(phyF), "matrix")
second = as.data.frame(samp1)

# Save as .csv files so you can use for other code
write.csv(main, 'main_Elwha_Nozeros.csv') 
write.csv(second, 'second_Elwha_Nozeros.csv') 

# Then, transpose the dataframe so OTUs are rows, and samples are columns: 
# Make sure that there are no OTUs with zeros? In excel, then Import the file: "iNext_Main_Elwha_Nozeros_flip.csv"

# Import the newly formatted iNEXT file in: "iNext_Elwha_update.csv" 
# - make sure header = Yes, and the first column = row names
# The code below might take a long time to run! Be prepared! 
main2 <- iNext_Elwha_update
out <- iNEXT(main2, datatype="abundance")
out
out$iNextEst$size_based
out$iNextEst$coverage_based
out$AsyEst

# This finally worked! Now, get it out of iNEXT:

size = as.data.frame(out$iNextEst$size_based)
cover = as.data.frame(out$iNextEst$coverage_based)
estimates = as.data.frame(out$AsyEst)

write.csv(size, 'ITS1_size-based.csv')
write.csv(cover, 'ITS1_coverage-based.csv')
write.csv(estimates, 'ITS1_estimates.csv')