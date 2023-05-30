# Carri J. LeRoy, 2-3-23

# Elwha fungal amplicon data from Argonne labs, 
# Clean up Fungal Amplicon Sequencing data
# Remove samples below sequencing depth cutoff
# Remove low abundance OTUs by prevalence
# Convert to proportional abundances

# Update R and RStudio: 
install.packages("installr")
library(installr)
updateR()
# To update RStudio, go to "Help" and click "Check for Updates"
# To update packages, go to "Tools" and "Check for Package Updates"

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
if (!require("magrittr", quietly = TRUE))
  install.packages("magrittr")
if (!require("vegan", quietly = TRUE))
  install.packages("vegan")
if(!require("BiocManager")){install.packages("BiocManager")}
BiocManager::install("phyloseq")
if (!require("microbiome", quietly = TRUE))
  install.packages("microbiome")

library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(microbiome)

# Import the three ITS files created from the Argonne QIIME output 
#(make sure both Column headers and row names are used for both OTU data and sample data)
# ITSflip.OTU.table.csv, ITSflip.taxonomy.table.csv, ITS.sample.data.csv
# then convert each into a phyloseq object


# Simplify each dataset name:
OTU1 <- ITSflip.OTU.table
Tax1 <- ITSflip.taxonomy.table
Samp1 <- ITS.sample.data


# Make some into matrices:
OTU2 = as.matrix(OTU1)
Tax2 = as.matrix(Tax1)

# Convert to phyloseq objects: 
OTU = otu_table(OTU2, taxa_are_rows = FALSE)
TAX = tax_table(Tax2)
Sample = sample_data(Samp1)
phylo = phyloseq(OTU, TAX, Sample)
phylo

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 326 taxa and 48 samples ]
#sample_data() Sample Data:       [ 48 samples by 45 sample variables ]
#tax_table()   Taxonomy Table:    [ 326 taxa by 7 taxonomic ranks ]

# Check names
table(rownames(Sample) %in% sample_names(OTU))
sample_names(OTU)
#
# Draw a plot - too big to show
plot_bar(phylo, fill = "Class")

#### Check meta data to confirm controls are removed
meta <- phylo %>% sample_data %>% data.frame()

minDepth <- 250
data.frame(SeqDepth=sort(sample_sums(phylo)), Region=sample_data(phylo)$region) %>%
  mutate(cutoff=SeqDepth>minDepth) %>%
  ggplot(aes(x=Region, y=SeqDepth)) +
  geom_violin() +
  geom_point(aes(color=cutoff),position=position_jitter(width=0.1)) +
  theme_classic()

#### Remove samples below sequencing depth cutoff
phylo1 <- phylo %<>% prune_samples(sample_sums(.)>minDepth,.)

#### Remove low abundance OTUs by prevalance
# Only keep OTUs present in at least 1% of samples
phylo2 <- phylo1 %>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE)

phylo2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 314 taxa and 48 samples ]
#sample_data() Sample Data:       [ 48 samples by 45 sample variables ]
#tax_table()   Taxonomy Table:    [ 314 taxa by 7 taxonomic ranks ]

saveRDS(phylo2, "phy_Elwha_Integer.rds") 

main = as(otu_table(phylo2), "matrix")
samp1 = as(sample_data(phylo2), "matrix")
second = as.data.frame(samp1)

write.csv(main, 'main_Elwha_Integer.csv')
write.csv(second, 'second_Elwha_Integer.csv')

## Do DiffAbundance test with Integer data before moving on to the next step! 

#### Convert to proportional abundance (maybe skip this step for Differential Abundance in Aldex2?)
phylo3 <- phylo2 %<>% transform_sample_counts(function(x){x*min(sample_sums(.)/sum(x))})
phylo3

saveRDS(phylo3, "phy_Elwha.rds")

phylo3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 314 taxa and 48 samples ]
#sample_data() Sample Data:       [ 48 samples by 45 sample variables ]
#tax_table()   Taxonomy Table:    [ 314 taxa by 7 taxonomic ranks ]

#### Remove samples below sequencing depth cutoff
phylo4 <- phylo3 %<>% prune_samples(sample_sums(.)>minDepth,.)

phylo4
#### Remove low abundance OTUs by prevalance
# Only keep OTUs present in at least 1% of samples
phylo5 <- phylo4 %>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE)
phylo5

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 314 taxa and 48 samples ]
#sample_data() Sample Data:       [ 48 samples by 45 sample variables ]
#tax_table()   Taxonomy Table:    [ 314 taxa by 7 taxonomic ranks ]

main2 = as(otu_table(phylo5), "matrix")
samp2 = as(sample_data(phylo5), "matrix")
second2 = as.data.frame(samp2)

write.csv(main2, 'main_Elwha_ALL.csv')
write.csv(second2, 'second_Elwha_ALL.csv')