# Carri J. LeRoy, 3-10-23
# FUNGI: ITS Elwha Decomp

#Differential abundance testing using the Aldex2 pipeline. 
#modular pipeline= aldex.clr --> aldex.ttest (for two groups)/ aldex.glm --> 
#aldex.effect --> aldex.plot

# Update R and RStudio: 
install.packages("installr")
library(installr)
updateR()
# To update RStudio, go to "Help" and click "Check for Updates"
# To update packages, go to "Tools" and "Check for Package Updates"

if (!require("magrittr")) install.packages("magrittr")
if (!require("ALDEx2")) install.packages("ALDEx2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("BiocManager")) install.packages("BiocManager")
library(BiocManager)
BiocManager::install("ALDEx2")
BiocManager::install("phyloseq")
#
#libraries
library(phyloseq)
library(magrittr)
library(ALDEx2)
library(reshape2)

###########
#read in phy object and filter (skip this step if coming in with phylo already in place)
phylo <- readRDS("phy_Elwha_Integer.rds")

phylo %<>% filter_taxa(., function(x) {sum(x>0) > 0.03*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.)
sampledata <- sample_data(phylo)%>% data.frame

#extract OTU table matrix from phyloseq object 
otu <- (as(otu_table(phylo, taxa_are_rows = FALSE), 'matrix'))
otu %<>% t()

#Make vectors of explanatory variables to test #this is for aldex.glm()
harvest <- sample_data(phylo)$harvest %>% as.character()
region <- sample_data(phylo)$region %>% as.character()
location <- sample_data(phylo)$location %>% as.character()
mills <- sample_data(phylo)$mills %>% as.character()
trib <- sample_data(phylo)$trib %>% as.character()
middle <- sample_data(phylo)$middle %>% as.character()
reservoir <- sample_data(phylo)$reservoir %>% as.character()
lower <- sample_data(phylo)$lower %>% as.character()
delta <- sample_data(phylo)$delta %>% as.character()

#If not all reads are integers round to 0 digits
#UNHASH# otuR <- round(otu,0) # if you need to round, change "otu" to "otuR" in the code below

otu.clr.harvest <- aldex.clr(otu, conds = harvest, mc.samples = 128)
otu.clr.region <- aldex.clr(otu, conds = region, mc.samples = 128)
otu.clr.location <- aldex.clr(otu, conds = location, mc.samples = 128)
otu.clr.mills <- aldex.clr(otu, conds = mills, mc.samples = 128)
otu.clr.trib <- aldex.clr(otu, conds = trib, mc.samples = 128)
otu.clr.middle <- aldex.clr(otu, conds = middle, mc.samples = 128)
otu.clr.reservoir <- aldex.clr(otu, conds = reservoir, mc.samples = 128)
otu.clr.lower <- aldex.clr(otu, conds = lower, mc.samples = 128)
otu.clr.delta <- aldex.clr(otu, conds = delta, mc.samples = 128)


###########
###Harvest
###########

ttest.harvest <- aldex.ttest(otu.clr.harvest)
head(ttest.harvest)
#gives you p-values and such in a table. will need to filter. 
#we.ep    we.eBH      wi.ep    wi.eBH
#OTU.1 0.02753902 0.3024985 0.02952823 0.2966020
#OTU.2 0.30692429 0.7639201 0.38186259 0.8321250
#OTU.3 0.48826925 0.8662602 0.51093766 0.8737292
#OTU.4 0.50736816 0.8593259 0.52140618 0.8852658
#OTU.6 0.01455829 0.3352627 0.01889812 0.3526463
#OTU.7 0.54671271 0.8694318 0.54315194 0.8869297


#effect size estimates
effect.harvest <-aldex.effect(otu.clr.harvest)
head(effect.harvest)

#rab.all  rab.win.1  rab.win.3   diff.btw diff.win      effect
#OTU.1  0.3577324  1.6996079 -0.4736825 -2.6346826 4.880413 -0.47646961
#OTU.2  1.3681244  1.6740030  1.0592874 -1.0574407 5.518494 -0.18743628
#OTU.3 -0.5839112 -0.7159638 -0.4520548  0.2778376 4.704403  0.05093705
#OTU.4 -0.4365343 -0.5457439 -0.2870349  0.1963230 4.998761  0.03278106
#OTU.6  8.0921997  6.9987448  8.8226358  1.8342069 3.789246  0.43869706
#OTU.7 -0.7239828 -0.6433250 -0.8070520 -0.2258273 4.605801 -0.04281104
#overlap
#OTU.1 0.2838542
#OTU.2 0.4144437
#OTU.3 0.4687500
#OTU.4 0.4791667
#OTU.6 0.2949219
#OTU.7 0.4762524

#combine ttest and effect sizes
ttest.effect.harvest <- data.frame(ttest.harvest, effect.harvest)

#make a table from the results
alpha = 0.05
sigHarvest = ttest.effect.harvest[which(ttest.effect.harvest$we.eBH < -1), ]
sigHarvest = cbind(as(sigHarvest, "data.frame"), as(tax_table(phylo)[rownames(sigHarvest), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent

## Maybe this is because there are no significant OTUs for Harvest?...
head(sigHarvest)
dim(sigHarvest)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigHarvest, "aldexHarvestA.Rda")
write.csv(sigHarvest, "AldexHarvestA.csv")

#Basically wrote an empty csv file - no significant OTUs for harvest. 

sigHarvest = ttest.effect.harvest[which(ttest.effect.harvest$we.eBH > 1), ]
sigHarvest = cbind(as(sigHarvest, "data.frame"), as(tax_table(phylo)[rownames(sigHarvest), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for Harvest...
head(sigHarvest)
dim(sigHarvest)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigHarvest, "aldexHarvestB.Rda")
write.csv(sigHarvest, "AldexSigHarvestB.csv")

sigHarvest = ttest.effect.harvest[which(ttest.effect.harvest$we.eBH < alpha), ] #this instead gets p-values
sigHarvest = cbind(as(sigHarvest, "data.frame"), as(tax_table(phylo)[rownames(sigHarvest), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for Harvest...
head(sigHarvest)
dim(sigHarvest)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigHarvest, "aldexHarvestC.Rda")
write.csv(sigHarvest, "AldexSigHarvestC.csv")

sigHarvest = ttest.effect.harvest[which(ttest.effect.harvest$wi.eBH < alpha), ] #this instead gets p-values
sigHarvest = cbind(as(sigHarvest, "data.frame"), as(tax_table(phylo)[rownames(sigHarvest), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for Harvest...
head(sigHarvest)
dim(sigHarvest)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigHarvest, "aldexHarvestD.Rda")
write.csv(sigHarvest, "AldexSigHarvestD.csv")

#four CSvs were exported, one for positive and one for negative effect sizes greater than one, two for p-values. those data frames were combined in excel and saved as "AldexSigFungi.csv"

###########
###Region and Location won't work because there are more than two categories! 
###Use dummy variables to compare each region with the rest? Add dummies to ITS.sample.data
### added: mills, trib, middle, reservoir, lower, delta as dummy variables to ITS.sample.data.csv
###########

###########
### Mills
###########

ttest.mills <- aldex.ttest(otu.clr.mills)
head(ttest.mills)
#gives you p-values and such in a table. will need to filter. 
#we.ep    we.eBH     wi.ep    wi.eBH
#OTU.1 0.2779724 0.7137356 0.3003207 0.7970735
#OTU.2 0.4208103 0.8388225 0.5001709 0.9144417
#OTU.3 0.5125539 0.8661104 0.5357434 0.9206955
#OTU.4 0.4474228 0.8534237 0.4891560 0.8964845
#OTU.6 0.8140679 0.9624061 0.9352506 0.9903976
#OTU.7 0.4260543 0.8002913 0.4634362 0.8800385

#effect size estimates
effect.mills <-aldex.effect(otu.clr.mills)
head(effect.mills)

#rab.all rab.win.Mills rab.win.other   diff.btw diff.win      effect
#OTU.1  0.3510784    -0.6329898     0.6210594  1.7560572 5.204090  0.30149658
#OTU.2  1.3898662     1.0903786     1.4926699  0.7152278 4.849120  0.16515799
#OTU.3 -0.6050403    -0.4355204    -0.6327085 -0.1831472 4.559220 -0.03034977
#OTU.4 -0.4454145    -0.9415452    -0.3046546  0.6099965 4.724284  0.12897659
#OTU.6  8.0716301     7.9686545     8.1207714  0.2301255 3.206669  0.05646918
#OTU.7 -0.6659325    -0.0610876    -0.7515510 -0.7805481 4.525089 -0.15349484
#overlap
#OTU.1 0.3574219
#OTU.2 0.4121094
#OTU.3 0.4873294
#OTU.4 0.4366472
#OTU.6 0.4707031
#OTU.7 0.4269006

#combine ttest and effect sizes
ttest.effect.mills <- data.frame(ttest.mills, effect.mills)

#make a table from the results
alpha = 0.05
sigmills = ttest.effect.mills[which(ttest.effect.mills$we.eBH < -1), ]
sigmills = cbind(as(sigmills, "data.frame"), as(tax_table(phylo)[rownames(sigmills), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent

## Maybe this is because there are no significant OTUs for mills?...
head(sigmills)
dim(sigmills)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmills, "aldexmillsA.Rda")
write.csv(sigmills, "AldexmillsA.csv")

#Basically wrote an empty csv file - likely because of the ERROR in dimnames. 

sigmills = ttest.effect.mills[which(ttest.effect.mills$we.eBH > 1), ]
sigmills = cbind(as(sigmills, "data.frame"), as(tax_table(phylo)[rownames(sigmills), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for mills...
head(sigmills)
dim(sigmills)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmills, "aldexmillsB.Rda")
write.csv(sigmills, "AldexSigmillsB.csv")

sigmills = ttest.effect.mills[which(ttest.effect.mills$we.eBH < alpha), ] #this instead gets p-values
sigmills = cbind(as(sigmills, "data.frame"), as(tax_table(phylo)[rownames(sigmills), ], "matrix"))
head(sigmills)
dim(sigmills)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmills, "aldexmillsC.Rda")
write.csv(sigmills, "AldexSigmillsC.csv")

sigmills = ttest.effect.mills[which(ttest.effect.mills$wi.eBH < alpha), ] #this instead gets p-values
sigmills = cbind(as(sigmills, "data.frame"), as(tax_table(phylo)[rownames(sigmills), ], "matrix"))
#Error in dimnames(x) <- dn : 
# length of 'dimnames' [1] not equal to array extent
head(sigmills)
dim(sigmills)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmills, "aldexmillsD.Rda")
write.csv(sigmills, "AldexSigmillsD.csv")

#four CSvs were exported, one for positive and one for negative effect sizes greater than one, two for p-values. those data frames were combined in excel and saved as "AldexSigFungi.csv"

###########
### Trib
###########

ttest.trib <- aldex.ttest(otu.clr.trib)
head(ttest.trib)
#gives you p-values and such in a table. will need to filter. 
#we.ep      we.eBH        wi.ep      wi.eBH
#OTU.1 1.270921e-01 0.432100955 7.397957e-02 0.356554519
#OTU.2 1.314349e-02 0.167098200 1.555271e-03 0.045151281
#OTU.3 5.371451e-01 0.812530502 5.751046e-01 0.833583576
#OTU.4 2.290855e-01 0.613576316 2.024254e-01 0.532751391
#OTU.6 1.374048e-06 0.000288492 1.163483e-05 0.002123895
#OTU.7 5.018335e-01 0.793713429 5.218122e-01 0.787682323

#effect size estimates
effect.trib <-aldex.effect(otu.clr.trib)
head(effect.trib)

#rab.all rab.win.other rab.win.Tributary   diff.btw diff.win
#OTU.1  0.3897636     0.0400879         3.2746408  2.9165177 4.745342
#OTU.2  1.4135386     0.9016680         6.0619423  5.3756451 4.645574
#OTU.3 -0.6169667    -0.5831648        -0.7326124 -0.2545742 4.415128
#OTU.4 -0.4547075    -0.6752536         1.0217102  2.0037133 5.601143
#OTU.6  8.0771243     7.5569536        11.6236824  3.8902110 2.601436
#OTU.7 -0.7121120    -0.6711382        -1.1075492 -0.3930460 4.486716
#effect    overlap
#OTU.1  0.53381643 0.27290457
#OTU.2  1.06529295 0.13085963
#OTU.3 -0.03532150 0.48538012
#OTU.4  0.33303605 0.33789068
#OTU.6  1.36649212 0.06237876
#OTU.7 -0.07969327 0.45614036

#combine ttest and effect sizes
ttest.effect.trib <- data.frame(ttest.trib, effect.trib)

#make a table from the results
alpha = 0.05
sigtrib = ttest.effect.trib[which(ttest.effect.trib$we.eBH < -1), ]
sigtrib = cbind(as(sigtrib, "data.frame"), as(tax_table(phylo)[rownames(sigtrib), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent

## Maybe this is because there are no significant OTUs for trib?...
head(sigtrib)
dim(sigtrib)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigtrib, "aldextribA.Rda")
write.csv(sigtrib, "AldextribA.csv")

#Basically wrote an empty csv file - likely because of the ERROR in dimnames. 

sigtrib = ttest.effect.trib[which(ttest.effect.trib$we.eBH > 1), ]
sigtrib = cbind(as(sigtrib, "data.frame"), as(tax_table(phylo)[rownames(sigtrib), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for trib...
head(sigtrib)
dim(sigtrib)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigtrib, "aldextribB.Rda")
write.csv(sigtrib, "AldexSigtribB.csv")

sigtrib = ttest.effect.trib[which(ttest.effect.trib$we.eBH < alpha), ] #this instead gets p-values
sigtrib = cbind(as(sigtrib, "data.frame"), as(tax_table(phylo)[rownames(sigtrib), ], "matrix"))
head(sigtrib)
dim(sigtrib)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigtrib, "aldextribC.Rda")
write.csv(sigtrib, "AldexSigtribC.csv")

sigtrib = ttest.effect.trib[which(ttest.effect.trib$wi.eBH < alpha), ] #this instead gets p-values
sigtrib = cbind(as(sigtrib, "data.frame"), as(tax_table(phylo)[rownames(sigtrib), ], "matrix"))
head(sigtrib)
dim(sigtrib)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigtrib, "aldextribD.Rda")
write.csv(sigtrib, "AldexSigtribD.csv")

#four CSvs were exported, one for positive and one for negative effect sizes greater than one, two for p-values. those data frames were combined in excel and saved as "AldexSigFungi.csv"

###########
### Middle
###########

ttest.middle <- aldex.ttest(otu.clr.middle)
head(ttest.middle)
#gives you p-values and such in a table. will need to filter. 
#we.ep    we.eBH     wi.ep    wi.eBH
#OTU.1 0.3475170 0.8590048 0.4178868 0.8917983
#OTU.2 0.3386818 0.8447126 0.3987976 0.8897373
#OTU.3 0.5003500 0.8858835 0.5126114 0.9053571
#OTU.4 0.5279126 0.8947900 0.5258660 0.9189863
#OTU.6 0.7981857 0.9649661 0.9146205 0.9869058
#OTU.7 0.5238419 0.8946209 0.5299697 0.9209898

#effect size estimates
effect.middle <-aldex.effect(otu.clr.middle)
head(effect.middle)

#rab.all rab.win.Middle rab.win.other    diff.btw diff.win
#OTU.1  0.4247476     -0.2028684     0.5723986  1.07901328 4.931989
#OTU.2  1.3436424      0.7105490     1.4894136  1.17803559 4.827669
#OTU.3 -0.5478190     -0.1051562    -0.6472132 -0.41433778 4.456525
#OTU.4 -0.4489936     -0.1955719    -0.5188859 -0.30270638 4.869269
#OTU.6  8.0669962      8.3775402     8.0669962  0.02048122 4.727873
#OTU.7 -0.7249516     -0.8004798    -0.7173013  0.07776128 4.415987
#effect   overlap
#OTU.1  0.197913459 0.4082032
#OTU.2  0.220277641 0.3879143
#OTU.3 -0.074992545 0.4550781
#OTU.4 -0.040473121 0.4658870
#OTU.6  0.002381826 0.4990253
#OTU.7  0.012050968 0.4912281

#combine ttest and effect sizes
ttest.effect.middle <- data.frame(ttest.middle, effect.middle)

#make a table from the results
alpha = 0.05
sigmiddle = ttest.effect.middle[which(ttest.effect.middle$we.eBH < -1), ]
sigmiddle = cbind(as(sigmiddle, "data.frame"), as(tax_table(phylo)[rownames(sigmiddle), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent

## Maybe this is because there are no significant OTUs for middle?...
head(sigmiddle)
dim(sigmiddle)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmiddle, "aldexmiddleA.Rda")
write.csv(sigmiddle, "AldexmiddleA.csv")

#Basically wrote an empty csv file - likely because of the ERROR in dimnames. 

sigmiddle = ttest.effect.middle[which(ttest.effect.middle$we.eBH > 1), ]
sigmiddle = cbind(as(sigmiddle, "data.frame"), as(tax_table(phylo)[rownames(sigmiddle), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for middle...
head(sigmiddle)
dim(sigmiddle)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmiddle, "aldexmiddleB.Rda")
write.csv(sigmiddle, "AldexSigmiddleB.csv")

sigmiddle = ttest.effect.middle[which(ttest.effect.middle$we.eBH < alpha), ] #this instead gets p-values
sigmiddle = cbind(as(sigmiddle, "data.frame"), as(tax_table(phylo)[rownames(sigmiddle), ], "matrix"))
#Error in dimnames(x) <- dn : 
#  length of 'dimnames' [1] not equal to array extent
head(sigmiddle)
dim(sigmiddle)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmiddle, "aldexmiddleC.Rda")
write.csv(sigmiddle, "AldexSigmiddleC.csv")

sigmiddle = ttest.effect.middle[which(ttest.effect.middle$wi.eBH < alpha), ] #this instead gets p-values
sigmiddle = cbind(as(sigmiddle, "data.frame"), as(tax_table(phylo)[rownames(sigmiddle), ], "matrix"))
#Error in dimnames(x) <- dn : 
#  length of 'dimnames' [1] not equal to array extent
head(sigmiddle)
dim(sigmiddle)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigmiddle, "aldexmiddleD.Rda")
write.csv(sigmiddle, "AldexSigmiddleD.csv")

#four CSvs were exported, one for positive and one for negative effect sizes greater than one, two for p-values. those data frames were combined in excel and saved as "AldexSigFungi.csv"

###########
### Reservoir
###########

ttest.reservoir <- aldex.ttest(otu.clr.reservoir)
head(ttest.reservoir)
#gives you p-values and such in a table. will need to filter. 
#we.ep    we.eBH      wi.ep    wi.eBH
#OTU.1 0.6184990 0.8917080 0.56146814 0.8828284
#OTU.2 0.3145613 0.7478950 0.32720490 0.7584169
#OTU.3 0.5072564 0.8157467 0.54183731 0.8573765
#OTU.4 0.4330363 0.7888621 0.45060152 0.8097128
#OTU.6 0.0536767 0.4976956 0.04295474 0.4961187
#OTU.7 0.5127284 0.8302034 0.53605764 0.8619018

#effect size estimates
effect.reservoir <-aldex.effect(otu.clr.reservoir)
head(effect.reservoir)

#rab.all rab.win.Aldwell rab.win.other   diff.btw diff.win
#OTU.1  0.4171245      -0.1192630     0.5371203  0.5920090 5.278936
#OTU.2  1.4524302       0.2222762     1.6440774  1.5164695 5.465095
#OTU.3 -0.6305538      -0.4837240    -0.6590755 -0.2244486 5.245427
#OTU.4 -0.3992919      -0.9072240    -0.2582541  0.7061878 4.467633
#OTU.6  8.0823399       6.5886880     8.2396394  2.4194626 3.951294
#OTU.7 -0.6712970      -0.5809298    -0.6871337 -0.3699275 4.653096
#effect   overlap
#OTU.1  0.09372871 0.4570313
#OTU.2  0.27711932 0.3801170
#OTU.3 -0.03361607 0.4834308
#OTU.4  0.11960668 0.4405458
#OTU.6  0.56817712 0.2904484
#OTU.7 -0.05948431 0.4628906

#combine ttest and effect sizes
ttest.effect.reservoir <- data.frame(ttest.reservoir, effect.reservoir)

#make a table from the results
alpha = 0.05
sigreservoir = ttest.effect.reservoir[which(ttest.effect.reservoir$we.eBH < -1), ]
sigreservoir = cbind(as(sigreservoir, "data.frame"), as(tax_table(phylo)[rownames(sigreservoir), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent

## Maybe this is because there are no significant OTUs for reservoir?...
head(sigreservoir)
dim(sigreservoir)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigreservoir, "aldexreservoirA.Rda")
write.csv(sigreservoir, "AldexreservoirA.csv")

#Basically wrote an empty csv file - likely because of the ERROR in dimnames. 

sigreservoir = ttest.effect.reservoir[which(ttest.effect.reservoir$we.eBH > 1), ]
sigreservoir = cbind(as(sigreservoir, "data.frame"), as(tax_table(phylo)[rownames(sigreservoir), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for reservoir...
head(sigreservoir)
dim(sigreservoir)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigreservoir, "aldexreservoirB.Rda")
write.csv(sigreservoir, "AldexSigreservoirB.csv")

sigreservoir = ttest.effect.reservoir[which(ttest.effect.reservoir$we.eBH < alpha), ] #this instead gets p-values
sigreservoir = cbind(as(sigreservoir, "data.frame"), as(tax_table(phylo)[rownames(sigreservoir), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for reservoir...
head(sigreservoir)
dim(sigreservoir)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigreservoir, "aldexreservoirC.Rda")
write.csv(sigreservoir, "AldexSigreservoirC.csv")

sigreservoir = ttest.effect.reservoir[which(ttest.effect.reservoir$wi.eBH < alpha), ] #this instead gets p-values
sigreservoir = cbind(as(sigreservoir, "data.frame"), as(tax_table(phylo)[rownames(sigreservoir), ], "matrix"))
head(sigreservoir)
dim(sigreservoir)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigreservoir, "aldexreservoirD.Rda")
write.csv(sigreservoir, "AldexSigreservoirD.csv")

#four CSvs were exported, one for positive and one for negative effect sizes greater than one, two for p-values. those data frames were combined in excel and saved as "AldexSigFungi.csv"

###########
### Lower
###########

ttest.lower <- aldex.ttest(otu.clr.lower)
head(ttest.lower)
#gives you p-values and such in a table. will need to filter. 
#we.ep    we.eBH      wi.ep    wi.eBH
#OTU.1 0.56921350 0.8705635 0.59062867 0.8893920
#OTU.2 0.68755376 0.9218902 0.67582563 0.9253251
#OTU.3 0.49682583 0.8381242 0.52357506 0.8667519
#OTU.4 0.52912752 0.8489609 0.56666497 0.8793782
#OTU.6 0.03001064 0.4090658 0.04180328 0.5252527
#OTU.7 0.44796965 0.7960363 0.48166724 0.8448030

#effect size estimates
effect.lower <-aldex.effect(otu.clr.lower)
head(effect.lower)

#rab.all rab.win.Lower rab.win.other    diff.btw diff.win      effect
#OTU.1  0.4352081     0.6029135     0.3895748 -0.36039143 5.705567 -0.05104772
#OTU.2  1.4057849     1.4736795     1.3933433 -0.19689744 5.441156 -0.03787294
#OTU.3 -0.6623315    -0.8066489    -0.6378294  0.12783991 4.577553  0.02303788
#OTU.4 -0.4195182    -0.5240394    -0.3767754  0.45043808 4.417676  0.07231169
#OTU.6  8.0653566     6.5286092     8.4992091  1.74658695 2.730545  0.55176412
#OTU.7 -0.6830480    -0.4801879    -0.7164246 -0.03328811 4.230552 -0.00783838
#overlap
#OTU.1 0.4736842
#OTU.2 0.4824219
#OTU.3 0.4804688
#OTU.4 0.4561404
#OTU.6 0.2695313
#OTU.7 0.4970760

#combine ttest and effect sizes
ttest.effect.lower <- data.frame(ttest.lower, effect.lower)

#make a table from the results
alpha = 0.05
siglower = ttest.effect.lower[which(ttest.effect.lower$we.eBH < -1), ]
siglower = cbind(as(siglower, "data.frame"), as(tax_table(phylo)[rownames(siglower), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent

## Maybe this is because there are no significant OTUs for lower?...
head(siglower)
dim(siglower)  #dimensions of the table. First number is how many significant OTUs

saveRDS(siglower, "aldexlowerA.Rda")
write.csv(siglower, "AldexlowerA.csv")

#Basically wrote an empty csv file - likely because of the ERROR in dimnames. 

siglower = ttest.effect.lower[which(ttest.effect.lower$we.eBH > 1), ]
siglower = cbind(as(siglower, "data.frame"), as(tax_table(phylo)[rownames(siglower), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for lower...
head(siglower)
dim(siglower)  #dimensions of the table. First number is how many significant OTUs

saveRDS(siglower, "aldexlowerB.Rda")
write.csv(siglower, "AldexSiglowerB.csv")

siglower = ttest.effect.lower[which(ttest.effect.lower$we.eBH < alpha), ] #this instead gets p-values
siglower = cbind(as(siglower, "data.frame"), as(tax_table(phylo)[rownames(siglower), ], "matrix"))
#Error in dimnames(x) <- dn : 
# length of 'dimnames' [1] not equal to array extent
head(siglower)
dim(siglower)  #dimensions of the table. First number is how many significant OTUs

saveRDS(siglower, "aldexlowerC.Rda")
write.csv(siglower, "AldexSiglowerC.csv")

siglower = ttest.effect.lower[which(ttest.effect.lower$wi.eBH < alpha), ] #this instead gets p-values
siglower = cbind(as(siglower, "data.frame"), as(tax_table(phylo)[rownames(siglower), ], "matrix"))
#Error in dimnames(x) <- dn : 
# length of 'dimnames' [1] not equal to array extent
head(siglower)
dim(siglower)  #dimensions of the table. First number is how many significant OTUs

saveRDS(siglower, "aldexlowerD.Rda")
write.csv(siglower, "AldexSiglowerD.csv")

#four CSvs were exported, one for positive and one for negative effect sizes greater than one, two for p-values. those data frames were combined in excel and saved as "AldexSigFungi.csv"

###########
### Delta
###########

ttest.delta <- aldex.ttest(otu.clr.delta)
head(ttest.delta)
#gives you p-values and such in a table. will need to filter. 
#we.ep    we.eBH     wi.ep    wi.eBH
#OTU.1 0.5700562 0.9092890 0.6103990 0.9445689
#OTU.2 0.4789355 0.8854984 0.5216931 0.9194938
#OTU.3 0.5273822 0.8837265 0.5126260 0.9088795
#OTU.4 0.5721913 0.8942827 0.5784301 0.9291848
#OTU.6 0.6285770 0.9341403 0.8377885 0.9806089
#OTU.7 0.5085694 0.8912251 0.5062730 0.9095207

#effect size estimates
effect.delta <-aldex.effect(otu.clr.delta)
head(effect.delta)

#rab.all rab.win.Delta rab.win.other   diff.btw diff.win      effect
#OTU.1  0.4049280     0.6918740     0.3597524 -0.3791213 4.984640 -0.06755874
#OTU.2  1.3771801     0.7064754     1.5063642  0.9250984 5.369610  0.15280595
#OTU.3 -0.5867877    -1.0284646    -0.4965963  0.6635642 4.444641  0.12321850
#OTU.4 -0.4457613    -0.5253356    -0.4309281  0.0951178 4.861740  0.01950724
#OTU.6  8.0603235     8.5511381     8.0278407 -0.3927671 3.151678 -0.10615386
#OTU.7 -0.7101479    -1.1022539    -0.6336935  0.4900750 4.329503  0.08852307
#overlap
#OTU.1 0.4639376
#OTU.2 0.4249513
#OTU.3 0.4307992
#OTU.4 0.4912281
#OTU.6 0.4600390
#OTU.7 0.4483431

#combine ttest and effect sizes
ttest.effect.delta <- data.frame(ttest.delta, effect.delta)

#make a table from the results
alpha = 0.05
sigdelta = ttest.effect.delta[which(ttest.effect.delta$we.eBH < -1), ]
sigdelta = cbind(as(sigdelta, "data.frame"), as(tax_table(phylo)[rownames(sigdelta), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent

## Maybe this is because there are no significant OTUs for delta?...
head(sigdelta)
dim(sigdelta)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigdelta, "aldexdeltaA.Rda")
write.csv(sigdelta, "AldexdeltaA.csv")

#Basically wrote an empty csv file - likely because of the ERROR in dimnames. 

sigdelta = ttest.effect.delta[which(ttest.effect.delta$we.eBH > 1), ]
sigdelta = cbind(as(sigdelta, "data.frame"), as(tax_table(phylo)[rownames(sigdelta), ], "matrix"))
## ERROR: Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent
## Maybe this is because there are no significant OTUs for delta...
head(sigdelta)
dim(sigdelta)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigdelta, "aldexdeltaB.Rda")
write.csv(sigdelta, "AldexSigdeltaB.csv")

sigdelta = ttest.effect.delta[which(ttest.effect.delta$we.eBH < alpha), ] #this instead gets p-values
sigdelta = cbind(as(sigdelta, "data.frame"), as(tax_table(phylo)[rownames(sigdelta), ], "matrix"))
#Error in dimnames(x) <- dn : 
#  length of 'dimnames' [1] not equal to array extenthead(sigdelta)
dim(sigdelta)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigdelta, "aldexdeltaC.Rda")
write.csv(sigdelta, "AldexSigdeltaC.csv")

sigdelta = ttest.effect.delta[which(ttest.effect.delta$wi.eBH < alpha), ] #this instead gets p-values
sigdelta = cbind(as(sigdelta, "data.frame"), as(tax_table(phylo)[rownames(sigdelta), ], "matrix"))
#Error in dimnames(x) <- dn : 
#  length of 'dimnames' [1] not equal to array extent
head(sigdelta)
dim(sigdelta)  #dimensions of the table. First number is how many significant OTUs

saveRDS(sigdelta, "aldexdeltaD.Rda")
write.csv(sigdelta, "AldexSigdeltaD.csv")

#four CSvs were exported, one for positive and one for negative effect sizes greater than one, two for p-values. those data frames were combined in excel and saved as "AldexSigFungi.csv"


