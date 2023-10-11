#' Carri J. LeRoy, PhD 2-5-23
#' 
#' Fungal (ITS) community data from 2016 Elwha study. Data in: main_Elwha_ALLFlip.csv and second_Elwha_ALL.csv
#' 
# Update R and RStudio: 
install.packages("installr")
library(installr)
updateR()
# To update RStudio, go to "Help" and click "Check for Updates"
# To update packages, go to "Tools" and "Check for Package Updates"
#' ---
#' ## Packages installed and loaded?
#' 
## ----load-libraries, eval=FALSE------------------------------------------
if(!require(vegan)){install.packages("vegan", dependencies = TRUE)}
if(!require(permute)){install.packages("permute")}
if(!require(MASS)){install.packages("MASS")}
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("gavinsimpson/ggvegan")
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(indicspecies)){install.packages("indicspecies")}
if(!requireNamespace("BiocManager")){install.packages("BiocManager")}
BiocManager::install("phyloseq")
a


library("vegan")
library("permute")
library("MASS")
library("ggvegan")
library("ggplot2")
library("tidyverse")
library("indicspecies")
library("phyloseq")
#' 
#' 
#' Set working directory under Session < Set working directory < choose directory
#' Load two matrices: main_ALL.csv and second_ALL.csv
#' 
#' Can't get the log(x+1) transformation to work in R, so doing it in Excel - all main matrices are already transformed: 
#' 
#' Important that R knows the row names are in column 1 and there are headers. 
#' Simplify your matrix names to "main" (taxa abundances) and "second" (environmental variables)
#' 
## ----loading matrices-----------------------------------------
main <- read.csv("main_Elwha_ALLFlip.csv", row.names = 1, header = TRUE)
second <- read.csv("second_Elwha_ALL.csv", row.names = 1, header = TRUE)

#------or use the matrices already loaded--------------------
#main <- main2
#second <- second2

#' Check headers and first three rows (sample names must be in same order in both matrices)
#' 
## ----head-BCI------------------------------------------------------------
str(main)
str(second)
#' 
head(main[,1:3], n = 3) 
head(second[,1:3], n = 3)
#' 
#' Make sure grouping variables in your second matrix are factors, not characters: 


second$harvest = factor(second$harvest,
                     levels=unique(second$harvest))
second$region = factor(second$region,
                  levels=unique(second$region))
second$location = factor(second$location,
                    levels=unique(second$location))


#' ## PERMANOVA using `adonis`
#' Analysis of variance using distance matrices and for fitting linear models to distance matrices
#' If you use permutations = 9999, you can get p-values down to 0.0001, if only = 999, only p-values down to 0.001
## ----adonis--------------------------------------------------------------
set.seed(1413)
adonis2(main ~ second$harvest*second$region, method = "bray", permutation = 9999)
#
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = main ~ second$harvest * second$region, permutations = 9999, method = "bray")
#                                 Df SumOfSqs      R2      F Pr(>F)    
#   second$harvest                1   0.3073 0.03339 2.1014 0.0600 .  
#   second$region                 5   2.5385 0.27588 3.4721 0.0001 ***
#   second$harvest:second$region  5   1.0917 0.11865 1.4933 0.0709 .  
#   Residual                     36   5.2639 0.57208                  
#   Total                        47   9.2014 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#'#####################
#'# Set your distance measure: "bray-curtis" "euclidean" etc.
dis <- vegdist(main, method = "bray")
# To determine the parsimonious axis number, we generate a scree plot. 
scree <- nmds(dis, mindim = 1, maxdim = 6, nits = 1)
#attributes(scree)
plot(scree$stress, xlab = "Number of Axes")
#'############ PermDisp? ####################
disp <- betadisper(dis, group = second$region)
disp

#Homogeneity of multivariate dispersions
#
#Call: betadisper(d = dis, group = second$region)
#
#No. of Positive Eigenvalues: 28
#No. of Negative Eigenvalues: 19
#
#Average distance to median:
#  Mills Reservoir         Tributary      Middle Elwha Aldwell Reservoir       Lower Elwha             Delta 
#0.09886           0.45927           0.43863           0.34409           0.40415           0.21744 
#
#Eigenvalues for PCoA axes:
#  (Showing 8 of 47 eigenvalues)
#PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#3.6707 1.6269 1.0315 0.7501 0.5661 0.5225 0.3880 0.2912 

permutest(disp)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#             Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#   Groups     5 0.80202 0.160404 7.7315    999  0.001 ***
#   Residuals 42 0.87136 0.020747                         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(disp)
#Analysis of Variance Table
#
#Response: Distances
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     5 0.80202 0.160404  7.7315 3.152e-05 ***
#  Residuals 42 0.87136 0.020747                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

plot(disp, hull=FALSE, ellipse=TRUE)


#'NMS Ordinations
summary(second)

#' ## Basic ordination and plotting
#' 
## ----NMDS-1, results='hide'----------------------------------------------
bray.ord <- metaMDS(main, distance = "bray", k = 3, trymax = 50, autotransform =FALSE)
# This code spits out the dimensions, stress, and information about the best solution: 
bray.ord

#Call:
#  metaMDS(comm = main, distance = "bray", k = 3, trymax = 50, autotransform = FALSE) 
#
#global Multidimensional Scaling using monoMDS
#
#Data:     main 
#Distance: bray 
#
#Dimensions: 3 
#Stress:     0.0914941 
#Stress type 1, weak ties
#Best solution was repeated 6 times in 20 tries
#The best solution was from try 0 (metric scaling or null solution)
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘main’ 
#' 
#' 
#' ## Basic ordination and plotting (using all defaults)
## --------------------------------------------------
plot(bray.ord, display = "sites")

################## This is working! Plotting axes 1 and 2 (changing this doesn't seem to change anything - "choices")
colvec <- c("red1","orange1", "yellow2", "green1", "blue1", "purple1")
shpvec <- c(16, 15, 16, 15, 16, 15)
plot(bray.ord, type = "n", choices=c(1, 2), cex.lab=1.5, cex.axis=1.3)
with(second, points(bray.ord, display = "sites", cex=1.5, col = colvec[region], pch = shpvec[region], bg = colvec[region]))
with(second, legend("topleft", legend = levels(region), bty = "n", col = colvec, pch = shpvec, pt.bg = colvec))
ordihull(bray.ord, groups = second$region, label = FALSE, col = colvec)

#Save as "NMS_Elwha_source.pdf"
#' 
#' ## Vectors in ordination space - vectors (correlations) with environmental variables (need to add location chems here!)
#' 
## -----------------------------------------------------------------
bray.ord.Sum <- envfit(bray.ord ~ Sum, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sum
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2  Pr(>r)  
      #Sum 0.10876 -0.99407 0.0596 0.3197
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Rich <- envfit(bray.ord ~ Rich, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Rich
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2   Pr(>r)   
      #Rich 0.083377 -0.996520 0.0819 0.2298
      #Permutation: free
      #Number of permutations: 1000

bray.ord.RichF <- envfit(bray.ord ~ RichF, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.RichF
#***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)
#RichF  0.024528 -0.999700 0.2078 0.006993 **
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Number of permutations: 1000

bray.ord.SimpF <- envfit(bray.ord ~ SimpF, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.SimpF
#.***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)
#SimpF  0.96692 -0.25507 0.5109 0.000999 ***
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.DivF <- envfit(bray.ord ~ DivF, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.DivF
#***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)
#DivF  0.93499 -0.35467 0.3917 0.000999 ***
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.k <- envfit(bray.ord ~ k, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.k
#n.s.
      #***VECTORS
      #
      #NMDS1     NMDS2     r2  Pr(>r)  
      #k 0.29637 -0.95507 0.0311 0.5554
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Depth <- envfit(bray.ord ~ Depth, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Depth
#***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)
#Depth -0.995190 -0.098012 0.2704 0.007992 **
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.width <- envfit(bray.ord ~ width, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.width
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #width  -0.81282 -0.58252 0.0992 0.1548
      #Permutation: free
      #Number of permutations: 1000
      
      #4 observations deleted due to missingness

bray.ord.CC <- envfit(bray.ord ~ CC, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.CC
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)   
#CC 0.96543 0.26064 0.3199 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Fines <- envfit(bray.ord ~ Fines, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Fines
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2   Pr(>r)    
      #Fines -0.51346  0.85812 0.0521 0.2867
      #Permutation: free
      #Number of permutations: 1000

bray.ord.D50 <- envfit(bray.ord ~ D50, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.D50
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2   Pr(>r)    
      #D50 0.70477 -0.70944 0.0059 0.8671
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Temp <- envfit(bray.ord ~ Temp, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Temp
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2  Pr(>r)  
      #Temp 0.998220 0.059699 0.0717 0.1788
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Cond <- envfit(bray.ord ~ Cond, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Cond
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#Cond 0.94268 0.33369 0.2522 0.003996 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.P <- envfit(bray.ord ~ P, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P
#***VECTORS
#
#NMDS1    NMDS2     r2  Pr(>r)  
#P 0.996110 0.088108 0.3624 0.001998 **
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000
#
#4 observations deleted due to missingness

bray.ord.n <- envfit(bray.ord ~ n, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.n
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #n 0.97446 -0.22458 0.1516 0.05295 .
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.po4 <- envfit(bray.ord ~ po4, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.po4
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2 Pr(>r)
      #po4 0.82721 0.56189 0.0651 0.2378
      #Permutation: free
      #Number of permutations: 1000

bray.ord.SiO4 <- envfit(bray.ord ~ SiO4, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.SiO4
#***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)
#SiO4 0.84199 0.53950 0.3799 0.001998 **
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000
#
#4 observations deleted due to missingness

bray.ord.No3 <- envfit(bray.ord ~ No3, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.No3
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2   Pr(>r)    
      #No3 0.93188 0.36278 0.1185 0.07692 .
      #Permutation: free
      #Number of permutations: 1000

      #4 observations deleted due to missingness

bray.ord.no2 <- envfit(bray.ord ~ no2, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.no2
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2  Pr(>r)  
      #no2 0.86178 0.50729 0.0589 0.3017
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.nh4 <- envfit(bray.ord ~ nh4, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.nh4
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #nh4 0.98432 -0.17639 0.0113 0.7832
      #Permutation: free
      #Number of permutations: 1000

bray.ord.no2no3 <- envfit(bray.ord ~ no2no3, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.no2no3
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #no2No3 1.000000 0.001229 0.0556 0.3147
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Chla <- envfit(bray.ord ~ Chla, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Chla
#***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)
#Chla 0.98529 0.17088 0.1668 0.02198 *
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.omRock <- envfit(bray.ord ~ omRock, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.omRock
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #omRock  -0.99123 -0.13213 0.0252 0.6284
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.imRock <- envfit(bray.ord ~ imRock, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.imRock
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #imRock  -0.54318 -0.83962 0.0447 0.4086
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.OM <- envfit(bray.ord ~ OM, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.OM
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2   Pr(>r)   
      #OM 0.98780 -0.15574 0.0368 0.4306
      #Permutation: free
      #Number of permutations: 1000

bray.ord.imSeston <- envfit(bray.ord ~ imSeston, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.imSeston
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2    r2 Pr(>r)
      #imSeston -0.10171 -0.99481 0.0591 0.2398
      #Permutation: free
      #Number of permutations: 1000

bray.ord.BI <- envfit(bray.ord ~ BI, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.BI
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #BI  -0.59111 -0.80659 0.0588 0.3267
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.SH <- envfit(bray.ord ~ SH, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.SH
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #SH  -0.03583  0.99936 0.06 0.3157
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.rkm <- envfit(bray.ord ~ rkm, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.rkm
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2   Pr(>r)    
      #rkm 0.0045716 0.9999900 0.0822 0.1449
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Elev <- envfit(bray.ord ~ Elev, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Elev
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2   Pr(>r)    
      #Elev -0.33151  0.94345 0.09 0.1039
      #Permutation: free
      #Number of permutations: 1000


bray.ord.Sed1 <- envfit(bray.ord ~ Sed1, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed1
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed1 -0.47480 -0.88009 0.0163 0.7093
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed2 <- envfit(bray.ord ~ Sed2, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed2
#***VECTORS
#
#NMDS1   NMDS2     r2 Pr(>r)
#Sed2 -0.71882  0.69520 0.2686 0.005994 **
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000
#
#4 observations deleted due to missingness

bray.ord.Sed3 <- envfit(bray.ord ~ Sed3, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed3
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed3 -0.84514  0.53454 0.149 0.05395 .
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed4 <- envfit(bray.ord ~ Sed4, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed4
#***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)
#Sed4 -0.96008  0.27974 0.1786 0.01998 *
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000
#
#4 observations deleted due to missingness

bray.ord.Sed5 <- envfit(bray.ord ~ Sed5, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed5
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2 Pr(>r)
      #Sed5 -0.83701  0.54719 0.1402 0.06394 .
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed6 <- envfit(bray.ord ~ Sed6, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed6
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed6  -0.83008  0.55764 0.1159 0.1099
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed7 <- envfit(bray.ord ~ Sed7, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed7
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed7  -0.44347  0.89629 0.0559 0.3437
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed8 <- envfit(bray.ord ~ Sed8, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed8
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed8 0.78407 0.62067 0.114 0.1029
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed9 <- envfit(bray.ord ~ Sed9, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed9
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed9 -0.36660  0.93038 0.0164 0.7353
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed10 <- envfit(bray.ord ~ Sed10, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed10
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed10 -0.93072  0.36573 0.0255 0.6264
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed11 <- envfit(bray.ord ~ Sed11, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed11
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed11 -0.999870  0.016233 7e-04   0.99
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed12 <- envfit(bray.ord ~ Sed12, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed12
#n.s.
      #***VECTORS
      #
      #NMDS1 NMDS2 r2 Pr(>r)
      #Sed12     0     0  0      1
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.SedTot <- envfit(bray.ord ~ SedTot, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.SedTot
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #SedTot -0.83793  0.54577 0.012 0.7952
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P1 <- envfit(bray.ord ~ P1, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P1
#n.s.
      #***VECTORS
      #
      #NMDS1     NMDS2    r2  Pr(>r)  
      #P1  0.057431 -0.998350 0.0444 0.4486
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P2 <- envfit(bray.ord ~ P2, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P2
#***VECTORS
#
#NMDS1   NMDS2     r2 Pr(>r)
#P2 -0.76074  0.64905 0.1712 0.02797 *
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000
#
#4 observations deleted due to missingness

bray.ord.P3 <- envfit(bray.ord ~ P3, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P3
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #P3  -0.97208  0.23465 0.0305 0.5524
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P4 <- envfit(bray.ord ~ P4, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P4
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #P4  -0.86790 -0.49674 0.0355 0.5265
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P5 <- envfit(bray.ord ~ P5, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P5
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2 Pr(>r)
      #P5 0.76859 0.63974 0.0234 0.6474
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P6 <- envfit(bray.ord ~ P6, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P6
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2 Pr(>r)
      #P6 0.92326 0.38417 0.0606 0.3297
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P7 <- envfit(bray.ord ~ P7, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P7
#***VECTORS
#
#NMDS1   NMDS2     r2 Pr(>r)
#P7 0.93057 0.36611 0.1568 0.04396 *
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000
#
#4 observations deleted due to missingness

bray.ord.P8 <- envfit(bray.ord ~ P8, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P8
#***VECTORS
#
#NMDS1     NMDS2     r2 Pr(>r)
#P8  0.98107 0.19366 0.2124 0.01499 *
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000
#
#4 observations deleted due to missingness

bray.ord.P9 <- envfit(bray.ord ~ P9, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P9
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2    r2 Pr(>r)
      #P9 0.51004 -0.86015 0.0545 0.3257
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P10 <- envfit(bray.ord ~ P10, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P10
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #P10 -0.84147  0.54031 0.047 0.4236
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P11 <- envfit(bray.ord ~ P11, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P11
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #P11 -0.98809 -0.15390 0.0477 0.4056
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.EvenM <- envfit(bray.ord ~ EvenM, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.EvenM
#***VECTORS
#
#NMDS1     NMDS2    r2   Pr(>r)    
#Even 0.47052 0.88239 0.2703 0.004995 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.DivM <- envfit(bray.ord ~ DivM, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.DivM
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#DivM 0.59336 0.80494 0.2608 0.003996 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.SimpsonM <- envfit(bray.ord ~ SimpsonM, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.SimpsonM
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)   
#SimpsonM 0.51542 0.85694 0.2657 0.003996 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Fine <- envfit(bray.ord ~ Fine, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Fine
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2  Pr(>r)  
      #Fine -0.98896  0.14817 0.0487 0.4096
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Gravel <- envfit(bray.ord ~ Gravel, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Gravel
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2   Pr(>r)   
      #Gravel 0.98317 -0.18271 0.1164 0.0979 .
      #Permutation: free
      #Number of permutations: 1000

bray.ord.Cobble <- envfit(bray.ord ~ Cobble, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Cobble
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Cobble -0.98298  0.18374 0.1071 0.1149
      #Number of permutations: 1000



# Add ALL significant vectors to plot
colvec <- c("red1","orange2", "gold", "green1", "blue1", "purple1")
shpvec <- c(16, 15, 16, 15, 16, 15)
plot(bray.ord, type = "n", choices=c(1,2), cex.lab=1.5, cex.axis=1.3)
with(second, points(bray.ord, display = "sites", cex=1.5, col = colvec[region], pch = shpvec[region], bg = colvec[region]))
with(second, legend("topright", legend = levels(region), bty = "n", col = colvec, pch = shpvec, pt.bg = colvec))
ordihull(bray.ord, groups = second$region, label = FALSE, col = colvec)
plot(bray.ord.Sum, add = TRUE, cex=1.5)
plot(bray.ord.Rich, add = TRUE, cex=1.5)
plot(bray.ord.DivM, add = TRUE, cex=1.5)
plot(bray.ord.Elev, add = TRUE, cex=1.5)
plot(bray.ord.Fines, add = TRUE, cex=1.5)
plot(bray.ord.D50, add = TRUE, cex=1.5)
plot(bray.ord.OM, add = TRUE, cex=1.5)
plot(bray.ord.Temp, add = TRUE, cex=1.5)
plot(bray.ord.Cond, add = TRUE, cex=1.5)
plot(bray.ord.No3, add = TRUE, cex=1.5)
plot(bray.ord.CC, add = TRUE, cex=1.5)

#save as pdf


################
# Extract X, Y, Z coordinates: 

nmds_coordinates <- bray.ord$points
plot(bray.ord$points)
write.csv(nmds_coordinates, 'Fungal_Elwha_points.csv')

# Species accumulation curve
spp.curve <- specaccum(comm = main, method = "random", permutations = 1000) 
plot(spp.curve, ci.type="bar", xlab="litterbags (random)", ylab="# species") 
plot(spp.curve, ci.type="line", ci.lty=2 , xlab="litterbags (random)", ylab="# species") 
