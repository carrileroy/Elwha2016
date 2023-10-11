#' Carri J. LeRoy, PhD 4-21-23
#' 
#' BMI community data from 2016 Elwha study. Data in: BMI_main.csv and BMI_second.csv
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
if(!require(ecodist)){install.packages("ecodist")}

library("vegan")
library("permute")
library("MASS")
library("ggvegan")
library("ggplot2")
library("tidyverse")
library("indicspecies")
library("ecodist")
library("agricolae")

#' 
#' 
#' Set working directory under Session < Set working directory < choose directory
#' Load two matrices: BMI_main.csv and BMI_second.csv
#' 
#' Can't get the log(x+1) transformation to work in R, so doing it in Excel - all main matrices are already transformed: 
#' 
#' Important that R knows the row names are in column 1 and there are headers. 
#' Simplify your matrix names to "main" (taxa abundances) and "second" (environmental variables)
#' 
## ----loading matrices-----------------------------------------
main <- read.csv("BMI_main.csv", row.names = 1, header = TRUE)
second <- read.csv("BMI_second.csv", row.names = 1, header = TRUE)

#------or use the matrices already loaded--------------------
#main <- BMI_main
#second <- BMI_second

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


second$region = factor(second$region,
                  levels=unique(second$region))
second$location = factor(second$location,
                    levels=unique(second$location))
second$mills = factor(second$mills,
                       levels=unique(second$mills))
second$trib = factor(second$trib,
                         levels=unique(second$trib))
second$middle = factor(second$middle,
                       levels=unique(second$middle))
second$reservoir = factor(second$reservoir,
                         levels=unique(second$reservoir))
second$lower = factor(second$lower,
                       levels=unique(second$lower))
second$delta = factor(second$delta,
                         levels=unique(second$delta))

# Set your distance measure: "bray-curtis" "euclidean" etc.
dis <- distance(main, method = "bray-curtis")
# To determine the parsimonious axis number, we generate a scree plot. 
scree <- nmds(dis, mindim = 1, maxdim = 6, nits = 1)
#attributes(scree)
plot(scree$stress, xlab = "Number of Axes")


#' ## PERMANOVA using `adonis`
#' Analysis of variance using distance matrices and for fitting linear models to distance matrices
#' If you use permutations = 9999, you can get p-values down to 0.0001, if only = 999, only p-values down to 0.001
## ----adonis--------------------------------------------------------------
set.seed(1413)
adonis2(main ~ second$region, method = "bray", permutation = 9999)
#
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999
#
#adonis2(formula = main ~ second$region, permutations = 9999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#second$region  4   1.5314 0.58128 5.2059  1e-04 ***
#  Residual      15   1.1031 0.41872                  
#Total         19   2.6345 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

############ PermDisp? ####################
disp <- betadisper(dis, group = second$region)
disp

#Homogeneity of multivariate dispersions
#
#Call: betadisper(d = dis, group = second$region)
#
#No. of Positive Eigenvalues: 12
#No. of Negative Eigenvalues: 6
#
#Average distance to median:
#  Tributary      Middle Elwha Aldwell Reservoir       Lower Elwha             Delta 
#0.3165            0.1802            0.1152            0.1811            0.2433 
#
#Eigenvalues for PCoA axes:
#  (Showing 8 of 18 eigenvalues)
#PCoA1   PCoA2   PCoA3   PCoA4   PCoA5   PCoA6   PCoA7   PCoA8 
#1.15193 0.51704 0.26991 0.18869 0.17552 0.12004 0.10869 0.07247 

permutest(disp)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     4 0.092509 0.023127 1.4729    999  0.269
#Residuals 15 0.235531 0.015702  

anova(disp)
#Analysis of Variance Table
#
#Response: Distances
#Df   Sum Sq  Mean Sq F value Pr(>F)
#Groups     4 0.092509 0.023127  1.4729 0.2596
#Residuals 15 0.235531 0.015702 

plot(disp, hull=FALSE, ellipse=TRUE)

#'#####################
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
#Stress:     0.05473148 
#Stress type 1, weak ties
#Best solution was repeated 5 times in 20 tries
#The best solution was from try 0 (metric scaling or null solution)
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘main’ 
#' 
#' 
#' ## Basic ordination and plotting (using all defaults)
## --------------------------------------------------
plot(bray.ord, display = "sites")

################## This is working! Plotting axes 1 and 2 (changing this doesn't seem to change anything - "choices")
colvec <- c("orange1", "yellow2", "green1", "blue1", "purple1")
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
#***VECTORS
#
#NMDS1    NMDS2     r2  Pr(>r)  
#Sum -0.14140 -0.98995 0.3349 0.03397 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Rich <- envfit(bray.ord ~ Rich, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Rich
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)   
#Rich -0.51032 -0.85998 0.4024 0.01499 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.RichF <- envfit(bray.ord ~ RichF, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.RichF
#n.s.
    #***VECTORS
    #
    #NMDS1    NMDS2     r2 Pr(>r)
    #RichF  0.80198 -0.59736 0.1926 0.1648
    #Permutation: free
    #Number of permutations: 1000

bray.ord.SimpF <- envfit(bray.ord ~ SimpF, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.SimpF
#n.s.
      #.***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #SimpF -0.85446 -0.51952 0.0773 0.5145
      #Permutation: free
      #Number of permutations: 1000

bray.ord.ShanF <- envfit(bray.ord ~ ShanF, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.ShanF
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #ShanF -0.53193 -0.84679 0.0663 0.5465
      #Permutation: free
      #Number of permutations: 1000

bray.ord.k <- envfit(bray.ord ~ k, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.k
#n.s.
      #***VECTORS
      #
      #NMDS1     NMDS2     r2  Pr(>r)  
      #k -0.99736  0.07264 0.255 0.0959 .
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      #Permutation: free
      #Number of permutations: 1000

bray.ord.depth <- envfit(bray.ord ~ depth, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.depth
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #depth  0.93902 -0.34387 0.0781 0.4755
      #Permutation: free
      #Number of permutations: 1000

bray.ord.width <- envfit(bray.ord ~ width, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.width
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #width  0.49401 -0.86946 0.1454 0.3297
      #Permutation: free
      #Number of permutations: 1000

      #4 observations deleted due to missingness

bray.ord.CC <- envfit(bray.ord ~ CC, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.CC
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)   
#CC -0.81193  0.58375 0.5311 0.002997 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Fines <- envfit(bray.ord ~ Fines, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Fines
#***VECTORS
#
#NMDS1   NMDS2     r2   Pr(>r)    
#Fines 0.59221 0.80579 0.6633 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.D50 <- envfit(bray.ord ~ D50, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.D50
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#D50 -0.90906 -0.41667 0.753 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Temp <- envfit(bray.ord ~ Temp, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Temp
#***VECTORS
#
#NMDS1    NMDS2     r2  Pr(>r)  
#Temp -0.79612 -0.60513 0.3575 0.03397 *
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Cond <- envfit(bray.ord ~ Cond, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Cond
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#Cond -0.92245  0.38611 0.7429 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.p <- envfit(bray.ord ~ p, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.p
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2  Pr(>r)  
      #p -0.82194  0.56957 0.3357 0.07293 .
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
      #n -0.96946 -0.24524 0.2052 0.2128
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
      #po4 0.025965 0.999660 0.1578 0.2308
      #Permutation: free
      #Number of permutations: 1000

bray.ord.sio4 <- envfit(bray.ord ~ sio4, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.sio4
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #sio4 -0.20194  0.97940 0.3363 0.06094 .
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.No3 <- envfit(bray.ord ~ No3, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.No3
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#No3 -0.94845  0.31693 0.466 0.003996 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

#4 observations deleted due to missingness

bray.ord.no2 <- envfit(bray.ord ~ no2, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.no2
#***VECTORS
#
#NMDS1    NMDS2     r2  Pr(>r)  
#no2 -0.14616  0.98926 0.35 0.04496 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
      #nh4 -0.35082  0.93644 0.0613 0.5674
      #Permutation: free
      #Number of permutations: 1000

bray.ord.no2no3 <- envfit(bray.ord ~ no2no3, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.no2no3
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #no2No3 -0.82614  0.56346 0.1767 0.1658
      #Permutation: free
      #Number of permutations: 1000

bray.ord.chla <- envfit(bray.ord ~ chla, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.chla
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #chla -0.83650 -0.54796 0.2203 0.1349
      #Permutation: free
      #Number of permutations: 1000

bray.ord.omRock <- envfit(bray.ord ~ omRock, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.omRock
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #omRock  0.93361 -0.35828 0.2179 0.1768
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
      #imRock  0.63562 -0.77200 0.1307 0.4326
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.OM <- envfit(bray.ord ~ OM, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.OM
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)   
#OM -0.82517  0.56489 0.5242 0.002997 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.imSeston <- envfit(bray.ord ~ imSeston, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.imSeston
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2    r2 Pr(>r)
      #imSeston -0.66223 -0.74930 0.2086 0.1449
      #Permutation: free
      #Number of permutations: 1000

bray.ord.BI <- envfit(bray.ord ~ BI, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.BI
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #BI  0.23805 -0.97125 0.192 0.2607
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
      #SH  0.86397 -0.50354 0.1818 0.2857
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.rkm <- envfit(bray.ord ~ rkm, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.rkm
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#rkm -0.74186  0.67055 0.6339 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Elev <- envfit(bray.ord ~ Elev, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Elev
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#Elev -0.70352  0.71068 0.6703 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Sed1 <- envfit(bray.ord ~ Sed1, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed1
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed1 -0.32667 -0.94514 0.2563 0.1469
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed2 <- envfit(bray.ord ~ Sed2, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed2
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2 Pr(>r)
      #Sed2 0.78160 0.62378 0.1815 0.2647
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
      #Sed3 -0.90806 -0.41883 0.0106 0.9141
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.Sed4 <- envfit(bray.ord ~ Sed4, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Sed4
#n.s.
      #***VECTORS
      #
      #NMDS1    NMDS2     r2 Pr(>r)
      #Sed4 -0.88140 -0.47236 0.0114 0.9411
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
      #Sed5 0.68955 0.72424 0.0016  0.991
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
      #Sed6  0.93345 -0.35872 0.0186 0.8771
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
      #Sed7  0.90150 -0.43279 0.1629 0.2997
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
      #Sed8 -0.84046 -0.54187 0.0982 0.5205
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
      #Sed9 -0.50899  0.86078 0.0587 0.6893
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
      #Sed10 -0.39870  0.91708 0.1138 0.4356
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
      #Sed11 -0.88194 -0.47136 0.0135 0.9231
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
      #SedTot -0.87700  0.48049 0.0415 0.7592
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
      #P1  0.01153 -0.99993 0.3322 0.07692 .
      #---
      #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P2 <- envfit(bray.ord ~ P2, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P2
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2 Pr(>r)
      #P2 0.98746 0.15788 0.2238 0.1968
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
      #P3  0.49727 -0.86760 0.0149 0.9021
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
      #P4  0.31074 -0.95049 0.0087 0.9381
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
      #P5 0.88791 0.46002 0.1362 0.3926
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
      #P6 0.94371 0.33078 0.2031 0.2258
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P7 <- envfit(bray.ord ~ P7, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P7
#n.s.
      #***VECTORS
      #
      #NMDS1   NMDS2     r2 Pr(>r)
      #P7 0.95853 0.28500 0.2419 0.1628
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.P8 <- envfit(bray.ord ~ P8, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.P8
#n.s.
      #***VECTORS
      #
      #NMDS1     NMDS2     r2 Pr(>r)
      #P8  0.9999800 0.0070179 0.0622 0.6494
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
      #P9 -0.33540 -0.94208 0.021 0.8861
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
      #P10 -0.83140  0.55568 0.1684 0.2987
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
      #P11 -0.16063 -0.98701 0.0104 0.9281
      #Permutation: free
      #Number of permutations: 1000
      #
      #4 observations deleted due to missingness

bray.ord.EvenM <- envfit(bray.ord ~ EvenM, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.EvenM
#***VECTORS
#
#NMDS1     NMDS2    r2   Pr(>r)    
#Even -0.058189  0.998310 0.6109 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.Div <- envfit(bray.ord ~ Div, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.Div
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)    
#Div -0.20906  0.97790 0.5462 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000

bray.ord.SimpsonM <- envfit(bray.ord ~ SimpsonM, data = second, permutations = 1000, na.rm = TRUE)
bray.ord.SimpsonM
#***VECTORS
#
#NMDS1    NMDS2     r2   Pr(>r)   
#SimpsonM -0.14641  0.98922 0.5464 0.000999 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 1000


# Add ALL significant vectors to plot
colvec <- c("orange2", "gold", "green1", "blue1", "purple1")
shpvec <- c(16, 15, 16, 15, 16, 15)
plot(bray.ord, type = "n", choices=c(1,2), cex.lab=1.5, cex.axis=1.3)
with(second, points(bray.ord, display = "sites", cex=1.5, col = colvec[region], pch = shpvec[region], bg = colvec[region]))
with(second, legend("topright", legend = levels(region), bty = "n", col = colvec, pch = shpvec, pt.bg = colvec))
ordihull(bray.ord, groups = second$region, label = FALSE, col = colvec)
plot(bray.ord.Sum, add = TRUE, cex=1.5)
plot(bray.ord.Rich, add = TRUE, cex=1.5)
plot(bray.ord.Div, add = TRUE, cex=1.5)
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

# Species accumulation curve
spp.curve <- specaccum(comm = main, method = "random", permutations = 1000) 
plot(spp.curve, ci.type="bar", xlab="litterbags (random)", ylab="# species") 
plot(spp.curve, ci.type="line", ci.lty=2 , xlab="litterbags (random)", ylab="# species") 

# Indicator Species Analysis:
if(!require(labdsv)){install.packages("labdsv")}
library(labdsv) 
ind.spp <- indval(main, second$region)
ind.spp
summary(ind.spp)

#cluster indicator_value probability
#BAETIS            2          0.3363       0.003
#OLIGOCH           4          0.4112       0.003
#ZAPADA            4          0.3992       0.025
#CHIRONO           4          0.2936       0.005
#RADIX.PHYSA       5          1.0000       0.003
