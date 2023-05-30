# Carri J. LeRoy, 1-24-23
#
# 2016 Elwha paper: Permutative Three-way ANCOVAs for decomp rate comparisons
# data in Elwha_AFDM.csv
#
# install required packages for blocked ANOVA
if(!require(agricolae)){install.packages("agricolae")}
if(!require(lmPerm)){install.packages("lmPerm")}
if(!require(psych)){install.packages("psych")}
if(!require(FSA)){install.packages("FSA")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(car)){install.packages("car")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(tidyr)){install.packages("tidyr")}
if(!require(dplyr)){install.packages("dplyr")}
if(!require(magrittr)){install.packages("magrittr")}
if(!require(gridExtra)){install.packages("gridExtra")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(readxl)){install.packages("readxl")}
if(!require(devtools)) install.packages("devtools")
if(!require(ggpubr)) install.packages("ggpubr")
#
# Call the packages
library(lmPerm)           # permutative analyses
library(psych)            # stats for psych research (useful)
library(FSA)              # Fisheries stock assessment methods
library(multcompView)     # Visualizations of paired comparisons
library(lsmeans)          # least-squared means
library(tidyr)          	# data re-shaping
library(ggplot2)        	# plotting & data
library(dplyr)          	# data manipulation
library(magrittr)       	# pipe operator
library(gridExtra)     	  # provides side-by-side plotting
library(car)     		      # companion to applied regression
library(tidyverse)        # tidyverse
library(readxl)           # reads Excel files into R
library(devtools)
library(ggpubr)
library(agricolae)
#
# Calculating 3-way ANCOVA for source*location*days:
#  Order factors by the order in data frame
#  Otherwise, R will alphabetize them
#
# Rename data set to "AFDMall"
AFDMall <- Elwha_AFDM
#
AFDMall$source = factor(AFDMall$source,
                      levels=unique(AFDMall$source))
AFDMall$location = factor(AFDMall$location,
                        levels=unique(AFDMall$location))
AFDMall$region = factor(AFDMall$region,
                   levels=unique(AFDMall$region))
#
#Check data structure: R needs to treat your grouping variables as "factors" not "numbers" or "chr" characters
str(AFDMall)
# Set a random seed
set.seed(1431)
#Run a permutative 3-way ANCOVA: with location
#
fit1 <- aovp(lnpAFDMR ~ source*location*days, data = AFDMall)
anova(fit1)
#
#Analysis of Variance Table
#
#Response: lnpAFDMR
#                         Df R Sum Sq R Mean Sq Iter Pr(Prob)    
#   source                 4    0.591     0.148 1301   0.5596    
#   location              19   72.625     3.822 5000   <2e-16 ***
#   source:location       76   20.527     0.270 5000   0.8582    
#   days                   1   99.155    99.155 5000   <2e-16 ***
#   source:days            4    0.442     0.111  353   0.9292    
#   location:days         19   67.717     3.564 5000   <2e-16 ***
#   source:location:days  76   20.848     0.274 5000   0.7714    
#   Residuals            194   65.807     0.339                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Run a permutative 3-way ANCOVA: with region
fit2 <- aovp(lnpAFDMR ~ source*region*days, data = AFDMall)
anova(fit2)
#
#Analysis of Variance Table
#
#Response: lnpAFDMR
#                       Df R Sum Sq R Mean Sq Iter Pr(Prob)    
#   source               4    0.509     0.127   51   1.0000    
#   region               4   20.488     5.122 5000   <2e-16 ***
#   source:region       16    5.443     0.340 1270   0.9283    
#   days                 1  107.922   107.922 5000   <2e-16 ***
#   source:days          4    0.395     0.099   81   1.0000    
#   region:days          4   19.629     4.907 5000   <2e-16 ***
#   source:region:days  16    5.366     0.335 1276   0.9718    
#   Residuals          344  200.316     0.582                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# If litter sources are not significantly different, run 2-way ANCOVAs instead:
#Run a permutative 2-way ANCOVA: with location
#
fit3 <- aovp(lnpAFDMR ~ location*days, data = AFDMall)
anova(fit3)
#
#Analysis of Variance Table
#
#Response: lnpAFDMR
#                   Df R Sum Sq R Mean Sq Iter  Pr(Prob)    
#   location       19   73.622     3.875 5000 < 2.2e-16 ***
#   days            1   99.974    99.974 5000 < 2.2e-16 ***
#   location:days  19   68.359     3.598 5000 < 2.2e-16 ***
#   Residuals     354  109.079     0.308                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#

#Run a permutative 2-way ANCOVA: with region
fit4 <- aovp(lnpAFDMR ~ region*days, data = AFDMall)
anova(fit4)
#
#Analysis of Variance Table
#
#Response: lnpAFDMR
#                 Df R Sum Sq R Mean Sq Iter  Pr(Prob)    F-ratio
#   region        4   20.489     5.122 5000 < 2.2e-16 *** 9.26
#   days          1  107.883   107.883 5000 < 2.2e-16 *** 195.08
#   region:days   4   19.630     4.908 5000 < 2.2e-16 *** 8.87
#   Residuals   384  212.272     0.553                   
#---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# OVerall F-ratio (9,384) = 213.22
#
# Tukey test is not working:
HSD.test(fit4, "region", group=TRUE)
out1 <- HSD.test(fit4, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#Study: fit4 ~ "region"
#
#HSD Test for lnpAFDMR 
#
#Mean Square Error:  0.5527922 
#
#region,  means
#
#lnpAFDMR       std  r      Min      Max
#Aldwell Reservoir 3.823752 1.2167013 75 0.000000 4.615121
#Elwha Middle      3.630329 1.2373331 80 0.000000 4.615121
#Estuary/Delta     4.312333 0.3050723 80 2.175484 4.615121
#Lower Elwha       4.031195 0.6677592 80 0.000000 4.615121
#Tributary         3.942140 0.9253570 79 0.000000 4.627886
#
#Alpha: 0.05 ; DF Error: 384 
#Critical Value of Studentized Range: 3.876054 
#
#Groups according to probability of means differences and alpha level( 0.05 )
#
#Treatments with the same letter are not significantly different.
#
#lnpAFDMR groups
#Estuary/Delta     4.312333      a
#Lower Elwha       4.031195     ab
#Tributary         3.942140     bc
#Aldwell Reservoir 3.823752     bc
#Elwha Middle      3.630329      c
#
fit5 <- aovp(lnpAFDMR ~ source*days, data = AFDMall)
anova(fit5)

#Analysis of Variance Table
#
#Response: lnpAFDMR
#                   Df R Sum Sq R Mean Sq Iter Pr(Prob)   F-ratio 
#     source        4    0.629     0.157   51        1    0.2404
#     days          1  108.396   108.396 5000   <2e-16 ***165.99
#     source:days   4    0.509     0.127   55        1    0.1944
#     Residuals   384  250.734     0.653                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



