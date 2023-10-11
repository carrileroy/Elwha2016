# Carri J. LeRoy, 5-2-23
#
# # For 2016 Elwha decomp paper: k-rate pANOVAs; and regressions with k and all other variables
# Datasets = kAll.csv
# Check to see if you need to install the lmPerm package and other packages needed for analysis and plotting
#
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
if(!require(devtools)){install.packages("devtools")}
if(!require(ggpubr)){install.packages("ggpubr")}
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
#
kAll <- read.csv("kAll.csv", header = TRUE)
#
# set region as a factor:
kAll$region = factor(kAll$region,
                     levels=unique(kAll$region))
#
# Set a random seed
set.seed(1431)
#
# Compare k-rates among regions:
res.aov3 <- aovp(k~region, data = kAll, perm="Prob")
anova(res.aov3)
#
#Analysis of Variance Table
#
#Response: k
#Df  R Sum Sq  R Mean Sq Iter Pr(Prob)
#region     4 0.0030303 0.00075758 1732   0.2552
#Residuals 15 0.0083998 0.00055998   
#
#


# Regressions between k and all other variables: use dataset kAll.csv 
kAll$region = factor(kAll$region,
                     levels=unique(kAll$region))
#
slr1 = lmp(depth~k, data=kAll)
summary(slr1)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    64.22   92    0.522
#
#Residual standard error: 7.635 on 18 degrees of freedom
#Multiple R-Squared: 0.04299,	Adjusted R-squared: -0.01018 
#F-statistic: 0.8085 on 1 and 18 DF,  p-value: 0.3804 
#
slr2 = lmp(elev~k, data=kAll)
summary(slr2)
#
#Coefficients:
#  Estimate Iter Pr(Prob)  
#k    675.1  997   0.0913 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 39.03 on 18 degrees of freedom
#Multiple R-Squared: 0.1596,	Adjusted R-squared: 0.1129 
#F-statistic: 3.419 on 1 and 18 DF,  p-value: 0.08094 
#
slr3 = lmp(chla~k, data=kAll)
summary(slr3)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -0.841   51     0.98
#
#Residual standard error: 0.6771 on 18 degrees of freedom
#Multiple R-Squared: 0.0009786,	Adjusted R-squared: -0.05452 
#F-statistic: 0.01763 on 1 and 18 DF,  p-value: 0.8958 
#
slr4 = lmp(omRock~k, data=kAll)
summary(slr4)
#
#Coefficients:
#  Estimate Iter Pr(Prob)  
#k   -1.754 1142   0.0806 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.0956 on 14 degrees of freedom
#Multiple R-Squared: 0.1848,	Adjusted R-squared: 0.1265 
#F-statistic: 3.173 on 1 and 14 DF,  p-value: 0.09656 

slr5 = lmp(imRock~k, data=kAll)
summary(slr5)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -6.226  389    0.206
#
#Residual standard error: 0.4253 on 14 degrees of freedom
#Multiple R-Squared: 0.1261,	Adjusted R-squared: 0.06368 
#F-statistic:  2.02 on 1 and 14 DF,  p-value: 0.1771 
#
slr6 = lmp(omSeston~k, data=kAll)
summary(slr6)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k -0.005151  136    0.426
#
#Residual standard error: 0.0005668 on 18 degrees of freedom
#Multiple R-Squared: 0.04983,	Adjusted R-squared: -0.002958 
#F-statistic: 0.944 on 1 and 18 DF,  p-value: 0.3441 
#
slr7 = lmp(imSeston~k, data=kAll)
summary(slr7)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k  0.03746   51    0.686
#
#Residual standard error: 0.005898 on 18 degrees of freedom
#Multiple R-Squared: 0.02497,	Adjusted R-squared: -0.0292 
#F-statistic: 0.461 on 1 and 18 DF,  p-value: 0.5058
#
slr8 = lmp(temp~k, data=kAll)
summary(slr8)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -1.035   51    0.882
#
#Residual standard error: 1.226 on 18 degrees of freedom
#Multiple R-Squared: 0.0004522,	Adjusted R-squared: -0.05508 
#F-statistic: 0.008143 on 1 and 18 DF,  p-value: 0.9291

slr9 = lmp(cond~k, data=kAll)
summary(slr9)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    623.3  248     0.29
#
#Residual standard error: 58.81 on 18 degrees of freedom
#Multiple R-Squared: 0.06659,	Adjusted R-squared: 0.01473 
#F-statistic: 1.284 on 1 and 18 DF,  p-value: 0.272 
#
slr10 = lmp(width~k, data=kAll)
summary(slr10)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -74.51   51    0.882
#
#Residual standard error: 25.28 on 14 degrees of freedom
#Multiple R-Squared: 0.005818,	Adjusted R-squared: -0.0652 
#F-statistic: 0.08192 on 1 and 14 DF,  p-value: 0.7789 
#
slr11 = lmp(p~k, data=kAll)
summary(slr11)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -89.13  138     0.42
#
#Residual standard error: 11.45 on 14 degrees of freedom
#Multiple R-Squared: 0.03922,	Adjusted R-squared: -0.02941 
#F-statistic: 0.5715 on 1 and 14 DF,  p-value: 0.4622
#
slr12 = lmp(n~k, data=kAll)
summary(slr12)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    -2334  569    0.151
#
#Residual standard error: 177.9 on 14 degrees of freedom
#Multiple R-Squared: 0.1039,	Adjusted R-squared: 0.03991 
#F-statistic: 1.624 on 1 and 14 DF,  p-value: 0.2233

slr13 = lmp(po4~k, data=kAll)
summary(slr13)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -6.838   53     0.66
#
#Residual standard error: 3.956 on 18 degrees of freedom
#Multiple R-Squared: 0.001894,	Adjusted R-squared: -0.05356 
#F-statistic: 0.03415 on 1 and 18 DF,  p-value: 0.8555 
#
slr14 = lmp(sio4~k, data=kAll)
summary(slr14)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -502.6   51     0.98
#
#Residual standard error: 2184 on 14 degrees of freedom
#Multiple R-Squared: 3.566e-05,	Adjusted R-squared: -0.07139 
#F-statistic: 0.0004992 on 1 and 14 DF,  p-value: 0.9825 
#
slr15 = lmp(no3~k, data=kAll)
summary(slr15)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k      110   51    0.941
#
#Residual standard error: 32 on 14 degrees of freedom
#Multiple R-Squared: 0.007902,	Adjusted R-squared: -0.06296 
#F-statistic: 0.1115 on 1 and 14 DF,  p-value: 0.7434
#
slr16 = lmp(no2~k, data=kAll)
summary(slr16)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k  -0.6822   51    0.863
#
#Residual standard error: 0.1744 on 14 degrees of freedom
#Multiple R-Squared: 0.0102,	Adjusted R-squared: -0.0605 
#F-statistic: 0.1443 on 1 and 14 DF,  p-value: 0.7097 

slr17 = lmp(nh4~k, data=kAll)
summary(slr17)
#
#Coefficients:
#  Estimate Iter Pr(Prob)  
#k   -608.9 1959    0.049 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 33.73 on 18 degrees of freedom
#Multiple R-Squared: 0.1715,	Adjusted R-squared: 0.1254 
#F-statistic: 3.725 on 1 and 18 DF,  p-value: 0.06952 
#
slr18 = lmp(no2no3~k, data=kAll)
summary(slr18)
#
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    135.1   81    0.556
#
#Residual standard error: 28.84 on 18 degrees of freedom
#Multiple R-Squared: 0.01374,	Adjusted R-squared: -0.04105 
#F-statistic: 0.2508 on 1 and 18 DF,  p-value: 0.6226 
#
slr19 = lmp(Rhytisma~k, data=kAll)
summary(slr19)
#
#Coefficients:
#  Estimate Iter Pr(Prob)  
#k    -76.3 5000   0.0114 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 2.963 on 18 degrees of freedom
#Multiple R-Squared: 0.2964,	Adjusted R-squared: 0.2573 
#F-statistic: 7.581 on 1 and 18 DF,  p-value: 0.01308

slr20 = lmp(fines~k, data=kAll)
summary(slr20)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -182.6  347    0.225
#
#Residual standard error: 14.75 on 18 degrees of freedom
#Multiple R-Squared: 0.08863,	Adjusted R-squared: 0.038 
#F-statistic: 1.751 on 1 and 18 DF,  p-value: 0.2024 

slr21 = lmp(D50~k, data=kAll)
summary(slr21)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    338.8  260    0.281
#
#Residual standard error: 33 on 18 degrees of freedom
#Multiple R-Squared: 0.06273,	Adjusted R-squared: 0.01066 
#F-statistic: 1.205 on 1 and 18 DF,  p-value: 0.2868 

slr22 = lmp(Biden~k, data=kAll)
summary(slr22)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k     2585   51    0.922
#
#Residual standard error: 2106 on 14 degrees of freedom
#Multiple R-Squared: 0.001013,	Adjusted R-squared: -0.07034 
#F-statistic: 0.0142 on 1 and 14 DF,  p-value: 0.9068 

slr23 = lmp(shred~k, data=kAll)
summary(slr23)
#Coefficients:
# Estimate Iter Pr(Prob)  
#k   -47.09 1499   0.0627 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 2.548 on 14 degrees of freedom
#Multiple R-Squared: 0.187,	Adjusted R-squared: 0.1289 
#F-statistic:  3.22 on 1 and 14 DF,  p-value: 0.09436

slr24 = lmp(CC~k, data=kAll)
summary(slr24)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    20.47   51        1
#
#Residual standard error: 36.12 on 14 degrees of freedom
#Multiple R-Squared: 0.0002162,	Adjusted R-squared: -0.0712 
#F-statistic: 0.003028 on 1 and 14 DF,  p-value: 0.9569

slr25 = lmp(rkm~k, data=kAll)
summary(slr25)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    101.4  804    0.111
#
#Residual standard error: 6.428 on 18 degrees of freedom
#Multiple R-Squared: 0.1365,	Adjusted R-squared: 0.08854 
#F-statistic: 2.846 on 1 and 18 DF,  p-value: 0.1089 

slr26 = lmp(Sed1~k, data=kAll)
summary(slr26)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   0.1033   51    0.941
#
#Residual standard error: 0.1586 on 14 degrees of freedom
#Multiple R-Squared: 0.0002853,	Adjusted R-squared: -0.07112 
#F-statistic: 0.003995 on 1 and 14 DF,  p-value: 0.9505

slr27 = lmp(Sed2~k, data=kAll)
summary(slr27)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   0.3269   51    0.824
#
#Residual standard error: 0.348 on 14 degrees of freedom
#Multiple R-Squared: 0.0005938,	Adjusted R-squared: -0.07079 
#F-statistic: 0.008319 on 1 and 14 DF,  p-value: 0.9286 

slr28 = lmp(Sed3~k, data=kAll)
summary(slr28)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -0.119   51     0.98
#
#Residual standard error: 0.5515 on 14 degrees of freedom
#Multiple R-Squared: 3.135e-05,	Adjusted R-squared: -0.07139 
#F-statistic: 0.000439 on 1 and 14 DF,  p-value: 0.9836 

slr29 = lmp(Sed4~k, data=kAll)
summary(slr29)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   0.6766   51        1
#
#Residual standard error: 0.4286 on 14 degrees of freedom
#Multiple R-Squared: 0.001675,	Adjusted R-squared: -0.06963 
#F-statistic: 0.0235 on 1 and 14 DF,  p-value: 0.8804 

slr30 = lmp(Sed5~k, data=kAll)
summary(slr30)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k     10.5  102      0.5
#
#Residual standard error: 1.302 on 14 degrees of freedom
#Multiple R-Squared: 0.04198,	Adjusted R-squared: -0.02645 
#F-statistic: 0.6135 on 1 and 14 DF,  p-value: 0.4465 

slr31 = lmp(Sed6~k, data=kAll)
summary(slr31)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    4.317  243    0.292
#
#Residual standard error: 0.5956 on 14 degrees of freedom
#Multiple R-Squared: 0.03417,	Adjusted R-squared: -0.03482 
#F-statistic: 0.4953 on 1 and 14 DF,  p-value: 0.4931 

slr32 = lmp(Sed7~k, data=kAll)
summary(slr32)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    1.876   51    0.863
#
#Residual standard error: 1.395 on 14 degrees of freedom
#Multiple R-Squared: 0.001217,	Adjusted R-squared: -0.07012 
#F-statistic: 0.01706 on 1 and 14 DF,  p-value: 0.8979 

slr33 = lmp(Sed8~k, data=kAll)
summary(slr33)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    38.53  183    0.355
#
#Residual standard error: 4.29 on 14 degrees of freedom
#Multiple R-Squared: 0.05152,	Adjusted R-squared: -0.01623 
#F-statistic: 0.7604 on 1 and 14 DF,  p-value: 0.3979 

slr34 = lmp(Sed9~k, data=kAll)
summary(slr34)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    140.5   68    0.603
#
#Residual standard error: 13.55 on 14 degrees of freedom
#Multiple R-Squared: 0.06751,	Adjusted R-squared: 0.0009057 
#F-statistic: 1.014 on 1 and 14 DF,  p-value: 0.3311 

slr35 = lmp(Sed10~k, data=kAll)
summary(slr35)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -402.4  138     0.42
#
#Residual standard error: 55.04 on 14 degrees of freedom
#Multiple R-Squared: 0.03476,	Adjusted R-squared: -0.03419 
#F-statistic: 0.5041 on 1 and 14 DF,  p-value: 0.4894 

slr36 = lmp(Sed11~k, data=kAll)
summary(slr36)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -818.2  205    0.332
#
#Residual standard error: 84.54 on 14 degrees of freedom
#Multiple R-Squared: 0.05933,	Adjusted R-squared: -0.007861 
#F-statistic: 0.883 on 1 and 14 DF,  p-value: 0.3633 

slr37 = lmp(SedTot~k, data=kAll)
summary(slr37)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    -1024  161    0.385
#
#Residual standard error: 123.5 on 14 degrees of freedom
#Multiple R-Squared: 0.04422,	Adjusted R-squared: -0.02405 
#F-statistic: 0.6477 on 1 and 14 DF,  p-value: 0.4344 

slr38 = lmp(P1~k, data=kAll)
summary(slr38)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k 0.008822   51    0.745
#
#Residual standard error: 0.001033 on 14 degrees of freedom
#Multiple R-Squared: 0.0468,	Adjusted R-squared: -0.02129 
#F-statistic: 0.6873 on 1 and 14 DF,  p-value: 0.421 

slr39 = lmp(P2~k, data=kAll)
summary(slr39)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k 0.000907   51        1
#
#Residual standard error: 0.003096 on 14 degrees of freedom
#Multiple R-Squared: 5.781e-05,	Adjusted R-squared: -0.07137 
#F-statistic: 0.0008093 on 1 and 14 DF,  p-value: 0.9777 

slr40 = lmp(P3~k, data=kAll)
summary(slr40)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   0.0472  111    0.477
#
#Residual standard error: 0.006703 on 14 degrees of freedom
#Multiple R-Squared: 0.03231,	Adjusted R-squared: -0.03681 
#F-statistic: 0.4675 on 1 and 14 DF,  p-value: 0.5053 

slr41 = lmp(P4~k, data=kAll)
summary(slr41)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k  0.05318  254    0.283
#
#Residual standard error: 0.005005 on 14 degrees of freedom
#Multiple R-Squared: 0.07065,	Adjusted R-squared: 0.004266 
#F-statistic: 1.064 on 1 and 14 DF,  p-value: 0.3197 

slr42 = lmp(P5~k, data=kAll)
summary(slr42)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   0.1373  180    0.361
#
#Residual standard error: 0.01606 on 14 degrees of freedom
#Multiple R-Squared: 0.04694,	Adjusted R-squared: -0.02113 
#F-statistic: 0.6896 on 1 and 14 DF,  p-value: 0.4203 

slr43 = lmp(P6~k, data=kAll)
summary(slr43)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k  0.04663   51    0.765
#
#Residual standard error: 0.008753 on 14 degrees of freedom
#Multiple R-Squared: 0.01875,	Adjusted R-squared: -0.05133 
#F-statistic: 0.2676 on 1 and 14 DF,  p-value: 0.613 

slr44 = lmp(P7~k, data=kAll)
summary(slr44)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k  -0.0848   51    0.961
#
#Residual standard error: 0.04609 on 14 degrees of freedom
#Multiple R-Squared: 0.002274,	Adjusted R-squared: -0.06899 
#F-statistic: 0.03191 on 1 and 14 DF,  p-value: 0.8608

slr45 = lmp(P8~k, data=kAll)
summary(slr45)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k  0.05343   51     0.98
#
#Residual standard error: 0.1384 on 14 degrees of freedom
#Multiple R-Squared: 0.0001004,	Adjusted R-squared: -0.07132 
#F-statistic: 0.001405 on 1 and 14 DF,  p-value: 0.9706

slr46 = lmp(P9~k, data=kAll)
summary(slr46)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    2.506  510    0.165
#
#Residual standard error: 0.152 on 14 degrees of freedom
#Multiple R-Squared: 0.1548,	Adjusted R-squared: 0.09438 
#F-statistic: 2.563 on 1 and 14 DF,  p-value: 0.1317 

slr47 = lmp(P10~k, data=kAll)
summary(slr47)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k  -0.3608   51    0.863
#
#Residual standard error: 0.2007 on 14 degrees of freedom
#Multiple R-Squared: 0.002171,	Adjusted R-squared: -0.0691 
#F-statistic: 0.03045 on 1 and 14 DF,  p-value: 0.864 

slr48 = lmp(P11~k, data=kAll)
summary(slr48)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -2.407  131    0.435
#
#Residual standard error: 0.2328 on 14 degrees of freedom
#Multiple R-Squared: 0.0672,	Adjusted R-squared: 0.0005686 
#F-statistic: 1.009 on 1 and 14 DF,  p-value: 0.3323 

slr49 = lmp(EvenM~k, data=kAll)
summary(slr49)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    2.051  202    0.332
#
#Residual standard error: 0.2266 on 18 degrees of freedom
#Multiple R-Squared: 0.04943,	Adjusted R-squared: -0.003378 
#F-statistic: 0.936 on 1 and 18 DF,  p-value: 0.3461 

slr50 = lmp(DivM~k, data=kAll)
summary(slr50)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    6.367  315    0.241
#
#Residual standard error: 0.4967 on 18 degrees of freedom
#Multiple R-Squared: 0.09449,	Adjusted R-squared: 0.04418 
#F-statistic: 1.878 on 1 and 18 DF,  p-value: 0.1874 

slr51 = lmp(SimpM~k, data=kAll)
summary(slr51)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    3.048  659    0.132
#
#Residual standard error: 0.2272 on 18 degrees of freedom
#Multiple R-Squared: 0.1025,	Adjusted R-squared: 0.05268 
#F-statistic: 2.056 on 1 and 18 DF,  p-value: 0.1687 

slr52 = lmp(SumM~k, data=kAll)
summary(slr52)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -27668   51    0.863
#
#Residual standard error: 4495 on 18 degrees of freedom
#Multiple R-Squared: 0.0235,	Adjusted R-squared: -0.03075 
#F-statistic: 0.4331 on 1 and 18 DF,  p-value: 0.5188 

slr53 = lmp(RichM~k, data=kAll)
summary(slr53)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    7.679   51    0.706
#
#Residual standard error: 2.88 on 18 degrees of freedom
#Multiple R-Squared: 0.004493,	Adjusted R-squared: -0.05081 
#F-statistic: 0.08124 on 1 and 18 DF,  p-value: 0.7789 

slr54 = lmp(RichF~k, data=kAll)
summary(slr54)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -371.5  360    0.219
#
#Residual standard error: 27.94 on 18 degrees of freedom
#Multiple R-Squared: 0.1009,	Adjusted R-squared: 0.05097 
#F-statistic:  2.02 on 1 and 18 DF,  p-value: 0.1723 

slr55 = lmp(SimpF~k, data=kAll)
summary(slr55)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k   -4.485   51    0.882
#
#Residual standard error: 1.986 on 18 degrees of freedom
#Multiple R-Squared: 0.00323,	Adjusted R-squared: -0.05215 
#F-statistic: 0.05832 on 1 and 18 DF,  p-value: 0.8119 

slr56 = lmp(ShanF~k, data=kAll)
summary(slr56)
#Coefficients:
#  Estimate Iter Pr(Prob)
#k    -4.37   51    0.882
#
#Residual standard error: 4.5 on 18 degrees of freedom
#Multiple R-Squared: 0.0005987,	Adjusted R-squared: -0.05492 
#F-statistic: 0.01078 on 1 and 18 DF,  p-value: 0.9184 
