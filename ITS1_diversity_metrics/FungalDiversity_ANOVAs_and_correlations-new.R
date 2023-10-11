# Carri J. LeRoy, 4-26-23
#
# Permutative ANOVAs + linear regressions with fungal diversity, k and environmental variables
#
# For 2016 Elwha paper: 
# Dataset = FungalDiversity.csv
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
if(!require(plyr)){install.packages("plyr")}
if(!require(dplyr)){install.packages("dplyr")}
if(!require(magrittr)){install.packages("magrittr")}
if(!require(gridExtra)){install.packages("gridExtra")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(readxl)){install.packages("readxl")}
if(!require(devtools)){install.packages("devtools")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(grid)){install.packages("grid")}

#
# Call the packages
library(lmPerm)           # permutative analyses
library(psych)            # stats for psych research (useful)
library(FSA)              # Fisheries stock assessment methods
library(multcompView)     # Visualizations of paired comparisons
library(lsmeans)          # least-squared means
library(tidyr)          	# data re-shaping
library(ggplot2)        	# plotting & data
library(plyr)
library(dplyr)          	# data manipulation
library(magrittr)       	# pipe operator
library(gridExtra)     	  # provides side-by-side plotting
library(car)     		      # companion to applied regression
library(tidyverse)        # tidyverse
library(readxl)           # reads Excel files into R
library(devtools)
library(ggpubr)
library(agricolae)
library(grid)

#
# set WD to folder
# then read the data in: 
fungal <- read.csv("FungalDiversity.csv", header = TRUE)
#
# Run 2-way ANOVAs for H, D, and S against harvest*region:
fungal$harvest = factor(fungal$harvest,
                     levels=unique(fungal$harvest))
fungal$region = factor(fungal$region,
                     levels=unique(fungal$region))
#
# Set a random seed
set.seed(1431)
#
res.aov1 <- aovp(H.est~region*harvest, data = fungal, perm="Prob")
anova(res.aov1)
#
#Analysis of Variance Table
#
#Response: H.est
#                     Df R Sum Sq R Mean Sq Iter Pr(Prob)  
#     region          5   153.13   30.6267 4677  0.0244 * F = 2.72
#     harvest         1     3.04    3.0387   69  0.7843   F = 0.27
#     region:harvest  5    40.34    8.0672  431  0.6186   F = 0.72
#     Residuals      36   404.65   11.2403                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
res.aov2 <- aovp(D.est~region*harvest, data = fungal, perm="Prob")
anova(res.aov2)
#Analysis of Variance Table
#
#Response: D.est
#Df R Sum Sq R Mean Sq Iter Pr(Prob)  
#region          5   47.686    9.5372 5000  0.0172 *
#  harvest         1    0.108    0.1085   51  1.000  
#region:harvest  5    7.394    1.4789  162  0.8788  
#Residuals      36  111.714    3.1032                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
res.aov3 <- aovp(S.est~region*harvest, data = fungal, perm="Prob")
anova(res.aov3)
#
#Analysis of Variance Table
#
#Response: S.est
#               Df R Sum Sq R Mean Sq Iter Pr(Prob)  
#region          5   3774.9    754.98 2001  0.30007  F = 1.44
#harvest         1   1546.6   1546.55 1790  0.07078  F = 2.95
#region:harvest  5   3030.1    606.02  595  0.3383   F = 1.16
#Residuals      36  18853.2    523.70                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Follow-up one-way permutative ANOVAs (merging harvests) to see which regions are different
#
fit1 <- aovp(H.est~region, data=fungal, perm="Prob")
anova(fit1)
#
#Analysis of Variance Table
#
#Response: H.est
#Df R Sum Sq R Mean Sq Iter Pr(Prob)  
#     region     5   153.13    30.627 5000   0.0192 * F = 2.8712
#     Residuals 42   448.02    10.667  

HSD.test(fit1, "region", group=TRUE)
out1 <- HSD.test(fit1, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#H.est groups
#Lower Elwha       6.921870      a
#Tributary         6.045359     ab
#Elwha Middle      5.418547     ab
#Aldwell Reservoir 3.321298     ab
#Delta             2.894550     ab
#Mills Reservoir   2.030153      b

fit2 <- aovp(D.est~region, data=fungal, perm="Prob")
anova(fit2)
#
#Analysis of Variance Table
#
#Response: D.est
#Df R Sum Sq R Mean Sq Iter Pr(Prob)  
#region     5   47.686    9.5372 5000   0.01 *
#  Residuals 42  119.217    2.8385                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

HSD.test(fit2, "region", group=TRUE)
out2 <- HSD.test(fit2, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#D.est groups
#Elwha Middle      3.841802      a
#Tributary         3.768699      a
#Lower Elwha       3.713548      a
#Aldwell Reservoir 2.415927      a
#Delta             1.865522      a
#Mills Reservoir   1.355692      a

fit3 <- aovp(S.est~region, data=fungal, perm="Prob")
anova(fit3)
#
#Analysis of Variance Table
#
#Response: S.est
#Df R Sum Sq R Mean Sq Iter Pr(Prob)
#region     5   3774.9    754.98  2233   0.166  F= 1.35
#Residuals 42  23429.9    557.85 

HSD.test(fit3, "region", group=TRUE)
out3 <- HSD.test(fit3, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#S.est groups
#Delta             83.73768      a
#Aldwell Reservoir 74.73649      a
#Mills Reservoir   68.89962      a
#Lower Elwha       63.66545      a
#Tributary         60.93789      a
#Elwha Middle      57.52549      a


# 2 bar charts - added following review: 
#Plot 1:
data_summary <- function(fungal, S.est, region){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]]))
  }
  data_sum<-ddply(fungal, region, .fun=summary_func,
                  S.est)
  data_sum <- plyr::rename(data_sum, c("mean" = S.est))
  return(data_sum)
}

fun1 <- data_summary(fungal, S.est="S.est", 
                     region=c("region"))
# Convert region to a factor variable
fun1$region=as.factor(fun1$region)
head(fun1)

plot1 <- p<- ggplot(fun1, aes(x=region, y=S.est, fill=region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=S.est-se, ymax=S.est+se), width=.2,
                position=position_dodge(.9)) +
  labs(x = element_blank(), y = "OTU Richness", tag = "A")+
  theme_classic() +
  scale_x_discrete(limits = c("Mills Reservoir","Tributary","Middle Elwha", "Aldwell Reservoir", "Lower Elwha", "Delta")) +
  scale_fill_manual(values=c("red1", "orange1", "yellow2", "green1","blue1", "purple2")) +
  ylim(0, 350) +
  theme(legend.position = c(0.3, 0.7)) +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 15))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  geom_text(size=4, x=1.0, y=100, label="a") +
  geom_text(size=4, x=2.0, y=100, label="a") +
  geom_text(size=4, x=3.0, y=100, label="a") +
  geom_text(size=4, x=4.0, y=100, label="a") +
  geom_text(size=4, x=5.0, y=100, label="a") +
  geom_text(size=4, x=6.0, y=120, label="a")
plot1

#plot 2:
data_summary <- function(fungal, H.est, region){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]]))
  }
  data_sum<-ddply(fungal, region, .fun=summary_func,
                  H.est)
  data_sum <- plyr::rename(data_sum, c("mean" = H.est))
  return(data_sum)
}

fun2 <- data_summary(fungal, H.est="H.est", 
                     region=c("region"))
# Convert region to a factor variable
fun2$region=as.factor(fun2$region)
head(fun2)

plot2 <- p<- ggplot(fun2, aes(x=region, y=H.est, fill=region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=H.est-se, ymax=H.est+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Elwha Section", y = "Shannon's Diversity Index", tag = "B")+
  theme_classic() +
  scale_x_discrete(limits = c("Mills Reservoir","Tributary","Middle Elwha", "Aldwell Reservoir", "Lower Elwha", "Delta")) +
  scale_fill_manual(values=c("red1", "orange1", "yellow2", "green1","blue1", "purple2")) +
  ylim(0, 10) +
  theme(legend.position ="none") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(text = element_text(size = 15))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  geom_text(size=4, x=1.0, y=3, label="b") +
  geom_text(size=4, x=2.0, y=8, label="ab") +
  geom_text(size=4, x=3.0, y=7.3, label="ab") +
  geom_text(size=4, x=4.0, y=4.7, label="ab") +
  geom_text(size=4, x=5.0, y=10, label="a") +
  geom_text(size=4, x=6.0, y=4.2, label="ab")
plot2

# Arrange plots A, B and C in one column:
g1 <- ggplotGrob(plot1)
g2 <- ggplotGrob(plot2)
g <- grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2)))

# Export pdf and jpg directly from viewing pane - can't get automatic saves to work...


# Permutative Simple Linear Regressions between Fungal Diversity metrics and 
# environmental values: use dataset FungalDiversity.csv 
#
slr1 = lmp(H.est~width, data=fungal)
summary(slr1)

slr2 = lmp(H.est~CC, data=fungal)
summary(slr2)

#Coefficients:
# Estimate Iter Pr(Prob)   
#CC  0.04207 5000   0.0044 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.326 on 46 degrees of freedom
#Multiple R-Squared: 0.1537,	Adjusted R-squared: 0.1353 
#F-statistic: 8.355 on 1 and 46 DF,  p-value: 0.005851

slr3 = lmp(H.est~fines, data=fungal)
summary(slr3)

slr4 = lmp(H.est~D50, data=fungal)
summary(slr4)

slr5 = lmp(H.est~temp, data=fungal)
summary(slr5)

slr6 = lmp(H.est~cond, data=fungal)
summary(slr6)

slr7 = lmp(H.est~p, data=fungal)
summary(slr7)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#p   0.2868 5000   0.0096 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.539 on 38 degrees of freedom
#Multiple R-Squared: 0.1456,	Adjusted R-squared: 0.1231 
#F-statistic: 6.477 on 1 and 38 DF,  p-value: 0.01511

slr8 = lmp(H.est~n, data=fungal)
summary(slr8)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#n  0.02823 4022   0.0244 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.586 on 38 degrees of freedom
#Multiple R-Squared: 0.1227,	Adjusted R-squared: 0.09956 
#F-statistic: 5.312 on 1 and 38 DF,  p-value: 0.02673 

slr9 = lmp(H.est~po4, data=fungal)
summary(slr9)

slr10 = lmp(H.est~sio4, data=fungal)
summary(slr10)

slr11 = lmp(H.est~no3, data=fungal)
summary(slr11)

slr12 = lmp(H.est~no2, data=fungal)
summary(slr12)

slr13 = lmp(H.est~nh4, data=fungal)
summary(slr13)


slr14 = lmp(H.est~chla, data=fungal)
summary(slr14)

slr15 = lmp(H.est~omRock, data=fungal)
summary(slr15)

slr16 = lmp(H.est~imRock, data=fungal)
summary(slr16)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#imRock   0.5326  958    0.095 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.623 on 38 degrees of freedom
#Multiple R-Squared: 0.1043,	Adjusted R-squared: 0.08073 
#F-statistic: 4.425 on 1 and 38 DF,  p-value: 0.04209 

slr17 = lmp(H.est~omSeston, data=fungal)
summary(slr17)

slr18 = lmp(H.est~imSeston, data=fungal)
summary(slr18)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#imSeston   0.5344 2869   0.0338 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.456 on 46 degrees of freedom
#Multiple R-Squared: 0.08595,	Adjusted R-squared: 0.06607 
#F-statistic: 4.325 on 1 and 46 DF,  p-value: 0.04315 

slr19 = lmp(H.est~BIden, data=fungal)
summary(slr19)

slr20 = lmp(H.est~shred, data=fungal)
summary(slr20)


slr20.1 = lmp(H.est~k, data=fungal)
summary(slr20.1)

slr20.2 = lmp(H.est~depth, data=fungal)
summary(slr20.2)

slr20.3 = lmp(H.est~elevation, data=fungal)
summary(slr20.3)

slr20.4 = lmp(H.est~no2no3, data=fungal)
summary(slr20.4)

slr20.5 = lmp(H.est~Sed1, data=fungal)
summary(slr20.5)

slr20.6 = lmp(H.est~Sed2, data=fungal)
summary(slr20.6)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#Sed2   -2.474 1650   0.0576 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.622 on 38 degrees of freedom
#Multiple R-Squared: 0.1049,	Adjusted R-squared: 0.08135 
#F-statistic: 4.454 on 1 and 38 DF,  p-value: 0.04147

slr20.7 = lmp(H.est~Sed3, data=fungal)
summary(slr20.7)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#Sed3   -2.712 2836   0.0342 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.572 on 38 degrees of freedom
#Multiple R-Squared: 0.1293,	Adjusted R-squared: 0.1064 
#F-statistic: 5.644 on 1 and 38 DF,  p-value: 0.02267 

slr20.8 = lmp(H.est~Sed4, data=fungal)
summary(slr20.8)

slr20.9 = lmp(H.est~Sed5, data=fungal)
summary(slr20.9)

slr20.10 = lmp(H.est~Sed6, data=fungal)
summary(slr20.10)

slr20.11 = lmp(H.est~Sed7, data=fungal)
summary(slr20.11)

slr20.12 = lmp(H.est~Sed8, data=fungal)
summary(slr20.12)

slr20.13 = lmp(H.est~Sed9, data=fungal)
summary(slr20.13)

slr20.14 = lmp(H.est~Sed10, data=fungal)
summary(slr20.14)

slr20.15 = lmp(H.est~Sed11, data=fungal)
summary(slr20.15)

slr20.16 = lmp(H.est~SedTot, data=fungal)
summary(slr20.16)
#Coefficients:
#  Estimate Iter Pr(Prob)
#SedTot -0.00929  397    0.202
#
#Residual standard error: 3.624 on 38 degrees of freedom
#Multiple R-Squared: 0.104,	Adjusted R-squared: 0.08046 
#F-statistic: 4.412 on 1 and 38 DF,  p-value: 0.04237


slr21 = lmp(H.est~P1, data=fungal)
summary(slr21)

slr22 = lmp(H.est~P2, data=fungal)
summary(slr22)

slr23 = lmp(H.est~P3, data=fungal)
summary(slr23)

slr24 = lmp(H.est~P4, data=fungal)
summary(slr24)

slr25 = lmp(H.est~P5, data=fungal)
summary(slr25)

slr26 = lmp(H.est~P6, data=fungal)
summary(slr26)

slr27 = lmp(H.est~P7, data=fungal)
summary(slr27)

slr28 = lmp(H.est~P8, data=fungal)
summary(slr28)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#P8     14.5 5000   0.0058 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.499 on 38 degrees of freedom
#Multiple R-Squared: 0.1645,	Adjusted R-squared: 0.1425 
#F-statistic: 7.482 on 1 and 38 DF,  p-value: 0.009418

slr29 = lmp(H.est~P9, data=fungal)
summary(slr29)


slr30 = lmp(H.est~P10, data=fungal)
summary(slr30)

slr31 = lmp(H.est~P11, data=fungal)
summary(slr31)

slr32 = lmp(H.est~Sum, data=fungal)
summary(slr32)

slr33 = lmp(H.est~Rich, data=fungal)
summary(slr33)

slr34 = lmp(H.est~EvenM, data=fungal)
summary(slr34)

slr35 = lmp(H.est~Div, data=fungal)
summary(slr35)

slr36 = lmp(H.est~SimpsonM, data=fungal)
summary(slr36)

######################################
slr41 = lmp(S.est~width, data=fungal)
summary(slr41)

slr42 = lmp(S.est~CC, data=fungal)
summary(slr42)

slr43 = lmp(S.est~fines, data=fungal)
summary(slr43)

slr44 = lmp(S.est~D50, data=fungal)
summary(slr44)

slr45 = lmp(S.est~temp, data=fungal)
summary(slr45)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#temp   -6.186 5000   0.0028 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 22.41 on 46 degrees of freedom
#Multiple R-Squared: 0.1509,	Adjusted R-squared: 0.1324 
#F-statistic: 8.173 on 1 and 46 DF,  p-value: 0.006369

slr46 = lmp(S.est~cond, data=fungal)
summary(slr46)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#cond  -0.1518 4070   0.0241 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 23.1 on 46 degrees of freedom
#Multiple R-Squared: 0.09769,	Adjusted R-squared: 0.07808 
#F-statistic:  4.98 on 1 and 46 DF,  p-value: 0.03055 

slr47 = lmp(S.est~p, data=fungal)
summary(slr47)

slr48 = lmp(S.est~n, data=fungal)
summary(slr48)

slr49 = lmp(S.est~po4, data=fungal)
summary(slr49)

slr50 = lmp(S.est~sio4, data=fungal)
summary(slr50)

slr51 = lmp(S.est~no3, data=fungal)
summary(slr51)

slr52 = lmp(S.est~no2, data=fungal)
summary(slr52)

slr53 = lmp(S.est~nh4, data=fungal)
summary(slr53)


slr54 = lmp(S.est~chla, data=fungal)
summary(slr54)
#Coefficients:
#  Estimate Iter Pr(Prob)
#chla   -11.72  804    0.111
#
#Residual standard error: 23.26 on 46 degrees of freedom
#Multiple R-Squared: 0.08482,	Adjusted R-squared: 0.06492 
#F-statistic: 4.263 on 1 and 46 DF,  p-value: 0.04461

slr55 = lmp(S.est~omRock, data=fungal)
summary(slr55)

slr56 = lmp(S.est~imRock, data=fungal)
summary(slr56)

slr57 = lmp(S.est~omSeston, data=fungal)
summary(slr57)

slr58 = lmp(S.est~imSeston, data=fungal)
summary(slr58)

slr59 = lmp(S.est~BIden, data=fungal)
summary(slr59)

slr60 = lmp(S.est~shred, data=fungal)
summary(slr60)


slr60.1 = lmp(S.est~k, data=fungal)
summary(slr60.1)

slr60.2 = lmp(S.est~depth, data=fungal)
summary(slr60.2)

slr60.3 = lmp(S.est~elevation, data=fungal)
summary(slr60.3)

slr60.4 = lmp(S.est~no2no3, data=fungal)
summary(slr60.4)

slr60.5 = lmp(S.est~Sed1, data=fungal)
summary(slr60.5)

slr60.6 = lmp(S.est~Sed2, data=fungal)
summary(slr60.6)

slr60.7 = lmp(S.est~Sed3, data=fungal)
summary(slr60.7)

slr60.8 = lmp(S.est~Sed4, data=fungal)
summary(slr60.8)

slr60.9 = lmp(S.est~Sed5, data=fungal)
summary(slr60.9)

slr60.10 = lmp(S.est~Sed6, data=fungal)
summary(slr60.10)

slr60.11 = lmp(S.est~Sed7, data=fungal)
summary(slr60.11)

slr60.12 = lmp(S.est~Sed8, data=fungal)
summary(slr60.12)

slr60.13 = lmp(S.est~Sed9, data=fungal)
summary(slr60.13)

slr60.14 = lmp(S.est~Sed10, data=fungal)
summary(slr60.14)

slr60.15 = lmp(S.est~Sed11, data=fungal)
summary(slr60.15)

slr60.16 = lmp(S.est~SedTot, data=fungal)
summary(slr60.16)


slr61 = lmp(S.est~P1, data=fungal)
summary(slr61)

slr62 = lmp(S.est~P2, data=fungal)
summary(slr62)

slr63 = lmp(S.est~P3, data=fungal)
summary(slr63)

slr64 = lmp(S.est~P4, data=fungal)
summary(slr64)

slr65 = lmp(S.est~P5, data=fungal)
summary(slr65)

slr66 = lmp(S.est~P6, data=fungal)
summary(slr66)

slr67 = lmp(S.est~P7, data=fungal)
summary(slr67)

slr68 = lmp(S.est~P8, data=fungal)
summary(slr68)

slr69 = lmp(S.est~P9, data=fungal)
summary(slr69)


slr70 = lmp(S.est~P10, data=fungal)
summary(slr70)

slr71 = lmp(S.est~P11, data=fungal)
summary(slr71)

slr72 = lmp(S.est~Sum, data=fungal)
summary(slr72)

slr73 = lmp(S.est~Rich, data=fungal)
summary(slr73)

slr74 = lmp(S.est~EvenM, data=fungal)
summary(slr74)

slr76 = lmp(S.est~SimpsonM, data=fungal)
summary(slr76)

slr77 = lmp(S.est~Div, data=fungal)
summary(slr77)


######################################