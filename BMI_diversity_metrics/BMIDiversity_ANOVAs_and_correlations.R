# Carri J. LeRoy, 4-26-23
#
# Permutative ANOVAs + linear regressions with BMI diversity, k and environmental variables
#
# For 2016 Elwha paper: 
# Dataset = BMIDiversity.csv
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
library(lmPerm)           # permutative analyses
library(grid)

#
# set WD to folder
# then read the data in: 
BMI <- read.csv("BMIDiversity.csv", row.names = 1, header = TRUE)
#BMI <- BMIDiv
#
# Run ANOVAs for H, D, and S against region:

BMI$region = factor(BMI$region,
                     levels=unique(BMI$region))

str(BMI)
#
# Set a random seed
set.seed(1431)

#permutative ANOVAs to see which regions are different:
#
#
fit1 <- aovp(Sum~region, data = BMI, perm="Prob")
anova(fit1)
#
#Analysis of Variance Table
#
#Response: Sum
#Df  R Sum Sq R Mean Sq Iter Pr(Prob)   
#region     4 164592990  41148247 5000    0.001 *** F= 2.97
#  Residuals 15 207769583  13851306                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
HSD.test(fit1, "region", group=TRUE)
out1 <- HSD.test(fit1, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#Sum groups
#Lower Elwha       7583.5      a
#Aldwell Reservoir 1101.5      a
#Elwha Middle       362.0      a
#Delta              278.0      a
#Tributary          107.5      a

fit2 <- aovp(Rich~region, data = BMI, perm="Prob")
anova(fit2)
#Analysis of Variance Table
#
#Response: Rich
#Df R Sum Sq R Mean Sq Iter Pr(Prob)  
#region     4       84      21.0 5000   0.0114 * F = 4.77
#  Residuals 15       66       4.4                
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#

HSD.test(fit2, "region", group=TRUE)
out1 <- HSD.test(fit2, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#Rich groups
#Aldwell Reservoir 12.5      a
#Lower Elwha       11.5      a
#Elwha Middle      10.0     ab
#Tributary          9.5     ab
#Delta              6.5      b

fit3 <- aovp(EvenM~region, data = BMI, perm="Prob")
anova(fit3)
#
#Analysis of Variance Table
#
#Response: EvenM
#Df R Sum Sq R Mean Sq Iter  Pr(Prob)    
#region     4  0.74928  0.187320 5000 < 2.2e-16 *** F= 12.60
#  Residuals 15  0.22295  0.014864                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#

HSD.test(fit3, "region", group=TRUE)
out1 <- HSD.test(fit3, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)
#EvenM groups
#Tributary         0.75575      a
#Elwha Middle      0.57700      a
#Delta             0.55725      a
#Aldwell Reservoir 0.26350      b
#Lower Elwha       0.25850      b

fit4 <- aovp(Div~region, data=BMI, perm="Prob")
anova(fit4)
#
#Analysis of Variance Table
#
#Response: Div
#Df R Sum Sq R Mean Sq Iter Pr(Prob)   
#region     4   3.3881   0.84702 5000    0.0026 ** F-ratio: 8.38
#  Residuals 15   1.5162   0.10108                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

HSD.test(fit4, "region", group=TRUE)
out1 <- HSD.test(fit4, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#Div groups
#Tributary         1.69575      a
#Elwha Middle      1.32775     ab
#Delta             1.01025    abc
#Aldwell Reservoir 0.64475     bc
#Lower Elwha       0.61700      c

fit5 <- aovp(SimpsonM~region, data=BMI, perm="Prob")
anova(fit5)
#
#Analysis of Variance Table
#
#Response: SimpsonM
#Df R Sum Sq R Mean Sq Iter Pr(Prob)   
#region     4  0.72597  0.181492 5000   0.0018 **
#  Residuals 15  0.30941  0.020627                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

HSD.test(fit5, "region", group=TRUE)
out2 <- HSD.test(fit5, "region", alpha = 0.05, group = TRUE, main = NULL, unbalanced = TRUE, console = TRUE)

#SimpsonM groups
#Tributary         0.748925      a
#Elwha Middle      0.656900      a
#Delta             0.507475     ab
#Lower Elwha       0.293600      b
#Aldwell Reservoir 0.271200      b

#ORIGINAL BAR CHART
#Now make a plot for Shannon's H (Div): The normal ggboxplots are weird because the medians and means 
#are pretty different, need bar charts with means and standard errors plotted:
data_summary <- function(BMI, Div, region){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]]))
  }
  data_sum<-ddply(BMI, region, .fun=summary_func,
                  Div)
  data_sum <- plyr::rename(data_sum, c("mean" = Div))
  return(data_sum)
}

fun2 <- data_summary(BMI, Div="Div", 
                     region=c("region"))
# Convert region to a factor variable
fun2$region=as.factor(fun2$region)
head(fun2)

plotF <- p<- ggplot(fun2, aes(x=region, y=Div, fill=region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Div-se, ymax=Div+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Elwha Section", y = "Shannon's Diversity Index")+
  theme_classic() +
  scale_x_discrete(limits = c("Tributary","Middle Elwha", "Aldwell Reservoir", "Lower Elwha", "Delta")) +
  scale_fill_manual(values=c("orange1", "yellow2", "green1", "blue1", "purple2")) +
  ylim(0, 5.0) +
  theme(legend.position = c(0.18, 0.8)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(text = element_text(size = 15))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  geom_text(size=4, x=1.0, y=2, label="a") +
  geom_text(size=4, x=2.0, y=1.8, label="ab") +
  geom_text(size=4, x=3.0, y=1.1, label="bc") +
  geom_text(size=4, x=4.0, y=1.0, label="c") +
  geom_text(size=4, x=5.0, y=1.5, label="abc")
plotF

ggsave(file="DivBMI.pdf", plotF, width=16.4, height=12.4, units="cm", dpi=800)
ggsave(plotF, file="DivBMI.eps", width=16.4, height=12.4, units="cm", dpi=800, device="eps")
ggsave(plotF, file="DivBMI.jpg", width=16.4, height=12.4, units="cm", dpi=800, device="jpg")

# 3 bar charts - added following review: 
#Plot A:
data_summary <- function(BMI, Sum, region){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]]))
  }
  data_sum<-ddply(BMI, region, .fun=summary_func,
                  Sum)
  data_sum <- plyr::rename(data_sum, c("mean" = Sum))
  return(data_sum)
}

funA <- data_summary(BMI, Sum="Sum", 
                     region=c("region"))
# Convert region to a factor variable
funA$region=as.factor(funA$region)
head(funA)

plotA <- p<- ggplot(funA, aes(x=region, y=Sum, fill=region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Sum-se, ymax=Sum+se), width=.2,
                position=position_dodge(.9)) +
  labs(x = element_blank(), y = "Total Abundance", tag = "A")+
  theme_classic() +
  scale_x_discrete(limits = c("Tributary","Middle Elwha", "Aldwell Reservoir", "Lower Elwha", "Delta")) +
  scale_fill_manual(values=c("orange1", "yellow2", "green1", "blue1", "purple2")) +
  ylim(0, 15000) +
  theme(legend.position = c(0.3, 0.7)) +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 15))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  geom_text(size=4, x=1.0, y=900, label="a") +
  geom_text(size=4, x=2.0, y=1500, label="ab") +
  geom_text(size=4, x=3.0, y=2500, label="bc") +
  geom_text(size=4, x=4.0, y=12500, label="c") +
  geom_text(size=4, x=5.0, y=1500, label="abc")
plotA

#plot B:
data_summary <- function(BMI, Rich, region){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]]))
  }
  data_sum<-ddply(BMI, region, .fun=summary_func,
                  Rich)
  data_sum <- plyr::rename(data_sum, c("mean" = Rich))
  return(data_sum)
}

funB <- data_summary(BMI, Rich="Rich", 
                     region=c("region"))
# Convert region to a factor variable
funB$region=as.factor(funB$region)
head(funB)

plotB <- p<- ggplot(funB, aes(x=region, y=Rich, fill=region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Rich-se, ymax=Rich+se), width=.2,
                position=position_dodge(.9)) +
  labs(x = element_blank(), y = "Taxa Richness", tag = "B")+
  theme_classic() +
  scale_x_discrete(limits = c("Tributary","Middle Elwha", "Aldwell Reservoir", "Lower Elwha", "Delta")) +
  scale_fill_manual(values=c("orange1", "yellow2", "green1", "blue1", "purple2")) +
  ylim(0, 17.0) +
  theme(legend.position ="none") +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 15))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  geom_text(size=4, x=1.0, y=12, label="ab") +
  geom_text(size=4, x=2.0, y=12, label="ab") +
  geom_text(size=4, x=3.0, y=15, label="a") +
  geom_text(size=4, x=4.0, y=14, label="a") +
  geom_text(size=4, x=5.0, y=9, label="b")
plotB

#Plot C:

data_summary <- function(BMI, Div, region){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]]))
  }
  data_sum<-ddply(BMI, region, .fun=summary_func,
                  Div)
  data_sum <- plyr::rename(data_sum, c("mean" = Div))
  return(data_sum)
}

funC <- data_summary(BMI, Div="Div", 
                     region=c("region"))
# Convert region to a factor variable
funC$region=as.factor(funC$region)
head(funC)

plotC <- p<- ggplot(funC, aes(x=region, y=Div, fill=region)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Div-se, ymax=Div+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Elwha Section", y = "Shannon's Diversity Index", tag = "C")+
  theme_classic() +
  scale_x_discrete(limits = c("Tributary","Middle Elwha", "Aldwell Reservoir", "Lower Elwha", "Delta")) +
  scale_fill_manual(values=c("orange1", "yellow2", "green1", "blue1", "purple2")) +
  ylim(0, 2.5) +
  theme(legend.position ="none") +
  theme(axis.text.x = element_text(angle=90)) +
  theme(text = element_text(size = 15))+
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.background = element_rect(fill="white", size=0.5)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  geom_text(size=4, x=1.0, y=2, label="a") +
  geom_text(size=4, x=2.0, y=1.8, label="ab") +
  geom_text(size=4, x=3.0, y=1.1, label="bc") +
  geom_text(size=4, x=4.0, y=1.0, label="c") +
  geom_text(size=4, x=5.0, y=1.5, label="abc")
plotC

# Arrange plots A, B and C in one column:
gA <- ggplotGrob(plotA)
gB <- ggplotGrob(plotB)
gC <- ggplotGrob(plotC)
g <- grid.draw(rbind(ggplotGrob(plotA), ggplotGrob(plotB), ggplotGrob(plotC)))

# Just export g from the viewing pane - for some reason the code below is not working...
#
ggsave(file="3-bars_BMI.pdf", g, width=12, height=20, units="cm", dpi=800)
ggsave(g, file="3-bars_BMI.jpg", width=12, height=20, units="cm", dpi=800, device="jpg")


# Regressions between BMI Diversity metrics and environmental values: use dataset BMIDiversity.csv 
#
slr1 = lmp(Div~width, data=BMI)
summary(slr1)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#width -0.01127 1277   0.0728 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4805 on 14 degrees of freedom
#Multiple R-Squared: 0.2636,	Adjusted R-squared: 0.211 
#F-statistic: 5.012 on 1 and 14 DF,  p-value: 0.04194 

slr2 = lmp(Div~CC, data=BMI)
summary(slr2)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#CC 0.008305 5000   0.0156 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
##Residual standard error: 0.4289 on 18 degrees of freedom
#Multiple R-Squared: 0.325,	Adjusted R-squared: 0.2875 
#F-statistic: 8.665 on 1 and 18 DF,  p-value: 0.008685 

slr3 = lmp(Div~Fines, data=BMI)
summary(slr3)

slr4 = lmp(Div~D50, data=BMI)
summary(slr4)

slr5 = lmp(Div~Temp, data=BMI)
summary(slr5)

slr6 = lmp(Div~Cond, data=BMI)
summary(slr6)

slr7 = lmp(Div~p, data=BMI)
summary(slr7)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#p  0.05345 2923   0.0332 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4855 on 14 degrees of freedom
##Multiple R-Squared: 0.248,	Adjusted R-squared: 0.1943 
#F-statistic: 4.616 on 1 and 14 DF,  p-value: 0.04965

slr8 = lmp(Div~n, data=BMI)
summary(slr8)

slr9 = lmp(Div~po4, data=BMI)
summary(slr9)

slr10 = lmp(Div~sio4, data=BMI)
summary(slr10)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#sio4 0.0002086 5000   0.0096 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.443 on 14 degrees of freedom
#Multiple R-Squared: 0.3741,	Adjusted R-squared: 0.3294 
#F-statistic: 8.367 on 1 and 14 DF,  p-value: 0.01182 

slr11 = lmp(Div~No3, data=BMI)
summary(slr11)

slr12 = lmp(Div~no2, data=BMI)
summary(slr12)

slr13 = lmp(Div~nh4, data=BMI)
summary(slr13)


slr14 = lmp(Div~chla, data=BMI)
summary(slr14)

slr15 = lmp(Div~omRock, data=BMI)
summary(slr15)

slr16 = lmp(Div~imRock, data=BMI)
summary(slr16)

slr17 = lmp(Div~OM, data=BMI)
summary(slr17)

slr18 = lmp(Div~imSeston, data=BMI)
summary(slr18)

slr19 = lmp(Div~BI, data=BMI)
summary(slr19)

slr20 = lmp(Div~SH, data=BMI)
summary(slr20)


slr20.1 = lmp(Div~k, data=BMI)
summary(slr20.1)

slr20.2 = lmp(Div~depth, data=BMI)
summary(slr20.2)

slr20.3 = lmp(Div~Elev, data=BMI)
summary(slr20.3)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#Elev 0.006644 5000   0.0148 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4533 on 18 degrees of freedom
#Multiple R-Squared: 0.2458,	Adjusted R-squared: 0.2039 
#F-statistic: 5.867 on 1 and 18 DF,  p-value: 0.0262

slr20.4 = lmp(Div~no2no3, data=BMI)
summary(slr20.4)

slr20.5 = lmp(Div~Sed1, data=BMI)
summary(slr20.5)

slr20.6 = lmp(Div~Sed2, data=BMI)
summary(slr20.6)

slr20.7 = lmp(Div~Sed3, data=BMI)
summary(slr20.7)

slr20.8 = lmp(Div~Sed4, data=BMI)
summary(slr20.8)

slr20.9 = lmp(Div~Sed5, data=BMI)
summary(slr20.9)

slr20.10 = lmp(Div~Sed6, data=BMI)
summary(slr20.10)

slr20.11 = lmp(Div~Sed7, data=BMI)
summary(slr20.11)

slr20.12 = lmp(Div~Sed8, data=BMI)
summary(slr20.12)

slr20.13 = lmp(Div~Sed9, data=BMI)
summary(slr20.13)

slr20.14 = lmp(Div~Sed10, data=BMI)
summary(slr20.14)

slr20.15 = lmp(Div~Sed11, data=BMI)
summary(slr20.15)

slr20.16 = lmp(Div~SedTot, data=BMI)
summary(slr20.16)


slr21 = lmp(Div~P1, data=BMI)
summary(slr21)

slr22 = lmp(Div~P2, data=BMI)
summary(slr22)

slr23 = lmp(Div~P3, data=BMI)
summary(slr23)

slr24 = lmp(Div~P4, data=BMI)
summary(slr24)

slr25 = lmp(Div~P5, data=BMI)
summary(slr25)

slr26 = lmp(Div~P6, data=BMI)
summary(slr26)

slr27 = lmp(Div~P7, data=BMI)
summary(slr27)

slr28 = lmp(Div~P8, data=BMI)
summary(slr28)

slr29 = lmp(Div~P9, data=BMI)
summary(slr29)


slr30 = lmp(Div~P10, data=BMI)
summary(slr30)

slr31 = lmp(Div~P11, data=BMI)
summary(slr31)

slr32 = lmp(Div~Sum, data=BMI)
summary(slr32)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#Sum -5.706e-05 5000    0.002 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4529 on 18 degrees of freedom
##Multiple R-Squared: 0.2472,	Adjusted R-squared: 0.2054 
#F-statistic: 5.911 on 1 and 18 DF,  p-value: 0.02572 

slr33 = lmp(Div~Rich, data=BMI)
summary(slr33)

slr34 = lmp(Div~EvenM, data=BMI)
summary(slr34)
#Coefficients:
#  Estimate Iter Pr(Prob)    
#EvenM    2.131 5000   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1647 on 18 degrees of freedom
##Multiple R-Squared: 0.9004,	Adjusted R-squared: 0.8949 
#F-statistic: 162.8 on 1 and 18 DF,  p-value: 1.87e-10 

slr36 = lmp(Div~SimpsonM, data=BMI)
summary(slr36)
#Coefficients:
#  Estimate Iter Pr(Prob)    
#SimpsonM    2.129 5000   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1088 on 18 degrees of freedom
#Multiple R-Squared: 0.9566,	Adjusted R-squared: 0.9541 
#F-statistic: 396.3 on 1 and 18 DF,  p-value: 1.043e-13

slr39 = lmp(Div~RichF, data=BMI)
summary(slr39)

slr40 = lmp(Div~SimpF, data=BMI)
summary(slr40)

slr40.1 = lmp(Div~ShanF, data=BMI)
summary(slr40.1)

######################################
slr41 = lmp(Rich~width, data=BMI)
summary(slr41)

slr42 = lmp(Rich~CC, data=BMI)
summary(slr42)

slr43 = lmp(Rich~Fines, data=BMI)
summary(slr43)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#Fines  -0.0937 3095   0.0313 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 2.51 on 18 degrees of freedom
#Multiple R-Squared: 0.2442,	Adjusted R-squared: 0.2023 
#F-statistic: 5.817 on 1 and 18 DF,  p-value: 0.02677 

slr44 = lmp(Rich~D50, data=BMI)
summary(slr44)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#D50  0.06007 5000   0.0044 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 2.172 on 18 degrees of freedom
#Multiple R-Squared: 0.4341,	Adjusted R-squared: 0.4027 
#F-statistic: 13.81 on 1 and 18 DF,  p-value: 0.001582 

slr45 = lmp(Rich~Temp, data=BMI)
summary(slr45)

slr46 = lmp(Rich~Cond, data=BMI)
summary(slr46)

slr47 = lmp(Rich~p, data=BMI)
summary(slr47)

slr48 = lmp(Rich~n, data=BMI)
summary(slr48)

slr49 = lmp(Rich~po4, data=BMI)
summary(slr49)

slr50 = lmp(Rich~sio4, data=BMI)
summary(slr50)

slr51 = lmp(Rich~No3, data=BMI)
summary(slr51)

slr52 = lmp(Rich~no2, data=BMI)
summary(slr52)

slr53 = lmp(Rich~nh4, data=BMI)
summary(slr53)


slr54 = lmp(Rich~chla, data=BMI)
summary(slr54)

slr55 = lmp(Rich~omRock, data=BMI)
summary(slr55)

slr56 = lmp(Rich~imRock, data=BMI)
summary(slr56)

slr57 = lmp(Rich~OM, data=BMI)
summary(slr57)

slr58 = lmp(Rich~imSeston, data=BMI)
summary(slr58)

slr59 = lmp(Rich~BI, data=BMI)
summary(slr59)

slr60 = lmp(Rich~SH, data=BMI)
summary(slr60)


slr60.1 = lmp(Rich~k, data=BMI)
summary(slr60.1)

slr60.2 = lmp(Rich~depth, data=BMI)
summary(slr60.2)

slr60.3 = lmp(Rich~Elev, data=BMI)
summary(slr60.3)

slr60.4 = lmp(Rich~no2no3, data=BMI)
summary(slr60.4)

slr60.5 = lmp(Rich~Sed1, data=BMI)
summary(slr60.5)

slr60.6 = lmp(Rich~Sed2, data=BMI)
summary(slr60.6)

slr60.7 = lmp(Rich~Sed3, data=BMI)
summary(slr60.7)

slr60.8 = lmp(Rich~Sed4, data=BMI)
summary(slr60.8)

slr60.9 = lmp(Rich~Sed5, data=BMI)
summary(slr60.9)

slr60.10 = lmp(Rich~Sed6, data=BMI)
summary(slr60.10)

slr60.11 = lmp(Rich~Sed7, data=BMI)
summary(slr60.11)

slr60.12 = lmp(Rich~Sed8, data=BMI)
summary(slr60.12)

slr60.13 = lmp(Rich~Sed9, data=BMI)
summary(slr60.13)

slr60.14 = lmp(Rich~Sed10, data=BMI)
summary(slr60.14)

slr60.15 = lmp(Rich~Sed11, data=BMI)
summary(slr60.15)

slr60.16 = lmp(Rich~SedTot, data=BMI)
summary(slr60.16)


slr61 = lmp(Rich~P1, data=BMI)
summary(slr61)

slr62 = lmp(Rich~P2, data=BMI)
summary(slr62)

slr63 = lmp(Rich~P3, data=BMI)
summary(slr63)

slr64 = lmp(Rich~P4, data=BMI)
summary(slr64)

slr65 = lmp(Rich~P5, data=BMI)
summary(slr65)

slr66 = lmp(Rich~P6, data=BMI)
summary(slr66)

slr67 = lmp(Rich~P7, data=BMI)
summary(slr67)

slr68 = lmp(Rich~P8, data=BMI)
summary(slr68)

slr69 = lmp(Rich~P9, data=BMI)
summary(slr69)


slr70 = lmp(Rich~P10, data=BMI)
summary(slr70)

slr71 = lmp(Rich~P11, data=BMI)
summary(slr71)

slr72 = lmp(Rich~Sum, data=BMI)
summary(slr72)


slr74 = lmp(Rich~EvenM, data=BMI)
summary(slr74)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#EvenM   -6.013 3279   0.0299 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 2.526 on 18 degrees of freedom
#Multiple R-Squared: 0.2343,	Adjusted R-squared: 0.1918 
#F-statistic: 5.509 on 1 and 18 DF,  p-value: 0.03056 

slr75 = lmp(Rich~Div, data=BMI)
summary(slr75)

slr76 = lmp(Rich~SimpsonM, data=BMI)
summary(slr76)


slr79 = lmp(Rich~RichF, data=BMI)
summary(slr79)

slr80 = lmp(Rich~SimpF, data=BMI)
summary(slr80)

slr80.1 = lmp(Rich~ShanF, data=BMI)
summary(slr80.1)

######################################

slr91 = lmp(Sum~width, data=BMI)
summary(slr91)

slr92 = lmp(Sum~CC, data=BMI)
summary(slr92)

slr93 = lmp(Sum~Fines, data=BMI)
summary(slr93)

slr94 = lmp(Sum~D50, data=BMI)
summary(slr94)


slr95 = lmp(Sum~Temp, data=BMI)
summary(slr95)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#Temp     1305 5000    0.019 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4003 on 18 degrees of freedom
##Multiple R-Squared: 0.2254,	Adjusted R-squared: 0.1824 
#F-statistic: 5.238 on 1 and 18 DF,  p-value: 0.03441 

slr96 = lmp(Sum~Cond, data=BMI)
summary(slr96)

slr97 = lmp(Sum~p, data=BMI)
summary(slr97)

slr98 = lmp(Sum~n, data=BMI)
summary(slr98)

slr99 = lmp(Sum~po4, data=BMI)
summary(slr99)

slr100 = lmp(Sum~sio4, data=BMI)
summary(slr100)

slr101 = lmp(Sum~No3, data=BMI)
summary(slr101)

slr102 = lmp(Sum~no2, data=BMI)
summary(slr102)

slr103 = lmp(Sum~nh4, data=BMI)
summary(slr103)


slr104 = lmp(Sum~chla, data=BMI)
summary(slr104)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#chla     3148 2374   0.0404 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
##Residual standard error: 4025 on 18 degrees of freedom
#Multiple R-Squared: 0.2167,	Adjusted R-squared: 0.1732 
#F-statistic: 4.979 on 1 and 18 DF,  p-value: 0.03861

slr105 = lmp(Sum~omRock, data=BMI)
summary(slr105)

slr106 = lmp(Sum~imRock, data=BMI)
summary(slr106)

slr107 = lmp(Sum~OM, data=BMI)
summary(slr107)

slr108 = lmp(Sum~imSeston, data=BMI)
summary(slr108)

slr109 = lmp(Sum~BI, data=BMI)
summary(slr109)

slr110 = lmp(Sum~SH, data=BMI)
summary(slr110)


slr6.1 = lmp(Sum~k, data=BMI)
summary(slr6.1)

slr6.2 = lmp(Sum~depth, data=BMI)
summary(slr6.2)

slr6.3 = lmp(Sum~Elev, data=BMI)
summary(slr6.3)

slr6.4 = lmp(Sum~no2no3, data=BMI)
summary(slr6.4)

slr6.5 = lmp(Sum~Sed1, data=BMI)
summary(slr6.5)

slr6.6 = lmp(Sum~Sed2, data=BMI)
summary(slr6.6)

slr6.7 = lmp(Sum~Sed3, data=BMI)
summary(slr6.7)

slr6.8 = lmp(Sum~Sed4, data=BMI)
summary(slr6.8)

slr6.9 = lmp(Sum~Sed5, data=BMI)
summary(slr6.9)

slr6.10 = lmp(Sum~Sed6, data=BMI)
summary(slr6.10)

slr6.11 = lmp(Sum~Sed7, data=BMI)
summary(slr6.11)

slr6.12 = lmp(Sum~Sed8, data=BMI)
summary(slr6.12)

slr6.13 = lmp(Sum~Sed9, data=BMI)
summary(slr6.13)

slr6.14 = lmp(Sum~Sed10, data=BMI)
summary(slr6.14)

slr6.15 = lmp(Sum~Sed11, data=BMI)
summary(slr6.15)

slr6.16 = lmp(Sum~SedTot, data=BMI)
summary(slr6.16)


slr121 = lmp(Sum~P1, data=BMI)
summary(slr121)

slr122 = lmp(Sum~P2, data=BMI)
summary(slr122)

slr123 = lmp(Sum~P3, data=BMI)
summary(slr123)

slr124 = lmp(Sum~P4, data=BMI)
summary(slr124)

slr125 = lmp(Sum~P5, data=BMI)
summary(slr125)

slr126 = lmp(Sum~P6, data=BMI)
summary(slr126)

slr127 = lmp(Sum~P7, data=BMI)
summary(slr127)

slr128 = lmp(Sum~P8, data=BMI)
summary(slr128)

slr129 = lmp(Sum~P9, data=BMI)
summary(slr129)


slr130 = lmp(Sum~P10, data=BMI)
summary(slr130)

slr131 = lmp(Sum~P11, data=BMI)
summary(slr131)

slr134 = lmp(Sum~EvenM, data=BMI)
summary(slr134)
#Coefficients:
#  Estimate Iter Pr(Prob)    
#EvenM   -10831 5000   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3788 on 18 degrees of freedom
#Multiple R-Squared: 0.3063,	Adjusted R-squared: 0.2678 
#F-statistic: 7.948 on 1 and 18 DF,  p-value: 0.01136 

slr135 = lmp(Sum~Div, data=BMI)
summary(slr135)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#Div    -4332 5000    0.003 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3946 on 18 degrees of freedom
#Multiple R-Squared: 0.2472,	Adjusted R-squared: 0.2054 
#F-statistic: 5.911 on 1 and 18 DF,  p-value: 0.02572 

slr136 = lmp(Sum~SimpsonM, data=BMI)
summary(slr136)
#Coefficients:
#  Estimate Iter Pr(Prob)    
#SimpsonM   -10174 5000    6e-04 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3838 on 18 degrees of freedom
#Multiple R-Squared: 0.2878,	Adjusted R-squared: 0.2483 
#F-statistic: 7.275 on 1 and 18 DF,  p-value: 0.01474 


slr139 = lmp(Sum~RichF, data=BMI)
summary(slr139)

slr140 = lmp(Sum~SimpF, data=BMI)
summary(slr140)

slr140.1 = lmp(Sum~ShanF, data=BMI)
summary(slr140.1)

###################################

slr140.2 = lmp(Div~Rhytisma, data=BMI)
summary(slr140.2)
#Coefficients:
#  Estimate Iter Pr(Prob)  
#Rhytisma  -0.9964 1504   0.0625 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4636 on 18 degrees of freedom
#Multiple R-Squared: 0.2112,	Adjusted R-squared: 0.1674 
#F-statistic:  4.82 on 1 and 18 DF,  p-value: 0.04148 

slr140.3 = lmp(Rich~Rhytisma, data=BMI)
summary(slr140.3)

slr140.4 = lmp(Sum~Rhytisma, data=BMI)
summary(slr140.4)

slr140.5 = lmp(k~Rhytisma, data=BMI)
summary(slr140.5)
#Coefficients:
#  Estimate Iter Pr(Prob)   
#Rhytisma -0.05786 5000   0.0056 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.021 on 18 degrees of freedom
#Multiple R-Squared: 0.3056,	Adjusted R-squared: 0.267 
#F-statistic: 7.921 on 1 and 18 DF,  p-value: 0.01147 

###################################