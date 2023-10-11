# Carri J. LeRoy, 9-20-23
#
# 2016 Elwha paper: ANCOVA figure. Data in Elwha_AFDM.csv
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
#

# Rename data to "AFDMall"
AFDMall <- read.csv("Elwha_AFDM.csv", header = TRUE)
#AFDMall <- Elwha_AFDM
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

# Final ANCOVA figure for Elwha paper:
#
plotLn <- ggplot(data=AFDMall, aes(x=days, y=lnpAFDMR, group=region))+
  geom_point(aes(color=region))+
  geom_smooth(aes(color=region), method="lm", se=FALSE)+
  scale_linetype_manual(values=c("twodash", "solid", "dashed", "dotdash", "dotted"))+
  scale_color_manual(values=c("orange1", "yellow2", "green1", "blue1", "purple1"))+
  theme_classic()+
  theme(axis.text = element_text(colour="black"))+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 15))+
  theme(legend.position = c(0.2,0.2))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=12))+
  labs(x = "Days in stream", y = "ln % AFDM Remaining")+
  annotate(geom="text", size=5, x=29, y=3.0, label="Middle Elwha") +
  annotate(geom="text", size=5, x=35, y=4.7, label="Delta") +
  annotate(geom="text", x = 2, y = 2.4,
           label = "F(9,384) = 213.22, P < 0.0001
     Section: P < 0.0001
     Days: P < 0.0001
     Days*Section: p < 0.0001", adj=0, size=4) +
  geom_text(size=4, x=47.0, y=4.1, label="a") +
  geom_text(size=4, x=47.0, y=3.5, label="ab") +
  geom_text(size=4, x=47.0, y=3.3, label="bc") +
  geom_text(size=4, x=47.0, y=2.9, label="bc") +
  geom_text(size=4, x=47.0, y=2.6, label="c")
plotLn

ggsave(file="ANCOVA.pdf", plotR, width=16.4, height=12.4, units="cm", dpi=800)
ggsave(plotR, file="ANCOVA.eps", width=16.4, height=12.4, units="cm", dpi=800, device="eps")
ggsave(plotR, file="ANCOVA.jpg", width=16.4, height=12.4, units="cm", dpi=800, device="jpg")

# Final 2-panel ANCOVA figure for Elwha paper:
#
plotLn <- ggplot(data=AFDMall, aes(x=days, y=lnpAFDMR, group=region))+
  geom_point(aes(color=region))+
  geom_smooth(aes(color=region), method="lm", se=FALSE)+
  scale_linetype_manual(values=c("twodash", "solid", "dashed", "dotdash", "dotted"))+
  scale_color_manual(values=c("orange1", "yellow2", "green1", "blue1", "purple1"))+
  theme_classic()+
  theme(axis.text = element_text(colour="black"))+
  theme(axis.text.x = element_blank()) +
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 15))+
  theme(legend.position ="none") +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=12))+
  labs(x = element_blank(), y = "ln % AFDM Remaining", tag = "A")+
  annotate(geom="text", size=5, x=29, y=3.0, label="Middle Elwha") +
  annotate(geom="text", size=5, x=35, y=4.7, label="Delta") +
  annotate(geom="text", x = 2, y = 1.4,
           label = "F(9,384) = 213.22, P < 0.0001
     Section: P < 0.0001
     Days: P < 0.0001
     Days*Section: p < 0.0001", adj=0, size=4) +
  geom_text(size=4, x=47.0, y=4.1, label="a") +
  geom_text(size=4, x=47.0, y=3.5, label="ab") +
  geom_text(size=4, x=47.0, y=3.3, label="bc") +
  geom_text(size=4, x=47.0, y=2.9, label="bc") +
  geom_text(size=4, x=47.0, y=2.6, label="c")
plotLn

plotM <- ggplot(data=AFDMall, aes(x=days, y=AFDMR, group=region))+
  geom_point(aes(color=region))+
  scale_color_manual(values=c("orange1", "yellow2", "green1", "blue1", "purple1"))+
  theme_classic()+
  theme(axis.text = element_text(colour="black"))+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 15))+
  theme(legend.position = c(0.25,0.3))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=12))+
  labs(x = "Days in stream", y = "% AFDM Remaining", tag = "B")
plotM

g1 <- ggplotGrob(plotLn)
g2 <- ggplotGrob(plotM)
g <- grid.draw(rbind(ggplotGrob(plotLn), ggplotGrob(plotM)))

# Export pdf and jpg directly from viewing pane - can't get automatic saves to work...

