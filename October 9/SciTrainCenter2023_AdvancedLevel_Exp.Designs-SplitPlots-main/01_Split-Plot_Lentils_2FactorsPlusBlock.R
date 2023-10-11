###############################################################################################
###                       Scientific Training for Plant Biotechnology                       ###
###                     Course in Modern Plant Breeding - Advanced Level                    ###
###                                14-25 March 2022, Skoltech                               ###        
###                                                                                         ### 
###                                       ----------                                        ###
###                                   Day1: Split Plot                                      ###
###                                                                                         ###
###             Case Study #1 : AUDPC for Ascochyta blight evaluation in lentil             ### 
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Dr Cécile Ben, Skoltech           ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################
## ascochyta blight due to Ascochyta lentis.
## 6 new strains of A. lentis assessed on 3 reference cultivar of Lentil  
## AUDPC values 

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the datafile and the associated R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed, no need to re-install them again when you close-open again RStudio.

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# use of the constraint 'sum-to-zero' for ANOVAs
options(contrasts=c('contr.sum','contr.poly'))

#############################################################################################################


###################
# Data import and modification in R
###################

rm(list = ls())
gc()

library(tidyverse)
library(dplyr)
#library(lme4)
library(lmerTest)
library(multcomp) ## for cld()

Lentils <- read.table("./01_Split-plot_Lentils.csv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
str(Lentils)
attach(Lentils)

# the physical whole-plots (the trays) are not identified in the raw dataset
# create them, using the physical design

Lentils <- Lentils %>%
    unite(Tray, Block, Strain,
          sep = ".", remove = FALSE) %>%
    mutate_if(is.character, as.factor) #%>%      ## I agree, not obvious to remember ....
#    select(Strain, Genotype, Block, AUDPC, Tray)

str(Lentils)
Lentils

## explore the design
##
ftable(xtabs(~ Strain + Block + Genotype, data = Lentils))
## looks like a perfectly balanced design


########################
# Graphic visualizations
########################

#Interaction plots to get a more comprehensive view of the effect of the interaction of 2 factors on an agronomical trait of interest.
x11()
interaction.plot(Strain,Genotype,AUDPC,
                 type="b",
                 xlab="Strains",ylab="AUDPC",trace.label="Genotypes",
                 pch=c(1:10),cex=2, lty=1,lwd=2,col=c(1:10))


## standard MLM, w/o p-values for F-tests
ana1 <- lme4::lmer( log10(AUDPC/100) ~ Strain * Genotype +
                        (1|Block) + (1|Tray), data = Lentils)
anova(ana1)
rand(ana1)

## MLM with F-test for fixed effects
ana2 <- lmerTest::lmer( log10(AUDPC/100) ~ Strain * Genotype +
                        (1|Block) + (1|Tray), data = Lentils)
anova(ana2)
lmerTest::rand(ana2)

## to get "old-school" decompositions of the sum-of-squares
anaFixe <- lm(log10(AUDPC/100) ~ Strain * Genotype + Block + Tray %in% Strain:Block, data = Lentils)

## WARNING : F-tests of the below table are WRONG. Use it only for source of variations and ddf
anova(anaFixe)  # attention les F-test de cette table sont faux



## OK, there no interaction, so we can reduce the model
## the lmerTest::step() is quit useful for that :
mystep <- lmerTest::step(ana2)
mystep
final <- get_model(mystep)
anova(final)

## evaluate each random effect :
lmerTest::rand(final)

#
# adujested means for AUDPC :
library(emmeans)

leastsquare1 <- lsmeans(ana2,
                       pairwise ~ Strain,
                       adjust = "tukey")  ### Tukey-adjusted comparisons
cld(leastsquare1,
    alpha = 0.05,
    Letters = letters,
    adjust = "tukey")
### Use lower-case letters for .group


leastsquare2 <- lsmeans(ana2,
                       pairwise ~ Genotype,
                       adjust = "tukey")  ### Tukey-adjusted comparisons
cld(leastsquare2,
    alpha = 0.05,
    Letters = letters,
    adjust = "tukey")
### Use lower-case letters for .group

## visualisation (not really meaningful, but why not)
x11()
plot(lmerTest::difflsmeans(ana2))
