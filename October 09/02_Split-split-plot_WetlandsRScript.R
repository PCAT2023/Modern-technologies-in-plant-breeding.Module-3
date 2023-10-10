###############################################################################################
###                       Scientific Training for Plant Biotechnology                       ###
###                     Course in Modern Plant Breeding - Advanced Level                    ###
###                                14-25 March 2022, Skoltech                               ###        
###                                                                                         ### 
###                                       ----------                                        ###
###                                                                                         ###
###                 - Day 1: Split-Split-Plot: Weed biomass in wetlands. -                  ### 
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Dr Cécile Ben, Skoltech           ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################
#Context: 
#An experiment studies the effect of nitrogen and weeds on plant growth in wetlands. 
#Effect of four levels of nitrogen, three weed treatments (no additional weeds, addition of weed species 1, 
#addition of weed species 2), and two herbivory treatments (clipping and no clipping) were investigated. 
#Experimental design includes eight trays; each tray holds three artificial wetlands consisting of 
#rectangular wire baskets containing wetland soil. The trays are full of water, so the artificial wetlands stay wet. 
#All of the artificial wetlands receive a standard set of seeds to start growth. Four of the trays are placed on a 
#table near the door of the greenhouse, and the other four trays are placed on a table in the center of the greenhouse. 
#On each table, we randomly assign one of the trays to each of the four nitrogen treatments. Within each tray, we randomly 
#assign the wetlands to the three weed treatments. Each wetland is split in half. One half is chosen at random and will be 
#clipped after 4 weeks, with the clippings removed; the other half is not clipped. After 8 weeks, we measure the fraction of biomass 
#in each wetland that is non-weed as our response. (Datafile: SplitSplitPlot.csv).

#Questions:
#1. Describe the experimental design implemented for this study.
#2. Do nitrogen levels, weed and herbivory treatments affect plant growth?

######### PREPARATION OF TE WORKING INTERFACE IN R ##########################################################
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

## Loading of the R packages needed for the analysis.
library(lmerTest)
library(lsmeans)#For cld function
library(lattice) # Needed for some graphs (e.g. bwplots)

###################
# Data import in R
###################


weeds <- read.table("05_Split-split-plot_WetLands.csv", sep = ";", header = TRUE)
weeds <- within(weeds, {
                nitrogen <- factor(nitrogen); weed <- factor(weed);
                clipping <- factor(clipping); table <- factor(table);
                tray <- factor(tray); wetland <- factor(wetland)}
                )
str(weeds)
attach(weeds)

########################
# Graphic visualizations
########################
##Boxplots to reveal the distribution and variance of the measured traits depending on the different factors of interest

x11()
bwplot(biomass ~ nitrogen|weed+clipping)

###############
out1a <- lmerTest::lmer(biomass ~ nitrogen*weed*clipping +
            (1|table) + (1|tray) + (1|wetland), data = weeds)

# summary(out1)
anova(out1a)
VarCorr(out1a)
rand(out1a)



# other coding for random effects
out1b <- lmerTest::lmer(biomass ~ nitrogen*weed*clipping +
            (1|table/tray/wetland), data = weeds)

# summary(out1)
anova(out1b)
VarCorr(out1b)
rand(out1b)


## in summary Nonweed biomass decreases as nitrogen increases, but the decrease is much larger for the weed-seeded treatments.



## playing with factors : let's imagine nitrogen level is a random factor  
out2 <- lmerTest::lmer(biomass ~ weed*clipping
                    + (1|nitrogen) + (1|nitrogen:weed) + (1|nitrogen:clipping)  ## interactions with a random factor are random
                    + (1|nitrogen:weed:clipping)
                    + (1|table/tray/wetland),  
                    data = weeds)

# summary(ou2a)
anova(out2)
# Extract the variance components
vc <- VarCorr(out2)
vc
# tests for variances
lmerTest::rand(out2)
