###############################################################################################
###                       Scientific Training for Plant Biotechnology                       ###
###                     Course in Modern Plant Breeding - Advanced Level                    ###
###                                14-25 March 2022, Skoltech                               ###        
###                                                                                         ### 
###                                       ----------                                        ###
###                                                                                         ###
### - Day 2: GXE interaction analysis: ANOVA & Joint (Linear) Regression Analysis (JRA) -   ###
###                                                                                         ###
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Pr Cécile Ben, Skoltech           ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################

#CASE STUDY 1. Analysis of the performance of different maize genotypes in different growing environments: Genotype X Environment (GXE) interaction analysis. (12 dots)

#Context:
#A seed company has developed 10 new F1 hybrid varieties of maize (Zea mays) 
#and wishes to know if some of these varieties are particularly well adapted to the growing conditions of 5 different production zones of the Corn Belt in the USA. 
#Which variety(ies) to develop and sell according to the environmental conditions of each production area? 
#This is the question the company wants to answer.

#To answer this question, a performance evaluation of the 10 maize varieties (namely Gen_1 to Gen_10) is conducted in 5 test stations (namely Env_1 to Env_5) characterized by different environmental conditions. 
#For each variety in each site, the grain yield (measured in quintals per hectare) was evaluated by measuring four samples taken randomly (Rep 1 to 4 for each environment).

#################################################################################################

rm(list = ls())
gc()

library(tidyverse)
library(lme4)
library(lmerTest)
library(pracma)
library(ggplot2)
library(multcomp)
library(lattice)

library(car)    # Levene test
library(agricolae)   # Newman-Keuls & Tukey test
library(ggplot2)
library(gridExtra)
library(dplyr)
library(MASS)#BoxCox

# to erase objects and reinitialize the workspace
rm(list = ls())
graphics.off()


# Constraints on ANOVA effects
options(contrasts=c("contr.sum", "contr.poly"))


### Import the data

Produc <- read.table("04_GxEwithBlocks.csv",sep = ";",dec = ".",header = TRUE, stringsAsFactors = TRUE)
str(Produc)

## modifications de code pour Repb, modifications de code pour Gen
Produc$Env <- factor(paste("Env", Produc$Env, sep = ''))
Produc$Gen <- factor(paste("Gen", Produc$Gen, sep = ''), levels = c('Gen1','Gen2','Gen3','Gen4','Gen5','Gen6','Gen7','Gen8','Gen9','Gen10')) 
#It is essential to indicate the order for Gen. 
#If not, discrepancies in the calculation of the slopes and coeff. of linear regression.
Produc$Rep <- factor(paste("Bloc", Produc$Rep, sep = ''))
Produc$Repb <- factor(paste("Block", Produc$Env, Produc$Rep, sep = '.'))  ## To indicate blocks nested within environments
str(Produc)

##Check for balanced dataset :
table(Produc$Env,Produc$Gen)
table(Produc$Gen,Produc$Env)
#Possible to test for effect of interactions between the 2 factors genotype and environment because for each combination "genotype x Env" data number =4 (>2)

############ Graphical visualisation

## 1. Basic graph: Boxplots considering only one factor of study
x11()
(graf1a <- ggplot(Produc, aes(x = Env, y = Prod, fill = Env)) +
    geom_boxplot())

x11()
(graf1b <- ggplot(Produc, aes(x = Gen, y = Prod, fill = Gen)) +
    geom_boxplot())

## 1. Basic graph: Boxplots considering both factors of study at the same time: much better and informative!
x11()
(graf1 <- ggplot(Produc, aes(x = Env, y = Prod, fill = Env)) +
    geom_boxplot() +
    facet_wrap( ~ Gen)
)


## 2. A variant: Violin plot
x11()
(graf1 <- ggplot(Produc, aes(x = Env, y = Prod, fill = Env)) +
    geom_violin() +
    facet_wrap( ~ Gen)
)


## Calculation  of means per site and per genotype, and measure numbers
Summaries <- Produc %>%
    group_by(Env, Gen) %>%
     summarise(avgProduc = mean(Prod, na.rm = TRUE),
               sdProduc = sd(Prod, na.rm = TRUE),
               nbData = n()
               ) %>%
     arrange(desc(avgProduc)) %>%
    print(n = Inf)  # to see all data


## Barplot classical but not very informative compared to boxplot
x11()
(graf2 <- ggplot( Summaries, aes( x = Env, y = avgProduc, fill = Env)) +
    geom_bar(stat = "identity") +
    facet_wrap( ~ Gen))


## un lollypop ?
x11()
ggplot(Summaries, aes(x = Env, y = avgProduc, color = Env)) + 
  geom_point(size = 5) + 
  geom_segment(aes(x = Env, 
                   xend = Env, 
                   y = 80,   ## astuce pour visus des diffÃ©rences
                   yend = avgProduc)) + 
  labs(title="Lollipop Chart", 
       subtitle="Average Production Vs Sites", 
       caption="source: personal data") + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6))  +
    facet_wrap( ~ Gen)


## classical for the study of the interaction between 2 factors
x11()
interaction.plot(Produc$Env,Produc$Gen,Produc$Prod,
                 type="b",
                 xlab="Lieux",ylab="Prod",trace.label="Genotypes",
                 pch=c(1:10),cex=2, lty=1,lwd=2)



###############  ANOVA 2 factors with interaction
###############  old school : GLM with aov function

## Model without blocks.
## all factors as fixed
modele1 <- aov( Prod ~ Gen * Env , data = Produc)
summary(modele1)

## Model with blocks nested in Environments
## all factors as fixed
modele2 <- aov( Prod ~ Rep %in% Env + Gen * Env , data = Produc)
summary(modele2)
## Why 15df for blocks within Env ?

## Model with blocks nested in Environments. Enumerator coding for the blocks
## all factors as fixed
modele3 <- aov( Prod ~ Gen * Env + Repb, data = Produc)
summary(modele3)
#Due to dependence between col. if I sum the Repb columns of the environment, we recreate the col of the environment because these columns are only in this environment.
#So R by calculating the design matrix systematically drops a col. Repb by enviro  . Therefore the df for Block within Env are the same in the model2 and modele3
#model.matrix(modele3)

#Test for ANOVA pre-requisites
#independence of the residuals
x11()
par(mfrow=c(2,2))
plot(modele3)

# Normality of the residuals
shapiro.test(modele3$residuals)  # in fact, three somewhat 'strong' values

# Homogeneity of variances
## attention, in this plan with nested blocks, it is necessary to "forget" a factor, otherwise no test possible
summary(aov(abs(modele3$residuals) ~ Produc$Gen:Produc$Env:Produc$Repb))

## so if one consider Rep not a pb
summary(aov(abs(modele3$residuals) ~ Produc$Gen:Produc$Env))   ## ouch, not so good!

############ Data transformation
bc <- boxcox(modele3)
x11()
boxcox(modele3)
trans <-bc$x[which.max(bc$y)]
trans

modele2trans <- aov( Prod^trans ~ Gen * Env + Repb, data = Produc)
summary(modele2trans)
# The Shapiro-Wilk test for normality is produced by:
shapiro.test(resid(modele2trans))
## Variances homogeneity for the residuals
leveneTest(Prod^trans ~ Gen * Env, data= Produc)

################# Multiple mean comparisons for the GXE interaction
##

# Assuming significant interaction between genotype and environments, we can use familiar code but 
# we may look at the output marked $gen:env for comparing ALL factor level combinations using Tukey:                                                
TukeyHSD(modele2trans, conf.level=0.95)

# The table shows the P-values for comparison among all factor level pairs.    
# These P-values below each t statistic are adjusted for the Tukey procedure.  
# conf.level=0.95 is actually the default level.  We could choose  
# another level if desired.  

######## a bit simpler, using the function agricolae::SNK.test (ou agricolae::HSD.test)
( SNK.test(modele2trans, c('Gen','Env') ) )
x11()
plot( SNK.test(modele2trans, c('Gen','Env') ) )

( HSD.test(modele2trans, c('Gen','Env') ) )
x11()
plot( HSD.test(modele2trans, c('Gen','Env') ) )

###############  ANOVA : 2 factors without interaction - to illustrate the matrix model of LMMs
###############  new school method: use of lmer function for Linear Mixed Models

## lme4::lmer allows to fit cross random effect

library(lme4)

ana1_MLM <- lme4::lmer( Prod ~ Gen * Env + (1 | Repb), data=Produc)  ## REML true par défaut: In statistics, the restricted (or residual, or reduced) maximum likelihood (REML) approach is a particular form of maximum likelihood estimation that does not base estimates on a maximum likelihood fit of all the information, but instead uses a likelihood function calculated from a transformed set of data, so that nuisance parameters have no effect.
summary(ana1_MLM)
anova(ana1_MLM)
#rand(ana1_MLM)

fixef(ana1_MLM)
ranef(ana1_MLM)

#Quality of prediction
x11()
plot(ana1_MLM)

##### Formal Tests of ANOVA Model Assumptions
# The Shapiro-Wilk test for normality of ANOVA residuals is produced by:
shapiro.test(resid(ana1_MLM))

###############################
library(lmerTest) ### To get the p-values

ana2_MLM <- lmerTest::lmer( Prod ~ Gen * Env + (1 | Repb), data=Produc)
summary(ana2_MLM)
anova(ana2_MLM)

fixef(ana2_MLM)  
ranef(ana2_MLM)  ## levels are assigned by alphabetic order, not in the order of the data. 

####Multiple mean comparison test
## Adjusted means (because the design is balanced, adjusted means = means from data)
library(lsmeans)
AdjustMoys1 <- lsmeans::lsmeans(ana1_MLM,   ## lmerTest is masking standard lsmeans. Unwanted !
                                pairwise ~ Gen:Env, 
                                adjust = "tukey")

## Multiple comparisons
(CompMoys1 <- cld(AdjustMoys1,  ## cld for compact letter display
                  alpha = 0.01 ,
                  Letters = letters,
                  adjust = "tukey"
)
)


################
## To go further ......
################ regression on environmental index: average of each location
##
## the idea is to no longer consider "Env" as a qualitative factor, but to transform it into a quantitative variable.
## we will use the environmental index, defined as being the averages of production among the 10 Genotypes of the different environments.

Index <- aggregate(Prod ~ Env , FUN = "mean", data = Produc)
str(Index)
Index 
Index[order(Index$Prod),]  ## by increasing environmental indices.

Produc$Index <- rep(Index$Prod, each= 4*10)  # pay attention to how data is entered
Produc

#the slope and the y-intercept can be different between the different lines
JointReg <- lm(Prod ~ Gen + Index + Gen:Index ,  data = Produc) 
# the production is a function of the Genotype (Gen), the general effect of the environment (Index) and the GxE interaction (Gen:Index)
summary(JointReg)

Analyses <- data.frame(Produc, Fitted = JointReg$fitted.values)  

# Find two predicted values by hand:
# genotype 1 in Env 1
# genotype 2 in Env 3


# pay attention to the interpretation of the tests on the coefficients.
# here, the test is "is the coeff different from zero?"
# BUT, for the Genotype x Env coefficients, this amounts to knowing whether or not the slope is equal to the "general" coefficient of the Env (here 1.000)
## *** WARNING ***: here "sum-to-zero" constraint (see beginning of the script) so the coefficient of the last genotype is - (sum of the coeffs of the other genotypes)
## (and therefore no edited test for this coeff)
# does the regression model describe the data well?

# slopes of lines
## from the coefficients of lm, we see that Gen10 was sum-to-zero (Sum-to-zero is used as constraints. See line 53: options() )
### *** BE CAUTIOUS *** and check which level of environment was set-to-zero. This sometimes depends on internal Windows configuration. ---> This may require to modify the below code
Pentes <- data.frame(Genotype=paste('Gen',seq(1:10),sep=""),
                     Pente=c( JointReg$coeff[11] + JointReg$coeff[12:20],JointReg$coeff[11]-sum(JointReg$coeff[12:20]))
)

Pentes[order(Pentes$Pente),]
# most stable? difference between stability and proportional response?


# colors
palette(topo.colors(10))

x11()
xyplot(Fitted ~ Index, group=Gen, data=Analyses,
       # what are we going to display? the "graphic panel"
       panel=function(x,y,...) {
         panel.xyplot(x, y,type="b",...)          # predicted data
         panel.abline(0,1,lty=2,col="darkgreen")  # bisector of slope 1
       },
       # graphic parameters for the panel
       par.settings = list(superpose.symbol=list(cex=1.5,pch=c(1:10),col=c(1:10)),   # dots
                           superpose.line = list(lwd = 2,col=c(1:10))                # lines
       ),
       # caption
       auto.key = list(text = levels(Analyses$Gen),
                       space = "right",
                       points = TRUE, lines=TRUE
       )
)


### or : 
x11()
xyplot(Fitted ~ Index | Gen, data=Analyses,
       # what are we going to display? the "graphic panel"
       panel=function(x,y,...) {
         panel.xyplot(x, y,type="b",...)          # predicted data
         panel.abline(0,1,lty=2,col="darkgreen")  # # bisector of slope 1
       },
       # graphic parameters for the panel
       par.settings = list(plot.symbol=list(cex=1.5,pch=19,col=c(1:10)),        # dots
                           plot.line = list(lwd = 2,col=c(1:10))                # lines
       )
)




### raw data
x11()
xyplot(Prod ~ Index, group=Gen, data=Analyses,
       type="p",
       
       par.settings = list(superpose.symbol=list(cex=1,pch=c(1:10),col=c(1:10)),
                           superpose.line = list(lwd = 2,col=c(1:10))
       ),
       
       auto.key = list(text = levels(Analyses$Gen),
                       space = "right",
                       points = TRUE, lines = TRUE
       )
)













