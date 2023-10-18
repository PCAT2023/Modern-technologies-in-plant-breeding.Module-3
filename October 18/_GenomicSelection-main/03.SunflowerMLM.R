#' ---
#' title:  |
#'   | A rapid introduction to MLM  
#'   | and getting started with prediction of genetic value
#' author:
#'   - Prof C. Ben, Project Center for Agro Technologies, Skoltech ^[c.ben@skoltech.ru]
#'   - Prof L. Gentzbittel, Project Center for Agro Technologies, Skoltech ^[l.gentzbittel@skoltech.ru]
#' date: "14-25 March 2022, Skoltech - Scientific Training for Plant Biotechnology, Course in Modern Plant Breeding - Advanced Level "
#' output: 
#'   rmdformats::readthedown:
#'     highlight: haddock
#'   pdf_document:
#'     keep_tex: true
#' use_bookdown: TRUE
#' latex_engine: xelatex
#' header-includes:
#'   - \usepackage{bbold}
#'   - \def\+#1{\mathbf{#1}}
#' geometry: left = 2cm, right = 1.5cm, top = 1.5cm, bottom = 1.5cm
#' ---
#'  
#' 
#+ echo = FALSE, message = FALSE, warning = FALSE
# just forgot the below lines. They are used to generate the printed version
knitr::opts_chunk$set(fig.width = 4.5, fig.height = 4.5, fig.align = "center",
                      warning = FALSE, message = FALSE,
                      tidy.opts = list(width.cutoff = 60), tidy = FALSE
)

# just forgot the below lines. They are used to generate the printed version OR the html version
# library(rmarkdown)
# render('03.SunflowerMLM.R',output_format = 'pdf_document')
# render('03.SunflowerMLM.R',output_format = 'html_document')

#+ echo = TRUE, message = TRUE, warning = FALSE 
###############################################################################################
###                       Scientific Training for Plant Biotechnology                       ###
###                     Course in Modern Plant Breeding - Advanced Level                    ###
###                                14-25 March 2022, Skoltech                               ###        
###                                                                                         ### 
###                                       ----------                                        ###
###                                                                                         ###
###     - Day 6: Prediction of genetic value of sunflower cultivars for yield -             ###
###                                                                                         ###
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Dr C�cile Ben, Skoltech           ###
###############################################################################################

######### CASE STUDY PRESENTATION #############################################################
### Prediction of genetic value of sunflower cultivars for yield 
### based on the phenotypic results of a Multi-Environment Trial (MET)
# 6 sunflower genotypes -"NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"- 
# grown in 5 environments near Toulouse (France) -"Odars", "Verfeil", "Villefranche", "Fourquevaux"- 
# and evaluated in a VERY UNBALANCED design.
# during 3 years -2003, 2005, 2007 -


######### PREPARATION OF THE WORKING INTERFACE IN R ######################################
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed, NO need to re-install them again 
# when you close and open again RStudio.

### III. Initialisation of the working space:clean workspace, clean memory 
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# garbage collector. To use to free memory
gc() 


## indicate here the folder where are located the data
## and subdirs where outputs will be stored
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(paste(main_dir, "./data3/", sep = ""))
# 
################################################################################

## Loading of the R packages needed for the analysis.
library(lme4) 
library(ggplot2)

###################
# Data import in R
###################
#########
## EX. 1: Computation of BLUEs & BLUPs to predict Yield using environnement 
## and Genotype
#########

sunflower <- read.table("SunflowerWeed.csv", sep = ";", header = TRUE, dec = ".", stringsAsFactors = TRUE)
sunflower #Dataset 'inspired' from Roumet et al., New Phytol, 2012
str(sunflower)
#6 genotypes of cultivated sunflowers, 5 environments, 3 years
sunflower$YEAR <- as.factor(sunflower$YEAR)
str(sunflower)

#' A very unbalanced dataset !
with(sunflower,
table(Geno, ENV, YEAR)
)

#' # Building a MLM 'by hand'

##Vector of observations--------------------------------------------------------
y <- sunflower$YLD
y #Dim (12,1)

#' ## Design matrix of the environments as fixed effects
#Let's create the X matrix by hand.

Xmanual <- matrix(c(1,1,0,0,0,   #µ + Effect of experimental field (i.e Odars) for Measure1 
                    1,1,0,0,0,	 #µ + Effect of experimental field (i.e Odars) for Measure2 
                    1,0,1,0,0,	 #µ + Effect of experimental field (i.e. Verfeil) for Measure3 
                    1,0,1,0,0,   #µ + Effect of experimental field (i.e. Verfeil) for Measure4 
                    1,0,0,1,0,   #µ + Effect of experimental field (i.e. Villefranche)for Measure5 
                    1,0,0,1,0,
                    1,0,0,1,0,	
                    1,0,0,0,1,
                    1,0,0,0,1,
                    1,0,0,0,1,
                    1,0,0,0,0,
                    1,0,0,0,0),
                  byrow = TRUE, nrow = 12, ncol = 5)
dimnames(Xmanual) <- list(c("Measure1","Measure2", "Measure3", "Measure4", "Measure5", "Measure6", "Measure7", "Measure8", "Measure9", "Measure10", "Measure11", "Measure12"),
                        c("Intercept", "Odars", "Verfeil", "Villefranche", "Fourquevaux"))
Xmanual

# another way to create the X matrix is :
X <- model.matrix( ~ sunflower$ENV)
# or, considering sites and years:
X_env_year <- model.matrix( ~ sunflower$ENV + sunflower$YEAR )

#Rq: In this dataset, no blocks, no replicates by enviro and year so no possibility to take
# into account the Site X Year interaction.
# If it was the case: X_envXyear <- model.matrix( ~ sunflower$ENV * sunflower$YEAR)

#' We new need to compute the 'effects' of the environments *i.e.* the numerical values of the effect of the 
#' different years and the different sites. It exists 2 solutions to solve this linear systems

# Solution 1: ANOVA Constraint "sum to 0" 
# alpha Env1 + alpha Env2 + alpha Env3 + ... = 0.

# Solution 2: ANOVA Constraint "set to 0"
# the effect of the 1st level by alphabetical order (here, `Balma`) is set to 0.
# In this case, the value of this environment = Intercept.
#' The X matrix created with the function `model.matrix` contains (Env-1) columns: Balma column is not present, 
#' The function `model.matrix` order the Factor level by alphabetical order and remove the 1st level (here, `Balma` )
X

EnvBetas <- solve(( t(X) %*% X)) %*% t(X) %*% y
EnvBetas
# Effect of Balma = 0 due to constraint "set-to-zero"
# Estimated mean = 1479.50000
# Balma estimated mean, with this model = 1479.5 + 0 = 1479.5

Env_yearBetas <- solve((t(X_env_year) %*% X_env_year)) %*% t(X_env_year) %*% y
Env_yearBetas
# Beta_Balma = 0 and Beta_2003 = 0 due to constraint "set-to-zero"
# Balma Estimated mean, with this model  = 1534.75 + 0 = 1534.75
# 2003 Estimated mean = intercept = 1534.75 + 0 = 1534.75

# If Site X Year interaction is considered in the model then intercept = Estimated mean for Balma_2003.

#' NB: Decomposition of the solving of this linear system:
# solve(t(X) %*% X)  #  Compute the inverse of ( t(X) %*% X )
# solve(t(X)%*%X) %*% (t(X)%*%X)  #  Identity Matrix proving what is computed
#t(X)%*%y

#' ## Mixed Linear Model with 'Genotypes' as random effects.
# Creation of Z matrix of random effects and of µ vector of BLUP
sunflower
Zmanual<- matrix(c(1,0,0,0,0,0,   #measure.1=u1
                   0,1,0,0,0,0,  	#measure.2=u2
                   1,0,0,0,0,0,		#measure.3=u1
                   0,1,0,0,0,0,		#measure.4=u2
                   0,1,0,0,0,0,
                   0,0,1,0,0,0, 
                   0,0,0,1,0,0,
                   1,0,0,0,0,0,	
                   0,0,0,0,1,0,
                   0,0,0,0,0,1,
                   0,1,0,0,0,0,
                   0,1,0,0,0,0
), byrow = TRUE, nrow = 12, ncol = 6)
dimnames(Zmanual) <- list(c("Measure1","Measure2", "Measure3", "Measure4", "Measure5", "Measure6", "Measure7", "Measure8", "Measure9", "Measure10", "Measure11", "Measure12"),
                        c("NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"))

Zmanual

## alternatively, one can compute this matrix using: --------------------------
Z <- model.matrix( ~ sunflower$Geno - 1)
# By default, model.matrix add an intercept. However, for random effect, there 
# are no intercept.  
# By default, Genotypes in alphabetical order in Z matrix 
Z

#' ## Calculation of the genetic values as BLUP from centered values of phenotypes 


#' We want to fit the following model for BLUP $y2 = mean(y) + Zu + E$.  
#' In this case, we do not account for the **fixed** environment effects.

#' Solving the system for random effects needs to take into account the `VarCov` matrix of the random effects (i.e. that will be later the pedigree matrix). 
#' In this first analysis, we do not have information about the relastionships between the genotypes. We will use the identity matrix. I = diag(6)

#' We will use $y2 = y - mean(y)$ as phenotypic values
#y2 <- (sunflower$YLD) - mean(sunflower$YLD) #Obs. Mean = 1502.333

# Solution of the linear model:
#uhatSimple <- solve(t(Z) %*% Z + 1*diag(6)) %*% t(Z) %*% y2
# 1*diag(6)=lambdaG with sigmaE/sigmaU arbitrarily set to 1
#uhatSimple


modellmer0 <- lmer(sunflower$YLD ~ 1 + (1|sunflower$Geno))
## Here are the genetic values with the MLM
ranef(modellmer0)

##Rq: Numerical diff. in the random effects (genetic values) computed manually and using lmer. This is because we compute the mean from observed data and lmer estimate the mean as a BLUE.

#' # Fitting a M.E.T. with a mixed model.

#' MultiEnvironment Trials: Estimation of genetic values but NOT of the breeding values (additive)
# ------------------------------------------------------------------------
# Mixed model: Fixed effect for Environnements & random effects for genotypes

## lme4::lmer()
modellmer <- lmer(sunflower$YLD ~ sunflower$ENV + sunflower$YEAR + (1|sunflower$Geno))
#' By default, lmer uses `G = SigmaU*I` and  `R = SigmaE*I`  
 
modellmer 
# Compare the BLUE for environments from this model and the effects of 
# environments computed with fixed effects only:
Env_yearBetas

## Here are the genetic values with the MLM
ranef(modellmer) 
# Compare the BLUP for genotypes from this model and the genetic values computed 
# with random effects only w/o accounting for environments:
#uhatSimple
ranef(modellmer0)

#' A very unbalanced dataset !
with(sunflower,
     table(Geno)
)

#To visualize the X design matrix and Z random effect matrix
getME(modellmer,'X')
t(as.matrix(getME(modellmer,'Zt')))#unexpected complex syntax...Sorry...


## Caution: If unbalanced design as here, the environment effects and genotype
## effects are not independent, so this does not provide an accurate 
## Broad sense heritability.
#' *Only if* the genotypes are pure lines and the design is balanced  we get narrow-sense heritability.

##estimate of broad-sense heritability = square of accuracy
as.data.frame(VarCorr(modellmer))#vcov for heritability calculation
# (i.e. VarG & Var E) are in column 4
(varianceestimates <- as.data.frame(VarCorr(modellmer))[,4])
h2hat <- varianceestimates[1] / sum(varianceestimates)
h2hat


#' # Estimating genetic effects, taking into account relationships between varieties

# Reminder: In the previous example, we do not know anything about the relationships between the varieties
# So we used the identity matrix I in the model.

#' ## Fitting a random model using pedigrees.
## Model P-BLUP (P for pedigree) based on a (fake) pedigree matrix

Amat <- matrix(c(1,0.5,0.25,0,0.25,0,
                0.5,1,0,0.25,0.5,0.25,
                0.25,0,1,0.5,0.5,0,
                0,0.25,0.5,1,0,0.25,
                0.25,0.5,0.5,0,1,0.5,
                0,0.25,0,0.25,0.5,1
), byrow = TRUE, nrow = 6, ncol = 6)
dimnames(Amat) <- list(c("NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"),
                     c("NK Countri", "PR64H41", "Atomic", "Santiago", "Boogy", "Aurasol"))
Amat

# Clone: 1; Full sibs: 0.5; Half sibs: 0.25; 1st cousins: 1/8; 
# Parent-offsprings: 0.5; GrandParents/Grandchild: 0.25; 
# Aunt(Uncle)-Nieces(nephew): 0.25

breedingvalue.withA <- solve(t(Z) %*% Z + 1*solve(Amat)) %*% t(Z) %*% y
breedingvalue.withA
# Compare the BLUP for genotypes from this model and the ones computed with 
# random effect model only w/o effect of environments
#uhatSimple
ranef(modellmer0)
# Compare the BLUP for genotypes from this model and the ones computed with
# the MLM accounting for environment BUT w/o knowing the genetic 
# relationships 
ranef(modellmer) 
# there are differences !


#' # The beauty of MLM: Predicting un-observed performances using full information.
#' Some of the above models allow to compute the performances of all varieties at all sites for all years. 
#' But the 'quality' of this prediction can be dramatically different depending on the model used.
#'
#' ## The prediction of observed conditions.  
#'  
#' fixed model w/o environment effect:
GLM0 <- aov( YLD ~ Geno, data = sunflower)
summary(GLM0)
#' fixed model with environment effect:
GLM <- aov( YLD ~ YEAR + ENV + Geno, data = sunflower)
summary(GLM)  ################  WARNING ! ! !


#' The MLM with full model for environment - but not accounting for genetic relationships - is computed above: `modellmer`.

sunflower$predictGLM0 <- predict(GLM0)
sunflower$predictGLM <- predict(GLM)
sunflower$predictMLM <- predict(modellmer)

#' Here are the predictions:
#' sunflower
#' 
#' #' a criteria of 'quality' for is the value of the residuals' sum-of-square. 
#' #' A `residual` is the difference between the **observed** value and the **predicted** value.
#' ( RSS.GLM0 <- t(residuals(GLM0)) %*%  residuals(GLM0) )
#' ( RSS.GLM <- t(residuals(GLM)) %*%  residuals(GLM) )
#' 
#' blirf <- sunflower$YLD - sunflower$predictMLM
#' ( RSS.MLM <- t(blirf) %*%  blirf )
#' 
#' ## The prediction of un-observed conditions.  
#' First create a design with ALL combinations
( ALLCOMB <- expand.grid( Site = levels(sunflower$ENV), 
                        Year = levels(sunflower$YEAR), 
                        Variety = levels(sunflower$Geno) )
)

#' Case 1: Let's use the model GLM0 to predict all combinations.

( DesignGLM0 <- model.matrix(GLM0) )
# estimation of coefficients
GLM0Betas <- solve( t(DesignGLM0) %*% DesignGLM0 ) %*% t(DesignGLM0) %*% y
GLM0Betas
# compute predictions for the observed conditions
( DesignGLM0 %*% GLM0Betas )
## compare with result of predict()function of GLM0 model, rounding number 
# to 6digits
all(round((DesignGLM0 %*% GLM0Betas), 6) == round(predict(GLM0), 6))


#' Now use this model to predict all cases. This is simply done by 
#' using some simple matrix algebra. 
( X.ALLCOMB.GLM0 <- model.matrix( ~ ALLCOMB$Variety ) )

ALLCOMB$PredictGLM0 <- X.ALLCOMB.GLM0 %*% GLM0Betas
ALLCOMB

#'
#' Case 2: Let's use the MLM to predict all combinations.
#' 
X.ALLCOMB <- model.matrix( ~ ALLCOMB$Site + ALLCOMB$Year )
Z.ALLCOMB <- model.matrix( ~ ALLCOMB$Variety - 1 )

#' We will use the coefficients of `modellmer` to estimates the components of the model
# Xb 
Fix.BLUEs <- X.ALLCOMB %*% fixef(modellmer)
# Fix.BLUEs
( cbind(ALLCOMB[, 1:3], Fix.BLUEs) )
# Zu
Ran.BLUPs <- Z.ALLCOMB %*% unlist(ranef(modellmer)) ## sorry for cryptic syntax!
# Ran.BLUPs
( cbind(ALLCOMB[, 1:3], Ran.BLUPs) )

## 'THE' full prediction:  Xb + Zu
ALLCOMB$PredictMLM <- Fix.BLUEs + Ran.BLUPs 
#' the 'Genetic Values" of the varieties - 'corrected' from environment effect
ALLCOMB$GeneticValues <- Ran.BLUPs + fixef(modellmer)[1]  

ALLCOMB
#Compare with the predicted values obtained for the observed conditions.
sunflower
#' 
#' Compare the genetic values with the performances computed using a naive approach:

x11()
( gr1 <- ggplot(ALLCOMB) + aes(y = GeneticValues, x = PredictGLM0, label = Variety) +
    geom_point() + geom_text(hjust = 0, vjust = 0) + 
    geom_abline(slope = 1, color = "red") +
    xlim(1460, 1580) + ylim(1460, 1580) +
    ggtitle("Genetic values and averaged phenotypic values")
  )


