## 
##
## scripts by Prof. Laurent Gentzbittel, Prof. Cecile Ben, Skoltech
##
## Visualising  population structure from SNP data
##
##




# clean workspace, clean memory 
rm(list = ls())
graphics.off()
gc()  # garbage collector. To use to free memory

# This is a trick to detect which folder contains the R script and the data
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)

library(adegenet)
library(ggplot2)


## -------------------  Using multivariate methods --- ------------------- ##

# data were formated using plink1.9 in plink's raw format
#
## This i a light dataset, for fast computations
SNPs <- read.PLINK("./data/KrAllSubset16K.raw", map.file = "./data/KrAllSubset16K.map", chunkSize = 50)
# chunkSize is the number of genomes to be read at a time.

SNPs


# additional information on accessions 
# shall be in the SAME order as in SNPs
Infos <- read.table("./data/AccessionsGeography.csv", sep=";", header=TRUE)

# geo-cordinates (lat/long) in the 'other' slot of SNPs
other(SNPs)$xy <- Infos[,c("Long", "Lat")]

SNPs
SNPs@other$xy
    

# faster to reload than generating the object again :
save(SNPs, file = "./data/SNPs.RData")
## # if needed reload SNPs object by :
## load("../data/SNPs.RData")


## ------------------ principal component Analysis ------------------ ##

### PCA , useful method to reduce the  dimension of dataset
ptm <- proc.time()

PCA1 <- glPca(SNPs, center = TRUE, scale = FALSE, nf = 200, loadings = TRUE, 
              alleleAsUnit = FALSE, useC = FALSE, parallel = TRUE,
              n.cores = 1, returnDotProd = FALSE, matDotProd = NULL)

(Temps.glPca <- proc.time() - ptm)
##

## save PCA analysis
save(PCA1, file = "./data/PCA1_glPca.RData")
## reload using
## load("../output/PCA1_glPca.RData")

## save ALL workspace before continuying, in case of problems
save.image(file = "./data/BioGeographyStep1.RData")


## ------------------ K-means ------------------ ##

## Kmeans :
grp1 <- find.clusters(SNPs, max.n.clust = 20, n.pca = 200,
                      choose.n.clust = FALSE, criterion = "min",
                      n.start = 120, glPca = PCA1)
grp1



## ------------------ Discriminant Analysis of principal components  DAPC ------------------ ##


## We want to assess the relationships between these groups using DAPC.
## Perform the DAPC and store the results in a new object called dapc:

dapc1 <- dapc(SNPs, pop = grp1$grp, scale = FALSE, n.pca = 200, n.da=8,
              glPca = PCA1)


save( dapc1, file= "./data/dapc1_ana.RData")
## load("../output/dapc1_ana.RData")
save.image(file = "./data/BioGeographyStep2.RData")
## load("../output/dapc1_ana.RData")


x11()
scatter(dapc1, scree.da = FALSE, bg="white",
        pch=20, cell=0, cstar=0, col=funky(5), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:5))


x11()
scatter(dapc1, xax=2, yax = 3,
        scree.da = FALSE, bg="white",
        pch=20, cell=0, cstar=0, col=funky(5), solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:5))


## Storing information about assignees to groups:
tbl <- dapc1$posterior
tbl <- as.data.frame(tbl)
# major group for each accession
maxi <- function(x){which(x==max(x))}
tbl$Major <- apply( tbl, 1, maxi)
tbl$Acc <- rownames(dapc1$tab)


save( tbl, file= "./data/dapc_tbl.RData")



## ------------------ Mapping the first DAPCs onto the geographic space  ------------------ ##


# subset SNPs genlight object for accessions that have a geo pos.
which(is.na(other(SNPs)$xy))
which(indNames(SNPs)=="HM288")  # remove "Germany" accession

ToRemove <- c(which(is.na(other(SNPs)$xy)), which(indNames(SNPs)=="HM288") )

## subsetting the SNP dataset :
SNPs2 <- SNPs[ - ToRemove, ]
# checking
which(is.na(other(SNPs2)$xy))

## libraries for drawing geographical maps:
library(mapplots)  # Data visualisation on maps
library(maptools)  # Tools for reading and handling spatial objects
library(rworldmap) # Mapping global data, vector and raster


## geographical map :
newmap <- getMap(resolution = "high")

x11(width=15, height=8)
## CairoPDF(file = "./DAPC_BassinMedit.pdf",
##          width = 20, height = 14, onefile = TRUE, family = "Helvetica",
##          title = "R Graphics Output", fonts = NULL, version = "1.1",
##          pointsize=18,paper="A4")
plot(newmap, axes=TRUE, ylim=c(29,45),xlim=c(-10,35), col="gray95")

points(other(SNPs2)$xy,
       cex = 3, pch = 20,
       col = funky(5)[tbl$Major[ - ToRemove ]])

title("Plot Mtr pop. structure using DAPC", col.main = "firebrick", line = -2, cex = 3)

##dev.off()


## ------------------ Using BioClimatic variables  ------------------ ##

library(raster)
library(rgdal)
library(RColorBrewer)

############## climates , WorldClim data
#
# slides
#


prec9 <-  raster("./data/prec9.bil")
prec10 <- raster("./data/prec10.bil")

ith <- raster("./data/bio3.bil")
#create a 'green' palette for rain intensities
rain <- colorRampPalette(brewer.pal(9,"Greens"))(100)


x11(width=14, height=8)
#png("./figs/Mahgreb_Precip.png",width=3000,height=1400)
par(mfrow=c(1,2))   ## TWO maps in the same figure

plot(prec9, 
     ylim=c(30,44),  ## latitude
     xlim=c(-11,12), ## longitude
     col = rain)
plot(newmap, ylim=c(32,44), xlim=c(-11,12), axes = TRUE,  add = TRUE)
points(other(SNPs2)$xy,
       cex = 3, pch = 20,
       col = funky(5)[tbl$Major[ - ToRemove ]])   ## colors function of group from dapc

plot(prec10, ylim=c(30,44),xlim=c(-11,12), col = rain)
plot(newmap, ylim=c(32,44), xlim=c(-11,12), axes=TRUE,  add=TRUE)
points(other(SNPs2)$xy,
       cex = 3, pch = 20,
       col = funky(5)[tbl$Major[ - ToRemove ]])
#dev.off()


## isothermality :
# I have interpreted Isothermality as "temperature even-ness" over the course of a year.
# BIO3 (Isothermality) is a quantification of how large the day-to-night temperature oscillation is in comparison to the summer-to-winter oscillation.
# A value of 100 would represent a site where the diurnal temperature range is equal to the annual temperature range.
# A value of 50 would indicate a location where the diurnal temperature range is half of the annual temperature range.
# autremendit dit : 100 pas de saisons contrastées, 0 TRES contrastées


iso <- colorRampPalette(brewer.pal(12,"Oranges"))(100)

x11(width=14, height=8)
#png("./figs/Mahgreb_Precip.png",width=3000,height=1400)

plot(ith, ylim=c(30,50),xlim=c(-10,40), col=iso)
plot(newmap, ylim=c(32,38), xlim=c(-11,12), axes=TRUE,  add=TRUE)
points(other(SNPs2)$xy,
       cex = 3, pch = 20,
       col = funky(5)[tbl$Major[ - ToRemove ]])


x11(width=14, height=9)
#png("./figs/Mahgreb_Precip.png",width=3000,height=1400)

plot(ith, ylim=c(30,40),xlim=c(-10,15), col=iso)
plot(newmap, ylim=c(32,38), xlim=c(-11,12), axes=TRUE,  add=TRUE)
points(other(SNPs2)$xy,
       cex = 3, pch = 20,
       col = funky(5)[tbl$Major[ - ToRemove ]])




## ------------------ WorldClim extracting data  ------------------ ##


Dir2 <-  "./data/"

predicteurs <- stack(
    paste(Dir2,"bio1.bil",sep=''),
    paste(Dir2,"bio2.bil",sep=''),
    paste(Dir2,"bio3.bil",sep=''),
    paste(Dir2,"bio4.bil",sep=''),
    paste(Dir2,"bio5.bil",sep=''),
    paste(Dir2,"bio6.bil",sep=''),
    paste(Dir2,"bio7.bil",sep=''),
    paste(Dir2,"bio8.bil",sep=''),
    paste(Dir2,"bio9.bil",sep=''),
    paste(Dir2,"bio10.bil",sep=''),
    paste(Dir2,"bio11.bil",sep=''),
    paste(Dir2,"bio12.bil",sep=''),
    paste(Dir2,"bio13.bil",sep=''),
    paste(Dir2,"bio14.bil",sep=''),
    paste(Dir2,"bio15.bil",sep=''),
    paste(Dir2,"bio16.bil",sep=''),
    paste(Dir2,"bio17.bil",sep=''),
    paste(Dir2,"bio18.bil",sep=''),
    paste(Dir2,"bio19.bil",sep='')
    )

## How information is stored in predictiers ?
names(predicteurs)
## for example for bioclimatic variable bio1 ?
predicteurs$bio1

## we need to extract the values of ALL bioclimatic variables at the locations of ALL the samples:
## this line is doing all the job for you :
records <- extract(predicteurs, other(SNPs2)$xy)
head(records)


## merge infos from different files
tblBioClim <- data.frame(indNames(SNPs2), tbl[tbl$Acc %in% indNames(SNPs2) , -6], other(SNPs2)$xy, records)
str(tblBioClim)
head(tblBioClim)


### Some correlations ?

library(psych)
library(corrplot)


CorrelTest <- corr.test(tblBioClim[ , c(10:28)], tblBioClim[ ,c(2:6)],use="complete.obs",method="pearson",alpha=0.01)  # de library(psych)
# p-values are listed within the object


dimnames(CorrelTest$r) <- list(c(" Annual Mean Temperature",
                                 " Mean Diurnal Range",
                                 " Isothermality",
                                 " Temperature Seasonality",
                                 " Max Temperature of Warmest Month",
                                 " Min Temperature of Coldest Month",
                                 " Temperature Annual Range",
                                 " Mean Temperature of Wettest Quarter",
                                 " Mean Temperature of Driest Quarter",
                                 " Mean Temperature of Warmest Quarter",
                                 " Mean Temperature of Coldest Quarter",
                                 " Annual Precipitation",
                                 " Precipitation of Wettest Month",
                                 " Precipitation of Driest Month",
                                 " Precipitation Seasonality",
                                 " Precipitation of Wettest Quarter",
                                 " Precipitation of Driest Quarter",
                                 " Precipitation of Warmest Quarter",
                                 " Precipitation of Coldest Quarter"
                                 ),
                               c("group 1", "group 2", "group 3", "group 4", "group 5"))



x11()
corrplot(CorrelTest$r, method = "color", addCoef.col = "black", order = "original", insig = "blank",
         number.digits = 1, number.cex = 0.8,

##         cl.cex = 3, ## legende  
         cl.pos = 'b', cl.offset = 4, ## position de la légende
         tl.offset = 1.5, tl.cex = 1, tl.col = "black", tl.srt = 45, #Rotation des etiquettes de textes
         p.mat = CorrelTest$p, sig.level = 0.05)

### conclusions ?
