library(jpeg)
library(dplyr)
library(imager)

## load image as an array
trogne <- readJPEG('Trogne.jpg')
ncol(trogne)
nrow(trogne)

## to vizualise
trogne2 <- load.image('Trogne.jpg')
x11()
plot(trogne2)

## SVD of image: U.Sigma.t(V)
trogne.svd <- svd(trogne)
str(trogne.svd)


## visualize information gathered by singular values
x11()
barplot(trogne.svd$d[1:30])

## reconstruct image 
for (i in seq.int(2, 30, by = 2 )) {
    trogne.compress <- trogne.svd$u[,1:i] %*% diag(trogne.svd$d[1:i]) %*% t(trogne.svd$v[,1:i])
    writeJPEG(trogne.compress, paste('Trogne_compressed_svd_rank_', round(i,0), '.jpg', sep=''))
  }

