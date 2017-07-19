Appendix S3. R scripts for the case-study using European longhorn beetles.
Data are available at www.cerambycidae.narod.ru (Danilevsky M.L. A check-list
                                                 of Longicorn Beetles (Coleoptera, Cerambycoidea) of Europe [Web page]).
European countries were divided into two groups according their mean latitude
(i.e. Northern countries have mean latitude higher than 48º).
## MULTIPLE-SITE DISSIMILARITIES
# Read the functions to compute de multiple-site (beta.SOR, beta.SIM and
# beta.NES)
source("beta-multi.R")
# Read the data. The Northern dataset has 19 countries, but the Southern
# dataset has only 15
ceram.s<-read.table(file="ceram-s.txt", header=T)
ceram.n<-read.table(file="ceram-n.txt", header=T)
# Create a matrix to save the results. In this case we want to compute
# three indices (beta.SOR, beta.SIM and beta.NES) and 100 samples
results.n<-matrix(nrow=100, ncol=3)
colnames(results.n)<-c("beta.SOR", "beta.SIM", "beta.NES")
# Loop on the selected number of samples
for(i in 1:100){
  # Take a sample of the Northern dataset with the same number of cases
  # as the Southern dataset
  position<-as.vector(1:nrow(ceram.n))
  sample.position<-sample(position, nrow(ceram.s))
  ceram.n.sample<-ceram.n[sample.position,1:ncol(ceram.n)]
  # Compute the three indices for this sample and save the results in the
  # results matrix
  results.n[i,1]<-beta.SOR(ceram.n.sample)
  results.n[i,2]<-beta.SIM(ceram.n.sample)
  results.n[i,3]<-beta.NES(ceram.n.sample)
}
# Compute the mean for each index
# These are the results for the Northern dataset
beta.SOR.n<-mean(results.n[,1])
beta.SOR.n<-mean(results.n[,2])
beta.NES.n<-mean(results.n[,3])
2
# Compute the three indices for the Southern dataset
beta.SOR.s<-beta.SOR(ceram.s)
beta.SIM.s<-beta.SIM(ceram.s)
beta.NES.s<-beta.NES(ceram.s)
## PAIR-WISE DISSIMILARITIES
# Read the functions to compute distance matrices using the three
# pair-wise dissimilarities (beta.sor, beta.sim and beta.nes)
source("beta-pairwise.R")
# Compute the distance matrices for the Northern and Southern datasets
dist.sor.n<-beta.sor(ceram.n)
dist.sim.n<-beta.sim(ceram.n)
dist.nes.n<-beta.nes(ceram.n)
dist.sor.s<-beta.sor(ceram.s)
dist.sim.s<-beta.sim(ceram.s)
dist.nes.s<-beta.nes(ceram.s)
# Read the geographic distance matrices for the Northern and Southern
# datasets
distgeo.s<-as.dist(read.table(file="distgeo-south.txt", header=TRUE))
distgeo.n<-as.dist(read.table(file="distgeo-north.txt", header=TRUE))
# Use the vegan package to conduct the Mantel tests
library(vegan)
# Conduct the Mantel tests
sor.geo.n<-mantel(distgeo.n, dist.sor.n, method="pearson",
                  permutations=1000)
sim.geo.n<-mantel(distgeo.n, dist.sim.n, method="pearson",
                  permutations=1000)
nes.geo.n<-mantel(distgeo.n, dist.nes.n, method="pearson",
                  permutations=1000)
sor.geo.s<-mantel(distgeo.s, dist.sor.s, method="pearson",
                  permutations=1000)
sim.geo.s<-mantel(distgeo.s, dist.sim.s, method="pearson",
                  permutations=1000)
nes.geo.s<-mantel(distgeo.s, dist.nes.s, method="pearson",
                  permutations=1000)
# Use the boot package to estimate the frequency distributions of
# parameters
library(boot)
# Unfold distance matrices into vectors
3
data.sor.n<-data.frame(as.vector(dist.sor.n), as.vector(distgeo.n))
colnames(data.sor.n)<-c("d", "distgeo")
data.sim.n<-data.frame(as.vector(dist.sim.n), as.vector(distgeo.n))
colnames(data.sim.n)<-c("d", "distgeo")
data.nes.n<-data.frame(as.vector(dist.nes.n), as.vector(distgeo.n))
colnames(data.nes.n)<-c("d", "distgeo")
data.sor.s<-data.frame(as.vector(dist.sor.s), as.vector(distgeo.s))
colnames(data.sor.s)<-c("d", "distgeo")
data.sim.s<-data.frame(as.vector(dist.sim.s), as.vector(distgeo.s))
colnames(data.sim.s)<-c("d", "distgeo")
data.nes.s<-data.frame(as.vector(dist.nes.s), as.vector(distgeo.s))
colnames(data.nes.s)<-c("d", "distgeo")
# This function returns the parameters (intercept and slope) of a linear
# model
boot.coefs<-function(data, indices){
  data<-data[indices,]
  mod<-lm(d~distgeo, data=data)
  coefficients(mod)
}
# Bootstrap the parameters
boot.sor.n<-boot(data.sor.n, boot.coefs, R=1000)
boot.sim.n<-boot(data.sim.n, boot.coefs, R=1000)
boot.nes.n<-boot(data.nes.n, boot.coefs, R=1000)
boot.sor.s<-boot(data.sor.s, boot.coefs, R=1000)
boot.sim.s<-boot(data.sim.s, boot.coefs, R=1000)
boot.nes.s<-boot(data.nes.s, boot.coefs, R=1000)
# Compute the empirical probabilities of a parameter being higher in the
# Northern dataset
# Slopes
n<-0
for (i in 1:1000){
  if (boot.sor.n$t[i,2]>boot.sor.s$t[i,2]){n<-n+1} else{n<-n+0}
}
p.slope.sor<-n/1000
n<-0
for (i in 1:1000){
  if (boot.sim.n$t[i,2]>boot.sim.s$t[i,2]){n<-n+1} else{n<-n+0}
}
p.slope.sim<-n/1000
n<-0
for (i in 1:1000){
  if (boot.nes.n$t[i,2]>boot.nes.s$t[i,2]){n<-n+1} else{n<-n+0}
}
p.slope.nes<-n/1000
# Intercepts
n<-0
4
for (i in 1:1000){
  if (boot.sor.n$t[i,1]>boot.sor.s$t[i,1]){n<-n+1} else{n<-n+0}
}
p.interc.sor<-n/1000
n<-0
for (i in 1:1000){
  if (boot.sim.n$t[i,1]>boot.sim.s$t[i,1]){n<-n+1} else{n<-n+0}
}
p.interc.sim<-n/1000
n<-0
for (i in 1:1000){
  if (boot.nes.n$t[i,1]>boot.nes.s$t[i,1]){n<-n+1} else{n<-n+0}
}
p.interc.nes<-n/1000