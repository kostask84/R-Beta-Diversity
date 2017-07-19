# Appendix S1 to: Partitioning abundance-based multiple-site dissimilarity into components: 
# balanced variation in abundance and abundance gradients

# Author: Andrés Baselga

# New functions to be included in R package betapart 1.4



###############################################

betapart.core.abund <- function(x){
  
  # test for a numeric matrix or data.frame
  if(! is.matrix(x)){
    x<-as.matrix(x)
  }
  
  if(! is.numeric(x))
    stop("The data in x is not numeric.",call.=TRUE)
  
  # simple test for positive data
  xvals <-  unique(as.vector(x))
  if (any(xvals<0)) 
    stop("The table contains negative values: data should be abundances.", call. = TRUE)
  
  pair.shared.abund<-matrix(nrow=nrow(x),ncol=nrow(x))
  rownames(pair.shared.abund)<-rownames(x)
  colnames(pair.shared.abund)<-rownames(x)
  
  pair.not.shared.abund<-matrix(nrow=nrow(x),ncol=nrow(x))
  rownames(pair.not.shared.abund)<-rownames(x)
  colnames(pair.not.shared.abund)<-rownames(x)
  
  for(i in 1:nrow(x)) {
    for(j in i:nrow(x)) {
      pair.shared.abund[j,i]<-sum(pmin(x[i,],x[j,]))
      pair.not.shared.abund[i,j]<-sum(x[i,])-sum(pmin(x[i,],x[j,]))
      pair.not.shared.abund[j,i]<-sum(x[j,])-sum(pmin(x[i,],x[j,]))
    }
  }
  
  pair.shared.abund<-as.dist(pair.shared.abund)
  pair.max.not.shared.abund<-pmax(as.dist(t(upper.tri(pair.not.shared.abund)*pair.not.shared.abund)), as.dist(pair.not.shared.abund))
  pair.min.not.shared.abund<-pmin(as.dist(t(upper.tri(pair.not.shared.abund)*pair.not.shared.abund)), as.dist(pair.not.shared.abund))
  
  multiple.shared.abund<-sum(x)-sum(apply(x,2,max))
  
  
  
  computations<-list(data=x, multiple.shared.abund=multiple.shared.abund, pair.shared.abund=pair.shared.abund, pair.min.not.shared.abund=pair.min.not.shared.abund, 
                     pair.max.not.shared.abund=pair.max.not.shared.abund, pair.sum.not.shared.abund=pair.min.not.shared.abund+pair.max.not.shared.abund)
  class(computations)<-"betapart.abund"
  
  return(computations)
} 


###############################################


###############################################

beta.pair.abund<-function (x, index.family = "bray") 
{
  index.family <- match.arg(index.family, c("bray", "ruzicka"))
  if (!inherits(x, "betapart.abund")) {
    x <- betapart.core.abund(x)
  }
  switch(index.family, bray = {
    beta.bray.bal <- x$pair.min.not.shared.abund/(x$pair.min.not.shared.abund + x$pair.shared.abund)
    beta.bray.gra <- ((x$pair.max.not.shared.abund - x$pair.min.not.shared.abund)/((2 * 
                                                                                      x$pair.shared.abund) + x$pair.sum.not.shared.abund)) * (x$pair.shared.abund/(x$pair.min.not.shared.abund + 
                                                                                                                                                                     x$pair.shared.abund))
    beta.bray <- x$pair.sum.not.shared.abund/(2 * x$pair.shared.abund + x$pair.sum.not.shared.abund)
    pairwise <- list(beta.bray.bal = as.dist(beta.bray.bal), beta.bray.gra = as.dist(beta.bray.gra), 
                     beta.bray = as.dist(beta.bray))
  }, ruzicka = {
    beta.ruz.bal <- (2 * x$pair.min.not.shared.abund)/((2 * x$pair.min.not.shared.abund) + 
                                                         x$pair.shared.abund)
    beta.ruz.gra <- ((x$pair.max.not.shared.abund - x$pair.min.not.shared.abund)/(x$pair.shared.abund + 
                                                                                    x$pair.sum.not.shared.abund)) * (x$pair.shared.abund/((2 * x$pair.min.not.shared.abund) + 
                                                                                                                                            x$pair.shared.abund))
    beta.ruz <- x$pair.sum.not.shared.abund/(x$pair.shared.abund + x$pair.sum.not.shared.abund)
    pairwise <- list(beta.ruz.bal = as.dist(beta.ruz.bal), beta.ruz.gra = as.dist(beta.ruz.gra), 
                     beta.ruz = as.dist(beta.ruz))
  })
  return(pairwise)
}

###############################################


###############################################

beta.multi.abund<-function (x, index.family = "bray") 
{
  index.family <- match.arg(index.family, c("bray", "ruzicka"))
  if (!inherits(x, "betapart.abund")) {
    x <- betapart.core.abund(x)
  }
  maxbibj <- sum(x$pair.max.not.shared.abund)
  minbibj <- sum(x$pair.min.not.shared.abund)
  switch(index.family, bray = {
    beta.bray.bal <- minbibj/(minbibj + x$multiple.shared.abund)
    beta.bray.gra <- (x$multiple.shared.abund/(minbibj + x$multiple.shared.abund)) * ((maxbibj - minbibj)/((2 * 
                                                                                                              x$multiple.shared.abund) + maxbibj + minbibj))
    beta.bray <- (minbibj + maxbibj)/(minbibj + maxbibj + 
                                        (2 * x$multiple.shared.abund))
    multi <- list(beta.BRAY.BAL = beta.bray.bal, beta.BRAY.GRA = beta.bray.gra, 
                  beta.BRAY = beta.bray)
  }, ruzicka = {
    beta.ruz.bal <- (2 * minbibj)/((2 * minbibj) + x$multiple.shared.abund)
    beta.ruz.gra <- (x$multiple.shared.abund/((2 * minbibj) + x$multiple.shared.abund)) * ((maxbibj - 
                                                                                              minbibj)/((x$multiple.shared.abund) + maxbibj + minbibj))
    beta.ruz <- (minbibj + maxbibj)/(minbibj + maxbibj + 
                                       x$multiple.shared.abund)
    multi <- list(beta.RUZ.BAL = beta.ruz.bal, beta.RUZ.GRA = beta.ruz.gra, 
                  beta.RUZ = beta.ruz)
  })
  return(multi)
}


###############################################


###############################################

beta.sample.abund<-function (x, index.family = "bray", sites = nrow(x), samples = 1) 
{
  index.family <- match.arg(index.family, c("bray", "ruzicka"))
  
  if (sites > nrow(x)) 
    stop("More sites requested for sample than are in the dataset")
  pb <- txtProgressBar(min = 0, max = samples, style = 3)
  results.n <- as.data.frame(matrix(nrow = samples, ncol = 3))
  
  for (i in 1:samples) {
    position <- as.vector(1:nrow(x))
    sample.position <- sample(position, sites)
    x.sample <- x[sample.position,]
    x.beta <- beta.multi.abund(x.sample, index.family)
    results.n[i, ] <- unlist(x.beta)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(results.n) <- names(x.beta)
  result <- list(sampled.values = results.n, mean.values = sapply(results.n, 
                                                                  mean), sd.values = sapply(results.n, sd))
  return(result)
}

###############################################


###############################################
multi.bray.minus1 <- function(x){
  
  # test for a numeric matrix or data.frame
  if(! is.matrix(x)){
    x<-as.matrix(x)
  }
  
  if(! is.numeric(x))
    stop("The data in x is not numeric.",call.=TRUE)
  
  # simple test for positive data
  xvals <-  unique(as.vector(x))
  if (any(xvals<0)) 
    stop("The table contains negative values: data should be abundances.", call. = TRUE)
  
  bray1<-(nrow(x)*sum(apply(x,2,max))-sum(x))/(sum(x)*(nrow(x)-1))
  return(bray1)
}
###############################################


###############################################
multi.bray.chao <- function(x){
  
  # test for a numeric matrix or data.frame
  if(! is.matrix(x)){
    x<-as.matrix(x)
  }
  
  if(! is.numeric(x))
    stop("The data in x is not numeric.",call.=TRUE)
  
  # simple test for positive data
  xvals <-  unique(as.vector(x))
  if (any(xvals<0)) 
    stop("The table contains negative values: data should be abundances.", call. = TRUE)
  
  species.means<-apply(x,2,sum)/nrow(x)
  
  difs<-abs(x-t(matrix(rep(species.means,nrow(x)), nrow=ncol(x))))
  
  delta.BCnorm<-sum(difs)/(2*(1-1/nrow(x))*sum(x))
  return(delta.BCnorm)
}
###############################################