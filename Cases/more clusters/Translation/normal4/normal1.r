library(BNPmix)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
library(MYpackage)
library(mvtnorm)
require(gplots)
require(ggpubr)


RAND = function(partition1, partition2){
  lung = nrow(partition2)
  vec = vector(length= lung)
  for(i in 1:lung){
    mat1=	1 * (sapply(partition1, 
                      function(x) x == partition1))
    mat2= 1 * (sapply(partition2[i,], 
                      function(x) x == partition2[i,]))
    vec[i]=(choose(length(partition1),2)-(sum(abs(mat1-mat2))/2))/ choose(length(partition1),2)
  }
  return( mean (vec) )
}


##### SIMULATION OF DATA TOY OF 2 MULTIVARIATE NORMAL, ONE OF THEM CENTERED IN THE ORIGIN, 
##### AND THE OTHER TRANSLATING FROM (-5, -5) TO (5, 5).

## Number of repetitions for each case 
set.seed(42)
iterate= 8
n= 15 # Sample Size.

## Means to test

means = list(c(-5, -5), c(-3, -3), c(-1, -1), c(1, 1), c(3, 3), c(5, 5))
i=1
print('Mean: ')
print(means[[i]])

## Structure for keeping the results:

Bindermin= matrix(0, nrow=iterate, ncol = n)
VImin =    matrix(0, nrow=iterate, ncol = n)
Hier1 = 	matrix(0, nrow=iterate, ncol = n)
Hier2 = 	matrix(0, nrow=iterate, ncol = n)
Hier3 = 	matrix(0, nrow=iterate, ncol = n)
Hier4 =		matrix(0, nrow=iterate, ncol = n)

for(it in 1:iterate){
  print (paste("Iteration Number ",it))
  data_toy <- #rbind(rmvnorm(5, means[[i]]), rmvnorm(10, c(0, 0)))
              rbind(rmvnorm(3, means[[i]]), rmvnorm(4, c(-2, 2)), rmvnorm(4, c(0, 0)), rmvnorm(4, c(2, -2)))
              # rbind(rmvnorm(2, means[[i]]), rmvnorm(2, c(0, 3)), rmvnorm(2, c(0, 0)), rmvnorm(2, c(0, -3),
              # rmvnorm(2, c(-3, 0)), rmvnorm(2, c(3, 0)),rmvnorm(2, c(-3, 3)),rmvnorm(1, c(3, -3)))
  
  ##### MCMC sampling for a DP mixture model
  model <- PYdensity(data_toy, 
                     mcmc = list(niter = 15000, nburn = 5000, 
                                 method = "ICS", model = "LS", hyper = F)) 
  
  
  Binder = B.loss.draws(model$clust)
  Bindermin[it,]=Binder$min
  
  hier=Hier(model)
  Hier1[it,] = hier$average.Binder
  Hier2[it,] = hier$complete.Binder
  Hier3[it,] =hier$average.VI
  Hier4[it,] =hier$complete.VI
  
  
  VI = VI.ineq.draws(model$clust)
  VImin[it,]= VI$min
}

truepart1 = c(rep(0,5), rep(1,10))
truepart2 = c(rep(0,3), rep(1,4),rep(2,4),rep(3,4))
truepart3 = c(rep(0,1), rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),rep(7,2))
truepart  = truepart2

print("RAND binder")
print(RAND(truepart,Bindermin[,])) 
print("RAND VI")
print(RAND(truepart,VImin[,]))
print("RAND hier average Binder")
print(RAND(truepart,Hier1[,]))
print("RAND hier complete Binder")
print(RAND(truepart,Hier2[,]))
print("RAND hier average VI")
print(RAND(truepart,Hier3[,]))
print("RAND hier complete VI")
print(RAND(truepart,Hier4[,]))

