library(BNPmix)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
library(MYpackage)
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


## Number of repetitions for each case 
set.seed(42)
iterate= 20
n=15 #sample size

##structure for keeping the results:

Bindermin= matrix(0, nrow=iterate, ncol = n)
VImin =    matrix(0, nrow=iterate, ncol = n)
Hier1 = 	matrix(0, nrow=iterate, ncol = n)
Hier2 = 	matrix(0, nrow=iterate, ncol = n)
Hier3 = 	matrix(0, nrow=iterate, ncol = n)
Hier4 =		matrix(0, nrow=iterate, ncol = n)

for(it in 1:iterate){

data_toy <- #c(rnorm(5, -2, 1), rnorm(10, 2, 1))
					 #c(rnorm(4, -2, 1), rnorm(4, 2, 1), rnorm(4, 6, 1), rnorm(3, -6, 1))
					c(rnorm(2, -2, 1), rnorm(2, 2, 1),rnorm(2, 6, 1),rnorm(2, -6, 1),rnorm(2, 10, 1),rnorm(2, -10, 1),rnorm(2, 14, 1),rnorm(1, -14, 1))

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
print (paste("iteration number ",it))
}

truepart1 = c(rep(0,5), rep(1,10))
truepart2 = c(rep(0,4), rep(1,4),rep(2,4),rep(3,3))
truepart3 = c(rep(0,2), rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),7)
truepart  = truepart3

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