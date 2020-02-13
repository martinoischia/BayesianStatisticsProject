library(coda)

library(BNPmix)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
library(MYpackage)
require(gplots)
require(ggpubr)
## VI method
VImethod = #greedy 
			#VI.loss.draws 
			VI.ineq.draws 

## RAND index

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

# RAND = function(partition1, partition2){

# mat1=	1 * (sapply(partition1, 
                              # function(x) x == partition1))
# mat2= 1 * (sapply(partition2, 
                              # function(x) x == partition2))

# return( (choose(length(partition1),2)-(sum(abs(mat1-mat2))/2))/ choose(length(partition1),2) )
# }



## Number of repetitions for each case 
set.seed(42)
iterate= 10
n=15 #sample size

##structure for keeping the results:

m=3 #number of models
Bindermin= array(0, c(iterate, n, m))
VImin =    array(0, c(iterate, n, m))
Hier1 = 	array(0, c(iterate, n, m))
Hier2 = 	array(0, c(iterate, n, m))
Hier3 = 	array(0, c(iterate, n, m))
Hier4 = 	array(0, c(iterate, n, m))

for(it in 1:iterate){
##### SIMULATION OF DATA TOY DA 2 GAUSSIANE UNIVARIATE
##### SIMULATE A MARKOV CHAIN BASED ON THE "CLASSIC" MIXTURE MODEL



data_toy <- rbind(c(rnorm(5, -2, 1), rnorm(10, 2, 1)),
					c(rnorm(4, -2, 1), rnorm(4, 2, 1), rnorm(4, 6, 1), rnorm(3, -6, 1)),
					c(rnorm(2, -2, 1), rnorm(2, 2, 1),rnorm(2, 6, 1),rnorm(2, -6, 1),rnorm(2, 10, 1),rnorm(2, -10, 1),rnorm(2, 14, 1),rnorm(1, -14, 1)))

##### MCMC sampling for a DP mixture model
for(mod in 1:3){
model <- PYdensity(data_toy[mod,], 
                       mcmc = list(niter = 15000, nburn = 5000, 
                                   method = "ICS", model = "LS", hyper = F)) 
								   
Binder = B.loss.draws(model$clust)
Bindermin[it,,mod]=Binder$min

hier=Hier(model)
Hier1[it,,mod] = hier$average.Binder
Hier2[it,,mod] = hier$complete.Binder
Hier3[it,,mod] =hier$average.VI
Hier4[it,,mod] =hier$complete.VI

VI = VImethod(model$clust)
VImin[it,,mod]=VI$min 
				#VI$cluster

}
print (paste("iteration number ",it))
}

truepart1 = c(rep(0,5), rep(1,10))
truepart2 = c(rep(0,4), rep(1,4),rep(2,4),rep(3,3))
truepart3 = c(rep(0,2), rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),7)
truepart  = rbind(truepart1, truepart2, truepart3)

for (model in 1:m){
cat("\n")
print(paste("Model", model))
cat("\n")
print("RAND binder")
print(RAND(truepart[model,],Bindermin[,,model])) 
print("RAND VI")
print(RAND(truepart[model,],VImin[,,model]))
print("RAND hier av Binder")
print(RAND(truepart[model,],Hier1[,,model]))
print("RAND hier co Binder")
print(RAND(truepart[model,],Hier2[,,model]))
print("RAND hier av VI")
print(RAND(truepart[model,],Hier3[,,model]))
print("RAND hier co VI")
print(RAND(truepart[model,],Hier4[,,model]))
}
