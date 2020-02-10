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
## RAND index

RAND = function(partition1, partition2){
vec = vector(length=nrow(partition2))
for(i in 1:nrow(partition2)){
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
iterate= 10
n=15 #sample size
#structure for keeping the results:
m=3 #number of models
Bindermin=list(matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n))
VImin = list(matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n))
Hier1 = list(matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n))
Hier2 = list(matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n))
Hier3 = list(matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n))
Hier4 = list(matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n),matrix(nrow=iterate, ncol=n))

for(it in 1:iterate){
##### SIMULATION OF DATA TOY DA 2 GAUSSIANE UNIVARIATE
##### SIMULATE A MARKOV CHAIN BASED ON THE "CLASSIC" MIXTURE MODEL



data_toy <- rbind(c(rnorm(5, -2, 1), rnorm(10, 2, 1)),
					c(rnorm(4, -2, 1), rnorm(4, 2, 1), rnorm(4, 6, 1), rnorm(3, -6, 1)),
					c(rnorm(2, -2, 1), rnorm(2, 2, 1),rnorm(2, 6, 1),rnorm(2, -6, 1),rnorm(2, 10, 1),rnorm(2, -10, 1),rnorm(2, 14, 1),rnorm(1, -14, 1)))

##### MCMC sampling for a DP mixture model
for(m in 1:3)
model <- PYdensity(data_toy[m,], 
                       mcmc = list(niter = 15000, nburn = 5000, 
                                   method = "SLI", model = "LS", hyper = F)) 
								   
								  
## Partition in the Markov chain minimizing the Binder Loss

Binder = B.loss.draws(model$clust)
Bindermin[[m]][it,]=Binder$min

hier=Hier(model)
Hier.list1[[m]][it,] = hier$average.Binder
Hier.list2[[m]][it,] = hier$complete.Binder
Hier.list3[[m]][it,] =hier$average.VI
Hier.list4[[m]][it,] =hier$complete.VI
## Partition in the Markov chain minimizing the Variation of Information Loss

VI = VI.ineq.draws(model$clust)
VImin[[m]][it,]= VI.ineq$min
}

# print("Binder")
# print(Bindermin)
# print("Hierarchical, 1st row average binder, 2nd complete, 3rd average VI, 4th complete")
# print(Hier.list)
# print("VI")
# print(VImin)
truepart1 = c(rep(0,5), rep(1,10))
truepart2 = c(rep(0,4), rep(1,4),rep(2,4),rep(3,3))
truepart3 = c(rep(0,2), rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),7)


print("Model 1")
print("RAND binder")
RAND(truepart1,Bindermin[[1]]) 
print("RAND VI")
RAND(truepart1,VImin[[1]]) 
print("RAND hier av Binder")
RAND(truepart1,Hier.list1[[1]])
print("RAND hier co Binder")
RAND(truepart1,Hier.list2[[1]])
print("RAND hier av VI")
RAND(truepart1,Hier.list3[[1]])
print("RAND hier co VI")
RAND(truepart1,Hier.list4[[1]])

print("Model 2")
print("Model 3")
