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

dimension = 2

# parameter = 1/8
# parameter = 1/4
# parameter = 1/2
# parameter = 1
# parameter = 2
# parameter = 4
# parameter = 8
parameter = 16

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

# distance = 6
distance = 5.5

means = rbind(
  c(distance/2,distance/2),
  c(distance/2,-distance/2),
  c(-distance/2,-distance/2),
  c(-distance/2,distance/2))
variance = var(means)[1,1]



##### SIMULATION OF DATA TOY OF 2 MULTIVARIATE NORMAL, ONE OF THEM CENTERED IN THE ORIGIN, 
##### AND THE OTHER TRANSLATING FROM (-5, -5) TO (5, 5).

## Number of repetitions for each case 
set.seed(42)
iterate= 20
n= 52 # must be divisible by 4!!!


## Structure for keeping the results:

data1 = list()
Bindermin= matrix(0, nrow=iterate, ncol = n)
VImin =    matrix(0, nrow=iterate, ncol = n)
Hier1 = 	matrix(0, nrow=iterate, ncol = n)
Hier2 = 	matrix(0, nrow=iterate, ncol = n)
Hier3 = 	matrix(0, nrow=iterate, ncol = n)
Hier4 =		matrix(0, nrow=iterate, ncol = n)

for(it in 1:iterate){
  print (paste("Iteration Number ",it))
  data_toy <- #rbind(rmvnorm(5, means[[i]]), rmvnorm(10, c(0, 0)))
    rbind(rmvnorm(n/4,means[1,] ), rmvnorm(n/4, means[2,]), rmvnorm(n/4, means[3,]), rmvnorm(n/4, means[4,]))
  # rbind(rmvnorm(2, means[[i]]), rmvnorm(2, c(0, 3)), rmvnorm(2, c(0, 0)), rmvnorm(2, c(0, -3),
  # rmvnorm(2, c(-3, 0)), rmvnorm(2, c(3, 0)),rmvnorm(2, c(-3, 3)),rmvnorm(1, c(3, -3)))
  
  ##### MCMC sampling for a DP mixture model
  model <- PYdensity(data_toy, 
                     mcmc = list(niter = 15000, nburn = 5000, 
                                 method = "ICS", model = "LS", hyper = F),prior=list(strength=parameter, k0=1/variance))
  
  
  Binder = B.loss.draws(model$clust)
  Bindermin[it,]=Binder$min
  data1[[it]]=data_toy
  hier=Hier(model)
  Hier1[it,] = hier$average.Binder
  Hier2[it,] = hier$complete.Binder
  Hier3[it,] =hier$average.VI
  Hier4[it,] =hier$complete.VI
  
  
  VI = VI.ineq.draws(model$clust)
  VImin[it,]= VI$min
  
}

truepart1 = c(rep(0,n/4), rep(1,10))
truepart2 = c(rep(0,n/4), rep(1,n/4),rep(2,n/4),rep(3,n/4))
truepart3 = c(rep(0,1), rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),rep(7,2))
truepart  = truepart2

print("PARAMETER")
print(parameter)
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

plot(model, show_points = T, show_hist = T, show_clust = T)

for (j in iterate:iterate) {
  # x11(); plot(data1[[j]],col=truepart+1, pch=16)
  # x11();plot(data1[[j]], col=Bindermin[j,]+1, pch=16  )
  x11();plot(data1[[j]], col =VImin[j,]+1,pch=16)
}


# Results:
# parameter= 1/8
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9783183
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9746229
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9783183
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9783183
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9728507
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9728507

# parameter=1/4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9782428
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9727376
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9786576
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9791101
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9727753
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9737179

# parameter=1/2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9851056
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9908748
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9855204
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.985181
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9886124
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9877451

# parameter=1
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9710784
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9785445
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9711161
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9695701
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9721342
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9725867

# parameter=2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9546003
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9758673
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9542986
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9553922
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9674208
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9685143

# parameter=4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9346531
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9750377
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9302413
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9310332
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9599925
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9625189

# parameter=8
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.8828431
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9410633
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.8780543
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.8766591
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9271493
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9316742

# parameter=16
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.8459276
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9087858
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.8395928
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.839178
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.8976621
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.8869532