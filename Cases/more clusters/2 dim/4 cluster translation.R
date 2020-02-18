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

parameter = c(1, 0)/5.5 # 45 degrees

# scalar = 0
# scalar = 1/4
# scalar = 1/2
# scalar = 1
# scalar = 2
# scalar = 4
# scalar = 8
scalar = 16

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
                                 method = "ICS", model = "LS", hyper = F),prior=list(strength=1.25, m0=distance*scalar*parameter, k0=1/variance))
  
  
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

print("SCALAR")
print(scalar)
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

# parameter=(1,1)/(sqrt(2) * 5.5)

# scalar = 0
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9695701
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9796003
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9689291
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9651584
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9742459
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9723982

# scalar = 1/4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9616893
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9745098
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"

# > print(RAND(truepart,Hier1[,]))
# [1] 0.9595023
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9593891
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9687406
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9707014

# scalar = 1/2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9649321
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9850679
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9631976
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9605204
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9807692
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9816742

# scalar = 1
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9644042
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9841629
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9635747
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.963537
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.980543
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9800528

# scalar = 2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9414404
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9698341
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9409879
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9396682
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9615762
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9615385

# scalar = 4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9353318
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9536576
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9341629
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9352187
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9453997
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9456637

# scalar = 8
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9204374
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9006787
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9206637
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9183635
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9106712
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9121418

# scalar = 16
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.6823906
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.6200603
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.6834087
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.6825038
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.6711538
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.6702489


# parameter = (1,0)/(sqrt(2) * 5.5)

# scalar = 0
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9695701
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9796003
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9689291
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9651584
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9742459
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9723982

# scalar = 1/4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9611237
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9769231
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9607466
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9601433
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9704751
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9713424

# scalar = 1/2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9622926
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9790347
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9616893
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.962368
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9711161
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9714178

# scalar = 1
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9608974
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9826923
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9588612
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9602187
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9745852
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9743213

# scalar = 2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9504902
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9674962
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9509804
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9488688
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9659879
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9645928

# scalar = 4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.954902
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9735294
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9536576
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9532428
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9702866
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.968552

# scalar = 8
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.9219457
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.9044495
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.9208522
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.9165913
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.9092006
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.9018854

# scalar = 16
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin[,])) 
# [1] 0.7395173
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin[,]))
# [1] 0.7140271
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1[,]))
# [1] 0.7403846
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2[,]))
# [1] 0.7403846
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3[,]))
# [1] 0.7403846
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4[,]))
# [1] 0.7403846