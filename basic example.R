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


##### SIMULATION OF DATA TOY DA 2 GAUSSIANE UNIVARIATE
##### SIMULATE A MARKOV CHAIN BASED ON THE "CLASSIC" MIXTURE MODEL


set.seed(42)
data_toy <- c(rnorm(5, -2, 1), rnorm(10, 2, 1))

##### MCMC sampling for a DP mixture model

model <- PYdensity(data_toy, 
                       mcmc = list(niter = 15000, nburn = 5000, 
                                   method = "SLI", model = "LS", hyper = F)) 
								   
								   
summary(model)
ENOC(1, length(data_toy))

## Plot of the estimated density and its bayesian credible intervals

data.frame(x = seq(min(model$grideval), max(model$grideval), by = 0.1),
            y = 0.33 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.1), -2, 1) +
                0.67 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.1), 2, 1)) %>%

ggplot()+
theme_minimal() + 
geom_line(aes(x = x, y = y,col = "True density")) + 
geom_line(data = data.frame(xx=model$grideval,yy=colMeans(model$density)), aes(x=xx, y=yy, col="Estimated density"))+
scale_color_manual(name = "Group",
values = c( "True density" = "blue", "Estimated density" = "red"))+
geom_point(data.frame(x = data_toy), 
           mapping =  aes(x = x, y = 0) )

## Coda diagnostic ( raw )

coda_model <- BNPdens2coda(model)
summary(coda_model)
x11()
plot(coda_model)
effectiveSize(coda_model)


## Partition in the Markov chain minimizing the Binder Loss

Binder = B.loss.draws(model$clust)
Binder$min
# Binder$exp.loss	
# Partitions found by eploiting different hierarchical clusterings
hier=Hier(model)
hier$average.Binder   
hier$complete.Binder
hier$average.VI
hier$complete.VI
hier$average.VI.ineq
hier$complete.VI.ineq
hier$average.PSM
hier$complete.PSM
hier$average.eu
hier$complete.eu

## Partition in the Markov chain minimizing the Variation of Information Loss

VI = VI.loss.draws(model$clust)
VI.ineq = VI.ineq.draws(model$clust)
VI$min
VI.ineq$min
# VI$exp.loss
# VI.ineq$VI
partitions.counter(model$clust)[1:5,] 


greedy (model$clust,l=	100*ncol(model$clust), maxiter= 10 ,distance="VI") 
greedy (model$clust,l=	100*ncol(model$clust), maxiter= 10 ,distance="VI",Jensen=FALSE) 
greedy (model$clust,l=	100*ncol(model$clust), maxiter= 20,distance="Binder")

# in this case all the estimates give the same result
# I'm a bit surprised that i had to put 100 above here, with fewer you would get stuck 
# in another partition, but in the paper was suggested to put 2



##### this is for the developing 
# library(devtools)
# setwd("./MYpackage")
# document()
# setwd("..")
# install("MYpackage")
# 3
# reload("MYpackage")
