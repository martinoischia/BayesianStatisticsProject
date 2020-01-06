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


##### SIMULATION OF DATA TOY OF 1 GAUSSIAN DIFFUSE DISTRIBUTION.

set.seed(42)
data_toy <- c(rnorm(15, 0, 50))

##### MCMC sampling for a DP mixture model.

model <- PYdensity(data_toy, 
                   mcmc = list(niter = 25000, nburn = 15000, 
                               method = "SLI", model = "LS", hyper = F)) 

summary(model)
ENOC(1, length(data_toy))

## Plot of the estimated density and its bayesian credible intervals

data.frame(x = seq(min(model$grideval), max(model$grideval), by = 0.01), 
           y = dnorm(seq(min(model$grideval), max(model$grideval), by = 0.01), 0, 50)) %>%
  ggplot()+
  theme_minimal() + 
  geom_line(aes(x = x, y = y,col = "True density")) + 
  geom_line(data = data.frame(xx=model$grideval,yy=colMeans(model$density)), aes(x=xx, y=yy, col="Estimated density"))+
  scale_color_manual(name = "Group",
                     values = c( "True density" = "blue", "Estimated density" = "red"))+
  geom_point(data.frame(x = data_toy), 
             mapping =  aes(x = x, y = 0) )

## Coda diagnostic (raw)

coda_model <- BNPdens2coda(model)
summary(coda_model)
x11()
plot(coda_model)
effectiveSize(coda_model)


## Partition in the Markov chain minimizing the Binder Loss

Binder = B.loss.draws(model$clust)
Binder$min
Binder$exp.loss

# Partitions found by exploiting different hierarchical clusterings
Hier(model)$average.Binder   
Hier(model)$complete.Binder
Hier(model)$average.PSM
Hier(model)$complete.PSM
Hier(model)$average.eu
Hier(model)$complete.eu

## Partition in the Markov chain minimizing the Variation of Information Loss

VI = VI.loss.draws(model$clust)
VI.ineq = VI.ineq.draws(model$clust)
VI$min
VI.ineq$min
VI$exp.loss
VI.ineq$VI
partitions.counter(model$clust)[1:5,] 


greedy (model$clust,l=	100*ncol(model$clust), maxiter= 10 ,distance="VI") 
greedy (model$clust,l=	100*ncol(model$clust), maxiter= 10 ,distance="VI",Jensen=FALSE) 
greedy (model$clust,l=	100*ncol(model$clust), maxiter= 20,distance="Binder")
