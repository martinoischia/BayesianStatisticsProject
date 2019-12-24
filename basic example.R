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
	# data.frame(x = seq(-5, 5, by = 0.1), 
           # y = 0.33 * dnorm(seq(-5, 5, by = 0.1), -2, 1) + 
             # 0.67 * dnorm(seq(-5, 5, by = 0.1), 2, 1)) %>%
  # ggplot() + 
  # theme_minimal() + 
  # geom_line(aes(x = x, y = y)) + 
  # geom_point(data.frame(x = data_toy, col = as.factor(c(rep(1, 5), rep(2, 10)))), 
             # mapping =  aes(x = x, y = 0, col = col)) +
  # theme(legend.position = "none")

##### MCMC sampling for a DP mixture model

model <- PYdensity(data_toy, 
                       mcmc = list(niter = 15000, nburn = 5000, 
                                   method = "SLI", model = "LS", hyper = F)) #output=list( out_type= "mean", mean_dens=TRUE))
								   
								   
summary(model)
ENOC(1, length(data_toy))

## Plot of the estimated density and its bayesian credible intervals

plot_model <- plot(model, show_points = T, show_hist = T, show_clust = T)
x11()
ggarrange(plot_model)

	##plot to be fixed
	# lines(x = seq(-5, 5, by = 0.1), 
	# y = 0.33 * dnorm(seq(-5, 5, by = 0.1), -2, 1) + 
	# 0.67 * dnorm(seq(-5, 5, by = 0.1), 2, 1), col= 'red')

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
# Partitions found by eploiting different hierarchical
# clusterings
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
# VI$exp.loss
# VI.ineq$VI
partitions.counter(model$clust)[1:5,] 


greedy (model$clust,l=	100*ncol(model$clust), maxiter= 10 ,distance="VI") 

greedy (model$clust,l=	100*ncol(model$clust), maxiter= 20,distance="Binder")

# in this case all the estimates give the same result
# I'm a bit surprised that i had to put 100 above here, with fewer you would get stuck 
# in another partition, but in the paper was suggested to put 2


