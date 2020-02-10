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


##### SIMULATION OF DATA TOY OF 1 GAMMA WITH DIFFERENT VARIANCES.

# Number of points
n = 10
# Alpha parameter of the distribution
alpha = 10
# Rate
rate_used = sqrt(10)
# Alphas to test
alphas = list(0.01, 0.5, 1, 3, 5)

for (i in alphas){
  set.seed(42)
  
  data_toy <- c(rgamma(n, alpha, rate=rate_used))
  
  ##### MCMC sampling for a DP mixture model.
  
  model <- PYdensity(data_toy, 
                     mcmc = list(niter = 25000, nburn = 15000, 
                                 method = "SLI", model = "LS", hyper = F),
                     prior=list(strength=i)) 
  
  summary(model)
  ENOC(1, length(data_toy))
  
  ## Plot of the estimated density and its bayesian credible intervals
  
  data.frame(x = seq(min(model$grideval), max(model$grideval), by = 0.01), 
             y = dgamma(seq(min(model$grideval), max(model$grideval), by = 0.01), alpha, rate=rate_used)) %>%
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
  print("CASE ALPHA:")
  print(i)
  
  Binder = B.loss.draws(model$clust)
  print("Binder´s Min.: ")
  print(Binder$min)
  
  VI = VI.loss.draws(model$clust)
  VI.ineq = VI.ineq.draws(model$clust)
  print("VI Min.:")
  print(VI$min)
  print("VI Ineq Min.:")
  print(VI.ineq$min)
  
  # Binder$exp.loss	
  # VI$exp.loss
  # VI.ineq$VI
  # Partitions found by exploiting different hierarchical clusterings
  hier=Hier(model)
  
  print("Average Binder:")
  print(hier$average.Binder)
  print("Complete Binder:")
  print(hier$complete.Binder)
  
  print("Average VI:")
  print(hier$average.VI)
  print("Complete VI")
  print(hier$complete.VI)
  
  print("Average VI Ineq:")
  print(hier$average.VI.ineq)
  print("Complete VI Ineq:")
  print(hier$complete.VI.ineq)
  
  print("Average PSM:")
  print(hier$average.PSM)
  print("Complete PSM:")
  hier$complete.PSM
  
  print("Average EU:")
  print(hier$average.eu)
  print("Complete EU:")
  print(hier$complete.eu)
  
  print("Partitions counter:")
  print(partitions.counter(model$clust)[1:5,])
  
  print("Greedy VI:")
  print(greedy (model$clust,l=	100*ncol(model$clust), maxiter= 10 ,distance="VI"))
  print("Greedy VI Ineq:")
  print(greedy (model$clust,l=	100*ncol(model$clust), maxiter= 10 ,distance="VI",Jensen=FALSE))
  print("Greedy Binder:")
  print(greedy (model$clust,l=	100*ncol(model$clust), maxiter= 20,distance="Binder"))
}
