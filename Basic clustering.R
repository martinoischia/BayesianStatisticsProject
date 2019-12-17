library(coda)
library(BNPmix)
require(gplots)
library(ggplot2)


set.seed(42)
data_toy <- c(rnorm(5, -2, 1), rnorm(10, 2, 1))
data.frame(x = seq(-5, 5, by = 0.1), 
           y = 0.4 * dnorm(seq(-5, 5, by = 0.1), -2, 1) + 
             0.6 * dnorm(seq(-5, 5, by = 0.1), 2, 1)) %>%
  ggplot() + 
  theme_minimal() + 
  geom_line(aes(x = x, y = y)) + 
  geom_point(data.frame(x = data_toy, col = as.factor(c(rep(1, 5), rep(2, 10)))), 
             mapping =  aes(x = x, y = 0, col = col)) +
  theme(legend.position = "none")

n <- length(data_toy)
model <- PYdensity(data_toy, 
                       mcmc = list(niter = 3000, nburn = 2000, 
                                   method = "SLI", model = "LS", hyper = F))
								   
								   
# compute the posterior similarity matrix
PSM <- matrix(0, ncol = n, nrow = n)
model$clust
for(j in 1:n){
  for(k in 1:n){
    PSM[j,k] <- mean(model$clust[,j] == model$clust[,k])
  }
}

# compute the distance for each partition to the PSM
d_abs <- c()
for(i in 1:nrow(model$clust)){
  d_abs[i] <- sum(abs(1 * (sapply(model$clust[i,], 
                              function(x) x == model$clust[i,])) - PSM))
  
}

# Get the optimal one
model$clust[which.min(d_abs),]								   