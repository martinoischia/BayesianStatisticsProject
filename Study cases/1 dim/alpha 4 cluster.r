library(BNPmix)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
library(MYpackage)
require(gplots)
require(ggpubr)

# parameter = 1/10
# parameter = 1/4
# parameter = 1/3
# parameter = 1/2
# parameter = 1/5
# parameter = 1
# parameter = 2
# parameter = 4
# parameter = 8
parameter = 16
# parameter = 20
# parameter = 50

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

distance = 4
#distance = 3.1
variance = var(c(-3/2*distance, distance/2, -distance/2,3/2*distance))

## Number of repetitions for each case 
set.seed(42)
iterate= 20
n=20 #sample size



##structure for keeping the results:
data1 = matrix(0, nrow=iterate, ncol = n)
Bindermin= matrix(0, nrow=iterate, ncol = n)
VImin =    matrix(-1, nrow=iterate, ncol = n)
Hier1 = 	matrix(-1, nrow=iterate, ncol = n)
Hier2 = 	matrix(0, nrow=iterate, ncol = n)
Hier3 = 	matrix(0, nrow=iterate, ncol = n)
Hier4 =		matrix(0, nrow=iterate, ncol = n)
for(it in 1:iterate){

data_toy <- #c(rnorm(5, -2, 1), rnorm(10, 2, 1))
					c(rnorm(5, distance/2, 1), rnorm(5, -distance/2, 1), rnorm(5, 3/2*distance, 1), rnorm(5, -3/2*distance, 1))
					#c(rnorm(2, -2, 1), rnorm(2, 2, 1),rnorm(2, 6, 1),rnorm(2, -6, 1),rnorm(2, 10, 1),rnorm(2, -10, 1),rnorm(2, 14, 1),rnorm(1, -14, 1))
data1[it,]=data_toy
##### MCMC sampling for a DP mixture model
model <- PYdensity(data_toy, 
                       mcmc = list(niter = 15000, nburn = 5000, 
                                   method = "ICS", model = "LS", hyper = F), prior=list(strength=parameter,a0=2,m0=0,b0=1, k0=1/variance)) 
								   

Binder = B.loss.draws(model$clust)
Bindermin[it,]=Binder$min

hier=Hier(model)
Hier1[it,] = hier$average.Binder
Hier2[it,] = hier$complete.Binder
Hier3[it,] =hier$average.VI
Hier4[it,] =hier$complete.VI


VI = VI.ineq.draws(model$clust)# greedy(model$clust, l=1000,distance="Binder",Jensen=FALSE)
VImin[it,]= VI$min #VI$cluster
print (paste("iteration number ",it))
}

data.frame(x = seq(min(model$grideval), max(model$grideval), by = 0.01),
            y = 0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.01), -distance/2,1) +
                0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.01), distance/2,1)+
				0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.01), 3/2*distance,1)  + 
				0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.01), -3/2*distance,1))%>%

ggplot()+
theme_minimal() + 
geom_line(aes(x = x, y = y,col = "True density")) + 
geom_line(data = data.frame(xx=model$grideval,yy=colMeans(model$density)), aes(x=xx, y=yy, col="Estimated density"))+
scale_color_manual(name = "Group",
values = c( "True density" = "blue", "Estimated density" = "red"))+
geom_point(data.frame(x = data_toy), 
           mapping =  aes(x = x, y = 0) )



truepart1 = c(rep(0,5), rep(1,10))
truepart2 = c(rep(0,5), rep(1,5),rep(2,5),rep(3,5))
truepart3 = c(rep(0,2), rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),7)
truepart  = truepart2

print("RAND binder")
print(RAND(truepart,Bindermin)) 
print("RAND VI")
print(RAND(truepart,VImin))
print("RAND hier average Binder")
print(RAND(truepart,Hier1))
print("RAND hier complete Binder")
print(RAND(truepart,Hier2))
print("RAND hier average VI")
print(RAND(truepart,Hier3))
print("RAND hier complete VI")
print(RAND(truepart,Hier4))
for (j in 1:1) {
x11(); plot(data1[j,], rep(0,n),ylim=c(0,2),col=truepart+1, pch=16, main= 'low true, middle Binder, above VI')
points(data1[j,], rep(1,n),col=Bindermin[j,]+1, pch=16  )
points(data1[j,], rep(2,n), col =VImin[j,]+1,pch=16)
}


print(parameter)


#graphics.off()

# plot(x = seq(min(model$grideval), max(model$grideval), by = 0.1),
            # y = 0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.1), 1.3,1) +
                # 0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.1), -1.3,1)+
				# 0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.1), -3.9,1)  + 
				# 0.25 * dnorm(seq(min(model$grideval), max(model$grideval), by = 0.1), 3.9,1))


# Results

# parameter=1/10
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.5934211
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.5089474
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.5934211
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.5934211
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.5615789
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.5615789

# parameter=1/5
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.7357895
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.6052632
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.7347368
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.7347368
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.7305263
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.7305263

# parameter=1/4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.8239474
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.7089474
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.8239474
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.8239474
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.7960526
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.7960526

# parameter=1/3
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.8726316
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.8055263
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.8721053
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.8726316
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.8528947
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.8594737

# parameter=1/2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.9247368
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.8555263
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.9244737
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.9263158
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.9002632
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.8844737

# parameter=1
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.9597368
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.9386842
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.9607895
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.9576316
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.9694737
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.9663158

# parameter=2
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.9286842
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.9473684
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.9276316
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.9276316
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.9371053
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.9392105

# parameter=4
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.8792105
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.9221053
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.8773684
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.8807895
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.9021053
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.9015789

# parameter=8
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.8481579
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.9142105
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.8415789
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.8397368
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.8655263
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.8663158

# parameter=16
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.8239474
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.8781579
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.8197368
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.8192105
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.8431579
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.84

# parameter=20
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.8263158
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.8721053
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.8142105
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.8152632
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.8460526
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.8452632

# parameter=50
# > print("RAND binder")
# [1] "RAND binder"
# > print(RAND(truepart,Bindermin)) 
# [1] 0.8178947
# > print("RAND VI")
# [1] "RAND VI"
# > print(RAND(truepart,VImin))
# [1] 0.8602632
# > print("RAND hier average Binder")
# [1] "RAND hier average Binder"
# > print(RAND(truepart,Hier1))
# [1] 0.8018421
# > print("RAND hier complete Binder")
# [1] "RAND hier complete Binder"
# > print(RAND(truepart,Hier2))
# [1] 0.8018421
# > print("RAND hier average VI")
# [1] "RAND hier average VI"
# > print(RAND(truepart,Hier3))
# [1] 0.8342105
# > print("RAND hier complete VI")
# [1] "RAND hier complete VI"
# > print(RAND(truepart,Hier4))
# [1] 0.8276316