## FSSA Example on Satelite Image Data
require(fda)
require(Rfssa)
require(Rssa)
## raw image data
NDVI=Jambi$NDVI
## Time vector
time <- Jambi$Date
NDVI[NDVI<0]<-NA
avgNDVI=sapply(1:448,function(x) mean(NDVI[,,x], na.rm=T))
time=Jambi$Date
plot(time,avgNDVI,lwd=2,cex.axis=1.2,cex.main=1.8,cex.lab=1.5,type='l',col='blue',xlab="Date",ylab="NDVI Value",main="Average NDVI Values",ylim=c(0,1))
grid(col='black',lwd=1.5)
## Decomposition stage of SSA
avgL=23
avgssa=ssa(x=avgNDVI,avgL)
plot(avgssa,type="values",numvalues = 15,main=list("Singular Values",cex=1.8),xlab=list("Index",cex=1.5),ylab=list("norms",cex=1.5),lwd=1)
plot(avgssa,type="wcor",main=list("W-Correlation",cex=1.8),groups=1:15,cex.axis=1.2)
plot(avgssa,type="vectors",idx=1:4,main=list("Left Singular Vectors",cex=1.8),par.strip.text=list(cex=1.6),lwd=3)
plot(avgssa,type="paired",idx = 1:4,idy=2:5,main=list("Paired-Plots",cex=1.8),par.strip.text=list(cex=1.6),lwd=3)
## Reconstruction stage of SSA
avgrecon=reconstruct(avgssa,groups=list(1,2:3,4:10))
plot(time,avgrecon[[1]],ylim=c(-0.5,0.5),lwd=2,type='l',lty=2,col='blue',xlab="Date",ylab="NDVI Value",main="Mean and Seasonal Components",cex.lab=1.5,cex.axis=1.2,cex.main=1.8)
points(time,avgrecon[[2]],lwd=2,col='blue',type='l',lty=1)
grid(col='black',lwd=2)
legend(x="topright",legend = c("Mean NDVI","Seasonal NDVI"),col = c("blue","blue"), lty=c(2,1),cex=1.5,lwd=c(2,2))

## kernel density estimation of pixel intensity
NDVI=Jambi$NDVI
D0_NDVI <- matrix(NA,nrow = 512, ncol = 448)
for(i in 1:448){
  D0_NDVI[,i] <- density(NDVI[,,i],from=0,to=1)$y
}
## define functional objects
d <- 11
basis <- create.bspline.basis(c(0,1),d)
u <- seq(0,1,length.out = 512)
y_NDVI <- smooth.basis(u,as.matrix(D0_NDVI),basis)$fd
## define functional time series
Y <- fts(y_NDVI,time=time)
plot(Y,xlab = "NDVI Value", ylab = "Relative Likelihood", main = "NDVI Kernel Density Functions")
plot(Y,type="heatmap",xlab = "NDVI Value", ylab = "Density",tlab="Date",main = "NDVI Density Functions")
L=45
## functional singular spectrum analysis
U=fssa(Y,L)
## original fts plots to choose from for upper left plot
## fssa plots
plot(U,d=15,type='values')
plot(U,d=15,type='wcor') # Have to build locally to update font size for w-plot
plot(U,d=4,type='vectors')
plot(U,d=4,type='lheats')
plot(U,d=4,type='lcurves')
# interesting reconstructions
recon <- freconstruct(U = U, group = list(c(1),c(2,3),c(4),c(1,4)))
plot(recon[[1]],type="heatmap",xlab = "NDVI Value", ylab = "Relative Likelihood",tlab="Date",main = "FSSA Mean Component")
plot(recon[[2]],type="heatmap",xlab = "NDVI Value", ylab = "Relative Likelihood",tlab="Date",main = "FSSA Yearly Component")
plot(recon[[3]],type="heatmap",xlab = "NDVI Value", ylab = "Relative Likelihood",tlab="Date",main = "FSSA Trend Component")
plot(recon[[4]],type="heatmap",xlab = "NDVI Value", ylab = "Relative Likelihood",tlab="Date",main = "FSSA Mean and Trend Components")
