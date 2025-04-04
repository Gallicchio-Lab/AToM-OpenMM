# cat r*/*.out > data
# R CMD BATCH uwham_analysis.R

#.libPaths("/home/emilio/R/x86_64-pc-linux-gnu-library/3.0/")
library("UWHAM")

bias.fcn <- function(epert, lam1, lam2, alpha, u0, w0){
# This is for the bias ilogistic potential
# (lambda2-lambda1) ln[1+exp(-alpha (u-u0))]/alpha + lambda2 u + w0
    ebias1 <- 0*epert
    if (alpha > 0) {
        ee <- 1 + exp(-alpha*(epert-u0))
        ebias1 <- (lam2 - lam1)*log(ee)/alpha
    }
    ebias1 + lam2*epert + w0
}

dbias.fcn <- function(epert, lam1, lam2, alpha, u0, w0){
# This is the derivative of the bias ilogistic potential
# (lambda2-lambda1)/(1+exp(-alpha (u-u0))) + lambda2
    ee <- 1. + exp(-alpha*(epert-u0))
    (lam2 - lam1)/ee + lam1
}


npot.fcn <- function(e0,epert, bet, lam1, lam2, alpha, u0, w0){ 
# This is the negative reduced energy 
# -beta*(U0+bias)
    -bet*(e0 + bias.fcn(epert, lam1, lam2, alpha, u0, w0))
}

uwham.r <- function(label,logQ,ufactormax,ufactormin=1){
  n <- dim(logQ)[1]
  m <- dim(logQ)[2]
  iniz <- array(0,dim=m) 
  uf <- ufactormax
  while(uf >= ufactormin & uf >= 1){
    mask <- seq(1,n,trunc(uf))
    out <- uwham(label=label[mask], logQ=neg.pot[mask,],init=iniz)
    show(uf)
    iniz <- out$ze
    uf <- uf/2
  }
  out$mask <- mask
  out
}

histw <-
function (x, w, xaxis, xmin, xmax, ymax, bar = TRUE, add = FALSE, 
            col = "black", dens = TRUE) 
{
  nbin <- length(xaxis)
  xbin <- cut(x, breaks = xaxis, include.lowest = T, labels = 1:(nbin -  1))
  y <- tapply(w, xbin, sum)
  y[is.na(y)] <- 0
  y <- y/sum(w)
  if (dens) 
    y <- y/(xaxis[-1] - xaxis[-nbin])
  if (!add) {
    plot.new()
    plot.window(xlim = c(xmin, xmax), ylim = c(0, ymax))
    axis(1, pos = 0)
    axis(2, pos = xmin)
  }
  if (bar == 1) {
    rect(xaxis[-nbin], 0, xaxis[-1], y)
  }
  else {
    xval <- as.vector(rbind(xaxis[-nbin], xaxis[-1]))
    yval <- as.vector(rbind(y, y))
    lines(c(min(xmin, xaxis[1]), xval, max(xmax, xaxis[length(xaxis)])), 
          c(0, yval, 0), lty = "11", lwd = 2, col = col)
  }
  invisible()
  list(y = y, breaks = xaxis)
}

args <- commandArgs(trailingOnly = F)
jobname <- sub("-","",args[length(args)-2])
mintimeid <- strtoi(sub("-","",args[length(args)-1]))
maxtimeid <- strtoi(sub("-","",args[length(args)  ]))

mintimeid
maxtimeid

#define states
tempt   <- c( 300 )
bet     <- 1.0/(0.001986209*tempt)
directn <-c(   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1 )
intermd <-c(   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 )
lambda1 <-c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 )
lambda2 <-c(0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00 )
alpha   <-c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 )
u0      <-c(110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110. )
w0      <-c(   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 )

nstates <- length(lambda1)
leg1istate <- which(intermd==1)[1]
leg2istate <- which(intermd==1)[2]

colnames <- c("stateid", "temperature", "direction", "lambda1", "lambda2", "alpha", "u0", "w0", "potE", "pertE") 
datafiles <- sprintf("r%d/%s.out",seq(0,length(lambda1)-1),jobname)
nfiles <- length(datafiles)
data <- read.table(datafiles[1])
colnames(data) <- colnames
data$timeid <- 1:length(data$stateid)
for ( i in 2:nfiles) {
    t <- read.table(datafiles[i])
    colnames(t) <- colnames
    t$timeid <- 1:length(t$stateid)
    data <- rbind(data,t)
}
data$bet <- 1.0/(0.001986209*data$temperature)
nsamples <- length(data$stateid)
samplesperreplica <- as.integer(nsamples/nstates)


#LEG1

data1 <- subset(data, stateid <= leg1istate - 1 & timeid >= mintimeid & timeid <= maxtimeid  )
mtempt <- length(bet)
leg1stateids <- 1:leg1istate
leg1stateids
mlam <- length(leg1stateids)
mlam
m <- mlam*mtempt
N <- length(data1$stateid)

#extract U0 values as U-bias
#this is relevant only if the states are at different temperatures
e0 <- data1$potE
for (i in 1:N) {
    e0[i] <- e0[i] - bias.fcn(data1$pertE[i],data1$lambda1[i],data1$lambda2[i],data1$alpha[i],data1$u0[i],data1$w0[i])
}

neg.pot <- matrix(0, N,m)
sid <- 1
# note the order of (be,te)
for (be in leg1stateids  ) {
     for (te in 1:mtempt) {
             neg.pot[,sid] <- npot.fcn(e0=e0,data1$pertE,bet[te],lambda1[be],lambda2[be],alpha[be],u0[be],w0[be])
             sid <- sid + 1
    }
}

#the alchemical state indexes start with 0, UWHAM's state labels start with 1
statelabels <- data1$stateid + 1

#runs UWHAM
out <- uwham.r(label=statelabels, logQ=neg.pot,ufactormax=1,ufactormin=1)
ze <- matrix(out$ze, nrow=mtempt, ncol=mlam)
-ze/bet
ve <- matrix(out$ve, nrow=mtempt, ncol=mlam)
sqrt(ve)/bet

dgbind1 <- (-ze[,mlam]/bet[]) - (-ze[,1]/bet[])
ddgbind1 <- sqrt(ve[,mlam]+ve[,1])/bet

dgbind1
ddgbind1


#estimate p0(u) lambda0 for leg 1
#uses a gaussian kernel estimator for each sample: g(u-uj)
#p0(u) = sum_j W0(j) g(u-uj)
#dp0(u)/du = (-1/sigma^2) sum_j W0(j) (u-uj) g(u-uj)
#and lambda0(u) = kbT d ln p0(u)/du = kbT (dp0(u)/du)/(p0(u))
ux <- seq(from = -20, to = 180, by = 5)
sigma <- 3.0
norm <- 1./sqrt(2.*pi*sigma**2)
sigma2 <- sigma**2
gg <- function(x, y, norm, s2){ norm*exp(-(x-y)**2/(2.*s2)) }
#matrix of M rows (u grid points) and N columns (samples)
gm <- outer(ux, data1$pertE, FUN = gg, norm, sigma2)
#this gives a column vector of p0(u) values for each u grid point
p0 <- gm %*% (out$W[,1]/N)
#do the same for derivative
ggd <- function(x, y, norm, s2){ ((y-x)/s2)*norm*exp(-(x-y)**2/(2.*s2)) }
gm <- outer(ux, data1$pertE, FUN = ggd, norm, sigma2)
dp0 <- gm %*% (out$W[,1]/N)
plot(ux, log(p0), type = "l")
plot(ux, (dp0/p0)/bet, type="l") 
dbf1 <- dbias.fcn(ux, 0., 0.4, 0.1, 105., 0.)
dbf2 <- dbias.fcn(ux, 0., 0.5, 0.1, 100., 0.)
lines(ux, dbf1)
lines(ux, dbf2)



#get plain be histograms at first temperature
umin <- min(data1$pertE)
umax <- max(data1$pertE)
hs <- hist(data1$pertE[ data1$stateid == mlam-1 ],plot=FALSE,breaks=10);
pmax = 1.2*max(hs$density)
plot(hs$mids,hs$density,type="l",xlim=c(umin,umax),ylim=c(0,pmax));
for ( i in 1:mlam ){ 
    hs <- hist(data1$pertE[ data1$stateid == i-1 ],plot=FALSE,breaks=10);
    lines(hs$mids,hs$density);
    outp <- cbind(hs$mids,hs$density);
    write(t(outp),file=sprintf("p-%d.dat",i-1),ncol=2)
}




#LEG2

data1 <- subset(data, stateid >= leg2istate - 1 & timeid >= mintimeid & timeid <= maxtimeid )
mtempt <- length(bet)
leg2stateids <- seq(from=nstates, to=leg2istate, by=-1)
leg2stateids
mlam <- length(leg2stateids )
mlam
m <- mlam*mtempt
N <- length(data1$stateid)

#extract U0 values as U-bias
#this is relevant only if the states are at different temperatures
e0 <- data1$potE
for (i in 1:N) {
    e0[i] <- e0[i] - bias.fcn(data1$pertE[i],data1$lambda1[i],data1$lambda2[i],data1$alpha[i],data1$u0[i],data1$w0[i])
}

neg.pot <- matrix(0, N,m)
sid <- 1
# note the order of (be,te)
for (be in leg2stateids ) {
     for (te in 1:mtempt) {
             neg.pot[,sid] <- npot.fcn(e0=e0,data1$pertE,bet[te],lambda1[be],lambda2[be],alpha[be],u0[be],w0[be])
             sid <- sid + 1
    }
}

#the alchemical state indexes in leg2 run backward
statelabels <- nstates - data1$stateid

#runs UWHAM
out <- uwham.r(label=statelabels, logQ=neg.pot,ufactormax=1,ufactormin=1)
ze <- matrix(out$ze, nrow=mtempt, ncol=mlam)
-ze/bet
ve <- matrix(out$ve, nrow=mtempt, ncol=mlam)
sqrt(ve)/bet


dgbind2 <- (-ze[,mlam]/bet[]) - (-ze[,1]/bet[])
ddgbind2 <- sqrt(ve[,mlam]+ve[,1])/bet

dgbind2
ddgbind2


#estimate p0(u) lambda0 for leg 2
ux <- seq(from = -20, to = 180, by = 5)
sigma <- 3.0
norm <- 1./sqrt(2.*pi*sigma**2)
sigma2 <- sigma**2
gg <- function(x, y, norm, s2){ norm*exp(-(x-y)**2/(2.*s2)) }
gm <- outer(ux, data1$pertE, FUN = gg, norm, sigma2)
p0 <- gm %*% (out$W[,1]/N)
ggd <- function(x, y, norm, s2){ ((y-x)/s2)*norm*exp(-(x-y)**2/(2.*s2)) }
gm <- outer(ux, data1$pertE, FUN = ggd, norm, sigma2)
dp0 <- gm %*% (out$W[,1]/N)
plot(ux, log(p0), type = "l")
plot(ux, (dp0/p0)/bet, type="l") 
dbf1 <- dbias.fcn(ux, 0., 0.4, 0.1, 105., 0.)
dbf2 <- dbias.fcn(ux, 0., 0.5, 0.1, 100., 0.)
lines(ux, dbf1)
lines(ux, dbf2)

#free energy difference
dgb <- dgbind1 - dgbind2
ddgb <- sqrt(ddgbind2*ddgbind2 + ddgbind1*ddgbind1)
maxsamples <- min(maxtimeid, samplesperreplica)
result <- sprintf("DDGb = %f +- %f range %d %d", dgb, ddgb, mintimeid, maxsamples)
write(result, "")
#noquote(result)

#get plain be histograms at first temperature
umin <- min(data1$pertE)
umax <- max(data1$pertE)
hs <- hist(data1$pertE[ data1$stateid == leg2istate - 1  ],plot=FALSE,breaks=10);
pmax = 1.2*max(hs$density)
plot(hs$mids,hs$density,type="l",xlim=c(umin,umax),ylim=c(0,pmax));
for ( i in nstates:leg2istate ){ 
    hs <- hist(data1$pertE[ data1$stateid == i-1 ],plot=FALSE,breaks=10);
    lines(hs$mids,hs$density);
    outp <- cbind(hs$mids,hs$density);
    write(t(outp),file=sprintf("p-%d.dat",i-1),ncol=2)
}
