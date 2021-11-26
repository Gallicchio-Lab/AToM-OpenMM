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



data.t <- read.table("repl.cycle.state.temp.lambda1.lambda2.alpha.u0.w0.epot.epert.dat")

states =   data.t$V3
temps =    data.t$V4
lambda1s = data.t$V5
lambda2s = data.t$V6
alphas =   data.t$V7
u0s =      data.t$V8
w0s =      data.t$V9
potEs =    data.t$V10
eperts =   data.t$V11


#parameters for ilogistic expected
tempt   <- c( 300 )
bet     <- 1.0/(0.001986209*tempt)
lambda1 <-c( 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
lambda2 <-c( 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
alpha   <-c( 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
u0      <-c( 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
w0coeff <-c( 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)


mtempt <- length(bet)
mlam <- length(lambda1)
m <- mlam*mtempt
N <- length(data.t$V1)

#extract U0 values as U-bias
#this is relevant only if the states are at different temperatures
e0 <- potEs
for (i in 1:N) {
    e0[i] <- e0[i] - bias.fcn(eperts[i],lambda1s[i],lambda2s[i],alphas[i],u0s[i],w0s[i])
}

neg.pot <- matrix(0, N,m)
sid <- 1
# note the order of (be,te)
for (be in 1:mlam) {
     for (te in 1:mtempt) {
             neg.pot[,sid] <- npot.fcn(e0=e0,eperts,bet[te],lambda1[be],lambda2[be],alpha[be],u0[be],w0coeff[be])
             sid <- sid + 1
    }
}

#the alchemical state indexes start with 0, UWHAM's state labels start with 1
statelabels <- states + 1

#runs UWHAM
out <- uwham.r(label=statelabels, logQ=neg.pot,ufactormax=1,ufactormin=1)
ze <- matrix(out$ze, nrow=mtempt, ncol=mlam)
-ze/bet
ve <- matrix(out$ve, nrow=mtempt, ncol=mlam)
sqrt(ve)/bet

dgbind <- (-ze[,mlam]/bet[]) - (-ze[,1]/bet[])
ddgbind <- sqrt(ve[,mlam]+ve[,1])/bet

#average energy
usl1 <- eperts[states == mlam - 1]
de <- mean(usl1)
sde <- sd(usl1)
dde <- sde/sqrt(length(usl1))
ub <- de + bet*sde*sde

#DGbind as a function of temperature
dgbind
sink("result.log")
cat("DGb = ", dgbind[1]," +- ",ddgbind[1]," DE = ", de, " +- ",dde,"\n")
sink()

#free energy profile at first temperature
dglambda <- cbind(-ze[1,]/bet[1],sqrt(ve[1,])/bet[1])
plot(-ze[1,]/bet[1], type="l")
write(t(dglambda),file="dglambda.dat",ncol=2)

#get plain be histograms at first temperature
umin <- min(eperts)
umax <- max(eperts)
hs <- hist(eperts[states == mlam-1 ],plot=FALSE,breaks=10);
pmax = 1.2*max(hs$density)
hs <- hist(eperts[states == 0 ],plot=FALSE,breaks=10);
plot(hs$mids,hs$density,type="l",xlim=c(umin,umax),ylim=c(0,pmax));
for ( i in 1:(mlam-1) ){ 
    hs <- hist(eperts[states == i ],plot=FALSE,breaks=10);
    lines(hs$mids,hs$density);
    outp <- cbind(hs$mids,hs$density);
    write(t(outp),file=sprintf("p-%d.dat",i),ncol=2)
}
