### Started 12 May 2024 ###
## By Lizzie ##
## Pulling code for simulation examples ... ###

### Below is snippets of code in decsensSims.R and decsensSimsFstar.R ###
# posted at https://github.com/temporalecologylab/labgit/tree/master/projects/decsenspost

## Simulation of the declining sensitivities ##
## This code shows two TYPES of simulations (each followed by plotting)
## (1) GDD model is constant by world warms (driving to grandma's house)
## (2) Required GDD increases (a simple version of the chilling model, but you see the prediction of declining sensitivities is wrong, they go UP)

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(113)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("ailene", getwd()))>0) { 
} else setwd("~/Documents/git/projects/others/fredi/SimulatingData/analyses")

##########################
# The below sets up data #
# for showing what happens #
# with GDD model and warming #
##########################

# Make some data ... note that this runs 100 times for 20 sites, via a loop 

# Step 1.1: Set up years, days per year, temperatures, required GDD (fstar)
daysperyr <- 60
yearz <- 30
sitez <- 20 # reps
simsnum.maintext <- 100
degreez.maintext <- seq(0, 2, length.out=simsnum.maintext)
degreez.forsupp <- c(0, 0.5, 1, 1.5, 2, 2.5, 4, 7)
sigma <- 4
basetemp <- 4 # alpha_0
dailytempchange <- 0.1 # alpha_1
fstar <- 150 # beta

# Step 1.2: Build the data and calculate sensitivities (note that alpha_1 -- spring temperature increase is set to 0)
df <- data.frame(degwarm=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric())

for (i in degreez.maintext){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(basetemp, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(daysperyr, basetemp + i, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange # add daily temp increase
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > fstar)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(degwarm=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2])
       df <- rbind(df, dfadd)
    }
}


dfsupp <- data.frame(degwarm=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric())

for (i in degreez.forsupp){
   for (j in 1:sitez){
       yearly_expected_temp <- rep(basetemp, yearz)
       daily_temp <- sapply(yearly_expected_temp, function(x) rnorm(daysperyr, basetemp + i, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange  # add daily temp increase
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > fstar)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfsuppadd <- data.frame(degwarm=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2])
       dfsupp <- rbind(dfsupp, dfsuppadd)
    }
}

plot(simplelm.trunc~degwarm, data=dfsupp, ylab="Sensitivity (days/C or log(days)/log(C)", xlab="Degree warming")
points(simplelm~degwarm, data=dfsupp, pch=16, col="gray")
points(loglm.trunc~degwarm, data=dfsupp, col="dodgerblue")
points(loglm~degwarm, data=dfsupp, col="dodgerblue", pch=16)
plot(perlm~degwarm, data=dfsupp, col="firebrick")

plot(abs(simplelm)~degwarm, data=dfsupp, col="lightgrey",
    ylab="Abs(Sensitivity (days/C or log(days)/log(C))", xlab="Degree warming")
dfsupp$degwarmJitter <- dfsupp$degwarm + 0.05
points(abs(loglm)~degwarmJitter, data=dfsupp, col="dodgerblue", cex=0.8)


##############################
## Plotting for basic model ##
##############################

# Summarize the sims
mean.sims <- aggregate(dfsupp[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], dfsupp["degwarm"], FUN=mean)
sd.sims <- aggregate(dfsupp[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], dfsupp["degwarm"], FUN=sd)

cexhere <- 0.95
pdf(file.path("figures/decsens_basicsims.pdf"), width = 6, height = 8)
par(mfrow=c(2,1))
par(xpd=FALSE)
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-6, -0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C to leafout)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="",
     bty="l", mgp=c(1.5,.5,0), tck=-.01)
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm.trunc[i]
  sdhere <- sd.sims$simplelm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm.trunc[i]
  sdhere <- sd.sims$loglm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Linear (untransformed)",
    "Non-linear (logged)"), cex=1, bty="n")
plot(x=NULL,y=NULL, xlim=c(-0.5, 8), ylim=c(-6, -0.1),
     ylab=expression(paste("Estimated sensitivity (days/", degree, "C over window)"), sep=""),
     xlab=expression(paste("Warming (", degree, "C)")), main="",
     bty="l", mgp=c(1.5,.5,0), tck=-.01)
# abline(h=0, lty=2, col="darkgrey")
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="darkblue")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="darkblue")
  }
for(i in 1:length(unique(mean.sims$degwarm))){
  pos.x <- mean.sims$degwarm[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col="salmon")
  points(pos.x, pos.y, cex=cexhere, pch=19, col="salmon")
  }
legend("bottomright", pch=c(19, 19), col=c("darkblue", "salmon"), legend=c("Linear (untransformed)", "Non-linear (logged)"),
   cex=1, bty="n")
dev.off()


##########################
# The below sets up data #
# for showing what happens #
# with GDD model when GDD #
# increases #
##########################


# Step 2.1: Set up additional items
simsnum <- 40
fstarsims <- seq(100, 300, length.out=simsnum)
sigma <- 4

# Step 2.2: Build the data and calculate sensitivities 
df <- data.frame(fstar=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric(), perlm=numeric(),
    simplelm.trunc=numeric(), loglm.trunc=numeric(), varx=numeric(), covarxy=numeric(), covarlogxy=numeric(),
    varx.trunc=numeric(), covarxy.trunc=numeric(), covarlogxy.trunc=numeric())

for (i in fstarsims){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) rnorm(daysperyr, basetemp, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > i)))
       yearly_temp <- colMeans(daily_temp)
       yearly_temp_trunc <- sapply(1:ncol(daily_temp), function(x) mean(daily_temp[1:leafout_date[x], x]))
       per_leafout_date <- leafout_date/mean(leafout_date)
       per_yearly_temp <- yearly_temp/mean(yearly_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(fstar=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2],
           perlm=coef(lm(per_leafout_date~per_yearly_temp))[2],
           simplelm.trunc=coef(lm(leafout_date~yearly_temp_trunc))[2],
           loglm.trunc=coef(lm(log(leafout_date)~log(yearly_temp_trunc)))[2],
           varx=var(yearly_temp), covarxy=cov(leafout_date, yearly_temp),
           covarlogxy=cov(log(leafout_date), log(yearly_temp)),
           varx.trunc=var(yearly_temp_trunc), covarxy.trunc=cov(leafout_date, yearly_temp_trunc),
           covarlogxy.trunc=cov(log(leafout_date), log(yearly_temp_trunc))
           )
       df <- rbind(df, dfadd)
    }
}

dfsm <- data.frame(fstar=numeric(), rep=numeric(), simplelm=numeric(), loglm=numeric())

for (i in fstarsims){
   for (j in 1:sitez){
       daily_temp <- sapply(rep(NA, yearz), function(x) rnorm(daysperyr, basetemp, sigma))
       daily_temp <- daily_temp + c(1:daysperyr)*dailytempchange 
       leafout_date <- sapply(1:ncol(daily_temp), function(x) min(which(cumsum(daily_temp[,x]) > i)))
       yearly_temp <- colMeans(daily_temp)
       plot(yearly_temp, leafout_date, pch=20)
       dfadd <- data.frame(fstar=i, rep=j, simplelm=coef(lm(leafout_date~yearly_temp))[2],
           loglm=coef(lm(log(leafout_date)~log(yearly_temp)))[2])
       dfsm <- rbind(dfsm, dfadd)
    }
}

mean.sims.sm <- aggregate(dfsm[c("simplelm", "loglm")], dfsm["fstar"], FUN=mean)


##############################################
## Plotting for increasing Fstar (GDD) sims ##
##############################################

pdf(file.path("figures/basicsims_fstaronly_varcov.pdf"), width = 10, height = 5)
par(mfrow=c(2,4))
plot(simplelm.trunc~fstar, data=df, xlab="thermal sum", ylab="lm sensitivity (temp until leafout)")
plot(varx.trunc~fstar, data=df, xlab="thermal sum", ylab="var(temp until leafout)")
plot(covarxy.trunc~fstar, data=df, xlab="thermal sum", ylab="covar(leafout day, temp until leafout)")
plot(covarlogxy.trunc~fstar, data=df, xlab="thermal sum", ylab="covar(log(leafout day), log(temp until leafout))")
plot(simplelm~fstar, data=df, xlab="thermal sum", ylab="lm sensitivity (temp over window)")
plot(varx~fstar, data=df, xlab="thermal sum", ylab="var(temp over window)")
plot(covarxy~fstar, data=df, xlab="thermal sum", ylab="covar(leafout day, temp over window)")
plot(covarlogxy~fstar, data=df, xlab="thermal sum", ylab="covar(log(leafout day), (log(temp over window))")
dev.off()

# Summarize the fstar sims
mean.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["fstar"], FUN=mean)
sd.sims <- aggregate(df[c("simplelm", "loglm", "perlm", "simplelm.trunc", "loglm.trunc")], df["fstar"], FUN=sd)

colz <- c("blue4", "salmon")
colzalpha <- adjustcolor(colz, alpha.f = 0.7)
cexhere <- 0.75
cexhere <- 1.2
cextext <- 0.75
jitterpep <- -0.04
pdf(file.path("figures/basicsims_fstaronly.pdf"), width = 7.5, height = 11)
par(xpd=FALSE)
par(mfrow=c(2,1))
par(mar=c(5,5,2,2))
plot(x=NULL,y=NULL, xlim=c(90, 310), ylim=c(-8, 0.1), yaxt="n",
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Thermal sum required (", degree, "C)")), main="", cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.5,0), tck=-.01)
axis(2,seq(-8,0,1),las=2)
tempsteps <- simsnum
tempdiffplot <- c(0,1)
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$simplelm.trunc[i]
  sdhere <- sd.sims$simplelm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[1])
}
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$loglm.trunc[i]
  sdhere <- sd.sims$loglm.trunc[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[2])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[2])
}
legend("topright", pch=c(19, 19), col=c(colzalpha[1], colzalpha[2]),legend=c("Linear (untransformed)", "Non-linear (logged)"),
   cex=1, bty="n")
plot(x=NULL,y=NULL, xlim=c(90, 310), ylim=c(-6.6, 0.1), yaxt="n",
     ylab=expression(paste("Estimated sensitivity"), sep=""),
     xlab=expression(paste("Thermal sum required (", degree, "C)")), main="", cex.lab=1.2,
     bty="l", mgp=c(1.5, 0.5, 0), tck=-.01)
axis(2,seq(-6,0,1),las=2)
tempsteps <- simsnum
tempdiffplot <- c(0,1)
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$simplelm[i]
  sdhere <- sd.sims$simplelm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[1])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[1])
}
for(i in 1:tempsteps){
  pos.x <- mean.sims$fstar[i]
  pos.y <- mean.sims$loglm[i]
  sdhere <- sd.sims$loglm[i]
  lines(x=rep(pos.x, 2), y=c(pos.y-sdhere, pos.y+sdhere), col=colzalpha[2])
  points(pos.x, pos.y, cex=cexhere, pch=19, col=colzalpha[2])
}
dev.off()
