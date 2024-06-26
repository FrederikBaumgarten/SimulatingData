## Started 12 May 2024 ##
## Code is by E Loken and A Gelman ##
## Lizzie pulled from below website today (and added some par(mfrow))
# https://statmodeling.stat.columbia.edu/wp-content/uploads/2023/11/graph-codes-to-share-for-science-paper-final.txt
## See also: https://statmodeling.stat.columbia.edu/2023/11/11/simulations-of-measurement-error-and-the-replication-crisis-update/

# First just the original two plots, high power N = 3000, low power N = 50, true slope = .15

r <- .15
sims<-array(0,c(1000,4))
xerror <- 0.5
yerror<-0.5

for (i in 1:1000) {
x <- rnorm(50,0,1)
y <- r*x + rnorm(50,0,1)
xx<-lm(y~x)
sims[i,1]<-summary(xx)$coefficients[2,1]
x<-x + rnorm(50,0,xerror)
y<-y + rnorm(50,0,yerror)
xx<-lm(y~x)
sims[i,2]<-summary(xx)$coefficients[2,1]

x <- rnorm(3000,0,1)
y <- r*x + rnorm(3000,0,1)
xx<-lm(y~x)
sims[i,3]<-summary(xx)$coefficients[2,1]
x<-x + rnorm(3000,0,xerror)
y<-y + rnorm(3000,0,yerror)
xx<-lm(y~x)
sims[i,4]<-summary(xx)$coefficients[2,1]

}

par(mfrow=c(1,3))

plot(sims[,2] ~ sims[,1],ylab="Observed with added error",xlab="Ideal Study")
abline(0,1,col="red")

plot(sims[,4] ~ sims[,3],ylab="Observed with added error",xlab="Ideal Study")
abline(0,1,col="red")


# third graph

# run 2000 regressions at points between N = 50 and N = 3050 

r <- .15

propor <-numeric(31)
powers<-seq(50,3050,100)

xerror<-0.5
yerror<-0.5

for (j in 1:31)  {

sims<-array(0,c(1000,4))
for (i in 1:1000) {
x <- rnorm(powers[j],0,1)
y <- r*x + rnorm(powers[j],0,1)
xx<-lm(y~x)
sims[i,1:2]<-summary(xx)$coefficients[2,1:2]
x<-x + rnorm(powers[j],0,xerror)
y<-y + rnorm(powers[j],0,yerror)
xx<-lm(y~x)
sims[i,3:4]<-summary(xx)$coefficients[2,1:2]
}

# find significant observations (t test > 2) and then check proportion

temp<-sims[abs(sims[,3]/sims[,4])> 2,]

propor[j] <- table((abs(temp[,3]/temp[,4])> abs(temp[,1]/temp[,2])))[2]/length(temp[,1])

print(j)
}

plot(powers,propor,type="l",xlab="Sample Size",ylab="Prop where error slope greater",col="blue")
