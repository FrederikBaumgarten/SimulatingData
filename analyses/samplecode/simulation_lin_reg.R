# Building up a mechanistic model
y = a+b*x + error #model formula with x being nitrogen concentration
error= rnorm(n , 0 , sigma ) #errordistribution


# specifying model parameters
a <- 30 #intercept
b <- 7 #slope (effect size)
sigma <- 5 #standard deviation of the errordistribution

n <- 50 #sample size
x<-rnorm(n , 10 , 4) #create data for nitrogen concentration
y<-a+b*x + rnorm(n , 0 , sigma ) #generate y values
fake<-data.frame (x, y)

plot(fake$x, fake$y, main="Fake data")

fit_1<-stan_glm(y~x , data=fake) # fit the model
print (fit_1, digits=2)
plot (fake$x, fake$y, main="Data and fitted regressionline")
a_hat<-coef(fit_1)[1]
b_hat<-coef(fit_1)[2]
abline(a_hat, b_hat)
