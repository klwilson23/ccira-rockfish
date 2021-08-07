nage <- 120
M <- -0.0345

survivorship <- rep(1,nage)
for(i in 2:nage)
{
  if(i==nage)
  {
    survivorship[i] <- survivorship[i-1]/(1-exp(M))
  }else{
    survivorship[i] <- survivorship[i-1]*exp(M)
  }
}

plot(1:nage,survivorship,ylab="proportions at age",type="l")
plot(1:nage,survivorship/sum(survivorship),ylab="proportions at age",type="l",ylim=c(0,1))

req <- 480000 # 480,000 age 1 recruits
num_age <- req*survivorship/sum(survivorship)

linf <- 65.6
vbk <- 0.04
t0 <- -9.23

fec_a <- 6.54E-06
fec_b <- 4.043
a50 <- 13.2
a95 <- 29.2 # sean anderson value
mat_age <- rep(0,nage)
mat_age[8:nage] <- (1+exp(log(19)*((8:nage-a50)/(a50-a95))))^-1
len_age <- linf*(1-exp(-vbk*(1:nage-t0)))
fec_age <- fec_a*len_age^fec_b

eggs_age <- num_age*fec_age*mat_age
plot(1:nage,eggs_age/sum(eggs_age),ylab="Proportional eggs produced per age",type="l")

plot(num_age)
plot(num_age/max(num_age))
lines(fec_age/max(fec_age),type="l",col="blue",lwd=2)
