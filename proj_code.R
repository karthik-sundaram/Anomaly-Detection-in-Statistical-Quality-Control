setwd("C:/Users/Karthik/Desktop/Sem 1/ISEN 614/project")
dat614=read.csv("CleanDat.csv",header=F,na.strings="?")

str(dat614)
#summary(dat614)
#dat614$V1=as.numeric(dat614$V1)
fix(dat614)
install.packages("MSQC")
library(MSQC)
mult.chart(dat614, type = "chi", alpha = 0.05)


#Phase I

mult.chart(type = "t2", dat614)
dat614=dat614[-c(26,256,274,341,514,836,839,969),] #1st iteration OC
dat614=dat614[-c(772),] #2nd iteration OC




install.packages(qcr)
library(qcr)
cusum_chart=mqcs.mcusum(dat614,k=1.5,h=20) #extrapolate the h/UCL value
plot(cusum_chart)
dat614=dat614[-c(764:884),] #k=2
dat614=dat614[-c(739:774),] #k=2
dat614=dat614[-c(164:364,764:924),] #k=1.5
dat614=dat614[-c(150:170,530:580),] #k=1.5


rm(list=ls())
