source("XIMCOR.R")
library(MASS)
library(independence)
library(TauStar)
library(ggplot2)
library(XICOR)
library(gridExtra)
library(dplyr)
library(MultiRNG)

# Table1 + Table2

# Revised Chatterjee’s
XIMpower<-function(n,M,rho,simtype){
  ntime<-1000
  B<-10000
  alpha<-0.05
  nrej<-0
  time.compute<-0
  coef.sim<-XIMsim(n,M,B)
  mu<-c(0,0)
  Sigma<-matrix(c(1,rho/sqrt(n),rho/sqrt(n),1),2,2)
  for (time in 1:ntime){
    if (simtype=="gaussian"){
      mysample<-mvrnorm(n,mu,Sigma)
    }
    if (simtype=="uniform"){
      mysample<-draw.d.variate.uniform(n,2,Sigma)
    }
    xvec <- mysample[,1]
    yvec <- mysample[,2]
    time.start<-Sys.time()
    coef.stat <- XIMstat(xvec,yvec,M)
    time.end<-Sys.time()
    time.compute<-time.compute+(time.end-time.start)
    nrej <- nrej + XIMtestT(coef.stat,coef.sim,alpha)
  }
  return(data.frame(n=n,M=M,rho=rho,rej=nrej/ntime,time=time.compute/ntime))
}


result<-NULL
for (n in c(1000,2000,5000)){
  for (M in c(1,20,100,200,n/2)){
    for (rho in c(0,1,2,5)){
      set.seed(123)
      result<-rbind(result,XIMpower(n,M,rho,"gaussian"))
      print(paste(n,M,rho,"complete"))
    }
  }
}


# Hoeffding's
myhoeffding<-function(n,rho,simtype){
  ntime<-1000
  alpha<-0.05
  nrej<-0
  time.compute<-0
  mu<-c(0,0)
  Sigma<-matrix(c(1,rho/sqrt(n),rho/sqrt(n),1),2,2)
  for (time in 1:ntime){
    if (simtype=="gaussian"){
      mysample<-mvrnorm(n,mu,Sigma)
    }
    if (simtype=="uniform"){
      mysample<-draw.d.variate.uniform(n,2,Sigma)
    }
    xvec <- mysample[,1]
    yvec <- mysample[,2]
    time.start<-Sys.time()
    ordering = relative.order(xvec, yvec) - 1
    coef.stat <- 36*(n-1)*.calc.hoeffding(ordering)
    time.end<-Sys.time()
    time.compute<-time.compute+(time.end-time.start)
    nrej <- nrej + (coef.stat>=qhoeffding)
  }
  return(data.frame(n=n,rho=rho,rej=nrej/ntime,time=time.compute/ntime))
}

alpha<-0.05
qhoeffding<-qHoeffInd(1-alpha)
result.hoeffding<-NULL
for (n in c(1000,2000,5000)){
  for (rho in c(0,1,2,5)){
    set.seed(123)
    result.hoeffding<-rbind(result.hoeffding,myhoeffding(n,rho,"gaussian"))
    print(paste(n,rho,"complete"))
  }
}


# Pearson's
mypearson<-function(n,rho,simtype){
  ntime<-1000
  alpha<-0.05
  nrej<-0
  time.compute<-0
  mu<-c(0,0)
  Sigma<-matrix(c(1,rho/sqrt(n),rho/sqrt(n),1),2,2)
  for (time in 1:ntime){
    if (simtype=="gaussian"){
      mysample<-mvrnorm(n,mu,Sigma)
    }
    if (simtype=="uniform"){
      mysample<-draw.d.variate.uniform(n,2,Sigma)
    }
    xvec <- mysample[,1]
    yvec <- mysample[,2]
    time.start<-Sys.time()
    coef.stat<-cor(xvec,yvec)
    time.end<-Sys.time()
    time.compute<-time.compute+(time.end-time.start)
    coef.stat<-cor.test(xvec,yvec)
    nrej <- nrej + (coef.stat$p.value <= alpha)
  }
  return(data.frame(n=n,rho=rho,rej=nrej/ntime,time=time.compute/ntime))
}


result.pearson<-NULL
for (n in c(1000,2000,5000)){
  for (rho in c(0,1,2,5)){
    set.seed(123)
    result.pearson<-rbind(result.pearson,mypearson(n,rho,"gaussian"))
    print(paste(n,rho,"complete"))
  }
}


# Deb's
debcoef<-function(xvec,yvec,M){
  n <- length(xvec)
  M0 <- floor(M/2)
  M1 <- ceiling(M/2)
  xrank <- rank(xvec, ties.method = "random")
  yrank <- rank(yvec, ties.method = "random")
  ord <- order(xrank)
  yrank <- yrank[ord]
  coef.sum <- 0
  for (i in (M0+1):(n-M1)){
    coef.sum.temp <- sum(pmin(yrank[i],yrank[(i-M0):(i+M1)])) - yrank[i]
    coef.sum <- coef.sum + coef.sum.temp
  }
  if (M0 > 0){
    for (i in 1:M0){
      coef.sum.temp <- sum(pmin(yrank[i],yrank[1:(M+1)])) - yrank[i]
      coef.sum <- coef.sum + coef.sum.temp
    }
  }
  for (i in ((n-M1+1):n)){
    coef.sum.temp <- sum(pmin(yrank[i],yrank[(n-M):n])) - yrank[i]
    coef.sum <- coef.sum + coef.sum.temp
  }
  coef.value <- coef.sum/(n*M)
  return(coef.value)
}

debcoef.sim<-function(n,M,B){
  coef.sim<-rep(0,B)
  M0 <- floor(M/2)
  M1 <- ceiling(M/2)
  for (b in 1:B){
    yrank <- sample(1:n,n)
    coef.sum <- 0
    for (i in (M0+1):(n-M1)){
      coef.sum.temp <- sum(pmin(yrank[i],yrank[(i-M0):(i+M1)])) - yrank[i]
      coef.sum <- coef.sum + coef.sum.temp
    }
    if (M0 > 0){
      for (i in 1:M0){
        coef.sum.temp <- sum(pmin(yrank[i],yrank[1:(M+1)])) - yrank[i]
        coef.sum <- coef.sum + coef.sum.temp
      }
    }
    for (i in ((n-M1+1):n)){
      coef.sum.temp <- sum(pmin(yrank[i],yrank[(n-M):n])) - yrank[i]
      coef.sum <- coef.sum + coef.sum.temp
    }
    coef.value <- coef.sum/(n*M)
    coef.sim[b]<-coef.value
  }
  return(coef.sim)
}

debpower<-function(n,M,rho,simtype){
  ntime<-1000
  B<-10000
  alpha<-0.05
  nrej<-0
  coef.sim<-debcoef.sim(n,M,B)
  mu<-c(0,0)
  Sigma<-matrix(c(1,rho/sqrt(n),rho/sqrt(n),1),2,2)
  for (time in 1:ntime){
    if (simtype=="gaussian"){
      mysample<-mvrnorm(n,mu,Sigma)
    }
    if (simtype=="uniform"){
      mysample<-draw.d.variate.uniform(n,2,Sigma)
    }
    xvec <- mysample[,1]
    yvec <- mysample[,2]
    coef.stat <- debcoef(xvec,yvec,M)
    nrej <- nrej + XIMtestT(coef.stat,coef.sim,alpha)
  }
  return(data.frame(n=n,M=M,rho=rho,rej=nrej/ntime))
}

result.deb<-NULL
for (n in c(1000,2000,5000)){
  for (M in c(1,20,100,200,n/2)){
    for (rho in c(0,1,2,5)){
      set.seed(123)
      result.deb<-rbind(result.deb,debpower(n,M,rho,"gaussian"))
      print(paste(n,M,rho,"complete"))
    }
  }
}


# Figure2
set.seed(123)
coef.true<-NULL
for (rho in seq(0,1,0.01)){
  n<-1000000
  mu<-c(0,0)
  Sigma<-matrix(c(1,rho,rho,1),2,2)
  mysample<-mvrnorm(n,mu,Sigma)
  xvec <- mysample[,1]
  yvec <- mysample[,2]
  XI<-calculateXI(xvec,yvec)
  coef.true<-rbind(coef.true,data.frame(rho=rho,xi=XI))
}

plot.as<-function(n,M){
  ntime<-1000
  result.as<-NULL
  for (rho in c(0,0.2,0.4,0.6,0.8,1)){
    for (time in 1:ntime){
      mu<-c(0,0)
      Sigma<-matrix(c(1,rho,rho,1),2,2)
      mysample<-mvrnorm(n,mu,Sigma)
      xvec <- mysample[,1]
      yvec <- mysample[,2]
      coef.stat <- XIMcalculate(xvec,yvec,M)
      result.as<-rbind(result.as,data.frame(rho=rho, coef = coef.stat))
    }
  }
  p<-ggplot()+
    geom_boxplot(data = result.as, aes(x=rho, y=coef, group = rho))+
    stat_summary(data = result.as, aes(x=rho, y=coef, group = rho),
                 fun = mean, colour="blue", geom="line", group= rho, linetype = "longdash")+
    geom_line(data = coef.true, aes(x=rho, y=xi), colour = "red")+
    scale_x_continuous(breaks=seq(0,1,0.2))+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(1,1,1,1))
  return(p)
}

p1<-plot.as(1000,1)
p2<-plot.as(2000,1)
p3<-plot.as(5000,1)
p4<-plot.as(1000,20)
p5<-plot.as(2000,20)
p6<-plot.as(5000,20)
p7<-plot.as(1000,100)
p8<-plot.as(2000,100)
p9<-plot.as(5000,100)
p10<-plot.as(1000,200)
p11<-plot.as(2000,200)
p12<-plot.as(5000,200)
ggsave("p1.pdf",plot = p1, dpi = 72)
ggsave("p2.pdf",plot = p2, dpi = 72)
ggsave("p3.pdf",plot = p3, dpi = 72)
ggsave("p4.pdf",plot = p4, dpi = 72)
ggsave("p5.pdf",plot = p5, dpi = 72)
ggsave("p6.pdf",plot = p6, dpi = 72)
ggsave("p7.pdf",plot = p7, dpi = 72)
ggsave("p8.pdf",plot = p8, dpi = 72)
ggsave("p9.pdf",plot = p9, dpi = 72)
ggsave("p10.pdf",plot = p10, dpi = 72)
ggsave("p11.pdf",plot = p11, dpi = 72)
ggsave("p12.pdf",plot = p12, dpi = 72)


# Figure1
gamma<-seq(0,1,0.001)
beta<-pmax(pmin(3/2*gamma-1/2,1/2*gamma),pmin(1/4+1/4*gamma,-1/4*gamma+1/2))
db<-data.frame(gamma=gamma,beta=beta)
ggplot()+
  geom_line(data=db%>%filter(gamma<=1/4), aes(x = gamma, y = beta), size = 2) +
  geom_line(data=db%>%filter(gamma>1/4), aes(x = gamma, y = beta), linetype = "dotted", size = 2) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"))
ggsave("db.pdf", dpi = 72)



