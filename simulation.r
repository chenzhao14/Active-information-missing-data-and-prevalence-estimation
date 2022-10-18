library(dplyr)
library(LaplacesDemon)
library(ggpubr)
library(ggplot2)

get_population_sympton_disease <- function(n,rho1,rho11,p0) {
  population <- data.frame(id = 1:n,symptom = rbinom(n,1,rho1),disease=0)
  rho01 <- p0-rho11
  rho0 <- 1-rho1
  rho10 <- rho1-rho11
  rho00 <- rho0-rho01
  p0_1 <- rho11/rho1
  p0_0 <- rho01/rho0
  sym_group <- filter(population, symptom==1)
  asym_group <- filter(population, symptom==0)
  n_1 <- floor(n*rho11)
  n_0 <- floor(n*rho01)
  sym_group[sample(nrow(sym_group),n_1),]$disease <-1
  asym_group[sample(nrow(asym_group),n_0),]$disease <-1
  population <- rbind(sym_group,asym_group)
  return(population)
}


test_group <- function(population,pi){
  a<-split(population,f = population[,c("symptom")])
  tested_group <- c()
  for(i in 1:length(pi)){
    tested_group <- rbind(tested_group,sample_frac(a[[i]],pi[i]))
  }
  return(tested_group)
}

get_v3 <- function(N,rho1,tested_group){
  N_T1 <- nrow(filter(tested_group,symptom == 1))
  N_T0 <- nrow(filter(tested_group,symptom == 0))
  N_T11 <- nrow(filter(tested_group,symptom == 1,disease==1))
  N_T01 <- nrow(filter(tested_group,symptom == 0,disease==1))
  rho0 <- 1-rho1
  a<-rho1*(1-(N_T1)/(N*rho1))/((N_T1)/(N*rho1))
  b<-rho0*(1-N_T0/(N*rho0))/(N_T0/(N*rho0))
  v31<- a*N_T11/N_T1*(1-N_T11/N_T1)
  v30<- b*N_T01/N_T0*(1-N_T01/N_T0)
  return(v31+v30)
}


get_test_prevalence <- function(n,tested_group){
  N <- n
  N_T <- nrow(tested_group)
  N_T1 <- nrow(filter(tested_group, symptom == 1))
  N_T0 <- nrow(filter(tested_group, symptom == 0))
  N_T11 <- nrow(filter(tested_group, symptom == 1, disease ==1))
  N_T_1 <- nrow(filter(tested_group, disease ==1))
  N_T01 <- nrow(filter(tested_group, symptom == 0, disease ==1))
  p_T1s <- N_T1/N_T
  rho_1_hat <- p_T1s/2*(N_T/N+1)
  p0_1 <- rho_1_hat*N_T11/N_T1
  p0_0 <- (1-rho_1_hat)*N_T01/N_T0
  sample <- N_T_1/N_T
  res <- p0_1+p0_0
  return(c(res,sample))
}


run_model_mar <- function(n,p0,rho1,pi,rho11){
  population <- get_population_sympton_disease(n,rho1,rho11,p0)
  tested_group <- test_group(population,pi)
  test_patient_prev <- get_test_prevalence(n,tested_group)
  v3 <- get_v3(n,rho1,tested_group)
  sigma <- sqrt(v3/n)
  moe <- qnorm(.975)*sigma/(test_patient_prev*(1-test_patient_prev))
  ci.upper <- invlogit(logit(test_patient_prev)+moe)
  ci.lower <- invlogit(logit(test_patient_prev)-moe)
  return(c(p0_hat = test_patient_prev, ci.upper = ci.upper, ci.lower = ci.lower))
}

repeat_function <- function(m,x){
  res <- matrix(0,ncol=3,nrow = m)
  i = 1
  repeat{
    if(i>m)
      break
    res[i,] <- run_model_mar(n=x,p0=0.2,rho1=0.2,pi=c(0.1,0.9),rho11=0.15)
    i = i+1
  }
  return(res)
}

m=500

set.seed(114)

outcome.ex3<-sapply(c(1000, 10000, 100000, 1000000), function(i) repeat_function(m=500,x=i))

outcome.ex3.df<-as.data.frame(outcome.ex3)



p0_hat.est <- outcome.ex3.df[1:m,]
p0_hat.sample <- outcome.ex3.df[(m+1):(2*m),]

colm <- colMeans(p0_hat.est)
colm2 <- colMeans(p0_hat.sample)
##
log(colMeans(p0_hat.est)/0.2)
log(colMeans(p0_hat.sample)/0.2)
log(colMeans(p0_hat.est)/colMeans(p0_hat.sample))
##
apply(p0_hat.sample,2,sd)

colMeans(log(p0_hat.sample/0.2))
colMeans(log(p0_hat.est/p0_hat.sample))

mse.est <- matrix(0,nrow = 500,ncol = 4)
for (i in 1:4) {
  mse.est[,i] <- (p0_hat.est[,i]-0.2)^2
}
                    
                    
p0_hat.value <- outcome.ex3.df[1:m,]
p0_hat.upper <- outcome.ex3.df[(m+1):(2*m),]
p0_hat.lower <- outcome.ex3.df[(2*m+1):(3*m),]

p0_hat.ci <- cbind(p0_hat.value,p0_hat.upper,p0_hat.lower)
colnames(p0_hat.ci) <- c("v1","v2","v3","v4","v1.u","v2.u","v3.u","v4.u","v1.l","v2.l","v3.l","v4.l")


poN_10000<-mean(p0_hat.ci$v1)
poN_100000<-mean(p0_hat.ci$v2)
poN_1000000<-mean(p0_hat.ci$v3)
poN_10000000<-mean(p0_hat.ci$v4)


g1<-ggplot(p0_hat.ci,aes(x=1:nrow(p0_hat.ci),y=v1))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v1.u, ymax=v1.l))+
  ggtitle("95% CI of MAR when population = 1000")+
  geom_hline(yintercept=poN_10000, linetype="dashed", color = "red")+
  ylim(0.000, 0.03)

g2<-ggplot(p0_hat.ci, aes(x=1:nrow(p0_hat.ci),y=v2))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v2.u, ymax=v2.l))+
  ggtitle("95% CI of MAR when population = 10000")+
  geom_hline(yintercept=poN_100000, linetype="dashed", color = "red")+
  ylim(0.007, 0.015)

g3<-ggplot(p0_hat.ci, aes(x=1:nrow(p0_hat.ci),y=v3))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v3.u, ymax=v3.l))+
  ggtitle("95% CI of MAR when population = 100000")+
  geom_hline(yintercept=poN_1000000, linetype="dashed", color = "red")+
  ylim(0.007, 0.015)

g4<-ggplot(p0_hat.ci, aes(x=1:nrow(p0_hat.ci),y=v4))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v4.u, ymax=v4.l))+
  ggtitle("95% CI of MAR when population = 1000000")+
  geom_hline(yintercept=poN_10000000, linetype="dashed", color = "red")+
  ylim(0.007, 0.015)

rho5.ci<-ggpubr::ggarrange(g1,g2,g3,g4,ncol=2,nrow=2)

rho5.ci



