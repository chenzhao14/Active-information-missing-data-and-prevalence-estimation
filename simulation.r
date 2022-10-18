library(dplyr)
library(LaplacesDemon)
library(ggpubr)
library(ggplot2)


## Generate population
get_population_sympton <- function(n,rho1) {
  population <- data.frame(id = 1:n, symptom = rbinom(n,1,rho1))
  return(population)
}


## People with Symptoms, rho1 = probability with symptom (rho10+rho11)
disease_group <- function(population,p0,rho11,rho1){
  sym_group <- filter(population, symptom==1)
  asym_group <- filter(population, symptom==0)
  ratio <- rho11/rho1
  sym_group <- sym_group %>% 
    mutate(disease = rbinom(nrow(sym_group),1,ratio))
  n_sysm_non_disease = nrow(population)*p0 - nrow(filter(sym_group,disease==1))
  row_name = sample(nrow(asym_group),n_sysm_non_disease)
  asym_group$disease=0
  asym_group[row_name,]$disease=1
  population_s_d = rbind(sym_group,asym_group)
  return(population_s_d)
}

## Get tested group

test_group <- function(population,pi){
  a<-split(population,f = population[,c("symptom")])
  tested_group <- c()
  for(i in 1:length(pi)){
    tested_group <- rbind(tested_group,sample_frac(a[[i]],pi[i]))
  }
  return(tested_group)
}



## Calculate Testing Prevalence
get_test_prevalence <- function(tested_group,rho1){
  N_T11 <- nrow(filter(tested_group,symptom==1, disease==1))
  N_T1 <- nrow(filter(tested_group,symptom==1))
  N_T01 <- nrow(filter(tested_group,symptom==0, disease==1))
  N_T0 <- nrow(filter(tested_group,symptom==0))  
  res <- (N_T11/N_T1)*rho1+(N_T01/N_T0)*(1-rho1)
  return(res)
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



## Run programming
run_model_mar <- function(n,p0,rho1,pi,rho11){
  population <- get_population_sympton(n,rho1)
  population.symptom <- disease_group(population,p0,rho11,rho1)
  tested_group <- test_group(population.symptom,pi)
  res <- get_test_prevalence(tested_group,rho1)
  v3 <- get_v3(n,rho1,tested_group)
  sigma <- sqrt(v3/n)
  moe <- qnorm(.975)*sigma/(res*(1-res))
  ci.upper <- invlogit(logit(res)+moe)
  ci.lower <- invlogit(logit(res)-moe)  
  #ci.upper <- res+qnorm(.975)*sqrt(res*(1-res)/n)
  #ci.lower <- res-qnorm(.975)*sqrt(res*(1-res)/n)
  return(c(p0_hat = res, ci.upper = ci.upper, ci.lower = ci.lower))
}

repeat_function <- function(m,x){
  res <- matrix(0,ncol=3,nrow = m)
  i = 1
  repeat{
    if(i>m)
      break
    res[i,] <- run_model_mar(n=x,p0=0.2,rho1 = 0.2,pi = c(0.1,0.9),rho11=0.15)
    i = i+1
  }
  return(res)
}


#### run 500 times

m=500

set.seed(114)

outcome.100<-sapply(c(1000, 10000, 100000, 1000000), function(i) repeat_function(m=500,x=i))
outcome.100.df<-as.data.frame(outcome.100)


p0_hat.value <- outcome.100.df[1:m,]
p0_hat.upper <- outcome.100.df[(m+1):(2*m),]
p0_hat.lower <- outcome.100.df[(2*m+1):(3*m),]

p0_hat.ci <- cbind(p0_hat.value,p0_hat.upper,p0_hat.lower)
colnames(p0_hat.ci) <- c("v1","v2","v3","v4","v1.u","v2.u","v3.u","v4.u","v1.l","v2.l","v3.l","v4.l")


poN_10000<-mean(p0_hat.ci$v1)
poN_100000<-mean(p0_hat.ci$v2)
poN_1000000<-mean(p0_hat.ci$v3)
poN_10000000<-mean(p0_hat.ci$v4)

p0_hat.ci <- p0_hat.ci %>%
  mutate(v1.i = as.factor(ifelse(v1.u>=poN_10000&v1.l<=poN_10000,1,0)),
         v2.i = as.factor(ifelse(v2.u>=poN_100000&v2.l<=poN_100000,1,0)),
         v3.i = as.factor(ifelse(v3.u>=poN_1000000&v3.l<=poN_1000000,1,0)),
         v4.i = as.factor(ifelse(v4.u>=poN_10000000&v4.l<=poN_10000000,1,0)))

g1<-ggplot(p0_hat.ci,aes(x=1:nrow(p0_hat.ci),y=v1,color = v1.i))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v1.u, ymax=v1.l))+
  ggtitle("Population = 1000")+
  geom_hline(yintercept=poN_10000,  color = "red")+
  xlab("Iteration")+ylab("Prevalence")+ scale_color_manual(values=c("#000000", "#999999"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+ylim(0.13, .25)

g2<-ggplot(p0_hat.ci,aes(x=1:nrow(p0_hat.ci),y=v2,color = v2.i))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v2.u, ymax=v2.l))+
  ggtitle("Population = 10000")+
  geom_hline(yintercept=poN_100000,  color = "red")+
  xlab("Iteration")+ylab("Prevalence")+ scale_color_manual(values=c("#000000", "#999999"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+ylim(0.13, .25)

g3<-ggplot(p0_hat.ci,aes(x=1:nrow(p0_hat.ci),y=v3,color = v3.i))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v3.u, ymax=v3.l))+
  ggtitle("Population = 100000")+
  geom_hline(yintercept=poN_1000000,  color = "red")+
  xlab("Iteration")+ylab("Prevalence")+ scale_color_manual(values=c("#000000", "#999999"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+ylim(0.13, .25)

g4<-ggplot(p0_hat.ci,aes(x=1:nrow(p0_hat.ci),y=v4,color = v4.i))+
  geom_point()+
  geom_errorbar(mapping=aes( ymin=v4.u, ymax=v4.l))+
  ggtitle("Population = 1000000")+
  geom_hline(yintercept=poN_10000000,  color = "red")+
  xlab("Iteration")+ylab("Prevalence")+ scale_color_manual(values=c("#000000", "#999999"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+ylim(0.13, .25)

rho5.ci<-ggpubr::ggarrange(g1,g2,g3,g4,ncol=2,nrow=2)

rho5.ci
                    
                    
p0_hat.est <- outcome.100.df[1:m,]
p0_hat.sample <- outcome.100.df[(m+1):(2*m),]

colm <- colMeans(p0_hat.est)
colm2 <- colMeans(p0_hat.sample)
## Active Information 
log(colMeans(p0_hat.est)/0.2)
log(colMeans(p0_hat.sample)/0.2)
log(colMeans(p0_hat.est)/colMeans(p0_hat.sample))
## Active Information 
apply(p0_hat.sample,2,sd)

colMeans(log(p0_hat.sample/0.2))
colMeans(log(p0_hat.est/p0_hat.sample))




