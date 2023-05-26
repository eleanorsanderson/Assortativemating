#Neil Davies 08/11/22
#Eleanor Sanderson 4/2/23

##simulations to show the bias induced/not by assortative mating 

##set up enviroment 
rm(list = ls())
set.seed(1245)

#install.packages("AER")
#install.packages("multiwayvcov")
#install.packages("plm")
#install.packages("ivpack")
#install.packages("data.table") 

library("AER")
library('data.table')
#library('plm')
#library('ivpack')
library('dplyr')
library('multiwayvcov')
library('lmtest')
library('sandwich')
#library('ivreg')

#parameters that could change
reps = 500
n=40000



#dataframe for output
results_ols_all<-data.frame(matrix(NA,nrow = 0,ncol = 5))
colnames(results_ols_all)<-c("am","m_beta","m_se","f_beta","f_se")
results_mr_all<-data.frame(matrix(NA,nrow = 0,ncol = 5))
colnames(results_mr_all)<-c("am","m_beta","m_se","f_beta","f_se")
results_mvmr_all<-data.frame(matrix(NA,nrow = 0,ncol = 5))
colnames(results_mvmr_all)<-c("am","m_beta","m_se","f_beta","f_se")
results_wfmr_all<-data.frame(matrix(NA,nrow = 0,ncol = 5))
colnames(results_wfmr_all)<-c("am","m_beta","m_se","f_beta","f_se")


for(j in 1:reps){

#set up results dataframes for this repetition
results_ols<-data.frame(matrix(NA,nrow = 10,ncol = 5))
colnames(results_ols)<-c("am","m_beta","m_se","f_beta","f_se")
results_mr<-data.frame(matrix(NA,nrow = 10,ncol = 5))
colnames(results_mr)<-c("am","m_beta","m_se","f_beta","f_se")
results_mvmr<-data.frame(matrix(NA,nrow = 10,ncol = 5))
colnames(results_mvmr)<-c("am","m_beta","m_se","f_beta","f_se")
results_wfmr<-data.frame(matrix(NA,nrow = 10,ncol = 5))
colnames(results_wfmr)<-c("am","m_beta","m_se","f_beta","f_se")

  
  #Generate genotype for fathers and mothers
  data_m<-data.frame(matrix(NA,nrow = n,ncol = 1))
data_f<-data.frame(matrix(NA,nrow = n,ncol = 1))

#Create ID variables
data_m$id<-seq.int(nrow(data_m)) 
data_f$id<-seq.int(nrow(data_f)) 

#Drop the first column (suggestions welcome for how to get rid of this line)
data_m<-subset(data_m,select="id")
data_f<-subset(data_f,select="id")

#Draw parental genotypes from the binomial distribution 
data_m$m1<-rbinom(n, 1, 0.6)
data_m$m2<-rbinom(n, 1, 0.6)
data_m$ma<-data_m$m1+data_m$m2
data_f$f1<-rbinom(n, 1, 0.6)
data_f$f2<-rbinom(n, 1, 0.6)
data_f$fa<-data_f$f1+data_f$f2

#Generate the variable to sort on based on the allele and a random error term
#am is the tuning parameter that determines how much of the assortment is chance and how much due to genotype.
#simulations loop over am 0 - 9 where 0 is the highest level of assortative mating and 9 is almost random mating. 

for (am_count in 0:9) {
  
  am_count_i=am_count+1
  am<-am_count*0.1
  
  data_m$am_index<-rnorm(n, 0, 1)*am+data_m$ma*(1-am)
  data_f$am_index<-rnorm(n, 0, 1)*am+data_f$fa*(1-am)
  
  #Sort 
  data_m<-data_m[order(data_m$am_index),] 
  data_f<-data_f[order(data_f$am_index),] 
  
  #Recreate ID variables
  data_m$id<-seq.int(nrow(data_m)) 
  data_f$id<-seq.int(nrow(data_f)) 
  
  
  #Merge
  data<-merge(data_m,data_f,by="id")
  
  summary(lm(formula = fa~ma, data = data))
  
  #Draw inherited allele for the child:
  data$mia_sib1 <- rbinom(n, 1, 0.5)
  data$fia_sib1 <- rbinom(n, 1, 0.5)
  

  #Define sib genotype, offspring defined as o1
 
  data$o11 <-data$m1*data$mia_sib1+data$m2*(1-data$mia_sib1)
  data$o12 <-data$f1*data$fia_sib1+data$f2*(1-data$fia_sib1)
  data$o1a<-data$o11+data$o12
  
  
  #Thus, we've got parents, and one child
  #Generate error terms for the exposure and the outcomes
  data$o1e_error<-rnorm(n,0,1)
  data$o1o_error<-rnorm(n,0,1)
  data$me_error<-rnorm(n,0,1)
  data$fe_error<-rnorm(n,0,1)
  
  #Generate parental BMI
  data$m_bmi<-data$ma+data$me_error
  data$f_bmi<-data$fa+data$fe_error
  
  #Generate exposure (BMI)
  data$o1_bmi<-data$o1a+data$o1e_error
  
  #Generate outcome
  #This is a function, offspring BMI, and parents' BMI
  data$o1_educ<-data$m_bmi+data$f_bmi+data$o1_bmi+data$o1o_error
 
  ####ESTIMATION#####
 
  #1. phenotypic education and parents bmi 
  est_ols<-lm(formula = o1_educ~f_bmi+m_bmi, data = data)
  beta<-as.data.frame(est_ols$coefficients)
  se<-as.data.frame(sqrt(diag(vcov(est_ols))))
  results_ols[am_count_i,1]<-am_count
  results_ols[am_count_i,2]<-beta["m_bmi",1]
  results_ols[am_count_i,3]<-se["m_bmi",1]
  results_ols[am_count_i,4]<-beta["f_bmi",1]
  results_ols[am_count_i,5]<-se["f_bmi",1]
  
  
  
  #MR of each parent in turn
  est_ma <- summary(ivreg(o1_educ~m_bmi|ma, data = data))
  beta<-as.data.frame(est_ma$coefficients)
  results_mr[am_count_i,1]<-am_count
  results_mr[am_count_i,2]<-beta["m_bmi","Estimate"]
  results_mr[am_count_i,3]<-beta["m_bmi","Std. Error"]
  
  est_fa <- summary(ivreg(o1_educ~f_bmi|fa, data = data))
  beta<-as.data.frame(est_fa$coefficients)
  results_mr[am_count_i,4]<-beta["f_bmi","Estimate"]
  results_mr[am_count_i,5]<-beta["f_bmi","Std. Error"]
  
  #Family MVMR  - including both parents but not the offspring. 
  est_mvmr <- summary(ivreg(o1_educ~m_bmi+f_bmi|fa+ma, data = data))
  beta<-as.data.frame(est_mvmr$coefficients)
  results_mvmr[am_count_i,1]<-am_count
  results_mvmr[am_count_i,2]<-beta["m_bmi","Estimate"]
  results_mvmr[am_count_i,3]<-beta["m_bmi","Std. Error"]
  results_mvmr[am_count_i,4]<-beta["f_bmi","Estimate"]
  results_mvmr[am_count_i,5]<-beta["f_bmi","Std. Error"]
  
  #Within family MR
  est_wfmr <- summary(ivreg(o1_educ~m_bmi+f_bmi+o1a|fa+ma+o1a, data = data))
  beta<-as.data.frame(est_wfmr$coefficients)
  results_wfmr[am_count_i,1]<-am_count
  results_wfmr[am_count_i,2]<-beta["m_bmi","Estimate"]
  results_wfmr[am_count_i,3]<-beta["m_bmi","Std. Error"]
  results_wfmr[am_count_i,4]<-beta["f_bmi","Estimate"]
  results_wfmr[am_count_i,5]<-beta["f_bmi","Std. Error"]
}  

##combine the results for the repetition with the whole results table
results_mr_all <- bind_rows(results_mr_all, results_mr)
results_mvmr_all <- bind_rows(results_mvmr_all, results_mvmr)
results_ols_all <- bind_rows(results_ols_all, results_ols)
results_wfmr_all <- bind_rows(results_wfmr_all, results_wfmr)

}

sim_out_ols <- results_ols_all %>% 
  group_by(am) %>% 
  summarise(across(m_beta:f_se, mean))

sim_out_mr <- results_mr_all %>% 
  group_by(am) %>% 
  summarise(across(m_beta:f_se, mean))

sim_out_mvmr <- results_mvmr_all %>% 
  group_by(am) %>% 
  summarise(across(m_beta:f_se, mean))

sim_out_wfmr <- results_wfmr_all %>% 
  group_by(am) %>% 
  summarise(across(m_beta:f_se, mean))

#Save results for plotting in Stata
library(foreign)
write.dta(sim_out_ols, "results_ols.dta")
write.dta(sim_out_mr, "results_mr.dta")
write.dta(sim_out_mvmr, "results_mvmr.dta")
write.dta(sim_out_wfmr, "results_wfmr.dta")