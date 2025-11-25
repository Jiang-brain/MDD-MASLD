#################################masld==0|1:joint reference group of deprivation:3#################################
new_data=subset(datas,deprivation==0|deprivation==3)
new_data$joint_group=NA
index_0=which(new_data$masld==0&new_data$deprivation==0)
index_1=which(new_data$masld==1&new_data$deprivation==0)
index_2=which(new_data$masld==0&new_data$deprivation==3)
index_3=which(new_data$masld==1&new_data$deprivation==3)
new_data$joint_group[index_0]=0
new_data$joint_group[index_1]=1
new_data$joint_group[index_2]=2
new_data$joint_group[index_3]=3
new_data$joint_group=factor(new_data$joint_group)
rm(index_0,index_1,index_2,index_3)

fit_joint<- coxph(Surv(time, status) ~ joint_group+age+gender+race+university_education+activity+smoking+alcohol+metabolic_syndrome,
                  data = new_data)
result<-summary(fit_joint)
rm(fit_joint)
results<-data.frame(id=seq(1:4),HR=NA,P_value=NA,lower=NA, upper=NA,CI=NA,N=NA,events=NA)
row.names(results)=c('Exposure0_Deprivation0','Exposure1_Deprivation0','Exposure0_Deprivation1','Exposure1_Deprivation1')
for (i in c(2:4)) {
  results$HR[i]=round(result$conf.int[i-1,'exp(coef)'],2)
  results$lower[i]=round(result$conf.int[i-1,'lower .95'],2)
  results$upper[i]=round(result$conf.int[i-1,'upper .95'],2)
  results$P_value[i]=signif(result$coefficients[i-1,'Pr(>|z|)'],3)
  results$CI[i]=paste('(',results$lower[i],'-',results$upper[i],')',sep ='')
  results$events[i]=length(which(new_data$joint_group==i-1&new_data$status==1))
}
results$N=table(new_data$joint_group)
results$events[1]=length(which(new_data$joint_group==0&new_data$status==1))

results$HR[1]=1
results$lower[1]=1
results$upper[1]=1
rm(new_data,i,result)
#write.csv(results,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_masld_deprivation_additive_interactive.csv')


###km curve
fit <- survfit(Surv(time, status) ~ joint_group, data = new_data)
ggsurvplot(
  fit,
  data = new_data,
  fun = "cumhaz",#cumhaz
  conf.int = TRUE,
  legend.title = "Group",
  censor = FALSE,
  size=1.0,
  palette = "Spectral",#Set1
  lty = c(2, 2,2,2),
  legend.labs = c("M0_D0","M1_D0","M0_D3","M1_D3"),
  #ggtheme = theme_bw(),
  xlab = "Follow-up days",
  ylab = "Cumulative hazard (%)",
  risk.table = TRUE,
)
#################################masld==0|1:joint reference group of gender#################################
new_data=datas
new_data$joint_group=NA
index_0=which(new_data$masld==0&new_data$gender==1)
index_1=which(new_data$masld==1&new_data$gender==1)
index_2=which(new_data$masld==0&new_data$gender==0)
index_3=which(new_data$masld==1&new_data$gender==0)
new_data$joint_group[index_0]=0
new_data$joint_group[index_1]=1
new_data$joint_group[index_2]=2
new_data$joint_group[index_3]=3
new_data$joint_group=factor(new_data$joint_group)
rm(index_0,index_1,index_2,index_3)

fit_joint<- coxph(Surv(time, status) ~ joint_group+age+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,
                  data = new_data)
result<-summary(fit_joint)
rm(fit_joint)
results<-data.frame(id=seq(1:4),HR=NA,P_value=NA,lower=NA, upper=NA,CI=NA,N=NA,events=NA)
row.names(results)=c('Exposure0_gender1','Exposure1_gender1','Exposure0_gender0','Exposure1_gender0')
for (i in c(2:4)) {
  results$HR[i]=round(result$conf.int[i-1,'exp(coef)'],2)
  results$lower[i]=round(result$conf.int[i-1,'lower .95'],2)
  results$upper[i]=round(result$conf.int[i-1,'upper .95'],2)
  results$P_value[i]=signif(result$coefficients[i-1,'Pr(>|z|)'],3)
  results$CI[i]=paste('(',results$lower[i],'-',results$upper[i],')',sep ='')
  results$events[i]=length(which(new_data$joint_group==i-1&new_data$status==1))
}
results$N=table(new_data$joint_group)
results$events[1]=length(which(new_data$joint_group==0&new_data$status==1))

results$HR[1]=1
results$lower[1]=1
results$upper[1]=1
rm(new_data,i,result)
#write.csv(results,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_masld_gender_additive_interactive.csv')


###km curve
fit <- survfit(Surv(time, status) ~ joint_group, data = new_data)
ggsurvplot(
  fit,
  data = new_data,
  fun = "cumhaz",#cumhaz
  conf.int = TRUE,
  legend.title = "Group",
  censor = FALSE,
  size=1.0,
  palette = "Spectral",#Set1
  lty = c(2, 2,2,2),
  legend.labs = c("M0_g1","M1_g1","M0_g0","M1_g0"),
  #ggtheme = theme_bw(),
  xlab = "Follow-up days",
  ylab = "Cumulative hazard (%)",
  risk.table = TRUE,
)

###################test addictive interaction: bootstrap, deprivation###################
result<-data.frame(id=seq(1:5000),reri=NA,ap=NA,s=NA)
colnames(result)=c('boot','reri','ap','s')

new_data=subset(datas,deprivation==0|deprivation==3)
new_data$joint_group=NA
index_0=which(new_data$masld==0&new_data$deprivation==0)
index_1=which(new_data$masld==1&new_data$deprivation==0)
index_2=which(new_data$masld==0&new_data$deprivation==3)
index_3=which(new_data$masld==1&new_data$deprivation==3)
new_data$joint_group[index_0]=0
new_data$joint_group[index_1]=1
new_data$joint_group[index_2]=2
new_data$joint_group[index_3]=3
new_data$joint_group=factor(new_data$joint_group)
rm(index_0,index_1,index_2,index_3)

##botstrapping 5000 times
for (boot_strap in 1:5000) {
  random_index=sample(1:nrow(new_data), nrow(new_data), replace = TRUE)
  boot_datas=new_data[random_index,]
  rm(random_index)
  model<- coxph(Surv(time, status) ~ joint_group+age+gender+race+university_education+activity+smoking+alcohol+metabolic_syndrome,
                data = boot_datas)
  coeffs<-summary(model)$coefficients
  rm(model)
  
  result$reri[boot_strap]=coeffs['joint_group3','exp(coef)']-coeffs['joint_group2','exp(coef)']-coeffs['joint_group1','exp(coef)']+1
  result$ap[boot_strap]=result$reri[boot_strap]/coeffs['joint_group3','exp(coef)']
  result$s[boot_strap]=(coeffs['joint_group3','exp(coef)']-1)/(coeffs['joint_group2','exp(coef)']-1+coeffs['joint_group1','exp(coef)']-1)
  print(c(boot_strap,result$reri[boot_strap],result$ap[boot_strap]))
}
write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_masld_deprivation_depression.csv')


###################test addictive interaction: bootstrap, gender###################
result<-data.frame(id=seq(1:5000),reri=NA,ap=NA,s=NA)
colnames(result)=c('boot','reri','ap','s')

new_data=datas
new_data$joint_group=NA
index_0=which(new_data$masld==0&new_data$gender==1)
index_1=which(new_data$masld==1&new_data$gender==1)
index_2=which(new_data$masld==0&new_data$gender==0)
index_3=which(new_data$masld==1&new_data$gender==0)
new_data$joint_group[index_0]=0
new_data$joint_group[index_1]=1
new_data$joint_group[index_2]=2
new_data$joint_group[index_3]=3
new_data$joint_group=factor(new_data$joint_group)
rm(index_0,index_1,index_2,index_3)

##botstrapping 5000 times
for (boot_strap in 1:5000) {
  random_index=sample(1:nrow(new_data), nrow(new_data), replace = TRUE)
  boot_datas=new_data[random_index,]
  rm(random_index)
  model<- coxph(Surv(time, status) ~ joint_group+age+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,
                data = boot_datas)
  coeffs<-summary(model)$coefficients
  rm(model)
  
  result$reri[boot_strap]=coeffs['joint_group3','exp(coef)']-coeffs['joint_group2','exp(coef)']-coeffs['joint_group1','exp(coef)']+1
  result$ap[boot_strap]=result$reri[boot_strap]/coeffs['joint_group3','exp(coef)']
  result$s[boot_strap]=(coeffs['joint_group3','exp(coef)']-1)/(coeffs['joint_group2','exp(coef)']-1+coeffs['joint_group1','exp(coef)']-1)
  print(c(boot_strap,result$reri[boot_strap],result$ap[boot_strap]))
}
write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_masld_gender_depression.csv')

