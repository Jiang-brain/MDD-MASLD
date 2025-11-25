############################################interaction########################################
fit0<- coxph(Surv(time, status) ~ depression+age+gender+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit1<- coxph(Surv(time, status) ~ depression*age+gender+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit2<- coxph(Surv(time, status) ~ depression*gender+age+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit3<- coxph(Surv(time, status) ~ depression*race+age+gender+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit4<- coxph(Surv(time, status) ~ depression*deprivation+age+gender+race+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit5<- coxph(Surv(time, status) ~ depression*university_education+age+gender+race+deprivation+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit6<- coxph(Surv(time, status) ~ depression*activity+age+gender+race+deprivation+university_education+smoking+alcohol+metabolic_syndrome,data = datas)
fit7<- coxph(Surv(time, status) ~ depression*smoking+age+gender+race+deprivation+university_education+activity+alcohol+metabolic_syndrome,data = datas)
fit8<- coxph(Surv(time, status) ~ depression*alcohol+age+gender+race+deprivation+university_education+activity+smoking+metabolic_syndrome,data = datas)
fit9<- coxph(Surv(time, status) ~ depression*metabolic_syndrome+age+gender+race+deprivation+university_education+activity+smoking+alcohol,data = datas)

fit_models=list(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9)
rm(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9)
P_interaction=data.frame(matrix(ncol = 1, nrow = 9))
row.names(P_interaction)=c('Age','Sex','Race','Deprivation','Education','Physical_activity','Smoking','Alcohol','Metabolic_syndrome')
colnames(P_interaction)=c('P_interaction')
for (i in c(1:9)) {
  P_interaction$P_interaction[i]=signif(anova(fit0,fit_models[[i]],test="Chisq")$`Pr(>|Chi|)`[2],3)
}
rm(i,fit0,fit_models)
P_interaction$interaction=-log(P_interaction$P_interaction)
#write.csv(P_interaction,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/interaction_depression_masld.csv')

##########################Additive interaction R#####################################
library(interactionR)
result<-data.frame(id=seq(1:11),reri=NA,ap=NA,s=NA,P_reri=NA,P_ap=NA,P_s=NA)
colnames(result)=c('index','reri','ap','s','P_reri','P_ap','P_s')
row.names(result)=c('Age','Sex','Race','Deprivation: Q2','Deprivation: Q3','Deprivation: Q4','Education','Physical_activity',
                    'Smoking','Alcohol','Metabolic_syndrome')
iteration=1

#age
# index_middle=which(datas$age==0)
# index_old=which(datas$age==1)
# new_datas=datas
# new_datas$age[index_middle]=1
# new_datas$age[index_old]=0
# rm(index_middle,index_old)
model<- coxph(Surv(time, status) ~ depression*age+
                gender+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','age'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model)
iteration=iteration+1

#gender
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*gender+
                age+race+deprivation+university_education+smoking+activity+alcohol+metabolic_syndrome,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','gender'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model,new_datas)
iteration=iteration+1

#race
index_white=which(datas$race==1)
index_nonwhite=which(datas$race==0)
new_datas=datas
new_datas$race[index_white]=0
new_datas$race[index_nonwhite]=1
rm(index_white,index_nonwhite)
model<- coxph(Surv(time, status) ~ depression*race+
                age+gender+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','race'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model,new_datas)
iteration=iteration+1

#deprivation
for (deprivation_index in 1:3) {
  new_datas=subset(datas,deprivation==0|deprivation==deprivation_index)
  new_datas$deprivation=factor(new_datas$deprivation)
  model<- coxph(Surv(time, status) ~ depression*deprivation+
                  age+gender+race+university_education+activity+smoking+alcohol+metabolic_syndrome,data = new_datas)
  summary(model)
  out<-interactionR(model, exposure_names = c('depression','deprivation'),ci.type = "delta",ci.level = 0.95,em=F)
  outcome<-out$dframe
  rm(out)
  result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
  result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
  result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
  result$P_reri[iteration]=signif(outcome[10,5],3)
  result$P_ap[iteration]=signif(outcome[11,5],3)
  result$P_s[iteration]=signif(outcome[12,5],3)
  rm(outcome,model,new_datas)
  iteration=iteration+1
}
rm(deprivation_index)

#university
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*university_education+
                age+gender+race+deprivation+activity+alcohol+smoking+metabolic_syndrome,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','university_education'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model,new_datas)
iteration=iteration+1

#activity
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*activity+
                age+gender+race+deprivation+university_education+smoking+alcohol+metabolic_syndrome,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','activity'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model,new_datas)
iteration=iteration+1

#smoking
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*smoking+
                age+gender+race+deprivation+university_education+activity+alcohol+metabolic_syndrome,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','smoking'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model,new_datas)
iteration=iteration+1

#alcohol
index_alcohol_never=which(datas$alcohol==1)
index_alcohol_ever=which(datas$alcohol==0)
new_datas=datas
new_datas$alcohol[index_alcohol_never]=0
new_datas$alcohol[index_alcohol_ever]=1
new_datas$alcohol=factor(new_datas$alcohol)
rm(index_alcohol_never,index_alcohol_ever)
model<- coxph(Surv(time, status) ~ depression*alcohol+
                age+gender+race+deprivation+university_education+activity+smoking+metabolic_syndrome,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','alcohol'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model,new_datas)
iteration=iteration+1

#metabolic_syndrome
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*metabolic_syndrome+
                age+gender+race+deprivation+university_education+activity+alcohol+smoking,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','metabolic_syndrome'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
rm(out)
result$reri[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri[iteration]=signif(outcome[10,5],3)
result$P_ap[iteration]=signif(outcome[11,5],3)
result$P_s[iteration]=signif(outcome[12,5],3)
rm(outcome,model,new_datas)
rm(iteration)
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/additive_depression_masld.csv')

#################################joint reference group of metabolic syndrome#################################
new_data=datas
new_data$joint_group=NA
index_0=which(new_data$depression==0&new_data$metabolic_syndrome==0)
index_1=which(new_data$depression==1&new_data$metabolic_syndrome==0)
index_2=which(new_data$depression==0&new_data$metabolic_syndrome==1)
index_3=which(new_data$depression==1&new_data$metabolic_syndrome==1)
new_data$joint_group[index_0]=0
new_data$joint_group[index_1]=1
new_data$joint_group[index_2]=2
new_data$joint_group[index_3]=3
new_data$joint_group=factor(new_data$joint_group)
rm(index_0,index_1,index_2,index_3)

fit_joint<- coxph(Surv(time, status) ~ joint_group+age+gender+race+deprivation+university_education+activity+smoking+alcohol,
                  data = new_data)
result<-summary(fit_joint)
rm(fit_joint)
results<-data.frame(id=seq(1:4),HR=NA,P_value=NA,lower=NA, upper=NA,CI=NA,N=NA,events=NA)
row.names(results)=c('Exposure0_Met0','Exposure1_Met0','Exposure0_Met1','Exposure1_Met1')
for (i in c(2:4)) {
  results$HR[i]=round(result$conf.int[i-1,'exp(coef)'],2)
  results$lower[i]=round(result$conf.int[i-1,'lower .95'],2)
  results$upper[i]=round(result$conf.int[i-1,'upper .95'],2)
  results$P_value[i]=signif(result$coefficients[i-1,'Pr(>|z|)'],3)
  results$CI[i]=paste(results$HR[i],'(',results$lower[i],'-',results$upper[i],')',sep ='')
  results$events[i]=length(which(new_data$joint_group==i-1&new_data$status==1))
}
results$N=table(new_data$joint_group)
results$events[1]=length(which(new_data$joint_group==0&new_data$status==1))

results$HR[1]=1
results$lower[1]=1
results$upper[1]=1
rm(new_data,i,result)
#write.csv(results,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_depression_metabolic_additive_interaction.csv')

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
  legend.labs = c("D0_M0","D1_M0","D0_M1","D1_M1"),
  #ggtheme = theme_bw(),
  xlab = "Follow-up days",
  ylab = "Cumulative hazard (%)",
  risk.table = TRUE,
)

###################test addictive interaction: bootstrap, metabolic syndrome###################
result<-data.frame(id=seq(1:5000),reri=NA,ap=NA,s=NA)
colnames(result)=c('boot','reri','ap','s')

new_data=datas
new_data$joint_group=NA
index_0=which(new_data$depression==0&new_data$metabolic_syndrome==0)
index_1=which(new_data$depression==1&new_data$metabolic_syndrome==0)
index_2=which(new_data$depression==0&new_data$metabolic_syndrome==1)
index_3=which(new_data$depression==1&new_data$metabolic_syndrome==1)
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
  model<- coxph(Surv(time, status) ~ joint_group+age+gender+race+deprivation+university_education+smoking+alcohol+activity,data = boot_datas)
  coeffs<-summary(model)$coefficients
  rm(model)
  
  result$reri[boot_strap]=coeffs['joint_group3','exp(coef)']-coeffs['joint_group2','exp(coef)']-coeffs['joint_group1','exp(coef)']+1
  result$ap[boot_strap]=result$reri[boot_strap]/coeffs['joint_group3','exp(coef)']
  result$s[boot_strap]=(coeffs['joint_group3','exp(coef)']-1)/(coeffs['joint_group2','exp(coef)']-1+coeffs['joint_group1','exp(coef)']-1)
  print(c(boot_strap,result$reri[boot_strap],result$ap[boot_strap]))
}
write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_depression_MetSyn_masld.csv')
