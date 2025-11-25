fit0<- coxph(Surv(time, status) ~ masld+age+gender+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit1<- coxph(Surv(time, status) ~ masld*age+gender+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit2<- coxph(Surv(time, status) ~ masld*gender+age+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit3<- coxph(Surv(time, status) ~ masld*race+age+gender+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit4<- coxph(Surv(time, status) ~ masld*deprivation+age+gender+race+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit5<- coxph(Surv(time, status) ~ masld*university_education+age+gender+race+deprivation+activity+smoking+alcohol+metabolic_syndrome,data = datas)
fit6<- coxph(Surv(time, status) ~ masld*activity+age+gender+race+deprivation+university_education+smoking+alcohol+metabolic_syndrome,data = datas)
fit7<- coxph(Surv(time, status) ~ masld*smoking+age+gender+race+deprivation+university_education+activity+alcohol+metabolic_syndrome,data = datas)
fit8<- coxph(Surv(time, status) ~ masld*alcohol+age+gender+race+deprivation+university_education+activity+smoking+metabolic_syndrome,data = datas)
fit9<- coxph(Surv(time, status) ~ masld*metabolic_syndrome+age+gender+race+deprivation+university_education+activity+smoking+alcohol,data = datas)

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
#write.csv(P_interaction,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/interaction_masld_depression.csv')

############################masld==0|1: subgroup by age, two models###################################
datas0 <- subset(datas, age==0)
datas1 <- subset(datas, age==1)
results<-data.frame(id=seq(1:4),HR=1,P_value=NA,lower=1, upper=1,CI=NA,N=0,events=0)
row.names(results)=c('nonmasld_middle','masld_middle','nonmasld_old','masld_old')
new_datas=list(datas0,datas1)
rm(datas0,datas1)

for (i in c(1:2)) {
  temp_datas=new_datas[[i]]#middle, old
  
  fit_middle <- coxph(Surv(time, status) ~ masld+gender+race+deprivation+university_education+activity+
                        smoking+alcohol+metabolic_syndrome,data = temp_datas)
  results$HR[i*2]=round((summary(fit_middle)$conf.int[1,'exp(coef)']),2)
  lower=round(summary(fit_middle)$conf.int[1,'lower .95'],2)
  upper=round(summary(fit_middle)$conf.int[1,'upper .95'],2)
  results$lower[i*2]=lower
  results$upper[i*2]=upper
  results$CI[i*2]=paste('(',lower,'-',upper,')',sep ='')
  results$N[i*2]=length(which(temp_datas$age==i-1&temp_datas$masld==1))
  results$events[i*2]=length(which(temp_datas$age==i-1&temp_datas$status==1&temp_datas$masld==1))
  results$P_value[i*2]=signif(summary(fit_middle)$coefficients[1,'Pr(>|z|)'],3)
  rm(lower,upper)
  
  results$CI[i*2-1]='1 (Reference)'
  results$N[i*2-1]=length(which(temp_datas$age==i-1&temp_datas$masld==0))
  results$events[i*2-1]=length(which(temp_datas$age==i-1&temp_datas$status==1&temp_datas$masld==0))
  rm(fit_middle,temp_datas)
}


###########################masld==0|1:interaction R#####################################
library(interactionR)
result<-data.frame(id=seq(1:11),
                   reri_masld=NA,ap_masld=NA,s_masld=NA,P_reri=NA,P_ap=NA,P_s=NA)
colnames(result)=c('index','reri_masld','ap_masld','s_masld','P_reri_masld','P_ap_masld','P_s_masld')
row.names(result)=c('Age','Sex','Race','Deprivation: Q2','Deprivation: Q3','Deprivation: Q4','Education','Physical_activity',
                    'Smoking','Alcohol','Metabolic_syndrome')
iteration=1

#age
datas1=datas
model<- coxph(Surv(time, status) ~ masld*age+
                gender+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','age'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)

rm(model,outcome,datas1,out)
iteration=iteration+1

#gender
datas1=datas
index_male=which(datas1$gender==1)
index_female=which(datas1$gender==0)
datas1$gender[index_male]=0
datas1$gender[index_female]=1
rm(index_male,index_female)
datas1$gender=factor(datas1$gender)

model<- coxph(Surv(time, status) ~ masld*gender+
                age+race+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','gender'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe

result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)
rm(model,outcome,datas1,out)
iteration=iteration+1

#race
datas1=datas
model<- coxph(Surv(time, status) ~ masld*race+
                age+gender+deprivation+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','race'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)
rm(model,outcome,datas1,out)
iteration=iteration+1

#deprivation
for (deprivation_index in 1:3) {
  index=which(datas$deprivation==0|datas$deprivation==deprivation_index)
  new_datas=datas[index,]
  new_datas$deprivation=factor(new_datas$deprivation)
  rm(index)
  datas1=new_datas
  rm(new_datas)
  model<- coxph(Surv(time, status) ~ masld*deprivation+
                  age+gender+race+university_education+activity+smoking+alcohol+metabolic_syndrome,data = datas1)
  summary(model)
  out<-interactionR(model, exposure_names = c('masld','deprivation'),ci.type = "delta",ci.level = 0.95,em=F)
  outcome<-out$dframe
  result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
  result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
  result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
  result$P_reri_masld[iteration]=signif(outcome[10,5],3)
  result$P_ap_masld[iteration]=signif(outcome[11,5],3)
  result$P_s_masld[iteration]=signif(outcome[12,5],3)
  rm(model,outcome,datas1,out)
  iteration=iteration+1
}

#education
datas1=datas
model<- coxph(Surv(time, status) ~ masld*university_education+
                age+gender+race+deprivation+activity+smoking+alcohol+metabolic_syndrome,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','university_education'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)
rm(model,outcome,datas1,out)
iteration=iteration+1

#activity
datas1=datas
model<- coxph(Surv(time, status) ~ masld*activity+
                age+gender+race+deprivation+university_education+smoking+alcohol+metabolic_syndrome,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','activity'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)
rm(model,outcome,datas1,out)
iteration=iteration+1

#smoking
datas1=datas
model<- coxph(Surv(time, status) ~ masld*smoking+
                age+gender+race+deprivation+university_education+activity+alcohol+metabolic_syndrome,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','smoking'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe

result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)
rm(model,outcome,datas1,out)
iteration=iteration+1

#alcohol
new_datas=datas
index_0=which(datas$alcohol==0)
index_1=which(datas$alcohol==1)
new_datas$alcohol[index_0]=1
new_datas$alcohol[index_1]=0
rm(index_0,index_1)
new_datas$alcohol=factor(new_datas$alcohol)
datas1=new_datas
rm(new_datas)
model<- coxph(Surv(time, status) ~ masld*alcohol+
                age+gender+deprivation+race+university_education+activity+smoking+metabolic_syndrome,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','alcohol'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe

result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)
rm(model,outcome,datas1,out)
iteration=iteration+1

#metabolic_syndrome
datas1=datas
model<- coxph(Surv(time, status) ~ masld*metabolic_syndrome+
                age+gender+race+deprivation+university_education+activity+smoking+alcohol,data = datas1)
summary(model)
out<-interactionR(model, exposure_names = c('masld','metabolic_syndrome'),ci.type = "delta",ci.level = 0.95,em=F)
outcome<-out$dframe
result$reri_masld[iteration]=paste(round(outcome[10,2],2),' (',round(outcome[10,3],2),', ',round(outcome[10,4],2),')',sep='')
result$ap_masld[iteration]=paste(round(outcome[11,2]*100,2),'% (',round(outcome[11,3]*100,2),'%, ',round(outcome[11,4]*100,2),'%)',sep='')
result$s_masld[iteration]=paste(round(outcome[12,2],2),' (',round(outcome[12,3],2),', ',round(outcome[12,4],2),')',sep='')
result$P_reri_masld[iteration]=signif(outcome[10,5],3)
result$P_ap_masld[iteration]=signif(outcome[11,5],3)
result$P_s_masld[iteration]=signif(outcome[12,5],3)
rm(model,outcome,datas1,out)
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/additive_masld_depression.csv')





