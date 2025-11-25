##################################################Cardiometabolic components#####################################################
cardiometabolic_factor=data.frame(cardiometabolic_factor)
data=data.frame(status=data_censoring$label,time=data_censoring$censoring,depression=depression,
                metabolic_syndrome=metSyn,
                diabetes=cardiometabolic_factor$diabetes,hdl=cardiometabolic_factor$hdl_cholesterol,
                hypertension=cardiometabolic_factor$hypertension,tg=cardiometabolic_factor$triglyceride,
                age=age,gender=sex,sites=sites,race=ethnicity,obesity=cardiometabolic_factor$obesity,
                university_education=university_education,deprivation=deprivation,
                smoking=smoking,alcohol=alcohol_status,activity=activity)
rm(cardiometabolic_factor)
rm(age, sex,alcohol_status,deprivation,ethnicity,income,sedentary,depression,activity,
   sites,smoking,university_education,pdff,obesity,
   data_censoring,metSyn)


data$race=factor(data$race)
data$university_education=factor(data$university_education)
data$gender=factor(data$gender)
#data$income=factor(data$income)
data$alcohol=factor(data$alcohol)
data$smoking=factor(data$smoking)
data$deprivation=factor(data$deprivation)
#data$sedentary=factor(data$sedentary)
data$activity=factor(data$activity)
data$metabolic_syndrome=factor(data$metabolic_syndrome)
data$diabetes=factor(data$diabetes)
data$hdl=factor(data$hdl)
data$hypertension=factor(data$hypertension)
data$tg=factor(data$tg)
data$obesity=factor(data$obesity)

index_young_age=which(data$age<=40)
data$age[index_young_age]=NA
rm(index_young_age)

## retain complete data
number_na<-rowSums(is.na(data))
index_nonNA<-which(number_na==0)
datas<-data[index_nonNA,]
metabolites=metabolites[index_nonNA,]
rm(number_na,index_nonNA)
datas$sites=factor(datas$sites)


# ##age group
# index_0=which(datas$age<45)
# index_1=which(datas$age>=45&datas$age<55)
# index_2=which(datas$age>=55&datas$age<65)
# index_3=which(datas$age>=65)
# datas$age[index_0]=0
# datas$age[index_1]=0
# datas$age[index_2]=2
# datas$age[index_3]=3
# rm(index_0,index_1,index_2,index_3)
# datas$age=factor(datas$age)

##age group
index_middle=which(datas$age<65)
index_old=which(datas$age>=65)
datas$age[index_middle]=0
datas$age[index_old]=1
datas$age=factor(datas$age)
rm(index_middle,index_old)

##depression
datas$depression=factor(datas$depression)
datas$depression = relevel(datas$depression, ref = "0")

##landmark
index_5year_lanmark=which(!(datas$time<365*5&datas$status==1))
datas=datas[index_5year_lanmark,]
metabolites=metabolites[index_5year_lanmark,]
rm(index_5year_lanmark)

########################################association results####################################################
res.cox <- coxph(Surv(time, status) ~ depression+age+gender+race+deprivation+university_education+activity+
                   smoking+alcohol+
                   diabetes+hdl+hypertension+tg+obesity,data = datas)
summary(res.cox)
coeff<-summary(res.cox)
rm(res.cox)
result=data.frame(HR=round(coeff$coefficients[,2],2),lower=round(coeff$conf.int[,3],2),upper=round(coeff$conf.int[,4],2),p_value=signif(coeff$coefficients[,5],3))
for (i in 1:nrow(result)){
  result$CI[i]=paste(result$HR[i],' (',result$lower[i],'-',result$upper[i],')',sep='')
}
result<-rbind(c(1, 1, 1, NA, '1.00 (Ref)'),result[1,],c(1, 1, 1, NA, '1.00 (Ref)'),result[2,],c(1, 1, 1, NA, '1.00 (Ref)'),result[3,],c(1, 1, 1, NA, '1.00 (Ref)'),result[4,],c(1, 1, 1, NA, '1.00 (Ref)'),
              result[5:7,],c(1, 1, 1, NA, '1.00 (Ref)'),result[8,],c(1, 1, 1, NA, '1.00 (Ref)'),result[9,],c(1, 1, 1, NA, '1.00 (Ref)'),result[10,],c(1, 1, 1, NA, '1.00 (Ref)'),result[11,],
              c(1, 1, 1, NA, '1.00 (Ref)'),result[12,],c(1, 1, 1, NA, '1.00 (Ref)'),result[13,],c(1, 1, 1, NA, '1.00 (Ref)'),result[14,],c(1, 1, 1, NA, '1.00 (Ref)'),result[15,],c(1, 1, 1, NA, '1.00 (Ref)'),result[16,])
result$N=NA
result$cases=NA
#depression
for (i in 1:2){
  result$N[i]=length(which(datas$depression==i-1))
  result$cases[i]=length(which(datas$depression==(i-1)&datas$status==1))
}
{ #age
  for (i in 1:2){
    result$N[i+2]=length(which(datas$age==i-1))
    result$cases[i+2]=length(which(datas$age==(i-1)&datas$status==1))
  }
  #gender
  for (i in 1:2){
    result$N[i+4]=length(which(datas$gender==i-1))
    result$cases[i+4]=length(which(datas$gender==(i-1)&datas$status==1))
  }
  #race
  for (i in 1:2){
    result$N[i+6]=length(which(datas$race==i-1))
    result$cases[i+6]=length(which(datas$race==(i-1)&datas$status==1))
  }
  #deprivation
  for (i in 1:4){
    result$N[i+8]=length(which(datas$deprivation==i-1))
    result$cases[i+8]=length(which(datas$deprivation==(i-1)&datas$status==1))
  }
  #education
  for (i in 1:2){
    result$N[i+12]=length(which(datas$university_education==i-1))
    result$cases[i+12]=length(which(datas$university_education==(i-1)&datas$status==1))
  }
  #activity
  for (i in 1:2){
    result$N[i+14]=length(which(datas$activity==i-1))
    result$cases[i+14]=length(which(datas$activity==(i-1)&datas$status==1))
  }
  #smoking
  for (i in 1:2){
    result$N[i+16]=length(which(datas$smoking==i-1))
    result$cases[i+16]=length(which(datas$smoking==(i-1)&datas$status==1))
  }
  #alcohol
  for (i in 1:2){
    result$N[i+18]=length(which(datas$alcohol==i-1))
    result$cases[i+18]=length(which(datas$alcohol==(i-1)&datas$status==1))
  }
  #diabetes
  for (i in 1:2){
    result$N[i+20]=length(which(datas$diabetes==i-1))
    result$cases[i+20]=length(which(datas$diabetes==(i-1)&datas$status==1))
  }
  #hdl
  for (i in 1:2){
    result$N[i+22]=length(which(datas$hdl==i-1))
    result$cases[i+22]=length(which(datas$hdl==(i-1)&datas$status==1))
  }
  #hypertension
  for (i in 1:2){
    result$N[i+24]=length(which(datas$hypertension==i-1))
    result$cases[i+24]=length(which(datas$hypertension==(i-1)&datas$status==1))
  }
  #TG
  for (i in 1:2){
    result$N[i+26]=length(which(datas$tg==i-1))
    result$cases[i+26]=length(which(datas$tg==(i-1)&datas$status==1))
  }
  #obesity
  for (i in 1:2){
    result$N[i+28]=length(which(datas$obesity==i-1))
    result$cases[i+28]=length(which(datas$obesity==(i-1)&datas$status==1))
  }
  rm(coeff,i)
}
result<-rbind(c(NA,NA,NA,NA,NA,NA,NA),result[1:2,],c(NA,NA,NA,NA,NA,NA,NA),result[3:4,],c(NA,NA,NA,NA,NA,NA,NA),result[5:6,],c(NA,NA,NA,NA,NA,NA,NA),result[7:8,],c(NA,NA,NA,NA,NA,NA,NA),
              result[9:12,],c(NA,NA,NA,NA,NA,NA,NA),result[13:14,],c(NA,NA,NA,NA,NA,NA,NA),result[15:16,],c(NA,NA,NA,NA,NA,NA,NA),result[17:18,],c(NA,NA,NA,NA,NA,NA,NA),result[19:20,],
              c(NA,NA,NA,NA,NA,NA,NA),result[21:22,],c(NA,NA,NA,NA,NA,NA,NA),result[23:24,],c(NA,NA,NA,NA,NA,NA,NA),result[25:26,],c(NA,NA,NA,NA,NA,NA,NA),result[27:28,],c(NA,NA,NA,NA,NA,NA,NA),result[29:30,])
row.names(result)=c('Depression status','No depression','Depression','Age group','Age<65 years','>=65 yeas',
                    'Sex','Females','Males','Race','Ethnic minorities','White',
                    'Deprivation','Quarter 1','Quarter 2','Quarter 3','Quarter 4',
                    'Education level','Above Colleage','Less than college',
                    'Activity','Regularly active','Inactive','Smoking','Never','Ever','Alcohol','Never_','Ever_',
                    'Diabetes','No1','Yes1','Low HDL','No2','Yes2','Hypertension','No3','Yes3','High triglycerides','No4','Yes4','Central obesity','No5','Yes5')
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_depression_masld_metabolic_components.csv')


##########################Additive interaction R#######################
library(interactionR)
result<-data.frame(id=seq(1:5),reri=NA,ap=NA,s=NA,P_reri=NA,P_ap=NA,P_s=NA)
colnames(result)=c('index','reri','ap','s','P_reri','P_ap','P_s')
row.names(result)=c('T2D','hdl','Hypertension','tg','Obesity')
iteration=1

#diabetes
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*diabetes+age+gender+race+deprivation+university_education+activity+
                smoking+alcohol+
                hdl+hypertension+tg+obesity,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','diabetes'),ci.type = "delta",ci.level = 0.95,em=F)
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

#hdl
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*hdl+age+gender+race+deprivation+university_education+activity+
                smoking+alcohol+
                diabetes+hypertension+tg+obesity,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','hdl'),ci.type = "delta",ci.level = 0.95,em=F)
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

#hypertension
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*hypertension+age+gender+race+deprivation+university_education+activity+
                smoking+alcohol+
                diabetes+hdl+tg+obesity,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','hypertension'),ci.type = "delta",ci.level = 0.95,em=F)
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

#tg
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*tg+age+gender+race+deprivation+university_education+activity+
                smoking+alcohol+
                diabetes+hdl+hypertension+obesity,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','tg'),ci.type = "delta",ci.level = 0.95,em=F)
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

#obesity
new_datas=datas
model<- coxph(Surv(time, status) ~ depression*obesity+age+gender+race+deprivation+university_education+activity+
                smoking+alcohol+
                diabetes+hdl+hypertension+tg,data = new_datas)
summary(model)
out<-interactionR(model, exposure_names = c('depression','obesity'),ci.type = "delta",ci.level = 0.95,em=F)
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
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/Additive_interaction_metabolic_components.csv')

#################################joint reference group of obesity#################################
new_data=datas
new_data$joint_group=NA
index_0=which(new_data$depression==0&new_data$obesity==0)
index_1=which(new_data$depression==1&new_data$obesity==0)
index_2=which(new_data$depression==0&new_data$obesity==1)
index_3=which(new_data$depression==1&new_data$obesity==1)
new_data$joint_group[index_0]=0
new_data$joint_group[index_1]=1
new_data$joint_group[index_2]=2
new_data$joint_group[index_3]=3
new_data$joint_group=factor(new_data$joint_group)
rm(index_0,index_1,index_2,index_3)

fit_joint<- coxph(Surv(time, status) ~ joint_group+age+gender+race+deprivation+university_education+activity+smoking+alcohol+
                    diabetes+hdl+hypertension+tg,
                  data = new_data)
result<-summary(fit_joint)
rm(fit_joint)
results<-data.frame(id=seq(1:4),HR=NA,P_value=NA,lower=NA, upper=NA,CI=NA,N=NA,events=NA)
row.names(results)=c('Exposure0_Obesity0','Exposure1_Obesity0','Exposure0_Obesity1','Exposure1_Obesity1')
for (i in c(2:4)) {
  results$HR[i]=round(result$conf.int[i-1,'exp(coef)'],2)
  results$lower[i]=round(result$conf.int[i-1,'lower .95'],2)
  results$upper[i]=round(result$conf.int[i-1,'upper .95'],2)
  results$P_value[i]=signif(result$coefficients[i-1,'Pr(>|z|)'],3)
  results$CI[i]=paste(results$HR[i],' (',results$lower[i],'-',results$upper[i],')',sep ='')
  results$events[i]=length(which(new_data$joint_group==i-1&new_data$status==1))
}
results$N=table(new_data$joint_group)
results$events[1]=length(which(new_data$joint_group==0&new_data$status==1))

results$HR[1]=1
results$lower[1]=1
results$upper[1]=1
rm(new_data,i,result)
#write.csv(results,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_depression_obesity_additive_interactive.csv')

##km curve
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
  lty = c(2,2,2,2),
  legend.labs = c("D0_Obe0","D1_Obe0","D0_Obe1","D1_Obe1"),
  #ggtheme = theme_bw(),
  xlab = "Follow-up days",
  ylab = "Cumulative hazard (%)",
  risk.table = TRUE,
)



###################test addictive interaction: bootstrap, obesity###################
result<-data.frame(id=seq(1:5000),reri=NA,ap=NA,s=NA)
colnames(result)=c('boot','reri','ap','s')

new_data=datas
new_data$joint_group=NA
index_0=which(new_data$depression==0&new_data$obesity==0)
index_1=which(new_data$depression==1&new_data$obesity==0)
index_2=which(new_data$depression==0&new_data$obesity==1)
index_3=which(new_data$depression==1&new_data$obesity==1)
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
  model<- coxph(Surv(time, status) ~ joint_group+age+gender+race+deprivation+university_education+activity+smoking+alcohol+
                  diabetes+hdl+hypertension+tg,data = boot_datas)
  coeffs<-summary(model)$coefficients
  rm(model)
  
  result$reri[boot_strap]=coeffs['joint_group3','exp(coef)']-coeffs['joint_group2','exp(coef)']-coeffs['joint_group1','exp(coef)']+1
  result$ap[boot_strap]=result$reri[boot_strap]/coeffs['joint_group3','exp(coef)']
  result$s[boot_strap]=(coeffs['joint_group3','exp(coef)']-1)/(coeffs['joint_group2','exp(coef)']-1+coeffs['joint_group1','exp(coef)']-1)
  print(c(boot_strap,result$reri[boot_strap],result$ap[boot_strap]))
}
write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_depression_obesity_masld.csv')
