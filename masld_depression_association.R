#####################################combine data#########################################################################
data=data.frame(status=data_censoring$label,time=data_censoring$censoring,fli=fli_score,masld=masld,
                metabolic_syndrome=metSyn,
                age=age,gender=sex,sites=sites,race=ethnicity,
                university_education=university_education,deprivation=deprivation,
                smoking=smoking,alcohol=alcohol_status,activity=activity)

# data=data.frame(status=data_censoring$label,time=data_censoring$censoring,fli=fli_score,masld=masld,
#                 metabolic_syndrome=metSyn,phone_head_side=phone_head_side,
#                 age=age,gender=sex,sites=sites,race=ethnicity,
#                 university_education=university_education,deprivation=deprivation,
#                 smoking=smoking,alcohol=alcohol_status,activity=activity)

rm(age, sex,alcohol_status,deprivation,ethnicity,income,sedentary,masld,activity,fli_score,
  sites,smoking,university_education,alcohol_amount,pdff,phone_head_side,
   data_censoring,metSyn,obesity,cancer,employement_status,sleep,illness)


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


## retain complete data
index_young_age=which(data$age<=40)
data$age[index_young_age]=NA
rm(index_young_age)

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

##masld
index_0=which(datas$masld==0|datas$masld==1)
index_1=which(datas$masld==2)
datas$masld[index_0]=0
datas$masld[index_1]=1
rm(index_0,index_1)
datas$masld=factor(datas$masld)
datas$masld = relevel(datas$masld, ref = "0")

##landmark
index_5year_lanmark=which(!(datas$time<365*5&datas$status==1))
datas=datas[index_5year_lanmark,]
metabolites=metabolites[index_5year_lanmark,]
rm(index_5year_lanmark)


##############################################################################################################
library("survival")
library("survminer")

res.cox <- coxph(Surv(time, status) ~ masld+age+gender+race+deprivation+university_education+activity+
                   smoking+alcohol+metabolic_syndrome,data = datas)
summary(res.cox)
test.zph <- cox.zph(res.cox,'rank')
test.zph


res.cox <- coxph(Surv(time, status) ~ scale(fli)+age+gender+race+deprivation+university_education+activity+
                   smoking+alcohol+metabolic_syndrome,data = datas)
summary(res.cox)
test.zph <- cox.zph(res.cox,'rank')
test.zph

###################################negative control exposure#################################################
library("survival")
library("survminer")
datas$phone_head_side=factor(datas$phone_head_side)
res.cox <- coxph(Surv(time, status) ~ phone_head_side+age+gender+race+deprivation+university_education+activity+
                   smoking+alcohol+metabolic_syndrome,data = datas)
summary(res.cox)


########################################association results####################################################
res.cox <- coxph(Surv(time, status) ~ masld+age+gender+race+deprivation+university_education+activity+
                   smoking+alcohol+metabolic_syndrome,data = datas)
summary(res.cox)
coeff<-summary(res.cox)
rm(res.cox)
result=data.frame(HR=round(coeff$coefficients[,2],2),lower=round(coeff$conf.int[,3],2),upper=round(coeff$conf.int[,4],2),p_value=signif(coeff$coefficients[,5],3))
for (i in 1:nrow(result)){
  result$CI[i]=paste(result$HR[i],' (',result$lower[i],'-',result$upper[i],')',sep='')
}
result<-rbind(c(1, 1, 1, NA, '1.00 (Ref)'),result[1,],c(1, 1, 1, NA, '1.00 (Ref)'),result[2,],c(1, 1, 1, NA, '1.00 (Ref)'),result[3,],c(1, 1, 1, NA, '1.00 (Ref)'),result[4,],c(1, 1, 1, NA, '1.00 (Ref)'),
              result[5:7,],c(1, 1, 1, NA, '1.00 (Ref)'),result[8,],c(1, 1, 1, NA, '1.00 (Ref)'),result[9,],c(1, 1, 1, NA, '1.00 (Ref)'),result[10,],c(1, 1, 1, NA, '1.00 (Ref)'),result[11,],
              c(1, 1, 1, NA, '1.00 (Ref)'),result[12,])
result$N=NA
result$cases=NA
#masld
for (i in 1:2){
  result$N[i]=length(which(datas$masld==i-1))
  result$cases[i]=length(which(datas$masld==(i-1)&datas$status==1))
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
  #metabolic_syndrome
  for (i in 1:2){
    result$N[i+20]=length(which(datas$metabolic_syndrome==i-1))
    result$cases[i+20]=length(which(datas$metabolic_syndrome==(i-1)&datas$status==1))
  }
  rm(coeff,i)
}
result<-rbind(c(NA,NA,NA,NA,NA,NA,NA),result[1:2,],c(NA,NA,NA,NA,NA,NA,NA),result[3:4,],c(NA,NA,NA,NA,NA,NA,NA),result[5:6,],c(NA,NA,NA,NA,NA,NA,NA),result[7:8,],c(NA,NA,NA,NA,NA,NA,NA),
              result[9:12,],c(NA,NA,NA,NA,NA,NA,NA),result[13:14,],c(NA,NA,NA,NA,NA,NA,NA),result[15:16,],c(NA,NA,NA,NA,NA,NA,NA),result[17:18,],c(NA,NA,NA,NA,NA,NA,NA),result[19:20,],
              c(NA,NA,NA,NA,NA,NA,NA),result[21:22,])
row.names(result)=c('MASLD status','No masld','MASLD','Age group','Age<65 years','>=65 yeas',
                    'Sex','Females','Males','Race','Ethnic minorities','White',
                    'Deprivation','Quarter 1','Quarter 2','Quarter 3','Quarter 4',
                    'Education level','Above Colleage','Less than college',
                    'Activity','Regularly active','Inactive','Smoking','Never','Ever','Alcohol','Never_','Ever_',
                    'Metabolic syndrome','No','Yes')
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_masld_depression.csv')

########################################subgroup analysis by all covariates#########################################################
result<-data.frame(HR=NA,P=NA,CI=NA,lower=1,upper=1,id=seq(1:29),N=NA, Cases=NA)
row.names(result)=c('Age group','Age<65 years','>=65 yeas',
                    'Sex','Females','Males','Race','Ethnic minorities','White',
                    'Deprivation','Quarter 1','Quarter 2','Quarter 3','Quarter 4',
                    'Education level','Above Colleage','Less than college',
                    'Activity','Regularly active','Inactive','Smoking','Never','Ever','Alcohol','Never_','Ever_',
                    'Metabolic syndrome','No','Yes')

temp_data <-datas
##age
fit1<- coxph(Surv(time, status) ~ masld:age+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+1]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+1]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+1]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+1]=lower
  result$upper[i+1]=upper
  result$N[i+1]=length(which(temp_data$age==i-1))
  result$Cases[i+1]=length(which(temp_data$age==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
##sex
fit1<- coxph(Surv(time, status) ~ masld:gender+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+4]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+4]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+4]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+4]=lower
  result$upper[i+4]=upper
  result$N[i+4]=length(which(temp_data$gender==i-1))
  result$Cases[i+4]=length(which(temp_data$gender==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
##race
fit1<- coxph(Surv(time, status) ~ masld:race+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+7]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+7]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+7]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+7]=lower
  result$upper[i+7]=upper
  result$N[i+7]=length(which(temp_data$race==i-1))
  result$Cases[i+7]=length(which(temp_data$race==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
##deprivation
fit1<- coxph(Surv(time, status) ~ masld:deprivation+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:4)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+10]=round((summary(fit1)$conf.int[nrow_summary-4+i,'exp(coef)']),2)
  result$P[i+10]=signif((summary(fit1)$coefficients[nrow_summary-4+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-4+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-4+i,'upper .95'],2)
  result$CI[i+10]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+10]=lower
  result$upper[i+10]=upper
  result$N[i+10]=length(which(temp_data$deprivation==i-1))
  result$Cases[i+10]=length(which(temp_data$deprivation==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
##education
fit1<- coxph(Surv(time, status) ~ masld:university_education+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+15]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+15]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+15]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+15]=lower
  result$upper[i+15]=upper
  result$N[i+15]=length(which(temp_data$university_education==i-1))
  result$Cases[i+15]=length(which(temp_data$university_education==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
#activity
fit1<- coxph(Surv(time, status) ~ masld:activity+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+18]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+18]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+18]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+18]=lower
  result$upper[i+18]=upper
  result$N[i+18]=length(which(temp_data$activity==i-1))
  result$Cases[i+18]=length(which(temp_data$activity==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
##smoking
fit1<- coxph(Surv(time, status) ~ masld:smoking+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+21]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+21]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+21]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+21]=lower
  result$upper[i+21]=upper
  result$N[i+21]=length(which(temp_data$smoking==i-1))
  result$Cases[i+21]=length(which(temp_data$smoking==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
##alcohol
fit1<- coxph(Surv(time, status) ~ masld:alcohol+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+24]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+24]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+24]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+24]=lower
  result$upper[i+24]=upper
  result$N[i+24]=length(which(temp_data$alcohol==i-1))
  result$Cases[i+24]=length(which(temp_data$alcohol==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
##metabolic syndrome
fit1<- coxph(Surv(time, status) ~ masld:metabolic_syndrome+age+gender+race+deprivation+university_education+activity+
               smoking+alcohol+metabolic_syndrome,data = temp_data)
for (i in c(1:2)) {
  nrow_summary<-nrow(summary(fit1)$conf.int)
  result$HR[i+27]=round((summary(fit1)$conf.int[nrow_summary-2+i,'exp(coef)']),2)
  result$P[i+27]=signif((summary(fit1)$coefficients[nrow_summary-2+i,'Pr(>|z|)']),3)
  lower=round(summary(fit1)$conf.int[nrow_summary-2+i,'lower .95'],2)
  upper=round(summary(fit1)$conf.int[nrow_summary-2+i,'upper .95'],2)
  result$CI[i+27]=paste('(',lower,'-',upper,')',sep ='')
  result$lower[i+27]=lower
  result$upper[i+27]=upper
  result$N[i+27]=length(which(temp_data$metabolic_syndrome==i-1))
  result$Cases[i+27]=length(which(temp_data$metabolic_syndrome==i-1 & temp_data$status==1))
  rm(lower,upper,nrow_summary)
}
rm(i,fit1)
rm(temp_data)
result$CI=paste(result$HR,result$CI,sep=' ')
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_subgroup_masld_depression.csv')
