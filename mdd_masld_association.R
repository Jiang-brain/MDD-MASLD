########################################association results####################################################
res.cox <- coxph(Surv(time, status) ~ depression+age+gender+race+deprivation+university_education+activity+
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
row.names(result)=c('Depression status','No depression','Depression','Age group','Age<65 years','>=65 yeas',
                    'Sex','Females','Males','Race','Ethnic minorities','White',
                    'Deprivation','Quarter 1','Quarter 2','Quarter 3','Quarter 4',
                    'Education level','Above Colleage','Less than college',
                    'Activity','Regularly active','Inactive','Smoking','Never','Ever','Alcohol','Never_','Ever_',
                    'Metabolic syndrome','No','Yes')
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_depression_masld.csv')


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
fit1<- coxph(Surv(time, status) ~ depression:age+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:gender+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:race+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:deprivation+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:university_education+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:activity+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:smoking+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:alcohol+age+gender+race+deprivation+university_education+activity+
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
fit1<- coxph(Surv(time, status) ~ depression:metabolic_syndrome+age+gender+race+deprivation+university_education+activity+
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
#write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_subgroup_depression_masld.csv')
#############################################################################################################

###################################negative control exposure#################################################
library("survival")
library("survminer")
datas$phone_head_side=factor(datas$phone_head_side)
res.cox <- coxph(Surv(time, status) ~ phone_head_side+age+gender+race+deprivation+university_education+activity+
                   smoking+alcohol+metabolic_syndrome,data = datas)
summary(res.cox)
