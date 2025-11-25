######################################PSM#####################################################################################
library(MatchIt)
index=which(datas$masld==0|datas$masld==1)
temp_datas=datas[index,]
rm(index)

row.names(temp_datas)=c(1:nrow(temp_datas))
#exact=c('race','alcohol','sedentary','age'),
m.out1 <- matchit(masld ~ age+gender+race+deprivation+university_education+activity+
                    smoking+alcohol+metabolic_syndrome,exact=c('age','deprivation','activity','smoking'),caliper = 0.05,
                  method = "nearest",distance = "glm", link = "logit",ratio=1,data = temp_datas)
match_index=m.out1$match.matrix
index_non=which(is.na(match_index)==0)
match_index=as.matrix(match_index[index_non,])
rm(index_non)
controlled_index=as.numeric(match_index[,1])
treated_index=as.numeric(row.names(match_index))
control_data=temp_datas[controlled_index,]
treat_data=temp_datas[treated_index,]
rm(controlled_index,treated_index,m.out1,match_index,temp_datas)

##cox proportional model
new_data=rbind(control_data,treat_data)
library("survival")
library("survminer")
model <- coxph(Surv(time, status) ~ masld,data = new_data)
model1=summary(model)
rm(model)

##############################demographic data for matched data#######################
summary_information=data.frame(matrix(ncol = 4, nrow =30))
colnames(summary_information)=c('Whole population','No MASLD','MASLD','P for difference')
row.names(summary_information)=c('Total N','Age group','Age<65 years','>=65 yeas',
                                 'Sex','Females','Males','Race','Ethnic minorities','White',
                                 'Deprivation','Quarter 1','Quarter 2','Quarter 3','Quarter 4',
                                 'Education level','Above Colleage','Less than college',
                                 'Activity','Regularly active','Inactive','Smoking','Never','Ever','Alcohol','Never_','Ever_',
                                 'Metabolic syndrome','No','Yes')
for (masld_status in c(10,0,1)) {
  if (masld_status==10){
    index_masld=c(1:nrow(new_data))#10=whole samples
    column_index=1
  }else{
    index_masld=which(new_data$masld==masld_status)
    column_index=masld_status+2
  }
  temp_datas=new_data[index_masld,]
  rm(index_masld)
  #N
  summary_information[1,column_index]=nrow(temp_datas)
  {
    #age
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$age==i-1))
      summary_information[i+2,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #gender
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$gender==i-1))
      summary_information[i+5,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #race
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$race==i-1))
      summary_information[i+8,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #deprivation
    for (i in c(1:4)) {
      n_temp=length(which(temp_datas$deprivation==i-1))
      summary_information[i+11,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #education
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$university_education==i-1))
      summary_information[i+16,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #activity
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$activity==i-1))
      summary_information[i+19,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #smoking
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$smoking==i-1))
      summary_information[i+22,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #alcohol
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$alcohol==i-1))
      summary_information[i+25,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp)
    #metabolic_syndrome
    for (i in c(1:2)) {
      n_temp=length(which(temp_datas$metabolic_syndrome==i-1))
      summary_information[i+28,column_index]=paste(n_temp,' (',round(n_temp/nrow(temp_datas)*100,2),'%)',sep='')
    }
    rm(n_temp,i)
    rm(temp_datas)
  }
}
rm(masld_status,column_index)
#P-value
{
  #age
  a<-t(matrix(c(tabulate(treat_data$age),tabulate(control_data$age)),nrow=2,ncol=2))
  summary_information[2,4]=signif(chisq.test(a)$p.value,3)
  #sex
  a<-t(matrix(c(tabulate(treat_data$gender),tabulate(control_data$gender)),nrow=2,ncol=2))
  summary_information[5,4]=signif(chisq.test(a)$p.value,3)
  #race
  a<-t(matrix(c(tabulate(treat_data$race),tabulate(control_data$race)),nrow=2,ncol=2))
  summary_information[8,4]=signif(chisq.test(a)$p.value,3)
  #deprivation
  a<-t(matrix(c(tabulate(treat_data$deprivation),tabulate(control_data$deprivation)),nrow=4,ncol=2))
  summary_information[11,4]=signif(chisq.test(a)$p.value,3)
  #education
  a<-t(matrix(c(tabulate(treat_data$university_education),tabulate(control_data$university_education)),nrow=2,ncol=2))
  summary_information[16,4]=signif(chisq.test(a)$p.value,3)
  #activity
  a<-t(matrix(c(tabulate(treat_data$activity),tabulate(control_data$activity)),nrow=2,ncol=2))
  summary_information[19,4]=signif(chisq.test(a)$p.value,3)
  #smoking
  a<-t(matrix(c(tabulate(treat_data$smoking),tabulate(control_data$smoking)),nrow=2,ncol=2))
  summary_information[22,4]=signif(chisq.test(a)$p.value,3)
  #alcohol
  a<-t(matrix(c(tabulate(treat_data$alcohol),tabulate(control_data$alcohol)),nrow=2,ncol=2))
  summary_information[25,4]=signif(chisq.test(a)$p.value,3)
  #metabolic_syndrome
  a<-t(matrix(c(tabulate(treat_data$metabolic_syndrome),tabulate(control_data$metabolic_syndrome)),nrow=2,ncol=2))
  summary_information[28,4]=signif(chisq.test(a)$p.value,3)
  rm(a,treat_data,control_data)
}
#write.csv(summary_information,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/demographic_masld_depression_psm.csv')





