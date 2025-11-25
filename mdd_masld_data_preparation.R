#########################################depression####################################
library(data.table)
depression_data<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv",
                       select=c('eid','53-0.0','130894-0.0','130895-0.0','130896-0.0','130897-0.0'),
                       header = T,sep=",")
depression_data=data.frame(depression_data)
diagnosis_depression<-data.frame(id=depression_data$eid,diagnosis_date=as.Date(NA),source=NA,time_interval=NA,attending_date=depression_data$X53.0.0)

diagnosis_depression$diagnosis_date=pmin(depression_data$X130894.0.0,depression_data$X130896.0.0,na.rm = TRUE)
diagnosis_depression$time_interval=diagnosis_depression$diagnosis_date-diagnosis_depression$attending_date
rm(depression_data)

depression=rep(0, nrow(diagnosis_depression))
index_depression=which(diagnosis_depression$time_interval<=0&diagnosis_depression$diagnosis_date>as.Date('1910-01-01'))###==0
index_na=which(diagnosis_depression$diagnosis_date<as.Date('1910-01-01'))
depression[index_depression]=1
#depression[index_depression]=diagnosis_depression$time_interval[index_depression]
#index_depression=which(diagnosis_depression$time_interval==0)###==0
#depression[index_depression]=-0.5
depression[index_na]=NA
rm(index_depression,index_na,diagnosis_depression)

###########################fli########################
datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv",header = T,sep=",")
datas=data.frame(datas)
tg=datas$X30870.0.0*88.545#from mmol/L to mg/dL
bmi=datas$X21001.0.0
ggt=datas$X30730.0.0
waist_circumference=datas$X48.0.0
fli_l=0.953*log(tg)+0.139*bmi+0.718*log(ggt)+0.053*waist_circumference-15.745
fli_score=exp(fli_l)/(1+exp(fli_l))*100
rm(fli_l,bmi,tg,ggt,waist_circumference)

##########five metabolic/cardiovascular disease risk factors######
##obesity
obesity=rep(0, nrow(datas))

index_obesity=which(datas$X21001.0.0>=25)
obesity[index_obesity]=1
rm(index_obesity)
index_asian=which((datas$X21000.0.0 %in% c(3001,3002,3,3003,3004,5)) & (datas$X21001.0.0>=23))
obesity[index_asian]=1
rm(index_asian)

waist_circumference=datas$X48.0.0
index_female=which(datas$X31.0.0==0&waist_circumference>=80)
index_male=which(datas$X31.0.0==1&waist_circumference>=94)
obesity[index_female]=1
obesity[index_male]=1

index_non=which(is.na(waist_circumference)&is.na(datas$X21001.0.0))
obesity[index_non]=NA
rm(index_female,index_male,index_non,waist_circumference)

#############hypertension
hypertension=rep(0, nrow(datas))
dbp=rowMeans(cbind(datas$X4079.0.0,datas$X4079.0.1),na.rm = TRUE)
sbp=rowMeans(cbind(datas$X4080.0.0,datas$X4080.0.1),na.rm = TRUE)
hypertension_diagnosis=pmax(datas$X6150.0.0,datas$X6150.0.1,datas$X6150.0.2,datas$X6150.0.3,na.rm = TRUE)
index_hypertension=which(hypertension_diagnosis==4)
index_high_dbp=which(dbp>=85)
index_high_sbp=which(sbp>=130)
hypertension[unique(c(index_hypertension,index_high_dbp,index_high_sbp))]=1
index_nan=which(hypertension_diagnosis==-3&is.na(dbp)==TRUE&is.na(sbp)==TRUE)
hypertension[index_nan]=NA
rm(dbp,sbp,index_high_dbp,index_high_sbp,index_nan,index_hypertension,hypertension_diagnosis)

#blood pressure medication
medications_female=cbind(datas$X6153.0.0,datas$X6153.0.1,datas$X6153.0.2,datas$X6153.0.3)
medications_female=medications_female==2
medications_female=rowSums(medications_female,na.rm = TRUE)
index_female=which(medications_female==1)
rm(medications_female)

medications_male=cbind(datas$X6177.0.0,datas$X6177.0.1,datas$X6177.0.2)
medications_male=medications_male==2
medications_male=rowSums(medications_male,na.rm = TRUE)
index_male=which(medications_male==1)
rm(medications_male)

index_nan=which((is.na(hypertension)==TRUE)&((datas$X6153.0.0==-1|-3)|(datas$X6177.0.0==-1|-3)))
hypertension[index_nan]=NA

hypertension[index_male]=1
hypertension[index_female]=1
rm(index_female,index_male,index_nan)

#############triglyceride
triglyceride=rep(0, nrow(datas))

index_high_tg=which(datas$X30870.0.0>=1.7)

#lipid-lowering
medications_female=cbind(datas$X6153.0.0,datas$X6153.0.1,datas$X6153.0.2,datas$X6153.0.3)
medications_female=medications_female==1
medications_female=rowSums(medications_female,na.rm = TRUE)
index_female=which(medications_female==1)
rm(medications_female)

medications_male=cbind(datas$X6177.0.0,datas$X6177.0.1,datas$X6177.0.2)
medications_male=medications_male==1
medications_male=rowSums(medications_male,na.rm = TRUE)
index_male=which(medications_male==1)
rm(medications_male)

index_non=which((is.na(datas$X30870.0.0)==TRUE)&((datas$X6153.0.0==-1|-3)|(datas$X6177.0.0==-1|-3)))

triglyceride[unique(c(index_high_tg,index_female,index_male))]=1
triglyceride[index_non]=NA
rm(index_female,index_male,index_high_tg,index_non)

#############hdl
hdl_cholesterol=rep(0,nrow(datas))
index_low_hdl_male=which(datas$X30760.0.0<=1.0&datas$X31.0.0==1)
index_low_hdl_female=which(datas$X30760.0.0<=1.3&datas$X31.0.0==0)
index_no=which(is.na(datas$X30760.0.0)==TRUE)
hdl_cholesterol[index_low_hdl_female]=1
hdl_cholesterol[index_low_hdl_male]=1
hdl_cholesterol[index_no]=NA
rm(index_low_hdl_female,index_low_hdl_male,index_no)

#lipid-lowering
medications_female=cbind(datas$X6153.0.0,datas$X6153.0.1,datas$X6153.0.2,datas$X6153.0.3)
medications_female=medications_female==1
medications_female=rowSums(medications_female,na.rm = TRUE)
index_female=which(medications_female==1)
rm(medications_female)

medications_male=cbind(datas$X6177.0.0,datas$X6177.0.1,datas$X6177.0.2)
medications_male=medications_male==1
medications_male=rowSums(medications_male,na.rm = TRUE)
index_male=which(medications_male==1)
rm(medications_male)

hdl_cholesterol[unique(c(index_female,index_male))]=1
rm(index_female,index_male)

############t2d
diabetes=rep(0,nrow(datas))

glucose=datas$X30740.0.0
index_glucose_higher=which(glucose>=5.6)

HbA1c=datas$X30750.0.0
index_HbA1c_higher=which(HbA1c>=39)

index_na=which((is.na(glucose)==TRUE)&(is.na(HbA1c)==TRUE))

diabetes[index_HbA1c_higher]=1
diabetes[index_glucose_higher]=1
diabetes[index_na]=NA##overlap to modify
rm(index_glucose_higher,index_HbA1c_higher,index_na,glucose,HbA1c)

##t2d diagnosis from icd
diagnosis_t2d=datas$X130708.0.0-datas$X53.0.0
index_t2d_diagnosis_baseline=which(diagnosis_t2d<=0)
diabetes[index_t2d_diagnosis_baseline]=1
rm(diagnosis_t2d,index_t2d_diagnosis_baseline)

##t2d from self-report and age diagnosed
diabetes_age=datas$X2976.0.0
index_diabetes_age=which(diabetes_age>=40)
diabetes[index_diabetes_age]=1
rm(diabetes_age,index_diabetes_age)

##calculate cardiometabolic factors
cardiometabolic_factor=cbind(diabetes,hdl_cholesterol,hypertension,obesity,triglyceride)
rm(diabetes,hdl_cholesterol,hypertension,obesity,triglyceride)
cardiometabolic_amount=rowSums(cardiometabolic_factor,na.rm = TRUE)
cardiometabolic_nan=rowSums(is.na(cardiometabolic_factor))

index_cardiometabolic=which(cardiometabolic_amount>=1)#at least one factor
index_cardiometabolic_nan=which(cardiometabolic_nan==5|(cardiometabolic_amount==0&cardiometabolic_nan>=1))##nan values
rm(cardiometabolic_factor)

cardio=rep(0,nrow(datas))
cardio[index_cardiometabolic]=1
cardio[index_cardiometabolic_nan]=NA
rm(index_cardiometabolic,index_cardiometabolic_nan)

########################################MetSyndrome
metSyn=rep(NA,nrow(datas))
index_met=which(cardiometabolic_amount>=3)
index_no_met=which(cardiometabolic_amount<3&cardiometabolic_nan==0)
metSyn[index_met]=1
metSyn[index_no_met]=0
rm(index_no_met,index_met)
#for some components with missing values
index_no=which(cardiometabolic_nan==1&cardiometabolic_amount<2)
index_no1=which(cardiometabolic_nan==2&cardiometabolic_amount==0)
metSyn[index_no]=0
metSyn[index_no1]=0
rm(cardiometabolic_amount,cardiometabolic_nan,index_no,index_no1)

#############################################################################

####################alcohol intake units,glass-g/d####
alcohol_frequency=datas$X1558.0.0
index_no_answer=which(alcohol_frequency==-3)
alcohol_frequency[index_no_answer]=NA
rm(index_no_answer)

red_wine_month=datas$X4407.0.0
index_no_answer=which(red_wine_month<0)
red_wine_month[index_no_answer]=NA
#red_wine_month=red_wine_month*16.8/4.3##convert to unit/wk
red_wine_month=red_wine_month*1.7*8/4.3##convert to unit/wk
rm(index_no_answer)

red_wine_week=datas$X1568.0.0
index_no_answer=which(red_wine_week<0)
red_wine_week[index_no_answer]=NA
#red_wine_week=red_wine_week*16.8
red_wine_week=red_wine_week*1.7*8
rm(index_no_answer)

champagne_month=datas$X4418.0.0#white wine
index_no_answer=which(champagne_month<0)
champagne_month[index_no_answer]=NA
#champagne_month=champagne_month*16.8/4.3
champagne_month=champagne_month*1.7*8/4.3
rm(index_no_answer)

champagne_week=datas$X1578.0.0
index_no_answer=which(champagne_week<0)
champagne_week[index_no_answer]=NA
#champagne_week=champagne_week*16.8
champagne_week=champagne_week*1.7*8
rm(index_no_answer)

beer_month=datas$X4429.0.0
index_no_answer=which(beer_month<0)
beer_month[index_no_answer]=NA
#beer_month=beer_month*16/4.3
beer_month=beer_month*2.4*8/4.3
rm(index_no_answer)

beer_week=datas$X1588.0.0
index_no_answer=which(beer_week<0)
beer_week[index_no_answer]=NA
#beer_week=beer_week*16
beer_week=beer_week*2.4*8
rm(index_no_answer)

spirits_month=datas$X4440.0.0
index_no_answer=which(spirits_month<0)
spirits_month[index_no_answer]=NA
spirits_month=spirits_month*8/4.3
rm(index_no_answer)

spirits_week=datas$X1598.0.0
index_no_answer=which(spirits_week<0)
spirits_week[index_no_answer]=NA
spirits_week=spirits_week*8
rm(index_no_answer)

fortified_month=datas$X4451.0.0
index_no_answer=which(fortified_month<0)
fortified_month[index_no_answer]=NA
#fortified_month=fortified_month*14.08/4.3
fortified_month=fortified_month*1.2*8/4.3
rm(index_no_answer)

fortified_week=datas$X1608.0.0
index_no_answer=which(fortified_week<0)
fortified_week[index_no_answer]=NA
#fortified_week=fortified_week*14.08
fortified_week=fortified_week*1.2*8
rm(index_no_answer)

others_month=datas$X4462.0.0
index_no_answer=which(others_month<0)
others_month[index_no_answer]=NA
others_month=others_month*12/4.3
rm(index_no_answer)

others_week=datas$X5364.0.0
index_no_answer=which(others_week<0)
others_week[index_no_answer]=NA
others_week=others_week*12
rm(index_no_answer)

alcohol=cbind(red_wine_month,red_wine_week,champagne_month,champagne_week,beer_month,beer_week,
              spirits_month,spirits_week,fortified_month,fortified_week,others_month,others_week)
alcohol_amount=rowSums(alcohol,na.rm = TRUE)
number_na<-rowSums(is.na(alcohol))
rm(red_wine_month,red_wine_week,champagne_month,champagne_week,beer_month,beer_week,
   spirits_month,spirits_week,fortified_month,fortified_week,others_month,others_week)
index_never=which(alcohol_frequency==6)
index_no_answer=which(is.na(alcohol_frequency)==TRUE)
index_nan=which(number_na==12&alcohol_frequency!=6)
index_special_occuasion=which(alcohol_frequency==5&number_na==12)
index_1_3_month=which(alcohol_frequency==4&number_na==12)
alcohol_amount[index_never]=0
alcohol_amount[index_no_answer]=NA
alcohol_amount[index_nan]=NA
alcohol_amount[index_special_occuasion]=0
alcohol_amount[index_1_3_month]=0
rm(index_never,index_no_answer,index_nan,index_special_occuasion,index_1_3_month)
rm(alcohol_frequency,alcohol,number_na)

sex=datas$X31.0.0
alcohol=rep(0,nrow(datas))
index_male_excessive=which(sex==1&alcohol_amount/7>30)
index_female_excessive=which(sex==0&alcohol_amount/7>20)
index_nan=which(is.na(alcohol_amount)==TRUE)
alcohol[index_male_excessive]=1
alcohol[index_female_excessive]=1
alcohol[index_nan]=NA
rm(index_female_excessive,index_male_excessive,index_nan,sex)
rm(alcohol_amount)
################################################remove participants with other liver diseases###########################################
icd_10<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/nafld_diagnosis.csv", select=c(seq(2,260)),
              header = T,sep=",")
icd_10_date<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/nafld_diagnosis.csv", select=c(seq(261,519)),
                   header = T,sep=",")
attending_date<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv",select=c('eid','53-0.0'),
                      header = T,sep=",")
attending_date=data.frame(attending_date)
icd_10=data.frame(icd_10)
icd_10_date=data.frame(icd_10_date)

diagnosis_other_disease<-data.frame(id=attending_date$eid,diagnosis_date=as.Date(NA),attending_dates=attending_date$X53.0.0)
for (diagnosis_col in 1:ncol(icd_10)){
  diagnosis<-which(icd_10[,diagnosis_col] %in% c('K700','K701','K702','K703','K704','K709','B160','B169','B170','B171','B172','B178','B179','B180','B181','B182','B188','B189','B199',
                                                 'K830','K743','K754','E831','E830','E880','I820','K765','K739','K732','K744','K745',
                                                 'F100','F101','F102','F103','F104','F105','F106','F107','F108','F109',
                                                 'F110','F111','F112','F113','F114','F115','F117','F119',
                                                 'F120','F121','F122','F123','F125','F128','F129','F130','F131','F132','F133','F134','F139',
                                                 'F140','F141','F142','F145','F149','F161','F162','F163','F165','F167','F168','F169',
                                                 'F181','F182','F183','F185','F189','F190','F191','F192','F193','F194','F195','F198','F199',
                                                 'X6509','X651','X652','X653','X654','X655','X656','X658','X659',
                                                 'E244','G621','I426','K292','G312','G721','K852','K860','T510','T519','Y573','Z502','Z714','Z721'))
  if (length(diagnosis)>0)
  {
    diagnosis_other_disease$diagnosis_date[diagnosis]=pmin(as.Date(icd_10_date[diagnosis,diagnosis_col],"%m/%d/%Y"),diagnosis_other_disease$diagnosis_date[diagnosis],na.rm = TRUE)
  }
  #print(c(diagnosis_col))
  rm(diagnosis)
}
rm(diagnosis_col,icd_10,icd_10_date,attending_date)

diagnosis_other_disease$time_interval=as.IDate(diagnosis_other_disease$diagnosis_date)-diagnosis_other_disease$attending_dates
########################################################################################################################################
#############################################MASLD baseline###########################
masld=rep(0,nrow(datas))
index_masld=which(fli_score>=60&alcohol==0&cardio==1)
index_excessive_alcohol_masld=which(fli_score>=60&alcohol==1& !is.na(cardio))#remove? sensitivity analysis
index_non=which(is.na(fli_score)|is.na(alcohol)|is.na(cardio))
index_masld_alcohol_no_cardio=which(fli_score>=60&alcohol==0&cardio==0)#remove or not remove
index_control=which(fli_score<60 & !is.na(alcohol) & !is.na(cardio))
masld[index_masld]=2
masld[index_excessive_alcohol_masld]=NA
masld[index_non]=NA
masld[index_masld_alcohol_no_cardio]=NA####?
masld[index_control]=0
rm(index_masld,index_excessive_alcohol_masld,index_masld_alcohol_no_cardio,index_non,index_control)

index_grade_1_masld=which(fli_score>=30&fli_score<60 & !is.na(alcohol) & !is.na(cardio))
masld[index_grade_1_masld]=1
rm(index_grade_1_masld)
rm(alcohol,cardio,fli_score)

########################################ICD-nafld: follow-up######################################################################
############ ICD-10 extraction of first diagnosis
library(data.table)
icd_10<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/nafld_diagnosis.csv", select=c(seq(2,260)),
              header = T,sep=",")
icd_10_date<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/nafld_diagnosis.csv", select=c(seq(261,519)),
                   header = T,sep=",")
attending_date<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv",select=c('eid','53-0.0'),
                      header = T,sep=",")

attending_date=data.frame(attending_date)
icd_10=data.frame(icd_10)
icd_10_date=data.frame(icd_10_date)


diagnosis_nafld<-data.frame(id=attending_date$eid,diagnosis_date=as.Date(NA),attending_dates=attending_date$X53.0.0)
for (diagnosis_col in 1:ncol(icd_10)){
  diagnosis<-which(icd_10[,diagnosis_col] %in% c('K758','K760'))#broad definition c('K758','K760','K768','K769','K740','K741','K742','K746')
  if (length(diagnosis)>0)
  {
    diagnosis_nafld$diagnosis_date[diagnosis]=pmin(as.Date(icd_10_date[diagnosis,diagnosis_col],"%m/%d/%Y"),diagnosis_nafld$diagnosis_date[diagnosis],na.rm = TRUE)#if there is a diagnosis before this column
  }
  #print(c(diagnosis_col))
}
rm(diagnosis)
rm(diagnosis_col)
rm(icd_10)
rm(icd_10_date)
rm(attending_date)
diagnosis_nafld$time_interval=as.IDate(diagnosis_nafld$diagnosis_date)-diagnosis_nafld$attending_dates
index_nonNan=which(is.na(diagnosis_nafld$diagnosis_date)==FALSE)
diagnosis_nafld=diagnosis_nafld[index_nonNan,]
rm(index_nonNan)


############################################# max date of icd record
library(data.table)
icd_10_date<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/nafld_diagnosis.csv", select=c(seq(261,519)),
                   header = T,sep=",")
icd_10_date=data.frame(icd_10_date)

for (temp_col in 1:ncol(icd_10_date)){
  temp_diagnosis<-as.Date(icd_10_date[,temp_col],"%m/%d/%Y")
  if (temp_col==1){
    result_date=max(temp_diagnosis,na.rm = TRUE)
  }else{
    result_date=c(result_date,max(temp_diagnosis,na.rm = TRUE))
  }
  rm(temp_diagnosis)
}
censoring_date=max(result_date,na.rm = TRUE)
rm(result_date)
rm(temp_col)
rm(icd_10_date)


############################################# extraction of death information
library(data.table)
datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv", select=c('eid','40000-0.0','40001-0.0','53-0.0'),
             header = T,sep=",")
death=data.frame(datas)
colnames(death)=c('id','death_date','death_reason','attending_date')
rm(datas)
index_death<-which(!is.na(death$death_date)==TRUE)
death_information<-death[index_death,]
rm(death)
rm(index_death)
death_k760_k758<-which(death_information$death_reason %in% c('K758','K760'))##broad definition c('K758','K760','K768','K769','K740','K741','K742','K746')

death_information$label=0;
death_information$label[death_k760_k758]=1
death_information$time_interval=death_information$death_date-death_information$attending_date

missing_reason=which(death_information$death_reason=='')
death_information$death_reason[missing_reason]=NaN
rm(death_k760_k758)
rm(missing_reason)


######## lost of follow-up
datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv", 
             select=c('eid','190-0.0','191-0.0','53-0.0'),header = T,sep=",")
lost_follow_up=data.frame(datas)
colnames(lost_follow_up)=c('id','lost_reason','lost_date','attending_date')
rm(datas)

index_noNaN<-which(is.na(lost_follow_up$lost_date)==FALSE)
lost_follow_up<-lost_follow_up[index_noNaN,]
lost_follow_up$lost_time_interval<-lost_follow_up$lost_date-lost_follow_up$attending_date
rm(index_noNaN)

################ censoring data
data_censoring<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv", select=c('eid','53-0.0'),header = T,sep=",")
data_censoring=data.frame(data_censoring)
colnames(data_censoring)=c('id','attending_date')
data_censoring$censoring=rep(as.IDate(censoring_date),nrow(data_censoring))-data_censoring$attending_date
rm(censoring_date)

data_censoring$death=NA
data_censoring$death_reason=NA
death_index=match(death_information$id,data_censoring$id)
data_censoring$death[death_index]=death_information$time_interval
data_censoring$death_reason[death_index]=death_information$label
rm(death_index)
rm(death_information)

data_censoring$lost=NA
data_censoring$lost_reason=NA
lost_index=match(lost_follow_up$id,data_censoring$id)
data_censoring$lost[lost_index]=lost_follow_up$lost_time_interval
data_censoring$lost_reason[lost_index]=lost_follow_up$lost_reason
rm(lost_index)
rm(lost_follow_up)

data_censoring$diagnosis=NA
diagnosis_index=match(diagnosis_nafld$id,data_censoring$id)
data_censoring$diagnosis[diagnosis_index]=diagnosis_nafld$time_interval
rm(diagnosis_index)
rm(diagnosis_nafld)


######calculate censoring days and label
#diagnosis
data_censoring$label=0

index_diagnosis=which(data_censoring$diagnosis>0)
data_censoring$label[index_diagnosis]=1
data_censoring$censoring[index_diagnosis]=data_censoring$diagnosis[index_diagnosis]
rm(index_diagnosis)
index_diagnosis=which(data_censoring$diagnosis<=0)# remove subjects with diagnosis before enrollment
data_censoring$label[index_diagnosis]=NA
rm(index_diagnosis)

#death
index_death_wo_diagnosis=which(data_censoring$death_reason==1 & data_censoring$label==0)
data_censoring$censoring[index_death_wo_diagnosis]=data_censoring$death[index_death_wo_diagnosis]
data_censoring$label[index_death_wo_diagnosis]=1
rm(index_death_wo_diagnosis)

index_death_no_diagnosis=which(data_censoring$death_reason==0 & data_censoring$label==0)
data_censoring$censoring[index_death_no_diagnosis]=data_censoring$death[index_death_no_diagnosis]
data_censoring$label[index_death_no_diagnosis]=0
rm(index_death_no_diagnosis)

#lost follow-up
index_lost_wo_diagnosis_enrollment=which(is.na(data_censoring$lost)==FALSE & is.na(data_censoring$label)==FALSE)
data_censoring$censoring[index_lost_wo_diagnosis_enrollment]=pmin(data_censoring$censoring[index_lost_wo_diagnosis_enrollment],data_censoring$lost[index_lost_wo_diagnosis_enrollment])
rm(index_lost_wo_diagnosis_enrollment)



########################################remove baseline masld####################################
index_masld_baseline=which(masld==2)
data_censoring$label[index_masld_baseline]=NA
index_nan=which(is.na(masld)==1)
data_censoring$label[index_nan]=NA
rm(index_masld_baseline,masld,index_nan)

index_other_liver_disease=which(diagnosis_other_disease$time_interval<=0)# remove subjects with other liver diseases before enrollment
data_censoring$label[index_other_liver_disease]=NA
rm(index_other_liver_disease,diagnosis_other_disease)

####################################################covariate variables###############################################
datas<-fread("/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/frailty_nafld.csv",header = T,sep=",")
datas=data.frame(datas)
age=datas$X21003.0.0
sites=datas$X54.0.0

#exclude dismatch between self-reported sex and genetic sex
sex=datas$X31.0.0
genetic_sex=datas$X22001.0.0
index_not_match=which((sex==0&genetic_sex==1)|(sex==1&genetic_sex==0))
sex[index_not_match]=NA
rm(genetic_sex,index_not_match)

#deprivation
deprivation=datas$X22189.0.0
quarter_point=quantile(deprivation, probs = seq(0, 1, 1/4),na.rm = TRUE)
index_q1=which(deprivation<quarter_point[2])
index_q2=which(deprivation>=quarter_point[2]&deprivation<=quarter_point[3])
index_q3=which(deprivation>quarter_point[3]&deprivation<quarter_point[4])
index_q4=which(deprivation>=quarter_point[4])
deprivation[index_q1]=0
deprivation[index_q2]=1
deprivation[index_q3]=2
deprivation[index_q4]=3
rm(quarter_point,index_q1,index_q2,index_q3,index_q4)


##race
ethnicity=cbind(datas$X21000.0.0,datas$X21000.1.0,datas$X21000.2.0)
index_nan=which(ethnicity[,1]==-1|ethnicity[,1]==-3)
ethnicity[index_nan,1]=pmax(ethnicity[index_nan,1],ethnicity[index_nan,2],ethnicity[index_nan,3],na.rm = TRUE)
race=ethnicity[,1]
rm(ethnicity,index_nan)
index_white=which((race %in% c(1,1001,1002,1003))==TRUE)
index_nonwhite=which((race %in% c(1,1001,1002,1003))==FALSE)
index_not_answer=which(race==-1|race==-3)
race[index_white]=1
race[index_nonwhite]=0
race[index_not_answer]=NA
rm(index_white,index_nonwhite,index_not_answer)
ethnicity=race
rm(race)


#income
income=datas$X738.0.0
index_not_answe=which(income==-3)
index_non=which(income==-1)
index_low=which(income<=3&income>=1)
index_medium=which(income==4)
index_high=which(income==5)
income[index_non]=3
income[index_low]=2
income[index_medium]=1
income[index_high]=0
income[index_not_answe]=3
rm(index_non,index_low,index_medium,index_high,index_not_answe)

#smoking
smoking=datas$X20116.0.0
index_non=which(is.na(smoking)==TRUE)
smoking[index_non]=NA
index_no_answer=which(smoking==-3)
smoking[index_no_answer]=NA
index_current=which(smoking==2)
smoking[index_current]=1
rm(index_non,index_no_answer,index_current)

#alcohol status
alcohol_status=datas$X20117.0.0
index_no_answer=which(alcohol_status<0)
alcohol_status[index_no_answer]=NA
index_current=which(alcohol_status==2)
alcohol_status[index_current]=1
rm(index_no_answer,index_current)

frequency=datas$X1558.0.0
index_never=which(is.na(alcohol_status)&frequency==6)
alcohol_status[index_never]=0##never alcohol
rm(frequency,index_never)

#alcohol frequency
#alcohol_amount=alcohol_amount

#education
education=cbind(datas$X6138.0.0,datas$X6138.0.1,datas$X6138.0.2,datas$X6138.0.3,datas$X6138.0.4,datas$X6138.0.5)
for (i in 1:6){
  index_no_answer=which(education[,i]==-3)
  education[index_no_answer,i]=NA
  index_none=which(education[,i]==-7)
  education[index_none,i]=7
  
  index_20=which(education[,i]==1)
  index_13=which(education[,i]==2)
  index_10=which(education[,i]==3|education[,i]==4)
  index_19=which(education[,i]==5)
  index_15=which(education[,i]==6)
  index_7=which(education[,i]==7)
  education[index_20,i]=20
  education[index_13,i]=13
  education[index_10,i]=10
  education[index_19,i]=19
  education[index_15,i]=15
  education[index_7,i]=7
  rm(index_20,index_13,index_10,index_19,index_15,index_7,index_none,index_no_answer)
}
rm(i)
university_education=pmax(education[,1],education[,2],education[,3],education[,4],education[,5],education[,6],na.rm = TRUE)
index_university=which(university_education==20)
index_no_university=which(university_education<20&university_education>7)
index_no_education=which(university_education==7)
university_education[index_no_education]=1
university_education[index_no_university]=1
university_education[index_university]=0
rm(index_no_education,index_no_university,index_university,education)

#watching tv
sedentary=datas$X1070.0.0
index_na=which(sedentary==-1|sedentary==-3)
index_low=which((sedentary<=2&sedentary>=0)|sedentary==-10)
index_medium=which(sedentary>2&sedentary<4)
index_high=which(sedentary>=4)
sedentary[index_na]=NA
sedentary[index_low]=0
sedentary[index_medium]=1
sedentary[index_high]=2
rm(index_na,index_low,index_medium,index_high)

##physical activity
{
  activity_type=cbind(datas$X6164.0.0,datas$X6164.0.1,datas$X6164.0.2,datas$X6164.0.3,datas$X6164.0.4)
  for (i in 1:5){
    index_light=which(activity_type[,i]==4)
    activity_type[index_light,i]=0
    
    index_not_answer=which(activity_type[,i]==-3)
    activity_type[index_not_answer,i]=NA
    rm(index_light,index_not_answer)
  }
  
  activity1=pmax(activity_type[,1],activity_type[,2],activity_type[,3],activity_type[,4],activity_type[,5],na.rm = TRUE)
  rm(i,activity_type)
  
  light_activity_frequency=datas$X1011.0.0
  index_none_activity=which(activity1==-7)
  index_light=which(activity1==0&light_activity_frequency<=3&light_activity_frequency>0)
  index_Nan=which(activity1==0&light_activity_frequency<0)
  
  activity1[index_light]=10
  activity1[index_none_activity]=10
  activity1[index_Nan]=NA
  
  index_high=which(activity1<10)
  activity1[index_high]=0
  activity1[index_none_activity]=1
  activity1[index_light]=1
  
  rm(index_none_activity,index_light,index_Nan,index_high,light_activity_frequency)
}
rm(activity1)

##physical activity
moderate_days=datas$X884.0.0
index_nan=which(moderate_days<0)
moderate_days[index_nan]=NA
rm(index_nan)

vigorous_days=datas$X904.0.0
index_nan=which(vigorous_days<0)
vigorous_days[index_nan]=NA
rm(index_nan)

moderate_duration=datas$X894.0.0
index_nan=which(moderate_duration<0)
moderate_duration[index_nan]=0###set to one
rm(index_nan)

vigorous_duration=datas$X914.0.0
index_nan=which(vigorous_duration<0)
vigorous_duration[index_nan]=0###set to one
rm(index_nan)

moderate_minutes=rep(NA,nrow(datas))
index_moderate=which(is.na(moderate_duration)==FALSE)
moderate_minutes[index_moderate]=moderate_days[index_moderate]*moderate_duration[index_moderate]
rm(index_moderate)
index_0=which(moderate_days==0)
moderate_minutes[index_0]=0
rm(index_0,moderate_duration)

vigorous_minutes=rep(NA,nrow(datas))
index_vigorous=which(is.na(vigorous_duration)==FALSE)
vigorous_minutes[index_vigorous]=vigorous_days[index_vigorous]*vigorous_duration[index_vigorous]
rm(index_vigorous)
index_0=which(vigorous_days==0)
vigorous_minutes[index_0]=0
rm(index_0,vigorous_duration)

total_activity=rowSums(cbind(moderate_minutes,2*vigorous_minutes),na.rm = TRUE)

activity=rep(1,nrow(datas))
#index_activity=which(total_activity>=150|moderate_days>=5|vigorous_days>=1)
index_activity=which(total_activity>=150)
index_nan=which((is.na(moderate_days)==TRUE)&(is.na(vigorous_days)==TRUE))
index_nan1=which(total_activity<90&(is.na(moderate_minutes)==TRUE | is.na(vigorous_minutes)==TRUE))
activity[index_nan]=NA
activity[index_nan1]=NA
activity[index_activity]=0
rm(moderate_days,vigorous_days,moderate_minutes,vigorous_minutes,index_activity,index_nan,total_activity)
rm(index_nan1)

#combine activity
# index_non_activity=which(is.na(activity)==TRUE)
# activity[index_non_activity]=activity1[index_non_activity]
# rm(index_non_activity,activity1)

#obesity
bmi=datas$X21001.0.0
obesity=rep(NA,nrow(datas))
index_underweight=which(bmi<18.5)
index_normal=which(bmi<25&bmi>=18.5)
index_overweight=which(bmi>=25&bmi<30)
index_obese1=which(bmi>=30)
#index_obese2=which(bmi>=35&bmi<40)
#index_obese3=which(bmi>=40)
obesity[index_underweight]=0
obesity[index_normal]=1
obesity[index_overweight]=2
obesity[index_obese1]=3
#obesity[index_obese2]=4
#obesity[index_obese3]=5
#rm(index_underweight,index_normal,index_overweight,index_obese1,index_obese2,index_obese3)
rm(index_underweight,index_normal,index_overweight,index_obese1,bmi)



##########new added covariates#######
#sleep
sleep=datas$X1160.0.0
index_no_answer=which(sleep<0)
index_excess=which(sleep>9)
index_lack=which(sleep<7)
index_normal=which(sleep>=7&sleep<=9)
sleep[index_no_answer]=NA
sleep[index_normal]=0
sleep[index_excess]=1
sleep[index_lack]=1
rm(index_no_answer,index_normal,index_excess,index_lack)

#cancer history
cancer=datas$X2453.0.0
index_not_answer=which(cancer<0)
cancer[index_not_answer]=NA
rm(index_not_answer)

#long-term illness
illness=datas$X2188.0.0
index_not_answer=which(illness<0)
illness[index_not_answer]=NA
rm(index_not_answer)

#employement
employement=cbind(datas$X6142.0.0,datas$X6142.0.1,datas$X6142.0.2,datas$X6142.0.3,datas$X6142.0.4,datas$X6142.0.5,datas$X6142.0.6)
for (i in 1:7){
  index_no_answer=which(employement[,i]<0)
  employement[index_no_answer,i]=NA
  rm(index_no_answer)
}
rm(i)
employement_status=pmin(employement[,1],employement[,2],employement[,3],employement[,4],employement[,5],employement[,6],employement[,7],na.rm = TRUE)
index_employed=which(employement_status==1)
index_not_employed=which(employement_status>1)
employement_status[index_employed]=1
employement_status[index_not_employed]=0
rm(employement,index_employed,index_not_employed)



##############################metabolites:!###############################################
metabolites=datas[,110:358]
pdff=datas$X24352.2.0

#########################################PHQ-4 scores###################################
phq=cbind(datas$X2050.0.0,datas$X2060.0.0,datas$X2070.0.0,datas$X2080.0.0)
for (i in 1:4) {
  index_no_answer=which(phq[,i]==-3|phq[,i]==-1)
  phq[index_no_answer,i]=NA
  rm(index_no_answer)
}
phq_4=rowSums(phq,na.rm = FALSE)
rm(i)

#########################################negative control measures########
phone_head_side=datas$X1150.0.0
index_nan=which(phone_head_side<0)
phone_head_side[index_nan]=NA
rm(index_nan)

#####################################combine data#########################################################################
data=data.frame(status=data_censoring$label,time=data_censoring$censoring,depression=depression,
                metabolic_syndrome=metSyn,phq_4=phq_4,
                age=age,gender=sex,sites=sites,race=ethnicity,
                university_education=university_education,deprivation=deprivation,
                smoking=smoking,alcohol=alcohol_status,activity=activity)

data=data.frame(status=data_censoring$label,time=data_censoring$censoring,depression=depression,
                metabolic_syndrome=metSyn,
                age=age,gender=sex,sites=sites,race=ethnicity,
                university_education=university_education,deprivation=deprivation,
                smoking=smoking,alcohol=alcohol_status,activity=activity)

# data=data.frame(status=data_censoring$label,time=data_censoring$censoring,depression=depression,
#                 metabolic_syndrome=metSyn,pdff=pdff,
#                 age=age,gender=sex,sites=sites,race=ethnicity,
#                 university_education=university_education,deprivation=deprivation,
#                 smoking=smoking,alcohol=alcohol_status,activity=activity)

# data=data.frame(status=data_censoring$label,time=data_censoring$censoring,depression=depression,
#                 metabolic_syndrome=metSyn,phone_head_side=phone_head_side,
#                 age=age,gender=sex,sites=sites,race=ethnicity,
#                 university_education=university_education,deprivation=deprivation,
#                 smoking=smoking,alcohol=alcohol_status,activity=activity)

# data=data.frame(status=data_censoring$label,time=data_censoring$censoring,depression=depression,
#                 metabolic_syndrome=metSyn,
#                 age=age,gender=sex,sites=sites,race=ethnicity,
#                 university_education=university_education,deprivation=deprivation,
#                 smoking=smoking,alcohol=alcohol_status,activity=activity,
#                 sleep=sleep, cancer=cancer,illness=illness,employement_status=employement_status)

rm(age, sex,alcohol_status,deprivation,ethnicity,income,sedentary,depression,activity,
   sites,smoking,university_education,pdff,obesity,phq_4,phq,phone_head_side,
   data_censoring,metSyn,sleep, cancer,illness,employement_status)


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

