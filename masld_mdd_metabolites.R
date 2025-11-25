###################################linear association between metabolites and fli#######################################
library(lmerTest)
result_metabolites=data.frame(matrix(nrow=ncol(metabolites),ncol=5))
row.names(result_metabolites)=colnames(metabolites)
colnames(result_metabolites)=c('coeff','se','T_values','P_value','N')
for (i in 1:ncol(metabolites)) {
  datas$metabolite=metabolites[,i]
  number_na<-rowSums(is.na(datas))
  index_nonNA<-which(number_na==0)
  new_data<-datas[index_nonNA,]
  rm(number_na,index_nonNA)
  
  new_data$metabolite=qnorm((rank(new_data$metabolite, ties.method="average") - 0.5) / length(new_data$metabolite))
  model <- lmer(scale(metabolite)~ scale(fli)+age+gender+race+deprivation+university_education+activity+
                  smoking+alcohol+metabolic_syndrome
                +(1|sites),data=new_data)
  coeff<-summary(model)$coefficients
  result_metabolites$coeff[i]=round(coeff[2,1],3)
  result_metabolites$se[i]=round(coeff[2,2],4)
  result_metabolites$T_values[i]=round(coeff['scale(fli)','t value'],2)
  result_metabolites$P_value[i]=signif(coeff[2,5],3)
  result_metabolites$N[i]=nrow(new_data)
  rm(model,coeff,new_data)
}
rm(i)
phenotype<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/datas/Metabolities_name.csv',header=T)[1:249,]
result_metabolites$meta_name=phenotype$Biomarker
result_metabolites$full_name=phenotype$Description
rm(phenotype)
#write.csv(result_metabolites,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv')

###################################masld-metabolome#######################################
library(lmerTest)
result_metabolites=data.frame(matrix(nrow=ncol(metabolites),ncol=5))
row.names(result_metabolites)=colnames(metabolites)
colnames(result_metabolites)=c('coeff','se','T_values','P_value','N')

for (i in 1:ncol(metabolites)) {
  datas$metabolite=metabolites[,i]
  number_na<-rowSums(is.na(datas))
  index_nonNA<-which(number_na==0)
  new_data<-datas[index_nonNA,]
  rm(number_na,index_nonNA)
  
  new_data$metabolite=qnorm((rank(new_data$metabolite, ties.method="average") - 0.5) / length(new_data$metabolite))
  new_data$masld=as.numeric(as.character(new_data$masld))
  model <- lmer(scale(metabolite)~ scale(masld)+age+gender+race+deprivation+university_education+activity+
                  smoking+alcohol+metabolic_syndrome
                +(1|sites),data=new_data)
  coeff<-summary(model)$coefficients
  result_metabolites$coeff[i]=round(coeff[2,1],3)
  result_metabolites$se[i]=round(coeff[2,2],4)
  result_metabolites$T_values[i]=round(coeff['scale(masld)','t value'],2)
  result_metabolites$P_value[i]=signif(coeff[2,5],3)
  result_metabolites$N[i]=nrow(new_data)
  rm(model,coeff,new_data)
}
rm(i)
phenotype<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/datas/Metabolities_name.csv',header=T)[1:249,]
result_metabolites$meta_name=phenotype$Biomarker
result_metabolites$full_name=phenotype$Description
rm(phenotype)
#write.csv(result_metabolites,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/masld_metabolites_irnt.csv')

###########################################metabolome-depression##############################
result_meta_nafld=data.frame(matrix(nrow=ncol(metabolites),ncol=8))
row.names(result_meta_nafld)=colnames(metabolites)
colnames(result_meta_nafld)=c('P_overall','P for nonlinear','P for linear','HR','lower','upper','CI','N')

for (index_metabolite in c(1:ncol(metabolites))) {
  new_data=datas
  new_data$metabolite=metabolites[,index_metabolite]
  number_na<-rowSums(is.na(new_data))
  index_nonNA<-which(number_na==0)
  new_data<-new_data[index_nonNA,]
  rm(number_na,index_nonNA)
  
  new_data$metabolite=qnorm((rank(new_data$metabolite, ties.method="average") - 0.5) / length(new_data$metabolite))
  ##linear model
  fit_linear<- coxph(Surv(time, status) ~ scale(metabolite) +age+gender+race+deprivation+university_education+activity+
                       smoking+alcohol+metabolic_syndrome, data=new_data)
  coeff<-summary(fit_linear)
  result_meta_nafld$HR[index_metabolite]=round(coeff$coefficients['scale(metabolite)','exp(coef)'],2)
  result_meta_nafld$`P for linear`[index_metabolite]=signif(coeff$coefficients['scale(metabolite)','Pr(>|z|)'],3)
  lower=round(coeff$conf.int['scale(metabolite)','lower .95'],2)
  upper=round(coeff$conf.int['scale(metabolite)','upper .95'],2)
  result_meta_nafld$lower[index_metabolite]=lower
  result_meta_nafld$upper[index_metabolite]=upper
  result_meta_nafld$CI[index_metabolite]=paste(result_meta_nafld$HR[index_metabolite],' (',lower,'-',upper,')',sep='')
  result_meta_nafld$N[index_metabolite]=nrow(new_data)
  rm(fit_linear,lower,upper,coeff,new_data)
}
rm(index_metabolite)
phenotype<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/datas/Metabolities_name.csv',header=T)[1:249,]
result_meta_nafld$meta_name=phenotype$Biomarker
result_meta_nafld$full_name=phenotype$Description
rm(phenotype)
#write.csv(result_meta_nafld,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_depression.csv')




###############################association between metabolites maps##################
##masld
masld<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
masld1<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/masld_metabolites_irnt.csv',header=T)

x=masld
y=masld1
cor.test(x$coeff,y$coeff)
ggplot(data=NULL,mapping=aes(x = x$coeff, y = y$coeff)) +
  geom_point(color = "lightblue3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "fli",y = "masld")+theme_bw()
length(which(x$P_value<0.05/249))
length(which(y$P_value<0.05/249))
length(which(x$P_vrm(eeealue<0.05/249 & y$P_value<0.05/249))


##masld-depression
masld<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
masld1<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/masld_metabolites_irnt.csv',header=T)
depression<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_depression.csv',header=T)

x=masld1
y=depression
cor.test(x$coeff,y$HR)
ggplot(data=NULL,mapping=aes(x = x$coeff, y = y$HR)) +
  geom_point(color = "lightblue3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "masld",y = "depression")+theme_bw()
length(which(x$P_value<0.05/249))
length(which(y$P.for.linear<0.05/249))
length(which(x$P_value<0.05/249 & y$P.for.linear<0.05/249))
