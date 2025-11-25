###########################################metabolome-MASLD##############################
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
#write.csv(result_meta_nafld,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv')


###################################depression-metabolome#######################################
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
  new_data$depression=as.numeric(as.character(new_data$depression))
  model <- lmer(scale(metabolite)~ scale(depression)+age+gender+race+deprivation+university_education+activity+
                  smoking+alcohol+metabolic_syndrome
                +(1|sites),data=new_data)
  coeff<-summary(model)$coefficients
  result_metabolites$coeff[i]=round(coeff[2,1],3)
  result_metabolites$se[i]=round(coeff[2,2],4)
  result_metabolites$T_values[i]=round(coeff['scale(depression)','t value'],2)
  result_metabolites$P_value[i]=signif(coeff[2,5],3)
  result_metabolites$N[i]=nrow(new_data)
  rm(model,coeff,new_data)
}
rm(i)
phenotype<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/datas/Metabolities_name.csv',header=T)[1:249,]
result_metabolites$meta_name=phenotype$Biomarker
result_metabolites$full_name=phenotype$Description
rm(phenotype)
#write.csv(result_metabolites,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/depression_metabolites_irnt.csv')
###################################rds-metabolome#######################################
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
  
  #new_data$metabolite=log(new_data$metabolite+1)
  new_data$metabolite=qnorm((rank(new_data$metabolite, ties.method="average") - 0.5) / length(new_data$metabolite))
  
  model <- lmer(scale(metabolite)~ scale(phq_4)+age+gender+race+deprivation+university_education+activity+
                  smoking+alcohol+metabolic_syndrome
                +(1|sites),data=new_data)
  coeff<-summary(model)$coefficients
  result_metabolites$coeff[i]=round(coeff[2,1],3)
  result_metabolites$se[i]=round(coeff[2,2],4)
  result_metabolites$T_values[i]=round(coeff['scale(phq_4)','t value'],2)
  result_metabolites$P_value[i]=signif(coeff[2,5],3)
  result_metabolites$N[i]=nrow(new_data)
  rm(model,coeff,new_data)
}
rm(i)
phenotype<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/datas/Metabolities_name.csv',header=T)[1:249,]
result_metabolites$meta_name=phenotype$Biomarker
result_metabolites$full_name=phenotype$Description
rm(phenotype)
#write.csv(result_metabolites,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv')
##################################################################################################

###############################association between metabolites maps##################
##depression
depression<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/depression_metabolites_irnt.csv',header=T)
depression1<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)


x=depression
y=depression1
cor.test(x$coeff,y$coeff)
ggplot(data=NULL,mapping=aes(x = x$coeff, y = y$coeff)) +
  geom_point(color = "lightblue3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "depression",y = "phq")+theme_bw()
length(which(x$P_value<0.05/249))
length(which(y$P_value<0.05/249))
length(which(x$P_value<0.05/249 & y$P_value<0.05/249))

##masld
masld<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_original_NAFLD.csv',header=T)
masld1<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv',header=T)

x=masld
y=masld1
cor.test(x$HR,y$HR)
ggplot(data=NULL,mapping=aes(x = x$HR, y = y$HR)) +
  geom_point(color = "lightblue3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "X Variable",y = "Y Variable")+theme_bw()
length(which(x$P.for.linear<0.05/249))
length(which(y$P.for.linear<0.05/249))
length(which(x$P.for.linear<0.05/249 & y$P.for.linear<0.05/249))

##phq_masld
depression<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)
masld<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv',header=T)

x=depression
y=masld
cor.test(x$coeff,y$HR)
ggplot(data=NULL,mapping=aes(x = x$coeff, y = y$HR)) +
  geom_point(color = "lightblue3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "depression",y = "depression")+theme_bw()
length(which(x$P_value<0.05/249))
length(which(y$P.for.linear<0.05/249))
length(which(x$P_value<0.05/249 & y$P.for.linear<0.05/249))

#######################################run mediation for top metabolites#################################
phq<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)
masld<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv',header=T)
depression<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_depression.csv',header=T)
fli<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
significance_index=which(phq$P_value<0.05/249 & masld$P.for.linear<0.05/249 & depression$P.for.linear<0.05/249 & fli$P_value<0.05/249
                         & fli$coeff*(depression$HR-1)>0 & phq$coeff*(masld$HR-1)>0)
rm(phq,masld,depression,fli)

library("survival")
library("survminer")
library("mediation")
result<-data.frame(matrix(nrow=ncol(metabolites),ncol=19))
row.names(result)<-colnames(metabolites)
colnames(result)<-c('a','P_a','se_a','b','P_b','se_b','c','P_c','se_c','c1','P_c1','se_c1','a*b','Proportion mediated','Proportion mediated_c','N','lower','upper','P_mediation')
for (i in c(significance_index)) {
  #retain complete data
  new_data=datas
  new_data$depression=new_data$phq_4
  new_data$metabolite=metabolites[,i]
  number_na<-rowSums(is.na(new_data))
  index_nonNA<-which(number_na==0)
  new_data<-new_data[index_nonNA,]
  rm(number_na,index_nonNA)
  
  new_data$metabolite=qnorm((rank(new_data$metabolite, ties.method="average") - 0.5) / length(new_data$metabolite))
  new_data$depression=scale(new_data$depression)[,1]
  new_data$metabolite=scale(new_data$metabolite)[,1]
  
  model_c <- survreg(Surv(time, status) ~ depression+
                       age+gender+race+deprivation+university_education+activity+
                       smoking+alcohol+metabolic_syndrome,data = new_data)
  model_c1 <- survreg(Surv(time, status) ~ depression+metabolite+
                        age+gender+race+deprivation+university_education+activity+
                        smoking+alcohol+metabolic_syndrome,data = new_data)
  model_a<-lm(metabolite~depression+
                age+gender+race+deprivation+university_education+activity+
                smoking+alcohol+metabolic_syndrome,data = new_data)
  Coeff_a<-model_a$coefficients['depression']
  Coeff_b<-model_c1$coefficients['metabolite']
  Coeff_c<-model_c$coefficients['depression']
  Coeff_c1<-model_c1$coefficients['depression']
  result$P_a[i]=signif(summary(model_a)$coefficients['depression','Pr(>|t|)'],3)
  result$P_c[i]=signif(summary(model_c)$table['depression','p'],3)
  result$P_b[i]=signif(summary(model_c1)$table['metabolite','p'],3)
  result$P_c1[i]=signif(summary(model_c1)$table['depression','p'],3)
  result$se_a[i]=summary(model_a)$coefficients['depression','Std. Error']
  result$se_b[i]=summary(model_c1)$table['metabolite','Std. Error']
  result$se_c[i]=summary(model_c)$table['depression','Std. Error']
  result$se_c1[i]=summary(model_c1)$table['depression','Std. Error']
  
  #mediation <- mediate(model_a, model_c1, sims=200,boot=TRUE, treat="depression", mediator="metabolite",outcome = "time")
  mediation <- mediate(model_a, model_c1, sims=200,boot=FALSE, treat="depression", mediator="metabolite")
  mediation_result<-summary(mediation)
  result$lower[i]=mediation_result$n.avg.ci['2.5%']*100
  result$upper[i]=mediation_result$n.avg.ci['97.5%']*100
  result$P_mediation[i]=mediation_result$n.avg.p
  rm(mediation,mediation_result)
  rm(model_c,model_c1,model_a)
  
  sd1=sqrt(var(new_data$depression)*Coeff_c^2+pi^2/3)
  sd2=sqrt(var(new_data$depression)*Coeff_c1^2+Coeff_b^2*var(new_data$metabolite)+2*Coeff_c1*Coeff_b*cov(new_data$depression,new_data$metabolite)+pi^2/3)
  result$b[i]=round(Coeff_b/sd2,3)
  result$c[i]=round(Coeff_c/sd1,3)
  result$c1[i]=round(Coeff_c1/sd2,3)
  result$a[i]=round(Coeff_a,3)
  result$`a*b`[i]=result$a[i]*result$b[i]
  result$se_a[i]=result$se_a[i]
  result$se_b[i]=result$se_b[i]/sd2
  result$se_c[i]=result$se_c[i]/sd1
  result$se_c1[i]=result$se_c1[i]/sd2
  rm(Coeff_a,Coeff_b,Coeff_c,Coeff_c1,sd1,sd2)
  
  result$`Proportion mediated`[i]=signif(result$`a*b`[i]/(result$`a*b`[i]+result$c1[i]),3)*100
  result$`Proportion mediated_c`[i]=signif(result$`a*b`[i]/result$c[i],3)*100
  result$N[i]=nrow(new_data)
  print(c(i,result$`Proportion mediated`[i],result$lower[i],result$upper[i]))
  rm(new_data)
}
phenotype<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2024-frailty_Nafld/datas/Metabolities_name.csv',header=T)[1:249,]
result$meta_name=phenotype$Biomarker
result$full_name=phenotype$Description
rm(phenotype)
write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_phq_masld_62.csv')

top_mediated_index=order(abs(result$`Proportion mediated`),decreasing = TRUE)
results=result[top_mediated_index,]

