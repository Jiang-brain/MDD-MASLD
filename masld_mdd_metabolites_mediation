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
  new_data$masld=new_data$fli
  new_data$metabolite=metabolites[,i]
  number_na<-rowSums(is.na(new_data))
  index_nonNA<-which(number_na==0)
  new_data<-new_data[index_nonNA,]
  rm(number_na,index_nonNA)
  
  new_data$metabolite=qnorm((rank(new_data$metabolite, ties.method="average") - 0.5) / length(new_data$metabolite))
  new_data$masld=scale(new_data$masld)[,1]
  new_data$metabolite=scale(new_data$metabolite)[,1]
  
  model_c <- survreg(Surv(time, status) ~ masld+
                       age+gender+race+deprivation+university_education+activity+
                       smoking+alcohol+metabolic_syndrome,data = new_data)
  model_c1 <- survreg(Surv(time, status) ~ masld+metabolite+
                        age+gender+race+deprivation+university_education+activity+
                        smoking+alcohol+metabolic_syndrome,data = new_data)
  model_a<-lm(metabolite~masld+
                age+gender+race+deprivation+university_education+activity+
                smoking+alcohol+metabolic_syndrome,data = new_data)
  Coeff_a<-model_a$coefficients['masld']
  Coeff_b<-model_c1$coefficients['metabolite']
  Coeff_c<-model_c$coefficients['masld']
  Coeff_c1<-model_c1$coefficients['masld']
  result$P_a[i]=signif(summary(model_a)$coefficients['masld','Pr(>|t|)'],3)
  result$P_c[i]=signif(summary(model_c)$table['masld','p'],3)
  result$P_b[i]=signif(summary(model_c1)$table['metabolite','p'],3)
  result$P_c1[i]=signif(summary(model_c1)$table['masld','p'],3)
  result$se_a[i]=summary(model_a)$coefficients['masld','Std. Error']
  result$se_b[i]=summary(model_c1)$table['metabolite','Std. Error']
  result$se_c[i]=summary(model_c)$table['masld','Std. Error']
  result$se_c1[i]=summary(model_c1)$table['masld','Std. Error']
  
  mediation <- mediate(model_a, model_c1, sims=5000,boot=TRUE, treat="masld", mediator="metabolite",outcome = "time")

  mediation_result<-summary(mediation)
  result$lower[i]=mediation_result$n.avg.ci['2.5%']*100
  result$upper[i]=mediation_result$n.avg.ci['97.5%']*100
  result$P_mediation[i]=mediation_result$n.avg.p
  rm(mediation,mediation_result)
  rm(model_c,model_c1,model_a)
  
  sd1=sqrt(var(new_data$masld)*Coeff_c^2+pi^2/3)
  sd2=sqrt(var(new_data$masld)*Coeff_c1^2+Coeff_b^2*var(new_data$metabolite)+2*Coeff_c1*Coeff_b*cov(new_data$masld,new_data$metabolite)+pi^2/3)
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
write.csv(result,'/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_fli_depression_62.csv')

top_mediated_index=order(abs(result$`Proportion mediated`),decreasing = TRUE)
results=result[top_mediated_index,]

