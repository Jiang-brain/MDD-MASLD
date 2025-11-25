####################################plot HR for subgroup######################################################################
library(forestplot)
library(dplyr)
library(ggthemes)

mydata<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_subgroup_depression_masld.csv',header=T)
mydata=rbind(mydata,mydata[1,])
mydata$HR[nrow(mydata)]=4
mydata$lower[nrow(mydata)]=4
mydata$upper[nrow(mydata)]=4


mydata$X[nrow(mydata)]='Null'
result1=mydata
result_HR <- tibble(HR  = result1$HR, subgroup = result1$X,lower = result1$lower,upper = result1$upper,N=result1$N,
                    P_value=result1$P,ci=result1$CI,Cases=result1$Cases)
#rm(result)
header <- tibble(subgroup='subgroup')
combined_data <- bind_rows(header,result_HR)
rm(header,result_HR)
combined_data %>% 
  forestplot(labeltext = c(subgroup,ci,N,Cases,P_value),
             mean=HR,
             graph.pos =2,boxsize=0.38,zero=1,
             is.summary = FALSE,clip = c(0.5,4), xlog = FALSE,graphwidth = unit(.38,"npc"),
             col = fpColors(box = "#3c5a66",line = "black",summary = "grey"),
             mar=unit(rep(1.0, times = 4), "cm"),
             txt_gp=fpTxtGp(label=gpar(cex=0.8), ticks=gpar(cex=0.6), xlab=gpar(cex = 0.8), title=gpar(cex = 0.8)),
             xlab = "HR (95% CI)",lwd.xaxis = 1.0)
########################################################################################################################

#####################################plot mediation correlation: phq-depression####################################################
library(ggplot2)
library(RColorBrewer)
require("ggrepel")

phq<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)
masld<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv',header=T)
depression<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_depression.csv',header=T)
fli<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
significance_index=which(phq$P_value<0.05/249 & masld$P.for.linear<0.05/249 & depression$P.for.linear<0.05/249 & fli$P_value<0.05/249
                         & fli$coeff*(depression$HR-1)>0 & phq$coeff*(masld$HR-1)>0)
rm(phq,masld,depression,fli)

mediation_phq_masld=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_phq_masld_62.csv',header=T)[significance_index,]
mediation_depression_masld=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_depression_masld_62.csv',header=T)[significance_index,]

result<- data.frame(phq=mediation_phq_masld$Proportion.mediated_c,depression=mediation_depression_masld$Proportion.mediated_c,
                    p_phq=mediation_phq_masld$P_mediation,p_depression=mediation_depression_masld$P_mediation,phenotype_name=mediation_depression_masld$meta_name,
                    category=1)
rm(mediation_phq_masld,mediation_depression_masld)
non_significant<-which(result$p_phq>0.05/62 | result$p_depression>0.06/62)
result$category[non_significant]=0
result$category=as.factor(result$category)

ggplot(result,aes(x=phq,y=depression))+
  geom_point(size=1.5,shape=19,aes(colour=category))+
  scale_color_manual(values=c("1"='#71a2c2',"0"="grey"))+
  geom_text_repel(
    data = subset(result, p_phq<0.05/62 & p_depression<0.06/62),
    aes(label = phenotype_name),size = 2.5,angle=0,max.overlaps = Inf,
    point.padding = unit(0.3, "lines"),segment.size = 0.2
  )+theme_classic()

#####################################plot mediation correlation: fli-masld####################################################
library(ggplot2)
library(RColorBrewer)
require("ggrepel")

phq<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)
masld<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv',header=T)
depression<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_depression.csv',header=T)
fli<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
significance_index=which(phq$P_value<0.05/249 & masld$P.for.linear<0.05/249 & depression$P.for.linear<0.05/249 & fli$P_value<0.05/249
                         & fli$coeff*(depression$HR-1)>0 & phq$coeff*(masld$HR-1)>0)
rm(phq,masld,depression,fli)

mediation_fli_depression=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_fli_depression_62.csv',header=T)[significance_index,]
mediation_masld_depression=read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_masld_depression_62.csv',header=T)[significance_index,]

result<- data.frame(fli=mediation_fli_depression$Proportion.mediated_c,masld=mediation_masld_depression$Proportion.mediated_c,
                    p_fli=mediation_fli_depression$P_mediation,p_masld=mediation_masld_depression$P_mediation,phenotype_name=mediation_masld_depression$meta_name,
                    category=1)
rm(mediation_fli_depression,mediation_masld_depression)
non_significant<-which(result$p_fli>0.05/62 | result$p_masld>0.06/62)
result$category[non_significant]=0
result$category=as.factor(result$category)

ggplot(result,aes(x=fli,y=masld))+
  geom_point(size=1.5,shape=19,aes(colour=category))+
  scale_color_manual(values=c("1"='#71a2c2',"0"="grey"))+
  geom_text_repel(
    data = subset(result, p_fli<0.05/62 & p_masld<0.06/62),
    aes(label = phenotype_name),size = 2.5,angle=0,max.overlaps = Inf,
    point.padding = unit(0.3, "lines"),segment.size = 0.2
  )+theme_classic()
##################################plot reri, ap, s######################################
library(ggplot2)
#sex
result<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_masld_gender_depression.csv',header=T)
ggplot(result, aes(x=reri))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='lightblue3')+coord_cartesian(xlim=c(-0.01,0.48))+
  geom_vline(xintercept = c(quantile(result$reri,0.025),quantile(result$reri,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=ap*100))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='lightblue3')+coord_cartesian(xlim=c(0,22))+
  geom_vline(xintercept = c(quantile(result$ap*100,0.025),quantile(result$ap*100,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=s))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='lightblue3')+coord_cartesian(xlim=c(0.98,1.75))+
  geom_vline(xintercept = c(quantile(result$s,0.025),quantile(result$s,0.975)),size = 0.2)+theme_classic()
#deprivation
result<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_masld_deprivation_depression.csv',header=T)
ggplot(result, aes(x=reri))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='lightblue3')+coord_cartesian(xlim=c(-0.01,0.62))+
  geom_vline(xintercept = c(quantile(result$reri,0.025),quantile(result$reri,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=ap*100))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='lightblue3')+coord_cartesian(xlim=c(0,30))+
  geom_vline(xintercept = c(quantile(result$ap*100,0.025),quantile(result$ap*100,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=s))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='lightblue3')+coord_cartesian(xlim=c(0.96,2.15))+
  geom_vline(xintercept = c(quantile(result$s,0.025),quantile(result$s,0.975)),size = 0.2)+theme_classic()

#metabolic syndrome
result<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_depression_MetSyn_masld.csv',header=T)
ggplot(result, aes(x=reri))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='darkseagreen')+coord_cartesian(xlim=c(-0.6,3.1))+
  geom_vline(xintercept = c(quantile(result$reri,0.025),quantile(result$reri,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=ap*100))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='darkseagreen')+coord_cartesian(xlim=c(-18,58))+
  geom_vline(xintercept = c(quantile(result$ap*100,0.025),quantile(result$ap*100,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=s))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='darkseagreen')+coord_cartesian(xlim=c(0.75,3.25))+
  geom_vline(xintercept = c(quantile(result$s,0.025),quantile(result$s,0.975)),size = 0.2)+theme_classic()

#obesity
result<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/reri_depression_obesity_masld.csv',header=T)
ggplot(result, aes(x=reri))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='darkseagreen')+coord_cartesian(xlim=c(-0.1,3.0))+
  geom_vline(xintercept = c(quantile(result$reri,0.025),quantile(result$reri,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=ap*100))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='darkseagreen')+coord_cartesian(xlim=c(0,70))+
  geom_vline(xintercept = c(quantile(result$ap*100,0.025),quantile(result$ap*100,0.975)),size = 0.2)+theme_classic()
ggplot(result, aes(x=s))+
  geom_density(alpha = 0.3,size=0.6,linetype=1,fill='darkseagreen')+coord_cartesian(xlim=c(0.96,6))+
  geom_vline(xintercept = c(quantile(result$s,0.025),quantile(result$s,0.975)),size = 0.2)+theme_classic()


#####################################plot phq-metabolites####################################################
library(ggplot2)
library(RColorBrewer)
require("ggrepel")

# mycolors<-(brewer.pal(9, "Blues"))[1:9]
# mycolors<-(brewer.pal(9, "RdYlBu"))[1:9]
#mycolors<-rev((brewer.pal(11, "Spectral"))[1:11])
mycolors<-rev((brewer.pal(11, "RdBu"))[1:11])

mydata   <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)
#mydata   <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/depression_metabolites_irnt.csv',header=T)

result<- data.frame(phenotype_name=mydata$meta_name,effect_size=mydata$coeff,
                    p_value=-log10(mydata$P_value))

ggplot(result,aes(x=effect_size,y=p_value,fill=(effect_size)))+
  geom_point(size=2,shape=21,alpha=1,stroke=0.35)+scale_fill_gradientn(colours=mycolors,limits=c(min(result$effect_size),max(result$effect_size)))+
  coord_cartesian(xlim=c(min(result$effect_size),max(result$effect_size)),ylim=c(0,165))+
  geom_text_repel(
    data = subset(result, abs(effect_size)>=sort(abs(result$effect_size),decreasing = TRUE)[10]),
    aes(label = phenotype_name),size = 2,angle=0,color='black',
    point.padding = 5,segment.size = 0.2, max.overlaps = 20,force=2.5)+
  geom_vline(xintercept = c(0,0),lty = 2,col = "#999999",lwd = 0.45)+
  geom_hline(yintercept = c(-log10(0.05/249),-log10(0.05/249)),lty = 2,col = "#999999",lwd = 0.45)+
  theme_classic()
#####################################plot metabolites-masld####################################################
library(ggplot2)
library(RColorBrewer)
require("ggrepel")

#mycolors<-rev((brewer.pal(11, "Spectral"))[1:11])
mycolors<-rev((brewer.pal(11, "RdBu"))[1:11])

mydata   <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv',header=T)

result<- data.frame(phenotype_name=mydata$meta_name,effect_size=mydata$HR,
                    p_value=-log10(mydata$P.for.linear))

ggplot(result,aes(x=effect_size,y=p_value,fill=(effect_size)))+
  geom_point(size=2,shape=21,alpha=1,stroke=0.35)+scale_fill_gradientn(colours=mycolors,limits=c(min(result$effect_size),max(result$effect_size)))+
  coord_cartesian(xlim=c(min(result$effect_size),max(result$effect_size)),ylim=c(0,20))+
  geom_text_repel(
    data = subset(result, abs(effect_size-1)>=sort(abs(result$effect_size-1),decreasing = TRUE)[10]),
    aes(label = phenotype_name),size = 2,angle=0,color='black',
    point.padding = 5,segment.size = 0.2, max.overlaps = 20,force=2.5)+
  geom_vline(xintercept = c(1,1),lty = 2,col = "#999999",lwd = 0.45)+
  geom_hline(yintercept = c(-log10(0.05/249),-log10(0.05/249)),lty = 2,col = "#999999",lwd = 0.45)+
  theme_classic()

###################################association between effect sizes####################################
result_1  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)
result_2  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_NAFLD.csv',header=T)

significance_group=rep(0,nrow(result_1))
index=which(result_1$P_value<0.05/249&result_2$P.for.linear<0.05/249)
significance_group[index]=1
rm(index)

result=data.frame(coeff=result_1$coeff,hr=result_2$HR,significance_group=significance_group,meta_name=result_1$meta_name)
result$significance_group=as.factor(result$significance_group)

ggplot(result,aes(x=coeff,y=hr))+
  geom_point(size=1.6,shape=19,aes(colour=significance_group))+
  scale_color_manual(values=c("1"='darkseagreen4',"0"="grey"))+theme_classic()

##phq_depression
result_1  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/phq_metabolites_irnt.csv',header=T)
result_2  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/depression_metabolites_irnt.csv',header=T)

significance_group=rep(0,nrow(result_1))
index=which(result_1$P_value<0.05/249 & result_2$P_value<0.05/249)
significance_group[index]=1
rm(index)

result=data.frame(coeff1=result_1$coeff,coeff2=result_2$coeff,significance_group=significance_group,meta_name=result_1$meta_name)
result$significance_group=as.factor(result$significance_group)

ggplot(result,aes(x=coeff2,y=coeff1))+
  geom_point(size=1.6,shape=19,aes(colour=significance_group))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values=c("1"='darkseagreen4',"0"="grey"))+theme_classic()
#########################plot mediated variance#######################################################
library(ggplot2)
mediation  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_phq_masld_62_sorted.csv',header=T)
#mediation  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_depression_masld_62_sorted.csv',header=T)

mediation=mediation[1:62,]
mediation$significance_group=1
index=which(mediation$P_mediation>0.05/62)
mediation$significance_group[index]=0
mediation$significance_group=factor(mediation$significance_group)
behavior_order=order((mediation$Proportion.mediated_c))
ggplot(mediation, aes(x=meta_name, y=Proportion.mediated_c)) +
  geom_pointrange(aes(ymin=lower, ymax=upper),colour='darkseagreen4')+
  geom_point(data=mediation, aes(x=meta_name, y=Proportion.mediated_c,colour=significance_group), size=2.6, shape=19)+
  scale_x_discrete(limits=mediation$meta_name[behavior_order])+
  scale_color_manual(values=c("1"='darkseagreen4',"0"="grey"))+ylim(-1,16)+
  theme_classic()+theme(axis.text.x = element_text(angle = 90,vjust = 0.9,hjust = 0.85))


#coord_flip()+
#########################plot phq distribution####################################
library(ggplot2)
datas=data.frame(phq=datas$phq_4)
ggplot(datas, aes(x=phq)) +
  geom_histogram( fill="darkseagreen", colour='darkseagreen4',alpha=0.9,binwidth = 1)+coord_cartesian(xlim=c(3,18))+theme_classic()
