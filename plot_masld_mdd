####################################plot HR for subgroup######################################################################
library(forestplot)
library(dplyr)
library(ggthemes)

mydata<- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/HR_subgroup_masld_depression.csv',header=T)
mydata=rbind(mydata,mydata[1,])
mydata$HR[nrow(mydata)]=2
mydata$lower[nrow(mydata)]=2
mydata$upper[nrow(mydata)]=2


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
             is.summary = FALSE,clip = c(0.,2), xlog = FALSE,graphwidth = unit(.38,"npc"),
             col = fpColors(box = "#3c5a66",line = "black",summary = "grey"),
             mar=unit(rep(1.0, times = 4), "cm"),
             txt_gp=fpTxtGp(label=gpar(cex=0.8), ticks=gpar(cex=0.6), xlab=gpar(cex = 0.8), title=gpar(cex = 0.8)),
             xlab = "HR (95% CI)",lwd.xaxis = 1.0)
########################################################################################################################

#####################################plot fli-metabolites####################################################
library(ggplot2)
library(RColorBrewer)
require("ggrepel")

# mycolors<-(brewer.pal(9, "Blues"))[1:9]
# mycolors<-(brewer.pal(9, "RdYlBu"))[1:9]
#mycolors<-rev((brewer.pal(11, "Spectral"))[1:11])
mycolors<-rev((brewer.pal(11, "RdBu"))[1:11])

mydata   <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
#mydata   <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/masld_metabolites_irnt.csv',header=T)

index_0=which(mydata$P_value==0)
mydata$P_value[index_0]=5e-324

result<- data.frame(phenotype_name=mydata$meta_name,effect_size=mydata$coeff,
                    p_value=-log10(mydata$P_value))

ggplot(result,aes(x=effect_size,y=p_value,fill=(effect_size)))+
  geom_point(size=2,shape=21,alpha=1,stroke=0.35)+scale_fill_gradientn(colours=mycolors,limits=c(min(result$effect_size),max(result$effect_size)))+
  coord_cartesian(xlim=c(min(result$effect_size),max(result$effect_size)),ylim=c(0,360))+
  geom_text_repel(
    data = subset(result, abs(effect_size)>=sort(abs(result$effect_size),decreasing = TRUE)[10]),
    aes(label = phenotype_name),size = 2,angle=0,color='black',
    point.padding = 5,segment.size = 0.2, max.overlaps = 20,force=2.5)+
  geom_vline(xintercept = c(0,0),lty = 2,col = "#999999",lwd = 0.45)+
  geom_hline(yintercept = c(-log10(0.05/249),-log10(0.05/249)),lty = 2,col = "#999999",lwd = 0.45)+
  theme_classic()
#####################################plot metabolites-depression####################################################
library(ggplot2)
library(RColorBrewer)
require("ggrepel")

#mycolors<-rev((brewer.pal(11, "Spectral"))[1:11])
mycolors<-rev((brewer.pal(11, "RdBu"))[1:11])

mydata   <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_depression.csv',header=T)

result<- data.frame(phenotype_name=mydata$meta_name,effect_size=mydata$HR,
                    p_value=-log10(mydata$P.for.linear))

ggplot(result,aes(x=effect_size,y=p_value,fill=(effect_size)))+
  geom_point(size=2,shape=21,alpha=1,stroke=0.35)+scale_fill_gradientn(colours=mycolors,limits=c(min(result$effect_size),max(result$effect_size)))+
  coord_cartesian(xlim=c(min(result$effect_size),max(result$effect_size)),ylim=c(0,23))+
  geom_text_repel(
    data = subset(result, abs(effect_size-1)>=sort(abs(result$effect_size-1),decreasing = TRUE)[10]),
    aes(label = phenotype_name),size = 2,angle=0,color='black',
    point.padding = 5,segment.size = 0.2, max.overlaps = 20,force=2.5)+
  geom_vline(xintercept = c(1,1),lty = 2,col = "#999999",lwd = 0.45)+
  geom_hline(yintercept = c(-log10(0.05/249),-log10(0.05/249)),lty = 2,col = "#999999",lwd = 0.45)+
  theme_classic()

###################################association between effect sizes####################################
result_1  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
result_2  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/metabolites_irnt_depression.csv',header=T)

significance_group=rep(0,nrow(result_1))
index=which(result_1$P_value<0.05/249&result_2$P.for.linear<0.05/249)
significance_group[index]=1
rm(index)

result=data.frame(coeff=result_1$coeff,hr=result_2$HR,significance_group=significance_group,meta_name=result_1$meta_name)
result$significance_group=as.factor(result$significance_group)

ggplot(result,aes(x=coeff,y=hr))+
  geom_point(size=1.6,shape=19,aes(colour=significance_group))+
  scale_color_manual(values=c("1"='lightblue4',"0"="grey"))+theme_classic()


##fli_masld
result_1  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/fli_metabolites_irnt.csv',header=T)
result_2  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/masld_metabolites_irnt.csv',header=T)

significance_group=rep(0,nrow(result_1))
index=which(result_1$P_value<0.05/249 & result_2$P_value<0.05/249)
significance_group[index]=1
rm(index)

result=data.frame(coeff1=result_1$coeff,coeff2=result_2$coeff,significance_group=significance_group,meta_name=result_1$meta_name)
result$significance_group=as.factor(result$significance_group)

ggplot(result,aes(x=coeff2,y=coeff1))+
  geom_point(size=1.6,shape=19,aes(colour=significance_group))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values=c("1"='lightblue4',"0"="grey"))+theme_classic()
#########################plot mediated variance#######################################################
library(ggplot2)
mediation  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_fli_depression_62_sorted.csv',header=T)
#mediation  <- read.csv('/Users/andyjiang/OneDrive/rtjiang/G/working_data/2025-UKB-depression/datas/mediation_masld_depression_62_sorted.csv',header=T)

mediation=mediation[1:62,]
mediation$significance_group=1
index=which(mediation$P_mediation>0.05/62)
mediation$significance_group[index]=0
mediation$significance_group=factor(mediation$significance_group)
behavior_order=order((mediation$Proportion.mediated_c))
ggplot(mediation, aes(x=meta_name, y=Proportion.mediated_c)) +
  geom_pointrange(aes(ymin=lower, ymax=upper,colour=significance_group))+
  geom_point(data=mediation, aes(x=meta_name, y=Proportion.mediated_c,colour=significance_group), size=2.6, shape=19)+
  scale_x_discrete(limits=mediation$meta_name[behavior_order])+
  scale_color_manual(values=c("1"='lightblue4',"0"="grey"))+
  theme_classic()+theme(axis.text.x = element_text(angle = 90,vjust = 0.9,hjust = 0.85))
##################################plot fli distribution###################################
library(ggplot2)
datas=data.frame(phq=datas$fli)
ggplot(datas, aes(x=phq)) +
  geom_histogram( fill="lightblue3", colour='lightblue4',alpha=0.9)+coord_cartesian(xlim=c(0,100))+theme_classic()
