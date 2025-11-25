#################################nonlinear#############################################
library(survival)
library(rms)
library(dplyr)
library(Gmisc)
library(splines)
library(Greg)
library(tidyverse)

dd <- datadist(datas) #为后续程序设定数据环境
options(datadist='dd') #为后续程序设定数据环境

fit<- cph(Surv(time, status) ~ rcs(fli,3)+age+gender+race+deprivation+university_education+activity+
            smoking+alcohol+metabolic_syndrome, data = datas)
anova(fit)#test the significance of nolinear regression

#####################plot nonlinear association##########
b=quantile(datas$fli,seq(0,1,0.0005))
datas$quantiled_fli=datas$fli
for (i in c(1:length(b)-1)) {
  index=which(datas$fli>b[i]&datas$fli<=b[i+1])
  datas$quantiled_fli[index]=b[i]
}
rm(i,index,b)
fit<- cph(Surv(time, status) ~ rcs(quantiled_fli,3)+age+gender+race+deprivation+university_education+activity+
            smoking+alcohol+metabolic_syndrome, data = datas)
anova(fit)
plotHR(fit, term="quantiled_fli", xlab="Fatty liver index",se=TRUE,
       polygon_ci=TRUE,
       col.term = "#0070b9",
       col.se = "#DEEBF7BB",
       lwd.term=2,
       lty.term = 1,
       ylog=FALSE,
       alpha=0.05,#alpha level, 95%
       cex=1,
       axes = TRUE,
       plot.bty="l", xlim=c(0,100),ylim=c(0.4,2.0),rug="density")
#Pre0 <-rms::Predict(fit,frailty,fun=exp,type="predictions",ref.zero=TRUE,conf.int = 0.95,digits=2)
