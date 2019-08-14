###### packages #####
library(ggplot2)
library(brms)
library(nlme)
library(lme4)
library(emmeans)
library(lmerTest)
library(plyr)
##### read data ########
OV.data<-read.csv("OV.data.csv")
gall<-read.csv("gall.summarize.data.csv")
### exclude hybrid data ###
OV.data<-OV.data[OV.data$HP.M==OV.data$HP.F,]
OV.data$Wasp<-OV.data$HP.M
gall<-gall[gall$HP.M==gall$HP.F,]
gall$Wasp<-gall$HP.M
### 1. HR.OV ################
OV.data$status<-"Resident"
OV.data$status[OV.data$Plant!=OV.data$Wasp]<-"Immigrant"
OV.data<-OV.data[-grep("fallen leaf",OV.data$Branch,value=FALSE),]
OV.data$Wasp<-factor(OV.data$Wasp,levels=c("Qv","Qg"))
OV.data$Plant<-factor(OV.data$Plant,levels=c("Qv","Qg"))
### log transform the HR data and do the linear regression ####
HR.OV.m1<-lmer(HR.OV.log~Plant*status+Treatment+(1|M.site)+(1|Tree.ID),data=OV.data)
summary(HR.OV.m1)
lsmeans(HR.OV.m1,list(pairwise~Plant*status))
HR.OV.bay<-brm(No..HR|trials(No..OV.scars)~Plant*status+Treatment+(1|Site)+(1|Tree.ID),data=OV.data,
family=binomial("logit"),control=list(adapt_delta=0.999,max_treedepth=18))
OV.pc1<-hypothesis(HR.OV.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
OV.pc2<-hypothesis(HR.OV.bay, c("interactionQg_Qv=0"),alpha=0.025)

### gall. OV ######
gall$gall.OV.log<-log(gall$Total.galls/gall$Total.OV+0.0001)
gall$status<-"Resident"
gall$status[gall$Plant!=gall$Wasp]<-"Immigrant"
gall$status<-factor(gall$status,levels=c("Resident","Immigrant"))
gall$Plant<-factor(gall$Plant,levels=c("Qv","Qg"))
gall.OV.m1<-lmer(gall.OV.log~Plant*status+Treatment+(1|M.site),data=gall)
summary(gall.OV.m1)
lsmeans(gall.OV.m1,list(pairwise~Plant*status))

gall.OV.bay<-brm(Total.gall|trials(No..OV.scars)~Plant*status+Treatment+(1|Site),data=data.sum,
                 family=binomial("logit"),control=list(adapt_delta=0.9999,max_treedepth=15))
summary(gall.OV.bay)

gallOV.pc1<-hypothesis(gall.OV.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
gallOV.pc2<-hypothesis(gall.OV.bay, c("interactionQg_Qv=interactionQv_Qv"),alpha=0.025)
### gall size ##############



### number of emergence ####
gall$Wasp.gall.log<-log(gall$No.wasps/gall$Total.galls+0.0001)
Wasp.gall.m1<-lmer(Wasp.gall.log~Plant*status+Treatment+(1|M.site),data=gall)
summary(Wasp.gall.m1)
lsmeans(Wasp.gall.m1,list(pairwise~Plant*status))

emerg.gall.bay<-brm(No.wasps|trials(Total.gall)~interaction+Treatment+(1|Site),data=data.sum.new,
                    family=binomial("logit"),control=list(adapt_delta=0.999,max_treedepth=15))
summary(emerg.gall.bay)
emerg.gall.bay.pc1<-hypothesis(emerg.gall.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
emerg.gall.bay.pc2<-hypothesis(emerg.gall.bay, c("interactionQv_Qv=interactionQg_Qv"),alpha=0.025)

### number of emergence per OV ####
emerg.OV.field.bay<-brm(No.wasps|trials(No..OV.scars)~interaction+(1|Site)+(1|Tree.ID),data=data.sum.new[data.sum.new$Treatment=="Field",],
                        family=binomial("logit"),control=list(adapt_delta=0.9999999999,max_treedepth=15))
emergOV.field.pc1<-hypothesis(emerg.OV.field.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
emergOV.field.pc2<-hypothesis(emerg.OV.field.bay, c("interactionQg_Qv=0"),alpha=0.025)
emerg.OV.cage.bay<-brm(No.wasps|trials(No..OV.scars)~interaction+(1|Site)+(1|Tree.ID)+(1|Branch),data=data.sum.new[data.sum.new$Treatment=="Cage",],
                       family=binomial("logit"),control=list(adapt_delta=0.99999,max_treedepth=15))
emergOV.cage.pc1<-hypothesis(emerg.OV.cage.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)

