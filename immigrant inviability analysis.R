###### packages #####
library(ggplot2)
library(brms)
library(nlme)
library(lme4)
#######################
OV.data<-read.csv("OV scars GFE.csv")
gall.size<-read.csv("Gall Size.csv")
small.gall<-read.csv("Gall size smaller than 2mm.csv")
### make the column and row names consistent across 3 datasets####
head(OV.data)
head(gall.size)
head(small.gall)
colnames(OV.data)[5]<-"Tree.ID"
colnames(gall.size)[5]<-"Tree.ID"
colnames(small.gall)[5]<-"Tree.ID"
OV.data$Tree.ID<-paste(OV.data$Plant,OV.data$Tree.ID,sep="_")
OV.data$Branch<-paste(OV.data$Tree.ID,OV.data$Branch,sep="_")
gall.size$Treatment<-as.character(gall.size$Treatment)
gall.size$Treatment[gall.size$Treatment=="cage"]<-"Cage"
gall.size$Treatment[gall.size$Treatment=="field"]<-"Field"
gall.size$Tree.ID<-paste(gall.size$Plant,gall.size$Tree.ID,sep="_")
gall.size$Branch<-paste(gall.size$Tree.ID,gall.size$Branch,sep="_")
small.gall$Treatment<-as.character(small.gall$Treatment)
small.gall$Treatment[small.gall$Treatment=="cage"]<-"Cage"
small.gall$Treatment[small.gall$Treatment=="field"]<-"Field"
small.gall$Tree.ID<-paste(small.gall$Plant,small.gall$Tree.ID,sep="_")
small.gall$Branch<-paste(small.gall$Tree.ID,small.gall$Branch,sep="_")
OV.data$No..OV.scars[OV.data$Enter.by=="rg"]<-OV.data$No..OV.scars[OV.data$Enter.by=="rg"]+OV.data$No..HR[OV.data$Enter.by=="rg"]
#### gall emergence data ####
gall.emerg<-read.csv("2017 GFE Bt Emergence.csv")
head(gall.emerg)
colnames(gall.emerg)[11]<-"Tree.ID"
colnames(gall.emerg)[10]<-"Plant"
gall.emerg$No.wasps<-gall.emerg$Number.of.Wasps..Alive.+gall.emerg$Number.of.Wasps..dead.
gall.emerg<-gall.emerg[gall.emerg$Treatment %in% c("Field","Cage"),]
gall.emerg<-gall.emerg[gall.emerg$Site %in% c("ABS","LL","Dick","Alva","KRE","Okee"),]
gall.emerg<-gall.emerg[gall.emerg$Gall.species=="B. treatae",]
gall.emerg<-gall.emerg[gall.emerg$Plant %in% c("Qv","Qg"),]
gall.emerg<-gall.emerg[gall.emerg$Tree.ID!="",]
gall.emerg$Gall.species<-factor(gall.emerg$Gall.species)
gall.emerg$Treatment<-factor(gall.emerg$Treatment)
gall.emerg$Site<-factor(gall.emerg$Site)
gall.emerg$Plant<-factor(gall.emerg$Plant,levels=c("Qv","Qg"))
gall.emerg$Tree.ID<-factor(gall.emerg$Tree.ID)
gall.emerg$Wasp[gall.emerg$Site %in% c("ABS","LL","Dick")]<-"Qg"
gall.emerg$Wasp[gall.emerg$Site %in% c("Alva","KRE","Okee")]<-"Qv"
gall.emerg$Branch.Number[gall.emerg$Treatment=="Field"]<-NA
gall.emerg$Tree.ID<-paste(gall.emerg$Plant,gall.emerg$Tree.ID,sep="_")
gall.emerg$Branch<-paste(gall.emerg$Tree.ID,gall.emerg$Branch.Number,sep="_")
emerg.sum<-aggregate(No.wasps~Site+Treatment+Plant+Tree.ID+Branch,sum,data=gall.emerg)
### summarize the number of OV and galls on each branch/tree ID ###
## OV sum ##
OV.sum<-aggregate(No..OV.scars~Site+Treatment+Plant+Tree.ID+Branch,sum,data=OV.data)
## gall sum ##
# galls from leaves ##
gall.leave<-aggregate(cbind(No..OV.scars,No..HR,No..galls...2.mm.)~Site+Treatment+Plant+Tree.ID+Branch,sum,data=OV.data)
# galls from gall data #
gall.size$count<-1
gall.sum<-aggregate(count~Site+Treatment+Plant+Tree.ID+Branch,sum,data=gall.size)
## remove the natual collection galls from the small gall data ###
small.gall<-small.gall[small.gall$Treatment!="Nat Coll",]
## merge four gall number data ###
data.sum<-merge(gall.leave,gall.sum,by=c("Site","Treatment","Plant","Tree.ID","Branch"),all.x=TRUE)
data.sum<-merge(data.sum,small.gall,by=c("Site","Treatment","Plant","Tree.ID","Branch"),all.x=TRUE)
data.sum<-data.sum[,-c(10,11,13)]
data.sum$Total.gall<-rowSums(cbind(data.sum$No..galls...2.mm.,data.sum$count,data.sum$Number.of.Galls..2mm),na.rm=TRUE)
data.sum<-data.sum[,-c(8,9,10)]
data.sum<-data.sum[-121,]
data.sum<-merge(data.sum,emerg.sum,by=c("Site","Treatment","Plant","Tree.ID","Branch"),all.x=TRUE)
data.sum$No.wasps[is.na(data.sum$No.wasps)]<-0
###################################################################################################
### 1. HR.OV ################
OV.data<-OV.data[OV.data$No..OV.scars!=0,]
OV.data$Wasp[OV.data$Site %in% c("ABS","Dick","LL")]<-"Qg"
OV.data$Wasp[OV.data$Site %in% c("Alva","KRE","Okee")]<-"Qv"
OV.data$interaction<-paste(OV.data$Wasp,OV.data$Plant,sep="_")
OV.data$native.not<-"native"
OV.data$native.not[OV.data$Plant!=OV.data$Wasp]<-"non-native"
OV.data<-OV.data[-grep("fallen leaf",OV.data$Branch,value=FALSE),]
OV.data$Wasp<-factor(OV.data$Wasp,levels=c("Qv","Qg"))
OV.data$Plant<-factor(OV.data$Plant,levels=c("Qv","Qg"))
### log transform the HR data and do the linear regression ####
OV.data$HR.OV.log<-log((OV.data$No..HR+0.0001)/OV.data$No..OV.scars)
OV.data$Tree.ID<-as.factor(OV.data$Tree.ID)

HR.OV.m0<-lm(HR.OV.log~Wasp*Plant+Treatment,data=OV.data)
summary(HR.OV.m0)

HR.OV.m01<-lm(HR.OV.log~interaction+Treatment,data=OV.data)
summary(HR.OV.m01)

HR.OV.m1<-lme(HR.OV.log~Wasp*Plant+Treatment,method="ML",random=~1+Plant|Tree.ID,
              data=OV.data)
summary(HR.OV.m1)

HR.OV.m11<-lme(HR.OV.log~interaction+Treatment,method="ML",random=~1+Plant|Tree.ID,
              data=OV.data)
summary(HR.OV.m11)

HR.OV.m2<-lme(HR.OV.log~Wasp*native.not+Treatment,method="ML",random=~1+Plant|Tree.ID,
              data=OV.data)
summary(HR.OV.m2)

HR.OV.m3<-lme(HR.OV.log~interaction+Treatment,method="ML",random=~1+Plant|Tree.ID,
              data=OV.data)
summary(HR.OV.m3)
###############################
HR.OV.bay<-brm(No..HR|trials(No..OV.scars)~interaction+Treatment+(1|Site)+(1|Tree.ID)+(1|Branch),data=OV.data,
               family=binomial("logit"),control=list(adapt_delta=0.999,max_treedepth=15))
OV.pc1<-hypothesis(HR.OV.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
OV.pc2<-hypothesis(HR.OV.bay, c("interactionQg_Qv=0"),alpha=0.025)
HR.OV.field.bay<-brm(No..HR|trials(No..OV.scars)~interaction+(1|Site)+(1|Tree.ID),data=OV.data[OV.data$Treatment=="Field",],
                     family=binomial("logit"),control=list(adapt_delta=0.99999,max_treedepth=15))
OV.field.pc1<-hypothesis(HR.OV.field.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
HR.OV.cage.bay<-brm(No..HR|trials(No..OV.scars)~interaction+Treatment+(1|Site)+(1|Tree.ID)+(1|Branch),data=OV.data[OV.data$Treatment=="Cage",],
                    family=binomial("logit"),control=list(adapt_delta=0.999,max_treedepth=15))
OV.cage.pc1<-hypothesis(HR.OV.cage.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
### gall. OV ######
data.sum$Wasp[data.sum$Site %in% c("ABS","Dick","LL")]<-"Qg"
data.sum$Wasp[data.sum$Site %in% c("Alva","KRE","Okee")]<-"Qv"
data.sum$interaction<-paste(data.sum$Wasp,data.sum$Plant,sep="_")
error.list<-data.sum[data.sum$No..OV.scars<data.sum$Total.gall,]
write.csv(error.list,"error.list.csv")

data.sum<-data.sum[!(data.sum$No..OV.scars<data.sum$Total.gall),]
data.sum<-data.sum[-grep("fallen leaf",data.sum$Branch,value=FALSE),]
gall.ov.m1<-lm()
gall.OV.bay<-brm(Total.gall|trials(No..OV.scars)~interaction+Treatment+(1|Site)+(1|Tree.ID),data=data.sum,
                 family=binomial("logit"),control=list(adapt_delta=0.999,max_treedepth=15))
gallOV.pc1<-hypothesis(gall.OV.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
gall.OV.bay.field<-brm(Total.gall|trials(No..OV.scars)~interaction+(1|Site)+(1|Tree.ID),data=data.sum[data.sum$Treatment=="Field",],
                       family=binomial("logit"),control=list(adapt_delta=0.9999,max_treedepth=15))
gallOV.field.pc1<-hypothesis(gall.OV.bay.field, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
gall.OV.bay.cage<-brm(Total.gall|trials(No..OV.scars)~interaction+(1|Site)+(1|Tree.ID)+(1|Branch),data=data.sum[data.sum$Treatment=="Cage",],
                      family=binomial("logit"),control=list(adapt_delta=0.9999,max_treedepth=15))
gallOV.cage.pc1<-hypothesis(gall.OV.bay.cage, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
### number of emergence ####
data.sum.new<-data.sum[data.sum$Total.gall>=data.sum$No.wasps,]
data.sum.new<-data.sum.new[data.sum.new$Total.gall!=0,]
emerg.gall.bay<-brm(No.wasps|trials(Total.gall)~interaction+Treatment+(1|Site)+(1|Tree.ID)+(1|Branch),data=data.sum.new,
                    family=binomial("logit"),control=list(adapt_delta=0.999,max_treedepth=15))
emerg.gall.field.bay<-brm(No.wasps|trials(Total.gall)~interaction+(1|Site)+(1|Tree.ID),data=data.sum.new[data.sum.new$Treatment=="Field",],
                          family=binomial("logit"),control=list(adapt_delta=0.999,max_treedepth=15))
emerg.gall.field.pc1<-hypothesis(emerg.gall.field.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
emerg.gall.cage.bay<-brm(No.wasps|trials(Total.gall)~interaction+(1|Site)+(1|Tree.ID),data=data.sum.new[data.sum.new$Treatment=="Cage",],
                         family=binomial("logit"),control=list(adapt_delta=0.9999,max_treedepth=15))
emerg.gall.cage.pc1<-hypothesis(emerg.gall.cage.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
### number of emergence per OV ####
emerg.OV.field.bay<-brm(No.wasps|trials(No..OV.scars)~interaction+(1|Site)+(1|Tree.ID),data=data.sum.new[data.sum.new$Treatment=="Field",],
                        family=binomial("logit"),control=list(adapt_delta=0.9999999999,max_treedepth=15))
emergOV.field.pc1<-hypothesis(emerg.OV.field.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)
emergOV.field.pc2<-hypothesis(emerg.OV.field.bay, c("interactionQg_Qv=0"),alpha=0.025)
emerg.OV.cage.bay<-brm(No.wasps|trials(No..OV.scars)~interaction+(1|Site)+(1|Tree.ID)+(1|Branch),data=data.sum.new[data.sum.new$Treatment=="Cage",],
                       family=binomial("logit"),control=list(adapt_delta=0.99999,max_treedepth=15))
emergOV.cage.pc1<-hypothesis(emerg.OV.cage.bay, c("interactionQv_Qv=interactionQv_Qg"),alpha=0.025)

