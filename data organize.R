
################
## leaf area regression formula ###
QV_leaf<-read.csv("Leaf size_PPS-QV.csv")
QG_leaf<-read.csv("Leaf size_PPS-QG.csv")
QV_la<-lm(Leaf.area~Leaf.length*Leaf.width,data=QV_leaf)
QG_la<-lm(Leaf.area~Leaf.length*Leaf.width,data=QG_leaf)
#### read hybrid data ####
Ov.Hyb<-read.csv("17 SP OV hybrid.csv")
Ov.Hyb$No..OV.scars[Ov.Hyb$Entered.by=="IC"]<-Ov.Hyb$No..OV.scars[Ov.Hyb$Entered.by=="IC"]+
  Ov.Hyb$No..HR[Ov.Hyb$Entered.by=="IC"]
head(Ov.Hyb)
unique(Ov.Hyb$Cross..M.)
Ov.Hyb<-Ov.Hyb[Ov.Hyb$Cross..M.!="unknown",]
Ov.Hyb[Ov.Hyb$Cross..M.=="n/a",]
Ov.Hyb<-Ov.Hyb[Ov.Hyb$Cross..M.!="n/a",]
####################################
unique(Ov.Hyb$Cross..M.)
unique(Ov.Hyb$Cross..F.)
Ov.Hyb$Cross..M.<-factor(Ov.Hyb$Cross..M.)
Ov.Hyb$Cross..F.<-factor(Ov.Hyb$Cross..F.)
Ov.Hyb$HP.M<-"Qv"
Ov.Hyb$HP.M[Ov.Hyb$Cross..M.=="ABS"|Ov.Hyb$Cross..M.=="Dick"|Ov.Hyb$Cross..M.=="LL"]<-"Qg"
Ov.Hyb$HP.F<-"Qv"
Ov.Hyb$HP.F[Ov.Hyb$Cross..F.=="ABS"|Ov.Hyb$Cross..F.=="Dick"|Ov.Hyb$Cross..F.=="LL"]<-"Qg"
Ov.Hyb$Tree.ID<-paste(Ov.Hyb$Plant,Ov.Hyb$Tree.No.,sep="_")
Ov.Hyb$Branch<-paste(Ov.Hyb$Tree.ID,Ov.Hyb$Branch,sep="_")
Ov.Hyb$Leaf.ID<-paste(Ov.Hyb$Branch,Ov.Hyb$Leaf.ID,sep="_")
colnames(Ov.Hyb)[12]<-"Leaf.length"
colnames(Ov.Hyb)[13]<-"Leaf.width"
Ov.Hyb<-Ov.Hyb[Ov.Hyb$No..OV.scars!=0,]
Ov.Hyb$Leaf.area[Ov.Hyb$Plant=="Qv"]<-predict(QV_la,newdata=Ov.Hyb[Ov.Hyb$Plant=="Qv",])
Ov.Hyb$Leaf.area[Ov.Hyb$Plant=="Qg"]<-predict(QG_la,newdata=Ov.Hyb[Ov.Hyb$Plant=="Qg",])
OV.Hyb1<-Ov.Hyb[,c(5,6,16,17,3,18,7,8,9,10,11,19)]
colnames(OV.Hyb1)[c(1,2,9,10,11)]<-c("M.site","F.site","No.OV","No.HR","No.galls")
OV.Hyb1$Treatment<-"Cage"
OV.Hyb1<-OV.Hyb1[,c(13,1:12)]

#### read OV scar data ###
OV<-read.csv("OV scars GFE.csv")
head(OV)
OV$Tree.ID<-paste(OV$Plant,OV$Tree..No,sep="_")
OV$Branch<-paste(OV$Tree.ID,OV$Branch,sep="_")
OV$Leaf.ID<-paste(OV$Branch,OV$Leaf.ID,sep="_")
OV$No..OV.scars[OV$Enter.by=="rg"]<-OV$No..OV.scars[OV$Enter.by=="rg"]+OV$No..HR[OV$Enter.by=="rg"]
OV<-OV[OV$No..OV.scars!=0,]
OV$HP.M<-"Qv"
OV$HP.M[OV$Site=="ABS"|OV$Site=="Dick"|OV$Site=="LL"]<-"Qg"
OV$HP.F<-OV$HP.M
OV$M.site<-OV$Site
OV$F.site<-OV$Site
colnames(OV)[11]<-"Leaf.length"
colnames(OV)[12]<-"Leaf.width"
OV$Leaf.area[OV$Plant=="Qv"]<-predict(QV_la,newdata=OV[OV$Plant=="Qv",])
OV$Leaf.area[OV$Plant=="Qg"]<-predict(QG_la,newdata=OV[OV$Plant=="Qg",])
###### rename the colnums of OV data 
OV<-OV[,c(3,20,21,18,19,4,17,6,7,8,9,10,15)]
colnames(OV)[10]<-"No.OV"
colnames(OV)[11]<-"No.HR"
colnames(OV)[12]<-"No.galls"
colnames(OV)[13]<-"Leaf.area"

### merge the OV scar dataset ###
OV.mix<-rbind(OV,OV.Hyb1)
OV.mix$Plant<-factor(OV.mix$Plant,levels=c("Qv","Qg"))
OV.mix$HP.M<-factor(OV.mix$HP.M,levels=c("Qv","Qg"))
OV.mix$HP.F<-factor(OV.mix$HP.F,levels=c("Qv","Qg"))
OV.mix$cross<-paste(OV.mix$HP.M,OV.mix$HP.F,sep=" x ")
OV.mix$cross<-factor(OV.mix$cross,levels=c("Qv x Qv","Qv x Qg","Qg x Qv","Qg x Qg"))
OV.mix$HR.OV<-OV.mix$No.HR/OV.mix$No.OV
OV.mix$HR.OV.log<-log(OV.mix$HR.OV+0.0001)
OV.mix$count<-1

## read gall data
gall.size<-read.csv("Gall Size.csv")
### make the column and row names consistent across 3 datasets####
gall.size$Tree.ID<-paste(gall.size$Plant,gall.size$Tree.Number,sep="_")
gall.size$Branch<-paste(gall.size$Tree.ID,gall.size$Branch,sep="_")
gall.size$HP.M<-"Qv"
gall.size$HP.M[gall.size$Site=="ABS"|gall.size$Site=="Dick"|gall.size$Site=="LL"]<-"Qg"
gall.size$HP.F<-gall.size$HP.M
gall.size$M.site<-gall.size$Site
gall.size$F.site<-gall.size$Site
colnames(gall.size)
gall.size<-gall.size[,c(3,21,22,19,20,4,18,6,9,10)]
gall.size$Treatment<-as.character(gall.size$Treatment)
gall.size$Treatment[gall.size$Treatment=="cage"]<-"Cage"
gall.size$Treatment[gall.size$Treatment=="field"]<-"Field"
### read hybrid gall size ###
gall.Hyb<-read.csv("Hybrid galls.csv")
gall.Hyb$count<-1
gall.Hyb$Tree.ID<-paste(gall.Hyb$Host.Plant,gall.Hyb$Host.Plant.ID..,sep="_")
gall.Hyb$Branch<-paste(gall.Hyb$Host.Plant,gall.Hyb$Host.Plant.ID..,gall.Hyb$Branch..,sep="_")
gall.Hyb$Treatment<-"Cage"
gall.Hyb<-gall.Hyb[,c(16,4,5,1,2,7,14,15,11)]
colnames(gall.Hyb)[c(2,3,4,5,6,9)]<-c("M.site","F.site","HP.M","HP.F","Plant","Gall.Size")
gall.Hyb$Hole.Present<-NA
#### merge two gall dataset ###
gall.size<-rbind(gall.size,gall.Hyb)
gall.size$count<-1
gall.size$count[gall.size$Gall.Size<2]<-0
gall.size$count.small[gall.size$count==1]<-0
gall.size$count.small[gall.size$count==0]<-1
gall.size$count2.62[gall.size$Gall.Size>=2.62]<-1
gall.size$count2.62[gall.size$Gall.Size<2.62]<-0
gall.size$count5.82[gall.size$Gall.Size>=5.82]<-1
gall.size$count5.82[gall.size$Gall.Size<5.82]<-0
head(gall.size)

## organize small.galls ##
small.gall<-read.csv("Gall size smaller than 2mm.csv")
head(small.gall)
small.gall<-small.gall[!is.na(small.gall$Tree.Number),]
small.gall$Tree.ID<-paste(small.gall$Plant,small.gall$Tree.Number,sep="_")
small.gall$Branch<-paste(small.gall$Tree.ID,small.gall$Branch,sep="_")

small.gall$Treatment<-as.character(small.gall$Treatment)
small.gall$Treatment[small.gall$Treatment=="cage"]<-"Cage"
small.gall$Treatment[small.gall$Treatment=="field"]<-"Field"

small.gall$HP.M<-"Qv"
small.gall$HP.M[small.gall$Site=="ABS"|small.gall$Site=="Dick"|small.gall$Site=="LL"]<-"Qg"
small.gall$HP.F<-small.gall$HP.M
small.gall$M.site<-small.gall$Site
small.gall$F.site<-small.gall$Site
head(small.gall)
small.gall<-small.gall[,c("M.site","F.site","HP.M","HP.F","Tree.ID","Branch","Number.of.Galls..2mm")]
colnames(small.gall)[7]<-"No.small.galls"
#######################################

### summarize the total number of OV scars per branch ###
OV.branch<-ddply(OV.mix,c("Treatment","M.site","F.site","HP.M","HP.F","Plant","Tree.ID","Branch"),
                 summarize,Total.OV=sum(No.OV),Total.HR=sum(No.HR),Total.small.galls=sum(No.galls))
head(OV.branch)
head(gall.size)
### summarize the total number of galls per branch ###
gall.branch<-ddply(gall.size,c("Treatment","M.site","F.site","HP.M","HP.F","Plant","Tree.ID","Branch"),
                   summarize,Total.galls=sum(count),Total.small.galls=sum(count.small),avg.gall.size=mean(Gall.Size),gall.2.62=sum(count2.62),gall.5.82=sum(count5.82))
head(gall.branch)
## merge galls on OV and gall data ###
gall.temp<-merge(OV.branch,gall.branch[,c(8,9,11,12,13,10)],by="Branch",all.x = TRUE)
head(gall.temp)
head(small.gall)
gall<-merge(gall.temp,small.gall[,c(6,7)],by="Branch",all.x = TRUE)
head(gall)
gall$Total.galls[is.na(gall$Total.galls)]<-0
gall$Total.small.galls.x[is.na(gall$Total.small.galls.x)]<-0
gall$Total.small.galls.y[is.na(gall$Total.small.galls.y)]<-0
gall$gall.2.62[is.na(gall$gall.2.62)]<-0
gall$gall.5.82[is.na(gall$gall.5.82)]<-0
gall$No.small.galls[is.na(gall$No.small.galls)]<-0
gall$Total.small.galls<-gall$Total.small.galls.x+gall$Total.small.galls.y+gall$No.small.galls
gall$Total.galls<-gall$Total.galls+gall$Total.small.galls
head(gall)
gall<-gall[,-c(11,16,17)]
gall<-gall[!(gall$Total.OV<gall$Total.galls),]
#############################################################
#### gall emergence data ####
gall.emerg<-read.csv("2017 GFE Bt Emergence.csv")
gall.emerg<-gall.emerg[gall.emerg$Tree.Number!="",]
gall.emerg$Branch.Number[gall.emerg$Treatment=="Field"]<-NA
gall.emerg$Tree.ID<-paste(gall.emerg$Host.Plant,gall.emerg$Tree.Number,sep="_")
gall.emerg$Branch<-paste(gall.emerg$Tree.ID,gall.emerg$Branch.Number,sep="_")
head(gall.emerg)
gall.emerg$No.wasps<-gall.emerg$Number.of.Wasps..Alive.+gall.emerg$Number.of.Wasps..dead.
gall.emerg<-gall.emerg[gall.emerg$Treatment %in% c("Field","Cage"),]
gall.emerg<-gall.emerg[gall.emerg$Site %in% c("ABS","LL","Dick","Alva","KRE","Okee"),]
gall.emerg<-gall.emerg[gall.emerg$Gall.species=="B. treatae",]
gall.emerg<-gall.emerg[gall.emerg$Host.Plant %in% c("Qv","Qg"),]
gall.emerg$Gall.species<-factor(gall.emerg$Gall.species)
gall.emerg$Treatment<-factor(gall.emerg$Treatment)
gall.emerg$Site<-factor(gall.emerg$Site)
gall.emerg$Host.Plant<-factor(gall.emerg$Host.Plant,levels=c("Qv","Qg"))
gall.emerg$Tree.ID<-factor(gall.emerg$Tree.ID)
head(gall.emerg)
colnames(gall.emerg)
gall.emerg<-gall.emerg[,c(8,9,10,24,25,26)]
gall.emerg$M.site<-gall.emerg$Site
gall.emerg$F.site<-gall.emerg$Site
gall.emerg$HP.M<-"Qv"
gall.emerg$HP.M[gall.emerg$Site %in% c("ABS","LL","Dick")]<-"Qg"
gall.emerg$HP.F<-gall.emerg$HP.M
head(gall.emerg)
gall.emerg<-gall.emerg[,c(2,7,8,9,10,3,4,5,6)]
### read hybrid emergence data ##
hyb.emerg<-read.csv("Hybrid Emergence.csv")
hyb.emerg$Treatment<-"Cage"
hyb.emerg$Tree.ID<-paste(hyb.emerg$Host.Plant,hyb.emerg$Tree.Number,sep="_")
hyb.emerg$Branch<-paste(hyb.emerg$Tree.ID,hyb.emerg$Branch.Number,sep="_")
hyb.emerg<-hyb.emerg[,c(16,5,6,7,17,18,12)]
head(hyb.emerg)
colnames(hyb.emerg)[c(2,3,7)]<-c("M.site","F.site","No.wasps")
hyb.emerg$HP.M<-"Qv"
hyb.emerg$HP.M[hyb.emerg$M.site %in% c("Dick","ABS","LL")]<-"Qg"
hyb.emerg$HP.F<-"Qv"
hyb.emerg$HP.F[hyb.emerg$F.site %in% c("Dick","ABS","LL")]<-"Qg"
hyb.emerg<-hyb.emerg[,c(1,2,3,8,9,4,5,6,7)]
emerg<-rbind(gall.emerg,hyb.emerg)
emerg.sum<-ddply(emerg,c("Tree.ID","Branch"),
                 summarize,No.wasps=sum(No.wasps))
### merge gall and emergence data ##
gall<-merge(gall,emerg[,c(8,9)],by="Branch",all.x = TRUE)
tm1<-gall[gall$Total.galls<gall$No.wasps,]
tm1<-tm1[!is.na(tm1$Total.OV),]
gall<-gall[!(gall$Total.galls<gall$No.wasps),]
gall<-gall[!is.na(gall$Total.OV),]
write.csv(gall,"gall.summarize.data.csv")
######################################################
### merge the OV scar dataset ###
OV.mix<-rbind(OV.cage1,OV.Hyb1)
OV.mix<-OV.mix[OV.mix$No..HR!=34.54,]
OV.mix$Plant<-factor(OV.mix$Plant,levels=c("Qv","Qg"))
OV.mix$HP.M<-factor(OV.mix$HP.M,levels=c("Qv","Qg"))
OV.mix$HP.F<-factor(OV.mix$HP.F,levels=c("Qv","Qg"))
OV.mix$cross<-paste(OV.mix$HP.M,OV.mix$HP.F,sep=" x ")
OV.mix$cross<-factor(OV.mix$cross,levels=c("Qv x Qv","Qv x Qg","Qg x Qv","Qg x Qg"))
OV.mix$HR.OV<-OV.mix$No..HR/OV.mix$No..OV.scars
OV.mix$HR.OV.log<-log(OV.mix$HR.OV+0.0001)
OV.mix$count<-1
write.csv(OV.mix,"OV.data.csv")