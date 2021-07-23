library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(rlist)
library("survAUC")
#------------necessary---------------
Suvival.type<-"OS"

#interest<-"ZHongZ"
interest<-"Beijing"

pre_chemo<-"ignore"
#pre_chemo<-"Ture"
#pre_chemo<-"False"

#-load data----------
setwd("/data/backup/Yuni/CRC_Liver/12_feature_calculation/")

clinical1<-  read_excel("20201031_ZZ_CRLM_final_Eng.xlsx")
clinical11<-  read_excel("20201031_ZZ_CRLM_final_Eng.xlsx",sheet="new_145cases")
clinical2<- read_excel("Beijing_modified_clinical_info_English.xlsx")
beijing_list <- read_excel("image_list.xlsx")
Temp<-strsplit(beijing_list$对应病理号,split='-')
TT<-NULL
for(t in Temp){
  TT<-append(TT,t[1])
}
beijing_list$...3<-TT

Temp<-strsplit(beijing_list$`0`,split='*.svs')
TT<-NULL
for(t in Temp){
  TT<-append(TT,t[1])
}
beijing_list$`0`<-TT

clinical2$Liver_number<-as.character(clinical2$Liver_number)
clinical2$image_number<-beijing_list$`0`[match(clinical2$Liver_number,beijing_list$...3)]

Exit<-beijing_list$`0`[which(beijing_list$...3%in%clinical2$Liver_number)]
table(is.na(Exit))

clinical2<-clinical2[which(is.na(clinical2$image_number)==FALSE),]
dim(clinical2)
clinical2$Liver_number<-clinical2$image_number
clinical2$image_number<-NULL

table(clinical2$Colon_side)
temp<-ifelse((clinical2$Colon_side=="1"|clinical2$Colon_side=="2"),"1","0")
clinical2$Colon_side<-temp

table(clinical2$T_stage)
clinical2$T_stage<-ifelse((clinical2$T_stage=="0"|clinical2$T_stage=="1"|clinical2$T_stage=="2"),"0","1")

table(clinical2$N_stage)
clinical2$N_stage<-ifelse((clinical2$N_stage%in%c("0")),"0",ifelse(clinical2$N_stage%in%c("1","1a","1b","1c","2","2a","2b"),"1",NA))


table(clinical2$Prechemotherapy)
clinical2$Prechemotherapy<-ifelse((clinical2$Prechemotherapy%in%c("1")),"0","1")

table(clinical2$Other_metastasis)
clinical2$Other_metastasis<-ifelse((clinical2$Other_metastasis%in%c("0")),"0","1")

table(clinical2$Liver_location)
clinical2$Liver_location<-ifelse((clinical2$Liver_location%in%c("0")),"1","2")

summary(clinical2$Metastasis_size)
clinical2$Metastasis_size<-clinical2$Metastasis_size/10.0

table(clinical2$Resection_margin)
clinical2$Resection_margin<-ifelse((clinical2$Resection_margin%in%c("0")),"0",ifelse(clinical2$Resection_margin%in%c("1","2"),"1",NA))

table(clinical2$CRS_score)
clinical2$CRS_group<-ifelse((clinical2$CRS_score%in%c("0","1","2")),"0","1")
table(clinical2$CRS_group)

table(clinical1$CRS_score)
clinical1$CRS_group<-ifelse((clinical1$CRS_score%in%c("0","1","2")),"0","1")
table(clinical1$CRS_group)
table(clinical11$CRS_score)
clinical11$CRS_group<-ifelse((clinical11$CRS_score%in%c("0","1","2")),"0","1")
table(clinical11$CRS_group)

summary(clinical2$RFS_months)
clinical2$RFS_months<-ifelse(clinical2$RFS_months>=0,clinical2$RFS_months/30,NA)

summary(clinical2$OS_months)
clinical2$OS_months<-ifelse(clinical2$OS_months>=0,clinical2$OS_months/30,NA)

clinical1$...24<-NULL
colnames(clinical1)==colnames(clinical2)
colnames(clinical1)==colnames(clinical11)
clinical<-rbind(clinical1,clinical11,clinical2)
dim(clinical)
colnames(clinical)<-c("Hospital_number","ID","Gender","Age","Original_location","Original_Tstage","Original_Nstage","RAS_mutation","pre_chemotherapy","Other_Metastasis","CEA","DFI","liver_M_distribution","liver_M_number","liver_M_size","RR","radio_ablation","CRS_score","post_chemotherapy","Rfevent","RFS_month","Osevent" ,"OS_month","CRS_group")


filelist<-c("01","02","03","04","05")
for (f in filelist){
  #f<-filelist[1]
  assign(paste("proportion",f,sep="_"),read.csv2(paste0(f,'/tum_interection.csv'),header = T,skipNul =T,sep=",",fileEncoding="utf-8"))
  assign(paste("tum_interection",f,sep="_"),read.csv2(paste0(f,'/ALL_Tum_intereaction_score.csv'),header = T,skipNul =T,sep=",",fileEncoding="utf-8"))
  assign(paste("tum_entropy",f,sep="_"),read.csv2(paste0(f,'/tum_Entropy.csv'),header = T,skipNul =T,sep=",",fileEncoding="utf-8"))
  assign(paste("CH_distance",f,sep="_"),read.csv2(paste0(f,'/within_cluster_disp_score.csv'),header = T,skipNul =T,sep=",",fileEncoding="utf-8"))
  assign(paste("distance_STR",f,sep="_"),read.csv2(paste0(f,'/distance_STR_score.csv'),header = T,skipNul =T,sep=",",fileEncoding="utf-8"))
  assign(paste("distance_LYM",f,sep="_"),read.csv2(paste0(f,'/distance_LYM_score.csv'),header = T,skipNul =T,sep=",",fileEncoding="utf-8"))
}
for (f in filelist){
  print(paste0("dimension of ",paste("proportion",f,sep="_")))
  print(dim(get(paste("proportion",f,sep="_"))))
  print(get(paste("proportion",f,sep="_"))[1:5,])
}

Temp<-strsplit(proportion_04$X,split='/')
TT<-NULL
for(t in Temp){
  TT<-append(TT,t[2])
}
proportion_04$X<-TT


variableList<-c("proportion","tum_interection","tum_entropy","CH_distance","distance_STR","distance_LYM")

for (v in variableList){
  assign(v,get(paste(v,filelist[1],sep="_")))
  for (f in filelist[2:length(filelist)]){
    a<-get(paste(v,f,sep="_"))
    assign(v,rbind(get(v),a))
  }
}
dim(get("proportion"))
print(get("proportion")[1:5,])

#clinical<-read.delim('clinical.txt', header = TRUE, sep = "\t")


#------data-preprocessing-------------
#change the name---
Temp<-strsplit(proportion$X,split='[-(\\s)_ ]')
TT<-NULL
for(t in Temp){
  TT<-append(TT,t[1])
}
proportion$X<-TT

Temp<-strsplit(clinical$ID,split='[-(\\s)_ ]')
TT<-NULL
for(t in Temp){
  TT<-append(TT,t[1])
}
clinical$ID<-TT

table(clinical$Hospital_number=='x')

clinical<-clinical[which(clinical$Hospital_number!='x'),]
Tumor_name<-intersect(clinical$ID,proportion$X)
length(Tumor_name)
setdiff(clinical$ID, proportion$X)

surdata<-data.frame(clinical[which(clinical$ID%in%Tumor_name),]$ID)
names(surdata)<-"Sample.ID"
head(surdata)

if (Suvival.type=="OS"){
  surdata$time<-clinical[match(surdata$Sample.ID,clinical$ID),]$OS_month
  surdata$time<-as.numeric(surdata$time)
  surdata$event<-clinical[match(surdata$Sample.ID,clinical$ID),]$Osevent
  surdata$event<-as.numeric(surdata$event)
}
if (Suvival.type=="RFS"){
  surdata$time<-clinical[match(surdata$Sample.ID,clinical$ID),]$RFS_month
  surdata$time<-as.numeric(surdata$time)
  surdata$event<-clinical[match(surdata$Sample.ID,clinical$ID),]$Rfevent
  surdata$event<-as.numeric(surdata$event)
}

surdata$Age<-clinical[match(surdata$Sample.ID,clinical$ID),]$Age
surdata$Gender<-clinical[match(surdata$Sample.ID,clinical$ID),]$Gender
surdata$CRC_location<-clinical[match(surdata$Sample.ID,clinical$ID),]$Original_location
surdata$Tstage<-clinical[match(surdata$Sample.ID,clinical$ID),]$Original_Tstage
surdata$Nstage<-clinical[match(surdata$Sample.ID,clinical$ID),]$Original_Nstage
surdata$RAS_mutation<-clinical[match(surdata$Sample.ID,clinical$ID),]$RAS_mutation
surdata$pre_chemotherapy<-clinical[match(surdata$Sample.ID,clinical$ID),]$pre_chemotherapy
surdata$Other_Metastasis<-clinical[match(surdata$Sample.ID,clinical$ID),]$Other_Metastasis
surdata$CEA<-clinical[match(surdata$Sample.ID,clinical$ID),]$CEA
surdata$DFI<-clinical[match(surdata$Sample.ID,clinical$ID),]$DFI
surdata$liver_M_distribution<-clinical[match(surdata$Sample.ID,clinical$ID),]$liver_M_distribution
surdata$liver_M_number<-clinical[match(surdata$Sample.ID,clinical$ID),]$liver_M_number
surdata$liver_M_size<-clinical[match(surdata$Sample.ID,clinical$ID),]$liver_M_size
surdata$RR<-clinical[match(surdata$Sample.ID,clinical$ID),]$RR
surdata$radio_ablation<-clinical[match(surdata$Sample.ID,clinical$ID),]$radio_ablation
surdata$CRS_score<-clinical[match(surdata$Sample.ID,clinical$ID),]$CRS_score
surdata$CRS_group<-clinical[match(surdata$Sample.ID,clinical$ID),]$CRS_group
surdata$post_chemotherapy<-clinical[match(surdata$Sample.ID,clinical$ID),]$post_chemotherapy

dim(surdata)
head(surdata)

surdata$RAS_mutation<-ifelse(surdata$RAS_mutation=="2",NA,surdata$RAS_mutation)
surdata$post_chemotherapy<-ifelse(surdata$post_chemotherapy=="2",NA,surdata$post_chemotherapy)


variableList
for (v in variableList){
  #strsplit(tum_interection$X[388],split='/')
  Temp<-strsplit(get(v)$X,split='/')
  TT<-NULL
  for(t in Temp){
    if (is.na(t[2])==TRUE){
      TT<-append(TT,t[1])
    }
    else{
      TT<-append(TT,t[2])
    }

  }
  
  Temp<-strsplit(TT,split='[-(\\s)_ ]')
  TT<-NULL
  for(t in Temp){
    TT<-append(TT,t[1])
  }
  assign(paste0(v,"_name"),TT)
  #assign(paste0(v,"[[X]]"),TT)
}
proportion$X<-proportion_name
tum_interection$X<-tum_interection_name
tum_entropy$X<-tum_entropy_name
CH_distance$X<-CH_distance_name
distance_STR$X<-distance_STR_name
distance_LYM$X<-distance_LYM_name


STR_Proportion<-proportion[match(surdata$Sample.ID,proportion$X),]$STR
surdata$STR_Proportion<-as.numeric(as.character(STR_Proportion))
Entropy_TUM<-tum_entropy[match(surdata$Sample.ID,tum_entropy$X),]$Entropy.of.slide
surdata$Entropy<-as.numeric(as.character(Entropy_TUM))
LYM_Proportion<-proportion[match(surdata$Sample.ID,proportion$X),]$LYM
surdata$LYM_Proportion<-as.numeric(as.character(LYM_Proportion))
STR_interaction<-tum_interection[match(surdata$Sample.ID,tum_interection$X),]$STR
surdata$STR_interaction<-as.numeric(as.character(STR_interaction))
LYM_interaction<-tum_interection[match(surdata$Sample.ID,tum_interection$X),]$LYM
surdata$LYM_interaction<-as.numeric(as.character(LYM_interaction))

Within_cluster_disp<-CH_distance[match(surdata$Sample.ID,CH_distance$X),]$Extra_disp
surdata$Within_cluster_disp<-as.numeric(as.character(Within_cluster_disp))

Between_cluster_disp<-CH_distance[match(surdata$Sample.ID,CH_distance$X),]$Intra_disp
surdata$Between_cluster_disp<-as.numeric(as.character(Between_cluster_disp))

CH_score<-CH_distance[match(surdata$Sample.ID,CH_distance$X),]$CH_score
surdata$CH_score<-as.numeric(as.character(CH_score))

far_STR<-distance_STR[match(surdata$Sample.ID,distance_STR$X),]$far_STR
surdata$far_STR<-as.numeric(as.character(far_STR))

around_STR<-distance_STR[match(surdata$Sample.ID,distance_STR$X),]$around_STR
surdata$around_STR<-as.numeric(as.character(around_STR))

inside_STR<-distance_STR[match(surdata$Sample.ID,distance_STR$X),]$inside_STR
surdata$inside_STR<-as.numeric(as.character(inside_STR))

far_LYM<-distance_LYM[match(surdata$Sample.ID,distance_LYM$X),]$far_LYM
surdata$far_LYM<-as.numeric(as.character(far_LYM))

around_LYM<-distance_LYM[match(surdata$Sample.ID,distance_LYM$X),]$around_LYM
surdata$around_LYM<-as.numeric(as.character(around_LYM))

inside_LYM<-distance_LYM[match(surdata$Sample.ID,distance_LYM$X),]$inside_LYM
surdata$inside_LYM<-as.numeric(as.character(inside_LYM))

surdata$infiltrated_LYM<-surdata$around_LYM+surdata$inside_LYM
surdata$infiltrated_STR<-surdata$around_STR+surdata$inside_STR

surdata$whole_slide_LYM<-surdata$around_LYM+surdata$inside_LYM+surdata$far_LYM
surdata$whole_slide_STR<-surdata$around_STR+surdata$inside_STR+surdata$far_STR

#surdata<-surdata[complete.cases(surdata),]

#---cohort---------
surdata<-surdata[which(surdata$Other_Metastasis=="0"),]
surdata<-surdata[which(surdata$time<=60),]#只看五年生存
if (pre_chemo=="ignore"){
  surdata<-surdata
}
if (pre_chemo=="Ture"){
  surdata<-surdata[which(surdata$pre_chemotherapy=="1"),]
}
if(pre_chemo=="False"){
  surdata<-surdata[which(surdata$pre_chemotherapy=="0"),]
}
#
Beijing_surdata<-surdata[which(surdata$Sample.ID%in%clinical2$Liver_number),]
table(Beijing_surdata$Other_Metastasis=="0")
ZhongZ_surdata<-surdata[which(surdata$Sample.ID%in%clinical2$Liver_number==FALSE),]
table(ZhongZ_surdata$Other_Metastasis=="0")

if (interest=="Beijing"){
  surdata_interest<-Beijing_surdata
}
if (interest=="ZHongZ"){
  surdata_interest<-ZhongZ_surdata
}

if (Suvival.type=="OS"){
  rootPath<-"/data/backup/Yuni/CRC_Liver/12_feature_calculation/OS_results"
}
if (Suvival.type=="RFS"){
  rootPath<-"/data/backup/Yuni/CRC_Liver/12_feature_calculation/RFS_results"
}
#------------necessary---------------
#=========single factor survival & group the factor for univarite or multivariate analysis============
p<-list()
for (factor in colnames(surdata_interest)[22:39]){
  #factor<-"Entropy"
  print(factor)
  sur.cut <- surv_cutpoint(surdata_interest, time = "time", event = "event",variables = factor,minprop = 0.45)
  #summary(sur.cut)#5.22
  sur.cut <- surv_categorize(sur.cut)
  #head(sur.cut)
  #sur.cut[[factor]]
  surdata_interest[[factor]]<-sur.cut[[factor]]
  surdata_interest[[factor]]<-ifelse(surdata_interest[[factor]]=="high",1,0)
  #surdata_interest[[factor]]<-as.factor(surdata_interest[[factor]])
  formu <- paste("Surv(time, event)", factor, sep = "~")#factor
  formu <- as.formula(formu)
  print(formu)
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit)
  p<-list.append(p,ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type," "),pval=TRUE,risk.table = FALSE, conf.int = FALSE,title=paste0(factor,""),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12))
}

m<-floor(length(p)/6)

if (interest=="ZHongZ"){
  PATH=paste0(rootPath,"/other/single_factor/ZhongZ")
}
if (interest=="Beijing"){
  PATH=paste0(rootPath,"/other/single_factor/Beijing")
}
dir.create(PATH,recursive = TRUE)
for (i in 0:m){
  #arrange_ggsurvplots(p[12:18], print = TRUE, ncol = 3, nrow = 2)
  print(i*6+1)
  print(i*6+6)
  a<-i*6+1
  b<-i*6+6
  if (b>length(p)){
    b<-length(p)
  }
  p_temp<-arrange_ggsurvplots(p[a:b], ncol = 3, nrow = 2,print = FALSE)
  ggsave(paste0(i,"_single_factor_survival.pdf"),path=PATH,p_temp,width = 6*m, height = 4*m)#
}

#clnical_info_single------------
p<-list()
for (factor in colnames(surdata_interest)[4:21]){
  #factor<-"Age"
  print(factor)
  if (factor%in%c("Age","CEA","liver_M_size","liver_M_number","DFI")){
    sur.cut<-data.frame(time = surdata_interest$time, surdata_interest$event,factor=surdata_interest[[factor]])
    colnames(sur.cut)<-c("time","event",factor)
    cutpoint<-median(sur.cut[[factor]])
    sur.cut<-sur.cut[which(is.na(sur.cut[[factor]])==FALSE),]
    sur.cut[[factor]]<-ifelse(sur.cut[[factor]]>cutpoint,"high","low")
  }
  else{
    sur.cut<-data.frame(time = surdata_interest$time, surdata_interest$event,factor=surdata_interest[[factor]])
    colnames(sur.cut)<-c("time","event",factor)
  }
  #surdata_interest[[factor]]<-sur.cut[[factor]]
  #surdata_interest[[factor]]<-ifelse(surdata_interest[[factor]]=="high",1,0)
  #surdata_interest[[factor]]<-as.factor(surdata_interest[[factor]])
  formu <- paste("Surv(time, event)", factor, sep = "~")#factor
  formu <- as.formula(formu)
  print(formu)
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit)
  p<-list.append(p,ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type," "),pval=TRUE,risk.table = FALSE, conf.int = FALSE,title=paste0(factor,""),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12))
}

number_of_display<-6
m<-floor(length(p)/number_of_display)

if (interest=="ZHongZ"){
  PATH=paste0(rootPath,"/other/single_factor/ZhongZ")
}
if (interest=="Beijing"){
  PATH=paste0(rootPath,"/other/single_factor/Beijing")
}
for (i in 0:m){
  #arrange_ggsurvplots(p[12:18], print = TRUE, ncol = 3, nrow = 2)
  print(i*number_of_display+1)
  print(i*number_of_display+number_of_display)
  a<-i*number_of_display+1
  b<-i*number_of_display+number_of_display
  if (b>length(p)){
    b<-length(p)
  }
  p_temp<-arrange_ggsurvplots(p[a:b], ncol = 3, nrow = 2,print = FALSE)
  ggsave(paste0(i,"_single_clincial_factor_survival.pdf"),path=PATH,p_temp,width = number_of_display*m, height = 4*m)#
}


#-----stratify analysis-------------
data_interest<-surdata_interest
if (interest=="ZHongZ"){
  PATH=paste0(rootPath,"/other/stratify_analysis/ZhongZ")
}
if (interest=="Beijing"){
  PATH=paste0(rootPath,"/other/stratify_analysis/Beijing")
}
dir.create(PATH,recursive = TRUE)
clinical_factor<-c("post_chemotherapy","Gender","RR","CRC_location","pre_chemotherapy","CRS_group","Nstage","Tstage","liver_M_distribution")
for (clinic in clinical_factor){
  p<-list()
  for (factor in colnames(data_interest)[22:39]){
    #factor<-"infiltrated_LYM"
    print(factor)
    sur.cut <- surv_cutpoint(data_interest, time = "time", event = "event",variables = factor,minprop = 0.40)
    summary(sur.cut)#5.22
    sur.cut <- surv_categorize(sur.cut)
    #colnames(sur.cut)<-c("time","event","group")
    head(sur.cut)
    #sur.cut[[factor]]<-ifelse(sur.cut[[factor]]=="high","post_che&high","post_che&low")
    #data_interest[[factor]]<-data_interest[[factor]]
    #data_interest[[factor]]<-as.factor(data_interest[[factor]])
    formu <- paste("Surv(time, event)", "group", sep = "~")#factor
    formu <- as.formula(formu)
    print(formu)

    if (clinic=="post_chemotherapy"){
      sur.cut$post_chemotherapy<-surdata_interest$post_chemotherapy
      sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$post_chemotherapy==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$post_chemotherapy==0),"0_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$post_chemotherapy==0),"0_low","1_low")))
    }
    
    if  (clinic=="Gender"){
    sur.cut$Gender<-surdata_interest$Gender
    sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$Gender==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$Gender==2),"2_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$Gender==2),"2_low","1_low")))
    }
    if  (clinic=="RR"){
      sur.cut$RR<-surdata_interest$RR
      sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$RR==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$RR==0),"0_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$RR==0),"0_low","1_low")))
    }
    
    if  (clinic=="CRC_location"){
    sur.cut$CRC_location<-surdata_interest$CRC_location
    sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$CRC_location==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$CRC_location==0),"0_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$CRC_location==0),"0_low","1_low")))
    }
    
    if  (clinic=="pre_chemotherapy"){
    sur.cut$pre_chemotherapy<-surdata_interest$pre_chemotherapy
    sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$pre_chemotherapy==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$pre_chemotherapy==0),"0_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$pre_chemotherapy==0),"0_low","1_low")))
    }
    
    if (clinic=="CRS_group"){
    sur.cut$CRS_group<-surdata_interest$CRS_group
    sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$CRS_group==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$CRS_group==0),"0_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$CRS_group==0),"0_low","1_low")))
    }
    
    if (clinic=="Nstage"){
    sur.cut$Nstage<-surdata_interest$Nstage
    sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$Nstage==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$Nstage==0),"0_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$Nstage==0),"0_low","1_low")))
    }
    
    if (clinic=="Tstage"){
    sur.cut$Tstage<-surdata_interest$Tstage
    sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$Tstage==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$Tstage==0),"0_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$Tstage==0),"0_low","1_low")))
    }
    
    if (clinic=="liver_M_distribution"){
    sur.cut$liver_M_distribution<-surdata_interest$liver_M_distribution
    sur.cut$group<-ifelse(((sur.cut[[factor]]=="high")&(sur.cut$liver_M_distribution==1)),"1_high",ifelse((sur.cut[[factor]]=="high"&sur.cut$liver_M_distribution==2),"2_high",ifelse((sur.cut[[factor]]=="low"&sur.cut$liver_M_distribution==2),"2_low","1_low")))
    }

    nrow(sur.cut)
    fit <- surv_fit(formu, data = sur.cut)
    # summary(fit)
    p<-list.append(p,ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = FALSE, conf.int = FALSE,title=paste(clinic,factor,sep="-"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12))
  }
  m<-floor(length(p)/12)
  for (i in 0:m){
    #arrange_ggsurvplots(p[12:18], print = TRUE, ncol = 3, nrow = 2)
    print(i*12+1)
    print(i*12+12)
    a<-i*12+1
    b<-i*12+12
    if (b>length(p)){
      b<-length(p)
    }
    p_temp<-arrange_ggsurvplots(p[a:b], ncol = 3, nrow = 4,print = FALSE)
    
    ggsave(paste(clinic,paste0(i,"_stratify_survival_5years.pdf"),sep="-"),path=PATH,p_temp,width = 18*m, height = 20*m)#
  }
}



#========as.factor-----
# surdata_interest$Gender<-as.factor(surdata_interest$Gender)
# surdata_interest$Tstage<-as.factor(surdata_interest$Tstage)
# surdata_interest$Nstage<-as.factor(surdata_interest$Nstage)
# 
# surdata_interest$CRC_location<-as.factor(surdata_interest$CRC_location)
# surdata_interest$RAS_mutation<-as.factor(surdata_interest$RAS_mutation)
# surdata_interest$pre_chemotherapy<-as.factor(surdata_interest$pre_chemotherapy)
# surdata_interest$Other_Metastasis<-as.factor(surdata_interest$Other_Metastasis)
# surdata_interest$liver_M_distribution<-as.factor(surdata_interest$liver_M_distribution)
# surdata_interest$RR<-as.factor(surdata_interest$RR)
# surdata_interest$radio_ablation<-as.factor(surdata_interest$radio_ablation)
# #surdata_interest$CRS_score<-as.factor(surdata_interest$CRS_score)
# surdata_interest$CRS_group<-as.factor(surdata_interest$CRS_group)
# surdata_interest$post_chemotherapy<-as.factor(surdata_interest$post_chemotherapy)
# colnames(surdata_interest)

for (factor in colnames(surdata_interest)[22:39]){
  surdata_interest[[factor]]<-scale(surdata_interest[[factor]], center=T,scale=T)
}
#=====univariate-------------
#group the fector and as.facot are needed for this part

data_interest<-surdata_interest
covariates <- c("Gender","Age","CRC_location","Tstage","Nstage","RAS_mutation","CEA","DFI","liver_M_distribution","liver_M_number","liver_M_size","RR","radio_ablation","CRS_score","CRS_group","post_chemotherapy","STR_Proportion","Entropy","LYM_Proportion","STR_interaction","LYM_interaction","Within_cluster_disp","Between_cluster_disp","CH_score","far_STR","around_STR","inside_STR","far_LYM","around_LYM","inside_LYM","infiltrated_LYM","infiltrated_STR","whole_slide_LYM","whole_slide_STR")#,,"pre_chemotherapy"
table(Beijing_surdata$post_chemotherapy)
## 得到一个列表分别为，对每个变量构建的生存对象公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, event)~', x)))
univ_formulas
## 巧妙的配合lapply                        
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data_interest)})
univ_models
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR2 <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR,HR2,HR.confint.lower, HR.confint.upper,wald.test, p.value)
                         names(res)<-c("beta", "HR","HR (95% CI for HR)","lower","upper","wald.test", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
class(univ_results)


#STR_Proportion(univ_results)
res <- t(as.data.frame(univ_results, check.names = FALSE))
res<-as.data.frame(res)
res$p.value<-as.numeric(as.character(res$p.value))
res$significance<-ifelse(res$p.value<0.05,ifelse(res$p.value<0.01,ifelse(res$p.value<0.001,"***","**"),"*"),"")
library("forestplot")

Data_str<-cbind(rownames(res),as.numeric(as.vector(res$beta)),as.vector(res$`HR (95% CI for HR)`),as.vector(res$p.value),as.vector(res$significance))

Data_str<-rbind(c("Factors","beta","HR (95% CI for HR)","P value","significance"),Data_str)
if (interest=="ZHongZ"){
  PATH=paste0(rootPath,"/other/univariate_multivariate")
  pdf(file=paste(PATH,paste(Suvival.type,"univariate_analysis_ZhongZ.pdf",sep="-"),sep="/"))}
if (interest=="Beijing"){
  PATH=paste0(rootPath,"/other/univariate_multivariate")
  pdf(file=paste(PATH,paste(Suvival.type,"univariate_analysis_Beijing.pdf",sep="-"),sep="/"))
}
dir.create(PATH,recursive = TRUE)

forestplot(Data_str,  #显示的文本
           graph.pos=4,
           c(NA,as.numeric(as.vector(res$HR))), 
           c(NA,as.numeric(as.vector(res$lower))), 
           c(NA,as.numeric(as.vector(res$upper))), 
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.3, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 3,   # 误差条的线的宽度
           hrzl_lines = list("1" = gpar(lwd=1, columns=1:6, col = "#000044"),
                             "2" = gpar(lwd=1, columns=1:6, col = "#000044"),
                             "18"=gpar(lty=2), 
                             "36" = gpar(lwd=1, columns=1:6, col = "#000044")
                             ),
           title=paste0("Univariate Analysis Hazard Ratio Plot - ",Suvival.type),
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5),cex = 0.7), #文本大小设置
           lineheight = "auto"
)
dev.off()
#=====multivariate=============
#set.seed(1000)
#a <- sample(1:nrow(surdata_interest),nrow(surdata_interest)*0.6)
#b <- setdiff(1:nrow(surdata_interest),a)
#training_surdata<-surdata_interest[a,]
#validation_surdata<-surdata_interest[b,]

che_surdata<-surdata_interest[which(surdata_interest$post_chemotherapy=="1"),]
nrow(che_surdata)
no_che_surdata<-surdata_interest[which(surdata_interest$post_chemotherapy=="0"),]
nrow(no_che_surdata)


data_interest<-surdata_interest


res.table<-list()
p<-list()
for (factor in colnames(data_interest)[22:39]){
  #factor<-"STR_Proportion"
  print(factor)
  #data_interest[[factor]]<-scale(data_interest[[factor]], center=T,scale=T)
  formu <- paste("Surv(time, event)",paste0("Age+Gender+CRC_location+Tstage+Nstage+RAS_mutation+pre_chemotherapy+CEA+DFI+liver_M_distribution+liver_M_number+liver_M_size+radio_ablation+",factor), sep = "~")#
  formu <- as.formula(formu)
  print(formu)
  res.cox <- coxph(formu, data =  data_interest)#+ Tstage +Nstage
  nrow(summary(res.cox)$coefficients)
  #validation_feature<-validation_surdata[c("Gender","Age","CRC_location","Tstage","Nstage","RAS_mutation","pre_chemotherapy","Other_Metastasis","CEA","DFI","liver_M_distribution","liver_M_number","liver_M_size","radio_ablation","CRS_score","post_chemotherapy",factor)]
  #surv.train = Surv(data_interest$time, data_interest$event)
  #surv.test = Surv(validation_surdata$time, validation_surdata$event)
  
  #pre.res<-predict(res.cox,newdata=validation_surdata)
  #times <- seq(10, 100, 10)
  #AUC_hc <- AUC.hc(surv.train, surv.test,  pre.res, times)
  #Cindex.pre<-survConcordance(surv.test ~ predict(res.cox, validation_surdata), data = validation_surdata)
  
  res.table<-list.append(res.table, c(factor,summary(res.cox)$coefficients[, 2][nrow(summary(res.cox)$coefficients)],summary(res.cox)$coefficients[, 5][nrow(summary(res.cox)$coefficients)]))#,summary(res.cox)$concordance[1],Cindex.pre$concordance
  p<-list.append(p,ggforest(res.cox, data = data_interest,main = paste0("Multivariate Hazard ratio - ",Suvival.type)))
  
}
unlist(res.table)
res.df <- data.frame(matrix(unlist(res.table), nrow=length(res.table), byrow=T))
colnames(res.df)<-c("Factors",paste0("HR(Multivariate)-",Suvival.type),"P value")#,"C-index","Validation.Cindex"
res.df$`P value`<-as.numeric(as.character(res.df$`P value`))
res.df$significance<-ifelse(res.df$`P value`<0.05,ifelse(res.df$`P value`<0.01,ifelse(res.df$`P value`<0.001,"***","**"),"*"),"")
res.df

if (interest=="ZHongZ"){
  PATH=paste0(rootPath,"/other/univariate_multivariate")
  pdf(file=paste(PATH,paste(Suvival.type,"multivariate_analysis_ZhongZ.pdf",sep="-"),sep="/"))}
if (interest=="Beijing"){
  PATH=paste0(rootPath,"/other/univariate_multivariate")
  pdf(file=paste(PATH,paste(Suvival.type,"multivariate_analysis_Beijing.pdf",sep="-"),sep="/"))
}
dir.create(PATH,recursive = TRUE)



tbody.style = tbody_style(color = "black",
                          fill = c("#e8f3de", "#d3e8bb"), just="left")
ggtexttable(res.df, rows = NULL,
            theme = ttheme(
              colnames.style = colnames_style(color = "white", fill = "#8cc257"),
              tbody.style = tbody.style
            )
)
dev.off()
#p



#--clinical.survival-------------
training_surdata<-surdata[which(surdata$Sample.ID%in%ZhongZ_surdata$Sample.ID),]
validation_surdata<-surdata[which(surdata$Sample.ID%in%Beijing_surdata
                                  $Sample.ID),]

colnames(training_surdata)
formu <- paste("Surv(time, event)","DFI+liver_M_number+liver_M_size+CEA+Nstage", sep = "~")#+RAS_mutation+Other_Metastasis
formu <- as.formula(formu)
res.cox <-coxph(formu, data =  training_surdata)
ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))

wi<-summary(res.cox)$coefficients[,1]
names(wi)<-names(unlist(res.cox[["assign"]]))
surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
surdata2<-apply(surdata2, 2, as.numeric)
wi<-wi[match(colnames(surdata2),names(wi))]
multiply.res<-sweep(surdata2,2,wi,"*")
validation_surdata$scores.clinic<-rowSums(multiply.res)

sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "scores.clinic",minprop = 0.5)
summary(sur.cut)#5.22
sur.cut <- surv_categorize(sur.cut)
head(sur.cut)
formu <- paste("Surv(time, event)", "scores.clinic", sep = "~")
formu <- as.formula(formu)
print(formu)

fit <- surv_fit(formu, data = sur.cut)
# summary(fit,times =5*12)
p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0("Clinic"," related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)

surv_median(fit)
info=surv_median(fit)

p$plot <- p$plot+ 
  ggplot2::annotate("text", 
                    x = info$median[1], y = 0.5, # x and y coordinates of the text
                    label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                           x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                           label = round(info$median[2]/12,digits=2), size = 4)
p

#--CRS.group.survival-----

surdata_interest<-ZhongZ_surdata
data_interest<-surdata_interest
CRS_data=data.frame(time=data_interest$time, event=data_interest$event,CRS.Group=data_interest$CRS_group)
CRS_data$CRS.Group<-ifelse(CRS_data$CRS.Group==1,"high","low")
fit <- surv_fit(Surv(time, event)~CRS.Group, data = CRS_data)
# summary(fit,times =5*12)
p1<-ggsurvplot(fit, data = CRS_data,ylab=paste0(Suvival.type,"  "),pval=TRUE,pval.size=4,pval.method.size=4,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0("CRS.Group","-Sun Yet-Sen"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
info<-surv_median(fit)
info$median[1]

p1$plot <- p1$plot+ 
  ggplot2::annotate("text", 
                    x = info$median[1], y = 0.5, # x and y coordinates of the text
                    label = round(info$median[1]/12,digits=2), size = 3)+ggplot2::annotate("text", 
                                                                                           x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                           label = round(info$median[2]/12,digits=2), size = 3)
CRS_data$CRS.Group<-relevel(as.factor(CRS_data$CRS.Group), ref = "low")
res_cox<-coxph(Surv(time, event)~CRS.Group, data=CRS_data)
ggforest(res_cox)
p1$plot<-p1$plot+ggplot2::annotate("text",x = 15, y = 0.13,
                                   label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 25, y = 0.06,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))

CRS_HR_ZhongZ<-round(summary(res_cox)$conf.int[1],2)

data_interest<-Beijing_surdata
CRS_data=data.frame(time=data_interest$time, event=data_interest$event,CRS.Group=data_interest$CRS_group)
CRS_data$CRS.Group<-ifelse(CRS_data$CRS.Group==1,"high","low")
fit <- surv_fit(Surv(time, event)~CRS.Group, data = CRS_data)
# summary(fit,times =5*12)
p2<-ggsurvplot(fit, data = CRS_data,ylab=paste0(Suvival.type,"  "),pval=TRUE,pval.size=4,pval.method.size=4,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0("CRS.Group","-Beijing"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
info<-surv_median(fit)
info$median[1]

p2$plot <- p2$plot+ 
  ggplot2::annotate("text", 
                    x = info$median[1], y = 0.5, # x and y coordinates of the text
                    label = round(info$median[1]/12,digits=2), size = 3)+ggplot2::annotate("text", 
                                                                                           x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                           label = round(info$median[2]/12,digits=2), size = 3)
CRS_data$CRS.Group<-relevel(as.factor(CRS_data$CRS.Group), ref = "low")
res_cox<-coxph(Surv(time, event)~CRS.Group, data=CRS_data)
ggforest(res_cox)
p2$plot<-p2$plot+ggplot2::annotate("text",x = 13, y = 0.13,
                                 label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 22, y = 0.06,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
CRS_HR_Beijing<-round(summary(res_cox)$conf.int[1],2)

ggsurvlist <- list(
  x = p1,
  y = p2
)
# Arrange multiple ggsurvplots and print the output
p<-arrange_ggsurvplots(ggsurvlist, ncol = 2, nrow = 1,print = FALSE)
PATH=paste0(rootPath,"/00_clinical_only")
dir.create(PATH,recursive = TRUE)
ggsave("CRS_group_survival.pdf",path=PATH,width = 10, height =6,p)


#---image+CRS survival------
#surdata<-surdata[which(surdata$post_chemotherapy=="1"),]
#set.seed(1000)
#a <- sample(1:nrow(surdata),nrow(surdata)*0.6)
#b <- setdiff(1:nrow(surdata),a)
#training_surdata<-surdata[a,]
#validation_surdata<-surdata[b,]

training_surdata<-surdata[which(surdata$Sample.ID%in%ZhongZ_surdata$Sample.ID),]
validation_surdata<-surdata[which(surdata$Sample.ID%in%Beijing_surdata
$Sample.ID),]

dim(validation_surdata)
dim(training_surdata)
colnames(training_surdata)
allp<-list()
effect_res<-list()
res.table<-list()
for (factor in colnames(training_surdata[,22:39])){
  #factor<-"Entropy"
  formu <- paste("Surv(time, event)",paste0("CRS_score+",factor), sep = "~")#+RAS_mutation+Other_Metastasis#Gender+Age+CRC_location+Tstage+pre_chemotherapy+liver_M_distribution+radio_ablation+###paste0("CRS_score+",factor)
  formu <- as.formula(formu)
  temp<-training_surdata[which(is.na(training_surdata[[factor]])==FALSE),]
  #dim(temp)
  #dim(training_surdata)
  res.cox <-coxph(formu, data =  temp)
  CRS.cox <-coxph(Surv(time, event) ~ CRS_score, data =  temp)
  anova.res<-anova(res.cox, CRS.cox)
  print(paste(factor,anova.res$`P(>|Chi|)`[2],sep='  '))
  ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))
  res.table<-list.append(res.table, c(factor,summary(res.cox)$coefficients[, 2][nrow(summary(res.cox)$coefficients)],summary(res.cox)$coefficients[, 5][nrow(summary(res.cox)$coefficients)]))
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-training_surdata[,which(colnames(training_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  training_surdata$risk.score<-rowSums(multiply.res)
  #print(factor)
  sur.cut <- surv_cutpoint(training_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.45)
  summary(sur.cut)#5.22
  cutpoint<-summary(sur.cut)[1]
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  #print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p1<-ggsurvplot(fit, font.main=10,data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,pval.size=4,pval.method.size=4,risk.table = FALSE,surv.median.line = "hv", conf.int = FALSE,title=paste0(factor,"+CRS_group - Sun Yet-Sen"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p1$plot <- p1$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 3)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 3)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  p1$plot<-p1$plot+ggplot2::annotate("text",x = 15, y = 0.13,
                                    label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 25, y = 0.06,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  
  HR1<-round(summary(res_cox)$conf.int[1],2)
  #validation..........................
  surdata3<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata3<-apply(surdata3, 2, as.numeric)
  multiply.res<-sweep(surdata3,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  #print(factor)
  #sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.45)
  #summary(sur.cut)#5.22
  #sur.cut <- surv_categorize(sur.cut)
  validation_surdata$risk.score<-ifelse(validation_surdata$risk.score>=as.numeric(cutpoint),"high","low")
  sur.cut<-data.frame(time=validation_surdata$time,event=validation_surdata$event,risk.score=validation_surdata$risk.score)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  #print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p2<-ggsurvplot(fit, font.main=10,data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,pval.size=4,pval.method.size=4,risk.table = FALSE,surv.median.line = "hv", conf.int = FALSE,title=paste0(factor,"+CRS_group - Beijing"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p2$plot <- p2$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 3)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 3)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  ggforest(res_cox)
  p2$plot<-p2$plot+ggplot2::annotate("text",x = 13, y = 0.13,
                                   label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 22, y = 0.06,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  
  HR2<-round(summary(res_cox)$conf.int[1],2)

  
  ggsurvlist <- list(
    x = p1,
    y = p2
  )
  # Arrange multiple ggsurvplots and print the output
  p<-arrange_ggsurvplots(ggsurvlist, ncol = 2, nrow = 1,print = FALSE)
  PATH="/data/backup/Yuni/CRC_Liver/12_feature_calculation/001_results/01_CRC_image"
  #ggsave(paste0(factor,"+CRS_survival.pdf"),path=PATH,width = 10, height = 5,p)
  allp<-list.append(allp,p)
  if (HR1>CRS_HR_ZhongZ&HR2>CRS_HR_Beijing){
    effect_res<-list.append(effect_res,factor)
  }
}
effect_res
#allp

unlist(res.table)
res.df <- data.frame(matrix(unlist(res.table), nrow=length(res.table), byrow=T))
colnames(res.df)<-c("Factors",paste0("HR(Multivariate)-",Suvival.type),"P value")#,"C-index","Validation.Cindex"
res.df$`P value`<-as.numeric(as.character(res.df$`P value`))
res.df$significance<-ifelse(res.df$`P value`<0.05,ifelse(res.df$`P value`<0.01,ifelse(res.df$`P value`<0.001,"***","**"),"*"),"")
res.df
tbody.style = tbody_style(color = "black",
                          fill = c("#e8f3de", "#d3e8bb"), just="left")
ggtexttable(res.df, rows = NULL,
            theme = ttheme(
              colnames.style = colnames_style(color = "white", fill = "#8cc257"),
              tbody.style = tbody.style
            )
)


#-----image+CRS+clinic survival-暂时不要-----
colnames(training_surdata)
allp<-list()
res.table<-list()
for (factor in colnames(training_surdata[,22:39])){
  #factor<-"whole_slide_LYM"
  formu <- paste("Surv(time, event)",paste0("pre_chemotherapy+Gender+CRS_group+",factor), sep = "~")#
  #with_che:X:radio, CRC_location,RAS_mutation,
  #with_che:ok:pre_chemotherapy,liver_M_distribution,Gender
  #CRC_location+Tstage+liver_M_distribution+radio_ablation+
  formu <- as.formula(formu)
  anova.res<-anova(res.cox, CRS.cox)
  print(paste(factor,anova.res$`P(>|Chi|)`[2],sep='  '))
  ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))
  res.table<-list.append(res.table, c(factor,summary(res.cox)$coefficients[, 2][nrow(summary(res.cox)$coefficients)],summary(res.cox)$coefficients[, 5][nrow(summary(res.cox)$coefficients)]))
  
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  
  
  sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.49)
  summary(sur.cut)#5.22
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0(factor,"+CRS+clinic Cox-Score related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p$plot <- p$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 4)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  ggforest(res_cox)
  p$plot<-p$plot+ggplot2::annotate("text",x = 15, y = 0.12,
                                   label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 15, y = 0.07,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  allp<-list.append(allp,p)
}
#allp


#------two+CRS------暂时不要------------
colnames(training_surdata)
allp<-list()
factors<-list(c("STR_Proportion","LYM_Proportion"),c("STR_interaction","LYM_interaction"),c("around_STR","around_LYM"),c("inside_STR","inside_LYM"),c("infiltrated_LYM","infiltrated_STR"),c("whole_slide_LYM","whole_slide_STR"),c("whole_slide_LYM","Between_cluster_disp"),c("whole_slide_LYM","CH_score"),c("whole_slide_STR","CH_score"),c("STR_Proportion","STR_interaction"),c("CH_score","LYM_Proportion"))
for (factor in factors){
  #print(factor)
  #factor<-"whole_slide_LYM"
  formu <- paste("Surv(time, event)",paste0("CRS_group+",paste(factor[1],factor[2],sep="+")), sep = "~")#
  #pre_chemotherapy+Gender+?
  formu <- as.formula(formu)
  res.cox <-coxph(formu, data =  training_surdata)
  Clinic.cox <-coxph(Surv(time, event) ~ CRS_group , data =  training_surdata)
  #anova.res<-anova(res.cox, Clinic.cox)
  #print(paste(paste(factor[1],factor[2],sep="+"),anova.res$`P(>|Chi|)`[2],sep='  '))
  
  #ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))
  
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  
  sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.49)
  summary(sur.cut)#5.22
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0(paste(factor[1],factor[2],sep="+"),"+CRS Cox-Score related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p$plot <- p$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 4)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  ggforest(res_cox)
  p$plot<-p$plot+ggplot2::annotate("text",x = 15, y = 0.12,
                                   label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 15, y = 0.07,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  allp<-list.append(allp,p)
}

allp

#-----image+clinic survival------
colnames(training_surdata)
f<-list()
allp<-list()
effect_res<-list()
res.table<-list()
for (factor in colnames(training_surdata[,22:39])){
  #factor<-"whole_slide_LYM"
  formu <- paste("Surv(time, event)",paste0("DFI+liver_M_number+CEA+",factor), sep = "~")
  
  #DFI+liver_M_number+CEA+(OK)
  #DFI+liver_M_number+CEA+Nstage+(OK)
  #DFI+liver_M_number+CEA+liver_M_distribution+(OK)
  #DFI+liver_M_number+CEA+liver_M_distribution+pre_chemotherapy+(pre加进来就会变得特别差！)
  #DFI+liver_M_number+CEA+liver_M_distribution+radio_ablation+(没有提升)
  #DFI+liver_M_number+CEA+liver_M_distribution+Gender+(Gender加进来就特别差)
  formu <- as.formula(formu)
  temp<-training_surdata[which(is.na(training_surdata[[factor]])==FALSE),]
  #dim(temp)
  #dim(training_surdata)
  res.cox <-coxph(formu, data =  temp)
  CRS.cox <-coxph(Surv(time, event) ~ DFI+liver_M_number+CEA, data =  temp)
  anova.res<-anova(res.cox, CRS.cox)
  print(paste(factor,anova.res$`P(>|Chi|)`[2],sep='  '))
  ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))
  res.table<-list.append(res.table, c(factor,summary(res.cox)$coefficients[, 2][nrow(summary(res.cox)$coefficients)],summary(res.cox)$coefficients[, 5][nrow(summary(res.cox)$coefficients)]))
  
  #Clinic.cox <-coxph(Surv(time, event) ~ DFI + liver_M_number + CEA+Nstage , data =  training_surdata)
  #anova.res<-anova(res.cox, Clinic.cox)
  #print(paste(factor,anova.res$`P(>|Chi|)`[2],sep='  '))
  
  f<-list.append(f,ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type)))
  
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-training_surdata[,which(colnames(training_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  training_surdata$risk.score<-rowSums(multiply.res)
  #print(factor)
  sur.cut <- surv_cutpoint(training_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.45)
  summary(sur.cut)#5.22
  cutpoint<-summary(sur.cut)[1]
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  #print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p1<-ggsurvplot(fit, font.main=10,data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,pval.size=4,pval.method.size=4,risk.table = FALSE,surv.median.line = "hv", conf.int = FALSE,title=paste0(factor,"+Clinic - Sun Yet-Sen"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p1$plot <- p1$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 3)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 3)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  p1$plot<-p1$plot+ggplot2::annotate("text",x = 15, y = 0.13,
                                     label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 25, y = 0.06,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  
  HR1<-round(summary(res_cox)$conf.int[1],2)
  #validation..........................
  surdata3<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata3<-apply(surdata3, 2, as.numeric)
  multiply.res<-sweep(surdata3,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  #print(factor)
  #sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.45)
  #summary(sur.cut)#5.22
  #sur.cut <- surv_categorize(sur.cut)
  validation_surdata$risk.score<-ifelse(validation_surdata$risk.score>=as.numeric(cutpoint),"high","low")
  sur.cut<-data.frame(time=validation_surdata$time,event=validation_surdata$event,risk.score=validation_surdata$risk.score)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  #print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p2<-ggsurvplot(fit, font.main=10,data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,pval.size=4,pval.method.size=4,risk.table = FALSE,surv.median.line = "hv", conf.int = FALSE,title=paste0(factor,"+Clinic - Beijing"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p2$plot <- p2$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 3)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 3)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  ggforest(res_cox)
  p2$plot<-p2$plot+ggplot2::annotate("text",x = 13, y = 0.13,
                                     label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 22, y = 0.06,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  
  HR2<-round(summary(res_cox)$conf.int[1],2)
  
  
  ggsurvlist <- list(
    x = p1,
    y = p2
  )
  # Arrange multiple ggsurvplots and print the output
  p<-arrange_ggsurvplots(ggsurvlist, ncol = 2, nrow = 1,print = FALSE)
  PATH="/data/backup/Yuni/CRC_Liver/12_feature_calculation/001_results/02_Clinic_image"
  #ggsave(paste0(factor,"+Clinic_survival.pdf"),path=PATH,width = 10, height = 5,p)
  allp<-list.append(allp,p)
  if (HR1>CRS_HR_ZhongZ&HR2>CRS_HR_Beijing){
    effect_res<-list.append(effect_res,factor)
  }
}
effect_res
#f

unlist(res.table)
res.df <- data.frame(matrix(unlist(res.table), nrow=length(res.table), byrow=T))
colnames(res.df)<-c("Factors",paste0("HR(Multivariate)-",Suvival.type),"P value")#,"C-index","Validation.Cindex"
res.df$`P value`<-as.numeric(as.character(res.df$`P value`))
res.df$significance<-ifelse(res.df$`P value`<0.05,ifelse(res.df$`P value`<0.01,ifelse(res.df$`P value`<0.001,"***","**"),"*"),"")
res.df
tbody.style = tbody_style(color = "black",
                          fill = c("#e8f3de", "#d3e8bb"), just="left")
ggtexttable(res.df, rows = NULL,
            theme = ttheme(
              colnames.style = colnames_style(color = "white", fill = "#8cc257"),
              tbody.style = tbody.style
            )
)
#allp
#------two+Clinic------
colnames(training_surdata)
allp<-list()
f<-list()
factors<-list(c("STR_Proportion","LYM_Proportion"),c("STR_interaction","LYM_interaction"),c("around_STR","around_LYM"),c("inside_STR","inside_LYM"),c("infiltrated_LYM","infiltrated_STR"),c("whole_slide_LYM","whole_slide_STR"),c("whole_slide_LYM","Between_cluster_disp"),c("whole_slide_LYM","CH_score"),c("whole_slide_STR","CH_score"),c("STR_Proportion","STR_interaction"),c("CH_score","LYM_Proportion"),c("whole_slide_STR","Between_cluster_disp"),c("Within_cluster_disp","infiltrated_STR"))
for (factor in factors){
  #print(factor)
  #factor<-"whole_slide_LYM"
  formu <- paste("Surv(time, event)",paste0("DFI+liver_M_number+CEA+",paste(factor[1],factor[2],sep="+")), sep = "~")
  #DFI+liver_M_number+CEA+(OK)
  #DFI+liver_M_number+CEA+Nstage+(OK)
  #DFI+liver_M_number+CEA+liver_M_distribution+(OK)
  formu <- as.formula(formu)
  res.cox <-coxph(formu, data =  training_surdata)
  
  #Clinic.cox <-coxph(Surv(time, event) ~ DFI + liver_M_number + CEA , data =  training_surdata)
  #anova.res<-anova(res.cox, Clinic.cox)
  #print(paste(paste(factor[1],factor[2],sep="+"),anova.res$`P(>|Chi|)`[2],sep='  '))
  
  
  f<-list.append(f,ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type)))
  
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  
  sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.5)
  summary(sur.cut)#5.22
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0(paste(factor[1],factor[2],sep="+"),"+clinic Cox-Score related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p$plot <- p$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 4)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  ggforest(res_cox)
  p$plot<-p$plot+ggplot2::annotate("text",x = 15, y = 0.12,
                                   label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 15, y = 0.07,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  allp<-list.append(allp,p)
}

allp
f
#-----image+clinic final survival------
colnames(training_surdata)
allp<-list()
formu <- paste("Surv(time, event)","CH_score+whole_slide_LYM+liver_M_number+Age+CEA", sep = "~")
#+liver_M_number+liver_M_size+liver_M_distribution+Nstage
formu <- as.formula(formu)
res.cox <-coxph(formu, data =  training_surdata)
ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))

wi<-summary(res.cox)$coefficients[,1]
names(wi)<-names(unlist(res.cox[["assign"]]))
surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
surdata2<-apply(surdata2, 2, as.numeric)
wi<-wi[match(colnames(surdata2),names(wi))]
multiply.res<-sweep(surdata2,2,wi,"*")
validation_surdata$risk.score<-rowSums(multiply.res)


sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.49)
summary(sur.cut)#5.22
sur.cut <- surv_categorize(sur.cut)
head(sur.cut)
formu <- paste("Surv(time, event)", "risk.score", sep = "~")
formu <- as.formula(formu)
print(formu)

fit <- surv_fit(formu, data = sur.cut)

# summary(fit,times =5*12)
p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0("image","+clinic Cox-Score related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)

info=surv_median(fit)

p$plot <- p$plot+ 
  ggplot2::annotate("text", 
                    x = info$median[1], y = 0.5, # x and y coordinates of the text
                    label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                           x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                           label = round(info$median[2]/12,digits=2), size = 4)
sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
res_cox<-coxph(formu, data=sur.cut)
ggforest(res_cox)
p$plot<-p$plot+ggplot2::annotate("text",x = 15, y = 0.12,
                                 label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 15, y = 0.07,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
p
summary(res_cox)$coef[2]

#--image.score.survival-----------
formu <- paste("Surv(time, event)","CH_score+whole_slide_LYM+inside_LYM", sep = "~")#STR_Proportion+STR_interaction+LYM_interaction+Within_cluster_disp+Between_cluster_disp+Entropy+infiltrated_STR+infiltrated_LYM+whole_slide_STR+whole_slide_LYMn+Other_Metastasis

#---image+clinic survival------
colnames(training_surdata)
allp<-list()
for (factor in colnames(training_surdata[,22:39])){
  #factor<-"whole_slide_LYM"
  formu <- paste("Surv(time, event)",paste0("CEA+liver_M_distribution+liver_M_number+Nstage+",factor), sep = "~")#+RAS_mutation+Other_Metastasis#Gender+Age+CRC_location+Tstage+pre_chemotherapy+liver_M_distribution+radio_ablation+###paste0("CRS_score+",factor)
  formu <- as.formula(formu)
  res.cox <-coxph(formu, data =  training_surdata)
  ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))
  
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  
  
  sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.49)
  summary(sur.cut)#5.22
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0(factor,"+Clinic Cox-Score related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p$plot <- p$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 4)
  allp<-list.append(allp,p)
}
allp


#---two.image.only------
colnames(training_surdata)
allp<-list()
factors<-list(c("STR_Proportion","LYM_Proportion"),c("STR_interaction","LYM_interaction"),c("around_STR","around_LYM"),c("inside_STR","inside_LYM"),c("infiltrated_LYM","infiltrated_STR"),c("whole_slide_LYM","whole_slide_STR"),c("whole_slide_LYM","Between_cluster_disp"),c("whole_slide_LYM","CH_score"),c("whole_slide_STR","CH_score"),c("STR_Proportion","STR_interaction"),c("CH_score","LYM_Proportion"),c("whole_slide_STR","Between_cluster_disp"),c("Within_cluster_disp","infiltrated_STR"))
for (factor in factors){
  print(factor)
  #factor<-c("STR_Proportion","LYM_Proportion")
  formu <- paste("Surv(time, event)",paste0("",paste(factor[1],factor[2],sep="+")), sep = "~")##Gender+Age+CRC_location+Tstage+pre_chemotherapy+liver_M_distribution+radio_ablation+###paste0("CRS_group+",paste(factor[1],factor[2],sep="+")), sep = "~")
  formu <- as.formula(formu)
  res.cox <-coxph(formu, data =  training_surdata)
  ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))
  
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  
  sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.40)
  summary(sur.cut)#5.22
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  #---
  # cutoff<-summary(validation_surdata$risk.score)[3]
  # validation_surdata$risk.score<-ifelse(validation_surdata$risk.score>=cutoff,"high","low")
  # sur.cut<-validation_surdata[c("time","event","risk.score")]
  #---
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0(paste(factor[1],factor[2],sep="+"),"+CRS Cox-Score related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p$plot <- p$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 4)
  sur.cut$risk.score<-relevel(as.factor(sur.cut$risk.score), ref = "low")
  res_cox<-coxph(formu, data=sur.cut)
  ggforest(res_cox)
  p$plot<-p$plot+ggplot2::annotate("text",x = 15, y = 0.12,
                                   label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 15, y = 0.07,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))
  allp<-list.append(allp,p)
}

allp
#---3.image.only------
colnames(training_surdata)
allp<-list()
factors<-list(c("STR_Proportion","LYM_Proportion","CH_score"),c("STR_interaction","LYM_interaction","CH_score"),c("around_STR","around_LYM","CH_score"),c("inside_STR","inside_LYM","CH_score"),c("infiltrated_LYM","infiltrated_STR","CH_score"),c("whole_slide_LYM","whole_slide_STR","CH_score"),c("whole_slide_LYM","Between_cluster_disp","Within_cluster_disp"),c("whole_slide_LYM","CH_score","Within_cluster_disp"),c("whole_slide_STR","CH_score","infiltrated_LYM"),c("STR_Proportion","STR_interaction","whole_slide_LYM"),c("CH_score","LYM_Proportion","infiltrated_STR"))
for (factor in factors){
  print(factor)
  #factor<-"whole_slide_LYM"
  formu <- paste("Surv(time, event)",paste0("",paste(factor[1],factor[2],factor[3],sep="+")), sep = "~")##Gender+Age+CRC_location+Tstage+pre_chemotherapy+liver_M_distribution+radio_ablation+###paste0("CRS_group+",paste(factor[1],factor[2],sep="+")), sep = "~")
  formu <- as.formula(formu)
  res.cox <-coxph(formu, data =  training_surdata)
  ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))
  
  wi<-summary(res.cox)$coefficients[,1]
  names(wi)<-names(unlist(res.cox[["assign"]]))
  surdata2<-validation_surdata[,which(colnames(validation_surdata)%in%names(unlist(res.cox[["assign"]])))]
  surdata2<-apply(surdata2, 2, as.numeric)
  wi<-wi[match(colnames(surdata2),names(wi))]
  multiply.res<-sweep(surdata2,2,wi,"*")
  validation_surdata$risk.score<-rowSums(multiply.res)
  
  sur.cut <- surv_cutpoint(validation_surdata, time = "time", event = "event",variables = "risk.score",minprop = 0.45)
  summary(sur.cut)#5.22
  sur.cut <- surv_categorize(sur.cut)
  head(sur.cut)
  formu <- paste("Surv(time, event)", "risk.score", sep = "~")
  formu <- as.formula(formu)
  print(formu)
  
  fit <- surv_fit(formu, data = sur.cut)
  # summary(fit,times =5*12)
  p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0(paste(factor[1],factor[2],sep="+"),"+CRS Cox-Score related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
  
  info=surv_median(fit)
  
  p$plot <- p$plot+ 
    ggplot2::annotate("text", 
                      x = info$median[1], y = 0.5, # x and y coordinates of the text
                      label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                             x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                             label = round(info$median[2]/12,digits=2), size = 4)
  allp<-list.append(allp,p)
}

allp
#----all.image.cox---------
a<-paste(colnames(training_surdata)[22:25],collapse ="+")
formu <- paste("Surv(time, event)",a, sep = "~")##Gender+Age+CRC_location+Tstage+pre_chemotherapy+liver_M_distribution+radio_ablation+###paste0("CRS_group+",paste(factor[1],factor[2],sep="+")), sep = "~")
formu <- as.formula(formu)
res.cox <-coxph(formu, data =  training_surdata)
ggforest(res.cox, data = training_surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))

#------traning.survival-----------
wi<-summary(res.cox)$coefficients[,1]
names(wi)<-names(unlist(res.cox[["assign"]]))
surdata2<-training_surdata[,which(colnames(training_surdata)%in%names(unlist(res.cox[["assign"]])))]
surdata2<-apply(surdata2, 2, as.numeric)
wi<-wi[match(colnames(surdata2),names(wi))]
multiply.res<-sweep(surdata2,2,wi,"*")
training_surdata$scores.image<-rowSums(multiply.res)

sur.cut <- surv_cutpoint(training_surdata, time = "time", event = "event",variables = "scores.image",minprop = 0.5)
summary(sur.cut)#5.22
sur.cut <- surv_categorize(sur.cut)
head(sur.cut)
formu <- paste("Surv(time, event)", "scores.image", sep = "~")
formu <- as.formula(formu)
print(formu)

fit <- surv_fit(formu, data = sur.cut)
# summary(fit,times =5*12)
ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0("scores.image"," related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)

surv_median(fit)



#========optimal cox===========
library("My.stepwise")
#readRDS("C:/Program Files/R/R-3.6.2/library/My.stepwise/R/My.stepwise.rdx")

my.data <- training_surdata#validation_surdata#training_surdata
colnames(my.data)
my.data$RR<-NULL
my.data$post_chemotherapy<-NULL
#my.data$RAS_mutation<-NULL
my.data$CRS_group<-NULL
#my.data$CRS_score<-NULL
my.data$Other_Metastasis<-NULL
my.data <- na.omit(my.data)
dim(my.data)
head(my.data)
colnames(my.data)
my.variable.list <-as.vector(colnames(my.data[,4:34]))
#my.variable.list<-c("CRS_score",my.variable.list)
My.stepwise.coxph(Time = "time", Status = "event", variable.list = my.variable.list, data = my.data, sle = 0.1, sls = 0.1,vif.threshold =5)#,in.variable = c("CRS_score")

My.stepwise.coxph(Time = "time", Status = "event", variable.list = my.variable.list,data = my.data, sle = 0.1, sls = 0.1)





#-------try-------
library(survival)
library(pec)
library(prodlim)
library(riskRegression)
library("pec")
f<-selectCox(Surv(time,event) ~ Age+Gender+CRC_location+Tstage+Nstage+CEA+DFI+liver_M_distribution+liver_M_size+radio_ablation+STR_Proportion+STR_interaction+LYM_interaction+Within_cluster_disp+Between_cluster_disp+Entropy+infiltrated_STR+infiltrated_LYM+whole_slide_STR+whole_slide_LYM,data=my.data,rule = "aic")
f


predictSurvProb(f,newdata=my.data[,4:34],times=10)




#------try2----------
library(rms)
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, TRUE))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
S <- Surv(dt,e)

f <- cph(S ~ age*sex, x=TRUE, y=TRUE)
# Validate full model fit
validate(f, B=10)               # normally B=150

# Validate a model with stratification.  Dxy is the only
# discrimination measure for such models, by Dxy requires
# one to choose a single time at which to predict S(t|X)
f <- cph(S ~ rcs(age)*strat(sex), 
         x=TRUE, y=TRUE, surv=TRUE, time.inc=2)
validate(f, B=40, bw=TRUE,u=2)


#=======lasso cox===========
library("glmnet")
library("survival")

surdata_numeric<-apply(surdata[2:ncol(surdata)], 2, as.numeric)

Z<-surdata_numeric[,3:ncol(surdata_numeric)]
Z<-Z[,-6]#remove RAS_mutation due to NA
colnames(Z)
Z<-Z[,-13]
Y<-surdata_numeric[,1:2]
colnames(Y)<-c('time','status')


table(is.na(Z))
index<-complete.cases(Z)
Z<-Z[index,]
Y<-Y[index,]

res1<-glmnet(Z, Y, family="cox")
#penalty parameter values:
res1$lambda

plot(res1)


#using 5-fold CV to select lambda:
res2<-cv.glmnet(Z, Y,family="cox", nfolds=10,type.measure = "C")
#Warning message:

plot(res2)
title("Lasso Cox regression")


res2$lambda.min

coef.min = coef(res2, s = "lambda.min")
active.min = which(coef.min != 0)
index.min = coef.min[active.min]
coef.min
summary(res2)
print(res2, digits = max(3, getOption("digits") - 3))
# final Lasso model:
res3<-glmnet(Z, Y, family="cox", lambda=res2$lambda.min)
coef(res3)

#predict(res3, newx, s = "lambda.min")
#Cindex(res3, Y, weights = rep(1, nrow(Y)))
#=======cox survival==========
library("glmnet")
#lasso 模型，交叉验证选最优---------------
surdata<-surdata[which(surdata$Other_Metastasis=="0"),]
#surdata<-surdata[which(surdata$post_chemotherapy=="1"),]

colnames(training_surdata)
my.data <- training_surdata#validation_surdata#training_surdata
colnames(my.data)
my.data$RR<-NULL
my.data$post_chemotherapy<-NULL
my.data$RAS_mutation<-NULL
my.data$CRS_group<-NULL
my.data$CRS_score<-NULL
my.data <- na.omit(my.data)
colnames(my.data)
x=as.matrix(my.data[17:34])
rownames(x)<-my.data$Sample.ID
y=data.matrix(Surv(my.data$time,my.data$event))
fit=glmnet(x, y, family = "cox", maxit = 1000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000,type.measure = "C")#,nfolds =5


coef=coef(fit, s = cvfit$lambda.1se)
index=which(coef != 0)
actCoef=coef[index]
lassoFeature=row.names(coef)[index]
FeatureCoef=cbind(Feature=lassoFeature,Coef=actCoef)
FeatureCoef

trainFinalGeneExp=training_surdata[,lassoFeature]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
training_surdata$risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))

valFinalGeneExp=validation_surdata[,lassoFeature]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
valScore=apply(valFinalGeneExp,1,myFun)
validation_surdata$risk=as.vector(ifelse(valScore>median(valScore),"high","low"))

sur.cut<-validation_surdata[c("time","event","risk")]

formu <- paste("Surv(time, event)", "risk", sep = "~")
formu <- as.formula(formu)
print(formu)

fit <- surv_fit(formu, data = sur.cut)
# summary(fit,times =5*12)
p<-ggsurvplot(fit, data = sur.cut,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0("","  "),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)

surv_median(fit)
info=surv_median(fit)

p$plot <- p$plot+ 
  ggplot2::annotate("text", 
                    x = info$median[1], y = 0.5, # x and y coordinates of the text
                    label = round(info$median[1]/12,digits=2), size = 4)+ggplot2::annotate("text", 
                                                                                           x = info$median[2], y = 0.5, # x and y coordinates of the text
                                                                                           label = round(info$median[2]/12,digits=2), size = 4)
p
#=================

#surdata[[factor]]<-scale(surdata[[factor]], center=T,scale=T)
factor<-"Within_cluster_disp"
formu <- paste("Surv(time, event)",paste0("Nstage+pre_chemotherapy+Other_Metastasis+liver_M_distribution+liver_M_number+liver_M_size+radio_ablation+",factor), sep = "~")
formu <- as.formula(formu)
print(formu)
res.cox <- coxph(formu, data =  surdata)#+ Tstage +Nstage
c(factor,summary(res.cox)$coefficients[, 2][15],summary(res.cox)$coefficients[, 5][15])
ggforest(res.cox, data = surdata,main = paste0("Multivariate Hazard ratio - ",Suvival.type))


#========ROC==============

Sys.setlocale('LC_ALL','C')
library(survivalROC)
head(validation_surdata)
colnames(validation_surdata)
nobs <- NROW(validation_surdata)
cutoff <- 12*1
ROCdata= survivalROC(Stime=validation_surdata$time,##生存时间
                     status=validation_surdata$event,## 终止事件    
                     marker = validation_surdata$scores.image.CRS, ## marker value    
                     predict.time = cutoff,## 预测时间截点
                     method="NNE",span = 0.25*nobs^(-0.20))##span,NNE法的namda#span = 0.25*nobs^(-0.20)

## 绘图
plot(ROCdata$FP, ROCdata$TP, ## x=FP,y=TP
     type="l",col="red", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(ROCdata$AUC,3)), ##连接
     ylab="TP",
     main=" Method = NNE \n  Year = 5")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色



#=====table=========
# install.packages("ggpubr")
# library(ggpubr)
# 
# tbody.style = tbody_style(color = "black",
#                           fill = c("#e8f3de", "#d3e8bb"), hjust=1, x=0.9)
# ggtexttable(res.df, rows = NULL,
#             theme = ttheme(
#               colnames.style = colnames_style(color = "white", fill = "#8cc257"),
#               tbody.style = tbody.style
#             )
# )


#===survival======
survfit(res.cox)
ggsurvplot(survfit(res.cox),data = surdata, palette = "#2E9FDF",
           ggtheme = theme_minimal())

# surdata$y<-Surv(surdata$time,surdata$event==1)
# head(surdata)
# survdiff(y~STR_Proportion+LYM_Proportion+strata(Age)+strata(Gender)+strata(Tstage),data = surdata)

summary(surdata$infiltrated_STR)
interest_df <- with(surdata,
                    data.frame(infiltrated_STR = c(-0.5, 0.05), 
                               Age = rep(mean(Age, na.rm = TRUE), 2),
                               Gender = as.factor(c(1, 1)),
                               CRC_location= as.factor(c(1, 1)),
                               Tstage=as.factor(c(1,1)),
                               Nstage=as.factor(c(1,1)),
                               RAS_mutation=as.factor(c(1,1)),
                               pre_chemotherapy=as.factor(c(1,1)),
                               Other_Metastasis=as.factor(c(1,1)),
                               CEA=rep(mean(CEA, na.rm = TRUE), 2),
                               DFI=rep(mean(DFI, na.rm = TRUE), 2),
                               liver_M_distribution=as.factor(c(1,1)),
                               liver_M_number=rep(mean(liver_M_number, na.rm = TRUE), 2),
                               liver_M_size=rep(mean(liver_M_size, na.rm = TRUE), 2),
                               radio_ablation=as.factor(c(1,1)),
                               LYM_Proportion=rep(mean(LYM_Proportion, na.rm = TRUE), 2)
                               
                    )
)
interest_df
fit <- survfit(res.cox, newdata = interest_df)
ggsurvplot(fit, conf.int = TRUE,data = surdata,ggtheme = theme_minimal())

#===============correlation==============
library(GGally)
colnames(surdata)
# Create data 
data<-surdata[4:39]
# Check correlations (as scatterplots), distribution and print corrleation coefficient 
ggpairs(data, title=paste0("correlogram - ",Suvival.type) ) 


library(GGally)

# Nice visualization of correlations
ggcorr(data, method = c("everything", "spearman"),label = TRUE,
       label_alpha = TRUE)

library(RColorBrewer)
library(corrplot)
head(mtcars)
M<-cor(data)
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(data)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
corrplot(M, method="number")

#=====draft=========

#boost cox
# Generate some survival data with 10 informative covariates
library("CoxBoost")
n <- 200; p <- 100
beta <- c(rep(1,10),rep(0,p-10))
x <- matrix(rnorm(n*p),n,p)
real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
cens.time <- rexp(n,rate=1/10)
status <- ifelse(real.time <= cens.time,1,0)
obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)
# 10-fold cross-validation
cv.res <- cv.CoxBoost(time=obs.time,status=status,x=x,maxstepno=500,K=10,type="verweij",penalty=100)
# examine mean partial log-likelihood in the course of the boosting steps
plot(cv.res$mean.logplik)
# Fit with optimal number of boosting steps
cbfit <- CoxBoost(time=obs.time,status=status,x=x,stepno=cv.res$optimal.step,penalty=100)
summary(cbfit)

coef.CoxBoost(cbfit)
p1 <- estimPVal(cbfit,x,permute.n=10)


#---------draft--------------
#---factor.ref===========================
# surdata$STR_Proportion<-as.factor(surdata$STR_Proportion)
# surdata$STR_interaction<-as.factor(surdata$STR_interaction)
# surdata$LYM_Proportion<-as.factor(surdata$LYM_Proportion)
# surdata$LYM_interaction<-as.factor(surdata$LYM_interaction)
# surdata$Entropy<-as.factor(surdata$Entropy)
# surdata$LYM_interaction<-relevel(surdata$LYM_interaction, ref = "low")
# surdata$LYM_Proportion<-relevel(surdata$LYM_Proportion, ref = "low")
# surdata$STR_Proportion<-relevel(surdata$STR_Proportion, ref = "low")
# surdata$STR_interaction<-relevel(surdata$STR_interaction, ref = "low")
# surdata$Entropy<-relevel(surdata$Entropy, ref = "low")


#----- with post chemotherapy and no_che ----------------
surdata_interest<-surdata_interest[which(surdata_interest$Other_Metastasis=="0"),]
#surdata_interest<-surdata_interest[!is.na(surdata_interest$post_chemotherapy),]
che_surdata<-surdata_interest[which(surdata_interest$post_chemotherapy=="1"),]
nrow(che_surdata)
no_che_surdata<-surdata_interest[which(surdata_interest$post_chemotherapy=="0"),]
nrow(no_che_surdata)

#====CRS.chemo.survival=========

data_interest<-surdata_interest
CRS_data=data.frame(time=data_interest$time, event=data_interest$event,Group=data_interest$CRS_group)
CRS_data$Group<-ifelse(CRS_data$Group==0,"chemotherapy.CRS_Low","chemotherapy.CRS_High")
noc<-no_che_surdata[c("time","event")]
noc$Group<-rep("no_chemotherapy",nrow(no_che_surdata))
CRS_data<-rbind(CRS_data,noc)
fit <- surv_fit(Surv(time, event)~Group, data = CRS_data)
# summary(fit,times =5*12)
ggsurvplot(fit, data = CRS_data,ylab=paste0(Suvival.type,"  "),pval=TRUE,risk.table = TRUE,surv.median.line = "hv", conf.int = FALSE,title=paste0("CRS.Group"," related Survival"),xscale = c("m_y"),xlab="Time in years",pval.method=TRUE,risk.table.y.text = FALSE,break.time.by = 5*12)
surv_median(fit)