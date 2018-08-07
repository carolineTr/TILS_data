# ---------------------------------------------
# Programme: CIT_study_OS
# Auteur CT
# Description: Determination of patient overall 
# and disease-free survival
# based on the tRNA TILS score for the CIT dataset
# --------------------------------------------------------

# On récupère les données cliniques
my_clin=read.table("E-MTAB365-expression-updated srdf for AE.csv",sep=";",h=T)
my_clin$Characteristics.Delay.Overall_Survival.event=as.numeric(gsub("ND","NA",my_clin$Characteristics.Delay.Overall_Survival.event))

# ---------------------------------------
# Overall Survival calculation
# ---------------------------------------

censor_temp=my_clin$Characteristics.Delay.Overall_Survival.event
TTS_temp=my_clin$Characteristics.Delay.Overall_Survival.months.
censor=TTS=vector(le=nrow(my_clin))
for (i in 1:nrow(my_clin)){
if(is.na(censor_temp[i])==T | is.na(TTS_temp[i])==T){
      censor[i]=NA
      TTS[i]=NA
}else if(censor_temp[i] == 0){
   if(TTS_temp[i]<120){
      TTS[i]=TTS_temp[i]
      censor[i]=0
      }else{
      TTS[i]=120
      censor[i]=0
      }  
   }else if(censor_temp[i] == 1){
   if(TTS_temp[i]<120){
      TTS[i]=TTS_temp[i]
      censor[i]=1
      }else{
      TTS[i]=120
      censor[i]=0
      }  
   }
} 
my_clin$OS_event=censor
my_clin$OS_time=TTS

# -------------------------------------------------------------
# Signature construction
# -------------------------------------------------------------

# -------------------------------------------------------------
# Computation of the bzsc_score
my_frma=read.table("CIT_data_E-MTAB-365.txt",h=T)
require(frma)
bczsc=barcode(as.matrix(my_frma),out="z-score",platform="GPL570")

combi_lympho=read.table("combi_lympho.txt", h=T)[,1]
combi_myelo=read.table("combi_myelo.txt", h=T)[,1]
combi_str=read.table("combi_str.txt", h=T)[,1]
combi_neg=read.table("combi_neg.txt", h=T)[,1]

tm_lympho=colMeans(bczsc[which(rownames(bczsc)%in%(combi_lympho)),])
tm_myelo=colMeans(bczsc[which(rownames(bczsc)%in%(combi_myelo)),])
tm_str=colMeans(bczsc[which(rownames(bczsc)%in%(combi_str)),])
tm_neg=colMeans(bczsc[which(rownames(bczsc)%in%(combi_neg)),])# On a q 4 marqueurs sur 10 
mon_indice= (tm_lympho+tm_myelo)/(tm_lympho+tm_myelo+tm_str+tm_neg)
my_var=scale(mon_indice,sc=T,ce=T )


# -------------------------------                                          
# Univariate Cox models using the
# tRNA TILS score  
# -------------------------------
# In the whole population
surv_var=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~my_var[-which(my_clin$pam=="Normal")])
exp(confint(surv_var))

# --------------------
# Without LumA
surv_var=coxph(Surv(my_clin_norm$OS_time[-which(my_clin_norm$pam=="LumA")],my_clin_norm$OS_event[-which(my_clin_norm$pam=="LumA")])~my_var[-which(my_clin$pam=="Normal")][-which(my_clin_norm$pam=="LumA")])
exp(confint(surv_var))

# --------------------
# In Luminal A patients
surv_lumA=coxph(Surv(my_clin[which(pam50=="LumA"),]$OS_time,my_clin[which(pam50=="LumA"),]$OS_event)~my_var[which(pam50=="LumA")]])
exp(confint(surv_lumA))

# --------------------
# In Luminal B patients
surv_lumB=coxph(Surv(my_clin[which(pam50=="LumB"),]$OS_time,my_clin[which(pam50=="LumB"),]$OS_event)~my_var[which(pam50=="LumB")]])
exp(confint(surv_lumB)

# --------------------
# In Her2 patients
surv_Her2=coxph(Surv(my_clin[which(pam50=="Her2"),]$OS_time,my_clin[which(pam50=="Her2"),]$OS_event)~my_var[which(pam50=="Her2")])
exp(confint(surv_Her2))

# --------------------
# In basal-like patients
surv_TN=coxph(Surv(my_clin[which(pam50=="Basal"),]$OS_time,my_clin[which(pam50=="Basal"),]$OS_event)~my_var[which(pam50=="Basal")])
exp(confint(surv_TN))

# ------------------------------------------
# Clinical variables
# ------------------------------------------

my_clin_norm=my_clin[-which(pam50=="Normal"),]
my_clin_norm$pam=factor(my_clin_norm$pam,levels=c("Basal",  "Her2",   "LumA",  "LumB"))

my_clin_norm$Lymph=as.numeric(gsub("D","NA",substr(as.factor(my_clin_norm$Characteristics..TNM..N.),2,2)))
surv_lymph=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~as.factor(my_clin_norm$Lymph))
exp(confint(surv_lymph))

my_clin_norm$Grade=my_clin_norm$Characteristics..Grade..Scarff.Bloom.Richardson.
grade=ifelse(my_clin_norm$Grade==3,2,1)
surv_grade=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~as.factor(grade))
exp(confint(surv_grade))

er=my_clin_norm$Characteristics..ESR1..Protein.expressio
er[which(as.factor(my_clin_norm$Characteristics..ESR1..Protein.expression.)=="ND")]=NA
my_clin_norm$er=factor(er,levels=c("N","Y"))
surv_er=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~as.factor(my_clin_norm$er))
exp(confint(surv_er))

my_clin_norm$Size=as.numeric(gsub("D","NA",substr(as.factor(my_clin_norm$Characteristics..TNM..T.),2,2)))
my_clin_norm$Size[which(my_clin_norm$Size>3)]=3
my_clin_norm$Size[which(my_clin_norm$Size==0)]=1
size=ifelse(my_clin_norm$Size==1,1,2)

surv_size=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~as.factor(size))
exp(confint(surv_size))

surv_pam=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~relevel(as.factor(my_clin_norm$pam),ref="LumA"))
exp(confint(surv_pam))


surv_age_fac=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~my_clin_norm$age_fac)
exp(confint(surv_age_fac))

erbb2_norm=factor(my_clin_norm$Characteristics..ERBB2.)
erbb2_norm=gsub("ND",NA,erbb2_norm)

# -------------------------------------
# Multivariate Cox model including PAM50
attach(my_clin_norm)

surv_all=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~relevel(as.factor(pam),ref="LumA")+
    as.factor(age_fac)+as.factor(size)+as.factor(grade)+as.factor(Lymph)+my_clin_norm$er+my_var[-which(my_clin$pam=="Normal")],data=my_clin_norm)
exp(confint(surv_all))


# -------------------------------------------
# Standard pathology classification
# -------------------------------------------
erbb22=my_clin_norm$Characteristics..ERBB2
erbb22[which(as.factor(my_clin_norm$Characteristics..ERBB2)=="ND")]=NA
my_clin_norm$erbb2=factor(erbb22,levels=c("N","Y"))
anapath=vector(le=nrow(my_clin_norm))
anapath[which(my_clin_norm$er=="N" & my_clin_norm$erbb2=="N")]="TN"
anapath[which(my_clin_norm$er=="Y" & my_clin_norm$erbb2=="N")]="Luminal"
anapath[which(my_clin_norm$erbb2=="Y" & (!is.na(my_clin_norm$er)))]="Her2+"
anapath=gsub(FALSE,NA,anapath)


surv_anapath=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~relevel(as.factor(anapath),ref="Luminal"))
exp(confint(surv_anapath))

# -------------------------------
# Multivariate model 
surv_all_score=coxph(Surv(my_clin_norm$OS_time,my_clin_norm$OS_event)~relevel(as.factor(anapath),ref="Luminal")+
    as.factor(age_fac)+as.factor(size)+as.factor(grade)+as.factor(Lymph)+my_var[-which(my_clin$pam=="Normal")],data=my_clin_norm)
exp(confint(surv_all_score))


# ----------------------
# Luminal patients
surv_lum=coxph(Surv(my_clin_norm$OS_time[which(anapath=="Luminal")],my_clin_norm$OS_event[which(anapath=="Luminal")])~my_var[-which(my_clin$pam=="Normal")][which(anapath=="Luminal")])
exp(confint(surv_lum))

# --------------------
# Her2 patients
surv_her2=coxph(Surv(my_clin_norm$OS_time[which(anapath=="Her2+")],my_clin_norm$OS_event[which(anapath=="Her2+")])~my_var[-which(my_clin$pam=="Normal")][which(anapath=="Her2+")])
exp(confint(surv_her2))

# --------------------
# Triple negative patients
surv_TN=coxph(Surv(my_clin_norm$OS_time[which(anapath=="TN")],my_clin_norm$OS_event[which(anapath=="TN")])~my_var[which(anapath=="TN")])
exp(confint(surv_TN))

# --------------------
# Without Luminal
surv_var=coxph(Surv(my_clin_norm$OS_time[-which(anapath=="Luminal")],my_clin_norm$OS_event[-which(anapath=="Luminal")])~my_var[-which(my_clin$pam=="Normal")][-which(anapath=="Luminal")])
exp(confint(surv_var))

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# The same program can be used for disease free survival by replacing 
# OS_time and OS_event by DFS_time and DFS_event. 
# Diseas-free survival was computed as follows: 

my_clin$Characteristics.Delay.Metastasis_Free_Survival.even=as.numeric(gsub("ND","NA",my_clin$Characteristics.Delay.Metastasis_Free_Survival.even))

# ---------------------------------------
# Disease free survival calculation
# ---------------------------------------

censor_temp=my_clin$Characteristics.Delay.Metastasis_Free_Survival.event
TTS_temp=my_clin$Characteristics.Delay.Metastasis_Free_Survival.months.
censor=TTS=vector(le=nrow(my_clin))
for (i in 1:nrow(my_clin)){
if(is.na(censor_temp[i])==T | is.na(TTS_temp[i])==T){
      censor[i]=NA
      TTS[i]=NA
}else if(censor_temp[i] == 0){
   if(TTS_temp[i]<120){
      TTS[i]=TTS_temp[i]
      censor[i]=0
      }else{
      TTS[i]=120
      censor[i]=0
      }  
   }else if(censor_temp[i] == 1){
   if(TTS_temp[i]<120){
      TTS[i]=TTS_temp[i]
      censor[i]=1
      }else{
      TTS[i]=120
      censor[i]=0
      }  
   }
} 
my_clin$DFS_event=censor
my_clin$DFS_time=TTS

  