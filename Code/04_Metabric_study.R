# ---------------------------------------------
# Programme: Metabric_study.R
# Auteur: CT
# Description: Determination of patient prognosis based 
# on the tRNA TILS score for the Metabric dataset    
# --------------------------------------------------------

my_clin=read.delim("data_clinical_supp_patient.txt",h=T)


# ---------------------------------------
# Survival calculation
# ---------------------------------------
censor_temp=ifelse(my_clin$VITAL_STATUS=="Died of Disease",1,0)
TTS_temp=my_clin$OS_MONTHS

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
# COMPLETE DATASET
# -------------------------------------------------------------


my_data_norm_disc=read.table("_ega-box-04_discovery_ExpressionMatrix.txt",h=T)
my_data_norm_valid=read.table("_ega-box-04_validation_ExpressionMatrix.txt",h=T)
my_data_norm=cbind(my_data_norm_disc,my_data_norm_valid)


# ------------------------------------------------------
# Correspondance of clinical and gene expression information
# ------------------------------------------------------
my_index_pat=vector(le=nrow(my_clin))
nom_temp=gsub("\\.","-",colnames(my_data_norm))
retrieve=NULL
for (i in 1:length(nom_temp)){
if (length(which(my_clin$PATIENT_ID==nom_temp[i]))>0){
my_index_pat[i]=which(my_clin$PATIENT_ID==nom_temp[i])
print(which(my_clin$PATIENT_ID==nom_temp[i]))
}else{
retrieve=c(retrieve,i)
}
}

my_clin_order=my_clin[my_index_pat,]
my_data_norm_temp=my_data_norm
my_data_norm=my_data_norm[,-retrieve]

# ---------------------------------------------------------
# Signature construction

liste_HGNC=read.delim("corres_HGU133p2_HGNC_selection.txt",sep="\t",h=T)
my_HGNC=liste_HGNC[,2]
my_index=NULL
for (i in my_HGNC){
my_index=c(my_index,which(rownames(my_data_norm)==i))
}

my_rnaseq_selec=t(my_data_norm)[,my_index]

my_probe=NULL
for (i in colnames(my_rnaseq_selec)){
my_probe=c(my_probe,as.character(liste_HGNC[which(liste_HGNC[,2]==i)[1],1]))
}

combi_lympho=read.table("combi_lympho.txt", h=T)[,1]
combi_myelo=read.table("combi_myelo.txt", h=T)[,1]
combi_str=read.table("combi_str.txt", h=T)[,1]
combi_neg=read.table("combi_neg.txt", h=T)[,1]
my_selec=c(as.character(combi_lympho),as.character(combi_myelo),as.character(combi_str),as.character(combi_neg))

tm_lympho=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_lympho))])
tm_myelo=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_myelo))])
tm_str=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_str))])
tm_neg=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_neg))])
tm_matrix=cbind(tm_lympho,tm_myelo,tm_str,tm_neg)
tm_matrix_old=tm_matrix
my_rnaseq_selec_noNorm=my_rnaseq_selec[-which(save_claudin=="Normal"),]
rown=gsub("\\.","-",rownames(my_rnaseq_selec_noNorm))
tm_matrix=tm_matrix[-which(save_claudin=="Normal"),]

mon_indice= (tm_matrix[,1]+tm_matrix[,2])/(tm_matrix[,1]+tm_matrix[,2]+tm_matrix[,3]+tm_matrix[,4])
my_var=scale(mon_indice,sc=T,cent=T)


# -------------------------------                                          
# Univariate Cox models using the
# tRNA TILS score  
# -------------------------------
# In the whole population
surv_Score=coxph(Surv(my_clin_order$OS_time, my_clin_order$OS_event)~my_var)
exp(confint(surv_Score))

# Without Luminal A
surv_all=coxph(Surv(my_clin_order$OS_time[-which(clin_pam50=="LumA")],my_clin_order$OS_event[-which(clin_pam50=="LumA")])~my_var[-which(clin_pam50=="LumA")])
exp(confint(surv_all))

# In Luminal A patients
surv_lumA=coxph(Surv(my_clin_order$OS_time[which(clin_pam50=="LumA")],my_clin_order$OS_event[which(clin_pam50=="LumA")])~my_var[which(clin_pam50=="LumA")])
exp(confint(surv_lumA))

# In Luminal B patients
surv_lumB=coxph(Surv(my_clin_order$OS_time[which(clin_pam50=="LumB")],my_clin_order$OS_event[which(clin_pam50=="LumB")])~my_var[which(clin_pam50=="LumB")])
exp(confint(surv_lumB))

# In Her2 patients
surv_her2=coxph(Surv(my_clin_order$OS_time[which(clin_pam50=="Her2")],my_clin_order$OS_event[which(clin_pam50=="Her2")])~my_var[which(clin_pam50=="Her2")])
exp(confint(surv_her2))

# In basal-like patients
surv_BL=coxph(Surv(my_clin_order$OS_time[which(clin_pam50=="TN")],my_clin_order$OS_event[which(clin_pam50=="TN")])~my_var[which(clin_pam50=="TN")])
exp(confint(surv_BL))

# ------------------------------------------
# Clinical variables
# ------------------------------------------
supp_sample=read.table("data_clinical_supp_sample.txt",h=T)
indic_supp=NULL
for (j in 1:nrow(my_clin_order)){
indic_supp=c(indic_supp,which(supp_sample$SAMPLE_ID==my_clin_order$PATIENT_ID[j]))
}
supp_sample_order=supp_sample[indic_supp,]

# Age
age_fac=my_clin_order$AGE>65
surv_age=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~as.factor(age_fac))
exp(confint(surv_age)) 

# PAM50
clin_pam50=factor(my_clin_order$CLAUDIN_SUBTYPE[-which(my_clin_order$CLAUDIN_SUBTYPE=="Normal" | my_clin_order$CLAUDIN_SUBTYPE=="NC")],levels=c("Basal", "claudin-low","Her2","LumA","LumB"))
clin_pam50=gsub("Basal","TN",clin_pam50)
clin_pam50=gsub("claudin-low","TN",clin_pam50)
clin_pam50=factor(clin_pam50,levels=c("TN","Her2","LumA","LumB"))

surv_pam50=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~relevel(as.factor(clin_pam50),ref="LumA"))
exp(confint(surv_pam50)) 

# Estrogen Receptor
surv_er=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~as.factor(supp_sample_order$ER_STATUS))
exp(confint(surv_er)) 

# SBR Grade
supp_sample_order$GRADE=factor(supp_sample_order$GRADE,levels=c( "1", "2", "3"))
grade=ifelse(supp_sample_order$GRADE==3,2,1)
surv_grade=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~as.factor(grade))
exp(confint(surv_grade))

# Tumor size
size_temp=(as.numeric(as.character(supp_sample_order$TUMOR_SIZE)))
size=ifelse(size_temp>20,2,1)
surv_size=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~as.factor(size))

# Nodal status
temp_stade=as.numeric(as.character(supp_sample_order$TUMOR_STAGE))
supp_sample_order$Node=as.factor(ifelse(temp_stade>=3,1,0))
surv_node=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~as.factor(supp_sample_order$Node))
exp(confint(surv_node))

# -------------------------------------
# Multivariate Cox model including PAM50

surv_all_score=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~relevel(as.factor(clin_pam50),ref="LumA")+
    as.factor(age_fac)+as.factor(size)+as.factor(grade)+as.factor(supp_sample_order$ER_STATUS)+supp_sample_order$Node+my_var,data=my_clin_order)
exp(confint(surv_all_score))

# -------------------------------------------
# Standard pathology classification
# -------------------------------------------

# Standard pathology classification status
anapath=vector(le=length(supp_sample_order$ER_STATUS))
anapath[which(supp_sample_order$ER_STATUS=="-" & supp_sample_order$HER2_STATUS=="-")]="TN"
anapath[which(supp_sample_order$ER_STATUS=="+" & supp_sample_order$HER2_STATUS=="-")]="Luminal"
anapath[which(supp_sample_order$HER2_STATUS=="+" & !is.na(supp_sample_order$ER_STATUS))]="Her2+"

surv_anapath=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~relevel(as.factor(anapath),ref="Luminal"))
exp(confint(surv_anapath))

# tRNA TILS score in the whole population and without Luminal A
surv_varanapath2=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~my_var)
surv_varanapath1=coxph(Surv(my_clin_order$OS_time[-which(anapath=="Luminal")],my_clin_order$OS_event[-which(anapath=="Luminal")])~my_var[-which(anapath=="Luminal")])

surv_lum=coxph(Surv(my_clin_order$OS_time[which(anapath=="Luminal")],my_clin_order$OS_event[which(anapath=="Luminal")])~my_var[which(anapath=="Luminal")])
exp(confint(surv_lum))

surv_her2=coxph(Surv(my_clin_order$OS_time[which(anapath=="Her2+")],my_clin_order$OS_event[which(anapath=="Her2+")])~my_var[which(anapath=="Her2+")])
exp(confint(surv_her2))

surv_TN=coxph(Surv(my_clin_order$OS_time[which(anapath=="TN")],my_clin_order$OS_event[which(anapath=="TN")])~my_var[which(anapath=="TN")])
exp(confint(surv_TN))

# --------------------------------------
# Multivariate model including 
# Standard pathology classification

surv_all_score=coxph(Surv(my_clin_order$OS_time,my_clin_order$OS_event)~relevel(as.factor(anapath),ref="Luminal")+
    as.factor(age_fac)+as.factor(size)+as.factor(grade)+supp_sample_order$Node+my_var,data=my_clin_order)
exp(confint(surv_all_score))



# --------------------------------------------------------------
# --------------------------------------------------------------
# VALIDATION OF THE ASSOCIATION BETWEEN THE tRNA TILS SCORE AND
# SURVIVAL FOR BASAL-LIKE PATIENTS
# --------------------------------------------------------------
# --------------------------------------------------------------

# ------------------------------------------------
# Estimation of the survival model using
# the discovery dataset
# -----------------------------------------------
my_data_norm_disc=read.table("_ega-box-04_discovery_ExpressionMatrix.txt",h=T)
my_data_norm=my_data_norm_disc

# ------------------------------------------------------
# Correspondance of clinical and gene expression information
my_index_pat=vector(le=nrow(my_clin))
nom_temp=gsub("\\.","-",colnames(my_data_norm))
retrieve=NULL
for (i in 1:length(nom_temp)){
if (length(which(my_clin$PATIENT_ID==nom_temp[i]))>0){
my_index_pat[i]=which(my_clin$PATIENT_ID==nom_temp[i])
print(which(my_clin$PATIENT_ID==nom_temp[i]))
}else{
retrieve=c(retrieve,i)
}
}
my_clin_order=my_clin[my_index_pat,]
my_data_norm_temp=my_data_norm

# ---------------------------------------------------------
# Signature construction
my_index=NULL
for (i in my_HGNC){
my_index=c(my_index,which(rownames(my_data_norm)==i))
}


my_rnaseq_selec=t(my_data_norm)[,my_index]
my_probe=NULL
for (i in colnames(my_rnaseq_selec)){
my_probe=c(my_probe,as.character(liste_HGNC[which(liste_HGNC[,2]==i)[1],1]))
}


combi_lympho=read.table("combi_lympho.txt", h=T)[,1]
combi_myelo=read.table("combi_myelo.txt", h=T)[,1]
combi_str=read.table("combi_str.txt", h=T)[,1]
combi_neg=read.table("combi_neg.txt", h=T)[,1]
my_selec=c(as.character(combi_lympho),as.character(combi_myelo),as.character(combi_str),as.character(combi_neg))

tm_lympho=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_lympho))])
tm_myelo=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_myelo))])
tm_str=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_str))])
tm_neg=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_neg))])
tm_matrix=cbind(tm_lympho,tm_myelo,tm_str,tm_neg)
mon_indice= (tm_lympho+tm_myelo)/(tm_lympho+tm_myelo+tm_str+tm_neg)
my_var=scale(mon_indice,sc=T,cent=T )


# -----------------------------------------------------
# Survival Cox model for basal-like patients 
# - discovery dataset

# Median cutoff: 
surv_TN=coxph(Surv(my_clin_order$OS_time[which(clin_pam50=="TN")],my_clin_order$OS_event[which(clin_pam50=="TN")])~as.factor(ifelse(my_var[which(clin_pam50=="TN")]>median(my_var[which(clin_pam50=="TN")]),2,1)))
median(my_var[which(clin_pam50=="TN")])
#0.7398921


# Kaplan-Meier Curves
source("KaplanMeier_ggplot.R")
km_TN <- survfit(Surv(my_clin_order$OS_time[which(clin_pam50=="TN")],my_clin_order$OS_event[which(clin_pam50=="TN")])~as.factor(ifelse(my_var[which(clin_pam50=="TN")]>median(my_var[which(clin_pam50=="TN")]),2,1)))
ggkm(km_TN, timeby=12, ystratalabs=c("Low tRNA TILS score", "High tRNA TILS score"),
main="Kaplan-Meier curves for Basal-like tumors  - discovery dataset",ystrataname="Groups",pval=T)

# ------------------------------------------------
# Validation of the survival model using
# the validation dataset
# -----------------------------------------------
my_data_norm_valid=read.table("_ega-box-04_validation_ExpressionMatrix.txt",h=T)
my_data_norm=my_data_norm_valid


my_index_pat=vector(le=nrow(my_clin))
nom_temp=gsub("\\.","-",colnames(my_data_norm))
retrieve=NULL
for (i in 1:length(nom_temp)){
if (length(which(my_clin$PATIENT_ID==nom_temp[i]))>0){
my_index_pat[i]=which(my_clin$PATIENT_ID==nom_temp[i])
print(which(my_clin$PATIENT_ID==nom_temp[i]))
}else{
retrieve=c(retrieve,i)
}
}

my_clin_order=my_clin[my_index_pat,]
my_data_norm_temp=my_data_norm
my_data_norm=my_data_norm[,-retrieve]

my_index=NULL
for (i in my_HGNC){
my_index=c(my_index,which(rownames(my_data_norm)==i))
}


my_rnaseq_selec=t(my_data_norm)[,my_index]
my_probe=NULL
for (i in colnames(my_rnaseq_selec)){
my_probe=c(my_probe,as.character(liste_HGNC[which(liste_HGNC[,2]==i)[1],1]))
}


combi_lympho=read.table("combi_lympho.txt", h=T)[,1]
combi_myelo=read.table("combi_myelo.txt", h=T)[,1]
combi_str=read.table("combi_str.txt", h=T)[,1]
combi_neg=read.table("combi_neg.txt", h=T)[,1]
my_selec=c(as.character(combi_lympho),as.character(combi_myelo),as.character(combi_str),as.character(combi_neg))

tm_lympho=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_lympho))])
tm_myelo=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_myelo))])
tm_str=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_str))])
tm_neg=rowMeans(my_rnaseq_selec[,which(my_probe%in%(combi_neg))])
tm_matrix=cbind(tm_lympho,tm_myelo,tm_str,tm_neg)
mon_indice= (tm_lympho+tm_myelo)/(tm_lympho+tm_myelo+tm_str+tm_neg)
my_var=scale(mon_indice,sc=T,cent=T )


# ---------------------------------------------
# Using the median computed on the
# discovery dataset (med=0.7398921)
surv_TN=coxph(Surv(my_clin_order$OS_time[which(clin_pam50=="TN")],my_clin_order$OS_event[which(clin_pam50=="TN")])~as.factor(ifelse(my_var[which(clin_pam50=="TN")]>0.7398921,2,1)))

# Kaplan Meier curves
source("KaplanMeier_ggplot.R")
km_TN <- survfit(Surv(my_clin_order$OS_time[which(clin_pam50=="TN")],my_clin_order$OS_event[which(clin_pam50=="TN")])~as.factor(ifelse(my_var[which(clin_pam50=="TN")]>0.7398921,2,1)))
ggkm(km_TN, timeby=12, ystratalabs=c("Low tRNA TILS score", "High tRNA TILS score"),
main="Kaplan-Meier curves for Basal-like tumors - validation dataset",ystrataname="Groups",pval=T)

