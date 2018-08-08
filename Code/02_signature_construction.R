# ---------------------------------------------
# Programme: 02_signature_construction 
# Auteur: CT
# Description: Construction of the tRNA TILS signature
# --------------------------------------------------------

# --------------------------------------
# Chargement des données et librairies
library(frma)
library(affy)
load("allGSM_frma.Rdata")   #my_frma
load("allGSM_bcbin.Rdata")  # bcbin
load("allGSM_bczsc.Rdata")  # bczsc

liste_gse=read.table("Add_File3.csv",h=T,sep=";")
load(file="info_clin.Rdata")

load(file="my_index.Rdata")   

# -------------------------------------------------------------
# Selection of probes to GPL570 et GPL96
# -------------------------------------------------------------
gpl570_annot=read.table("GPL570_annot.csv",h=T,sep=";") 
gpl96_annot=read.table("GPL96_annot.csv",h=T,sep=";")
whichat=grep('[1234567890]_at',gpl570_annot$ID )
temp=gpl570_annot[whichat,]
garde=temp$ID[which(temp$Gene!="" & temp$Gene%in%gpl96_annot$Gene.symbol)]

exprs_bczsc=bczsc[which(rownames(bczsc)%in%garde),]
exprs_bcbin=bcbin[which(rownames(bcbin)%in%garde),]
write.table(garde,"GPL570_GPL96_communs.txt") 

# --------------------------
# Use of binary information
cpte_0_classif=apply(exprs_bcbin,1,function(x) tapply(x,my_classif,sum,na.rm=T))

lympho=which(cpte_0_classif[1,]>0.8 & cpte_0_classif[2,]<0.2 & cpte_0_classif[3,]<0.2 & cpte_0_classif[4,]<0.2)
myelo=which(cpte_0_classif[1,]<0.2 & cpte_0_classif[2,]>0.8 & cpte_0_classif[3,]<0.2 & cpte_0_classif[4,]<0.2)
negctrl=which(cpte_0_classif[1,]<0.2 & cpte_0_classif[2,]<0.2 & cpte_0_classif[3,]>0.8 & cpte_0_classif[4,]<0.2)
strom=which(cpte_0_classif[1,]<0.2 & cpte_0_classif[2,]<0.2 & cpte_0_classif[3,]<0.2 & cpte_0_classif[4,]>0.8)


# -------------------------------------------------------------------------------
# Limma differential test
# ----------------------------------------------

require(limma)
f <- my_classif
design <- model.matrix(~0+f)
colnames(design) <- c("lymphoids", "myeloid",  "cancer",  "stroma")

fit <- lmFit(exprs_bczsc, design)

# Version 2, 1 condition vs les autres
contrast.matrix <- makeContrasts("LvsOther"=lymphoids-(myeloid+cancer+stroma)/3,
                                 "MvsOther"=myeloid-(lymphoids+cancer+stroma)/3,
                                 "NvsOther"=cancer-(lymphoids+myeloid+stroma)/3,
                                 "SvsOther"=stroma-(lymphoids+myeloid+cancer)/3,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


top_limma=topTableF(fit2, number=30, adjust="BH")
top_lympho=topTable(fit2, coef="LvsOther",number=nrow(exprs_bczsc), adjust="BH")
top_myelo=topTable(fit2, coef="MvsOther",number=nrow(exprs_bczsc), adjust="BH")
top_neg=topTable(fit2, coef="NvsOther",number=nrow(exprs_bczsc), adjust="BH")
top_str=topTable(fit2, coef="SvsOther",number=nrow(exprs_bczsc), adjust="BH")


combi_lympho=top_lympho[which(rownames(top_lympho)%in%names(lympho)),][1:10,]
combi_myelo=top_myelo[which(rownames(top_myelo)%in%names(myelo)),][1:10,]
combi_neg=top_neg[which(rownames(top_neg)%in%names(negctrl)),][1:10,]
combi_str=top_str[which(rownames(top_str)%in%names(strom)),][1:10,] 

liste_combi=c(rownames(combi_lympho),rownames(combi_myelo),rownames(combi_neg),rownames(combi_str))#,rownames(combi_svsneg) )
combi_exprs=exprs_bczsc[which(rownames(exprs_bczsc)%in%liste_combi),]
res.pca <- PCA(t(combi_exprs), graph=FALSE)
hc2 <- HCPC(res.pca, nb.clust=4,method="ward",graph=F)
table(hc2$data.clust[,ncol(hc2$data.clust)], as.factor(my_classif))



library(ade4)
my_pca=dudi.pca(t(combi_exprs))
s.class(my_pca$li,fac=as.factor(my_classif_papier),col=rainbow(4))
s.class(my_pca$li,fac=as.factor(my_classif_papier))


write.table(rownames(combi_lympho),"combi_lympho.txt", row=F)
write.table(rownames(combi_myelo),"combi_myelo.txt", row=F)
write.table(rownames(combi_str),"combi_str.txt", row=F)
write.table(rownames(combi_neg),"combi_neg.txt", row=F)
write.table(rownames(combi_svsneg),"combi_svsneg.txt", row=F)

