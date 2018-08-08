# ---------------------------------------------
# Programme: 01_normalization.R
# Auteur CT
# Description: Normalization of GSM files using frma
# The example code is given for the GSM files
# used to construct the tRNA TILS score signature
# --------------------------------------------------------

library(frma)
library(affy)
library("hgu133plus2frmavecs")
library("hgu133plus2hsensgprobe")


# ----------------------------------------------------------------
# Cel files normalization uzing frma
# ----------------------------------------------------------------
celfiles <- list.files(patt=".CEL", full = TRUE,rec=T)

for (i in celfiles){
abatch <- ReadAffy(filenames=i,compress=T)  
save(abatch,file=paste(dirname(i),"/",basename(dirname(i)),"_affy.Rdata",sep=""))
my_frma=frma(abatch, background="rma", normalize="quantile",target="probeset",
    input.vecs=NULL, output.param=NULL, verbose=FALSE) 
save(my_frma,file=paste(dirname(i),"/",basename(dirname(i)),"_frma.Rdata",sep=""))   
}

# ----------------------------------------------------------
# Construction of a unique data file
# ----------------------------------------------------------

setwd("/fhgfs/data/work/ptbc/ctruntze/TILS/Database")
frmafiles <- list.files(patt="frma", full = TRUE,rec=T)
load(frmafiles[1])
frma_all=my_frma
for (i in frmafiles[-1]){
load(file=i)   
frma_all=combine(frma_all,my_frma)
rm(my_frma)
}
save(frma_all,file="allGSM_frma.Rdata")   

# ----------------------------------------------------------
# Generation of the barcode information
# ----------------------------------------------------------
load(file="allGSM_frma.Rdata")   
bcbin=barcode(frma_all)
bczsc=barcode(frma_all,out="z-score")
save(bcbin,file="allGSM_bcbin.Rdata")   
save(bczsc,file="allGSM_bczsc.Rdata")
   
# --------------------------------------------------
# Construction of the immune cell types information 
# --------------------------------------------------

load(file="allGSM_frma.Rdata")   #frma_all
liste_gse=read.table("Add_File3.csv",h=T,sep=";")

my_index=NULL
for (i in colnames(frma_all)){
 my_gsm_temp=gsub("\\.","_",i)
 my_gsm=strsplit(my_gsm_temp,"_")[[1]][1]
 my_index=c(my_index,which(liste_gse[,1]==my_gsm))
}

save(my_index,file="my_index.Rdata")   
info_clin=liste_gse[my_index,]
save(info_clin,file="info_clin.Rdata")


