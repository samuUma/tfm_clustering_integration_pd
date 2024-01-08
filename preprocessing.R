# set wd
setwd("/home/eidrian/AMP_PD_Data_v3_2022")

# load vroom library to read files very fast
library(vroom)
library(FactoMineR)

clinic=vroom("./AMP PD Data v3 2022/clinical_metadata/releases_2022_v3release_1115_amp_pd_case_control.csv")
# patients with idiopathic parkinson
pacients=clinic[which(clinic$diagnosis_at_baseline=="Idiopathic PD" | clinic$diagnosis_at_baseline=="Parkinson's Disease" ),]
#--------------------------------------------------------------------------------
# TRANSCRIPTOME
# read counts 
wb=vroom("./AMP PD Data v3 2022/transcriptomics/rnaseq-WB-RWTS/subread_feature-counts/releases_2022_v3release_1115_rnaseq-WB-RWTS_subread_feature-counts_matrix.featureCounts.tsv")
# create dataframe
wb=as.data.frame(wb)
# rename rows
rownames(wb)=wb[,1]
# remove gene ids column
wb=wb[,-1]
# get baseline measures
wb=wb[,grep("BLM0T1",colnames(wb))]
# PROTEOME 
# PLA D01 SAMPLES
#PLA_PPEA_D01_neurology
PLA_PPEA_D01_neurology=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D01/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D01_olink-explore_protein-expression_PLA-PPEA-D01_matrix_neurology.csv")
#PLA_PPEA_D01_neurology=matrix(PLA_PPEA_D01_neurology)
PLA_PPEA_D01_neurology=as.data.frame(PLA_PPEA_D01_neurology)
# rename rows
rownames(PLA_PPEA_D01_neurology)=PLA_PPEA_D01_neurology[,1]
# remove uniprot id column
PLA_PPEA_D01_neurology=PLA_PPEA_D01_neurology[,-1]
PLA_PPEA_D01_neurology=PLA_PPEA_D01_neurology[,grep("BLM0T1",colnames(PLA_PPEA_D01_neurology))]
#------------------------------------------------------------------------------------
#PLA_PPEA_D01_cardiometabolic
PLA_PPEA_D01_cardiometabolic=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D01/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D01_olink-explore_protein-expression_PLA-PPEA-D01_matrix_cardiometabolic.csv")
#PLA_PPEA_D01_cardiometabolic=matrix(PLA_PPEA_D01_cardiometabolic)
PLA_PPEA_D01_cardiometabolic=as.data.frame(PLA_PPEA_D01_cardiometabolic)
# rename rows
rownames(PLA_PPEA_D01_cardiometabolic)=PLA_PPEA_D01_cardiometabolic[,1]
# remove uniprot id column
PLA_PPEA_D01_cardiometabolic=PLA_PPEA_D01_cardiometabolic[,-1]
PLA_PPEA_D01_cardiometabolic=PLA_PPEA_D01_cardiometabolic[,grep("BLM0T1",colnames(PLA_PPEA_D01_cardiometabolic))]
#------------------------------------------------------------------------------------
#PLA_PPEA_D01_inflammation
PLA_PPEA_D01_inflammation=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D01/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D01_olink-explore_protein-expression_PLA-PPEA-D01_matrix_inflammation.csv")
#PLA_PPEA_D01_inflammation=matrix(PLA_PPEA_D01_inflammation)
PLA_PPEA_D01_inflammation=as.data.frame(PLA_PPEA_D01_inflammation)
# rename rows
rownames(PLA_PPEA_D01_inflammation)=PLA_PPEA_D01_inflammation[,1]
# remove uniprot id column
PLA_PPEA_D01_inflammation=PLA_PPEA_D01_inflammation[,-1]
PLA_PPEA_D01_inflammation=PLA_PPEA_D01_inflammation[,grep("BLM0T1",colnames(PLA_PPEA_D01_inflammation))]
#------------------------------------------------------------------------------------
#PLA_PPEA_D01_oncology
PLA_PPEA_D01_oncology=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D01/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D01_olink-explore_protein-expression_PLA-PPEA-D01_matrix_oncology.csv")
#PLA_PPEA_D01_oncology=matrix(PLA_PPEA_D01_oncology)
PLA_PPEA_D01_oncology=as.data.frame(PLA_PPEA_D01_oncology)
# rename rows
rownames(PLA_PPEA_D01_oncology)=PLA_PPEA_D01_oncology[,1]
# remove uniprot id column
PLA_PPEA_D01_oncology=PLA_PPEA_D01_oncology[,-1]
PLA_PPEA_D01_oncology=PLA_PPEA_D01_oncology[,grep("BLM0T1",colnames(PLA_PPEA_D01_oncology))]
#--------------------------------------------------------------------------------------
# PLA D02 SAMPLES
#PLA_PPEA_D02_neurology
PLA_PPEA_D02_neurology=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D02/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D02_olink-explore_protein-expression_PLA-PPEA-D02_matrix_neurology.csv")
#PLA_PPEA_D02_neurology=matrix(PLA_PPEA_D02_neurology)
PLA_PPEA_D02_neurology=as.data.frame(PLA_PPEA_D02_neurology)
# rename rows
rownames(PLA_PPEA_D02_neurology)=PLA_PPEA_D02_neurology[,1]
# remove uniprot id column
PLA_PPEA_D02_neurology=PLA_PPEA_D02_neurology[,-1]
PLA_PPEA_D02_neurology=PLA_PPEA_D02_neurology[,grep("BLM0T1",colnames(PLA_PPEA_D02_neurology))]
#------------------------------------------------------------------------------------
#PLA_PPEA_D02_cardiometabolic
PLA_PPEA_D02_cardiometabolic=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D02/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D02_olink-explore_protein-expression_PLA-PPEA-D02_matrix_cardiometabolic.csv")
#PLA_PPEA_D02_cardiometabolic=matrix(PLA_PPEA_D02_cardiometabolic)
PLA_PPEA_D02_cardiometabolic=as.data.frame(PLA_PPEA_D02_cardiometabolic)
# rename rows
rownames(PLA_PPEA_D02_cardiometabolic)=PLA_PPEA_D02_cardiometabolic[,1]
# remove uniprot id column
PLA_PPEA_D02_cardiometabolic=PLA_PPEA_D02_cardiometabolic[,-1]
PLA_PPEA_D02_cardiometabolic=PLA_PPEA_D02_cardiometabolic[,grep("BLM0T1",colnames(PLA_PPEA_D02_cardiometabolic))]
#------------------------------------------------------------------------------------
#PLA_PPEA_D02_inflammation
PLA_PPEA_D02_inflammation=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D02/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D02_olink-explore_protein-expression_PLA-PPEA-D02_matrix_inflammation.csv")
#PLA_PPEA_D02_inflammation=matrix(PLA_PPEA_D02_inflammation)
PLA_PPEA_D02_inflammation=as.data.frame(PLA_PPEA_D02_inflammation)
# rename rows
rownames(PLA_PPEA_D02_inflammation)=PLA_PPEA_D02_inflammation[,1]
# remove uniprot id column
PLA_PPEA_D02_inflammation=PLA_PPEA_D02_inflammation[,-1]
PLA_PPEA_D02_inflammation=PLA_PPEA_D02_inflammation[,grep("BLM0T1",colnames(PLA_PPEA_D02_inflammation))]
#------------------------------------------------------------------------------------
#PLA_PPEA_D02_oncology
PLA_PPEA_D02_oncology=vroom("./AMP PD Data v3 2022/proteomics/proteomics-PLA-PPEA-D02/olink-explore_protein-expression/datasets_proteomics-PLA-PPEA-D02_olink-explore_protein-expression_PLA-PPEA-D02_matrix_oncology.csv")
#PLA_PPEA_D02_oncology=matrix(PLA_PPEA_D02_oncology)
PLA_PPEA_D02_oncology=as.data.frame(PLA_PPEA_D02_oncology)
# rename rows
rownames(PLA_PPEA_D02_oncology)=PLA_PPEA_D02_oncology[,1]
# remove uniprot id column
PLA_PPEA_D02_oncology=PLA_PPEA_D02_oncology[,-1]
PLA_PPEA_D02_oncology=PLA_PPEA_D02_oncology[,grep("BLM0T1",colnames(PLA_PPEA_D02_oncology))]

# remove repeated proteins: tCSF_PPEA_D01_inflammation=tCSF_PPEA_D01_inflammation[,!(colnames(tCSF_PPEA_D01_inflammation) %in% c("P01375","P05231","P10145"))]
# pla D01
tPLA_PPEA_D01_cardiometabolic=as.data.frame(t(PLA_PPEA_D01_cardiometabolic))
tPLA_PPEA_D01_cardiometabolic=tPLA_PPEA_D01_cardiometabolic[,!(colnames(tPLA_PPEA_D01_cardiometabolic) %in% c("P01375","P05231","P10145"))]
#colnames(tPLA_PPEA_D01_cardiometabolic)=paste(colnames(tPLA_PPEA_D01_cardiometabolic),"_card",sep = "")
tPLA_PPEA_D01_neurology=as.data.frame(t(PLA_PPEA_D01_neurology))
tPLA_PPEA_D01_neurology=tPLA_PPEA_D01_neurology[,!(colnames(tPLA_PPEA_D01_neurology) %in% c("P01375","P05231","P10145"))]
#colnames(tPLA_PPEA_D01_neurology)=paste(colnames(tPLA_PPEA_D01_neurology),"_neu",sep = "")
tPLA_PPEA_D01_oncology=as.data.frame(t(PLA_PPEA_D01_oncology))
tPLA_PPEA_D01_oncology=tPLA_PPEA_D01_oncology[,!(colnames(tPLA_PPEA_D01_oncology) %in% c("P01375","P05231","P10145","P42331","Q15797"))]
#colnames(tPLA_PPEA_D01_oncology)=paste(colnames(tPLA_PPEA_D01_oncology),"_onc",sep = "")
tPLA_PPEA_D01_inflammation=as.data.frame(t(PLA_PPEA_D01_inflammation))
tPLA_PPEA_D01_inflammation=tPLA_PPEA_D01_inflammation[,!(colnames(tPLA_PPEA_D01_inflammation) %in% c("P01375","P05231","P10145"))]
#colnames(tPLA_PPEA_D01_inflammation)=paste(colnames(tPLA_PPEA_D01_inflammation),"_inf",sep = "")


# merge in tpla_d01
tPLA_D01 = merge(tPLA_PPEA_D01_cardiometabolic, tPLA_PPEA_D01_neurology, by = "row.names", all = TRUE) 
rownames(tPLA_D01)=tPLA_D01[,"Row.names"]
tPLA_D01=tPLA_D01[,-1]
tPLA_D01 = merge(tPLA_D01, tPLA_PPEA_D01_oncology, by = "row.names", all = TRUE) 
rownames(tPLA_D01)=tPLA_D01[,"Row.names"]
tPLA_D01=tPLA_D01[,-1]
tPLA_D01 = merge(tPLA_D01, tPLA_PPEA_D01_inflammation, by = "row.names", all = TRUE) 
rownames(tPLA_D01)=tPLA_D01[,"Row.names"]
tPLA_D01=tPLA_D01[,-1]
# pla D02
tPLA_PPEA_D02_cardiometabolic=as.data.frame(t(PLA_PPEA_D02_cardiometabolic))
tPLA_PPEA_D02_cardiometabolic=tPLA_PPEA_D02_cardiometabolic[,!(colnames(tPLA_PPEA_D02_cardiometabolic) %in% c("P01375","P05231","P10145"))]
#colnames(tPLA_PPEA_D02_cardiometabolic)=paste(colnames(tPLA_PPEA_D02_cardiometabolic),"_card",sep = "")
tPLA_PPEA_D02_neurology=as.data.frame(t(PLA_PPEA_D02_neurology))
tPLA_PPEA_D02_neurology=tPLA_PPEA_D02_neurology[,!(colnames(tPLA_PPEA_D02_neurology) %in% c("P01375","P05231","P10145"))]
#colnames(tPLA_PPEA_D02_neurology)=paste(colnames(tPLA_PPEA_D02_neurology),"_neu",sep = "")
tPLA_PPEA_D02_oncology=as.data.frame(t(PLA_PPEA_D02_oncology))
tPLA_PPEA_D02_oncology=tPLA_PPEA_D02_oncology[,!(colnames(tPLA_PPEA_D02_oncology) %in% c("P01375","P05231","P10145"))]
#colnames(tPLA_PPEA_D02_oncology)=paste(colnames(tPLA_PPEA_D02_oncology),"_onc",sep = "")
tPLA_PPEA_D02_inflammation=as.data.frame(t(PLA_PPEA_D02_inflammation))
tPLA_PPEA_D02_inflammation=tPLA_PPEA_D02_inflammation[,!(colnames(tPLA_PPEA_D02_inflammation) %in% c("P01375","P05231","P10145"))]
#colnames(tPLA_PPEA_D02_inflammation)=paste(colnames(tPLA_PPEA_D02_inflammation),"_inf",sep = "")
# merge in tpla_d02
tPLA_D02 = merge(tPLA_PPEA_D02_cardiometabolic, tPLA_PPEA_D02_neurology, by = "row.names", all = TRUE) 
rownames(tPLA_D02)=tPLA_D02[,"Row.names"]
tPLA_D02=tPLA_D02[,-1]
tPLA_D02 = merge(tPLA_D02, tPLA_PPEA_D02_oncology, by = "row.names", all = TRUE) 
rownames(tPLA_D02)=tPLA_D02[,"Row.names"]
tPLA_D02=tPLA_D02[,-1]
tPLA_D02 = merge(tPLA_D02, tPLA_PPEA_D02_inflammation, by = "row.names", all = TRUE) 
rownames(tPLA_D02)=tPLA_D02[,"Row.names"]
tPLA_D02=tPLA_D02[,-1]

# change patient names
filas=c()
for(i in rownames(tPLA_D01)){
  fila=substr(i,start = 0,stop = gregexpr("-",i)[[1]][3]-1)
  filas=c(filas,fila)
}
rownames(tPLA_D01)=filas

filas=c()
for(i in rownames(tPLA_D02)){
  fila=substr(i,start = 0,stop = gregexpr("-",i)[[1]][3]-1)
  filas=c(filas,fila)
}
rownames(tPLA_D02)=filas


tPLA_D01$Dataset="D01"
tPLA_D02$Dataset="D02"

PLA=rbind(tPLA_D01,tPLA_D02)

