# plot clinical

# read metadata
demo <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Demographics.csv")
family <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Family_History_PD.csv")
smoke_alcohol <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Smoking_and_alcohol_history.csv")
updrs <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_UPDRS.csv")
medical_history <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_PD_Medical_History.csv")
mds_updrs1 <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MDS_UPDRS_Part_I.csv")

#tcpms2=data.frame(tcpms)
#tClustMx2=data.frame(tClustMx)
# preprocessing
smoke_alcohol=smoke_alcohol[smoke_alcohol$visit_name=="M0",]
smoke_alcohol=smoke_alcohol[grepl("PD-",smoke_alcohol$participant_id)|grepl("PP-",smoke_alcohol$participant_id),]

family=family[grepl("PD-",family$participant_id)|grepl("PP-",family$participant_id),]
family=family[family$visit_name=="M0",]

demo=demo[grepl("PD-",demo$participant_id)|grepl("PP-",demo$participant_id),]
demo=demo[demo$visit_name=="M0",]

updrs=updrs[grepl("PD-",updrs$participant_id)|grepl("PP-",updrs$participant_id),]
updrs=updrs[updrs$visit_name=="M0",]

mds_updrs1=mds_updrs1[grepl("PD-",mds_updrs1$participant_id)|grepl("PP-",mds_updrs1$participant_id),]
mds_updrs1=mds_updrs1[mds_updrs1$visit_name=="M0",]

medical_history=medical_history[grepl("PD-",medical_history$participant_id)|grepl("PP-",medical_history$participant_id),]
medical_history=medical_history[medical_history$visit_name=="M0",]

# rename columns
filas=c()
for(i in rownames(genesprot)){
  fila=substr(i,start = 0,stop = gregexpr("-",i)[[1]][2]-1)
  filas=c(filas,fila)
}
genesprot2=genesprot
rownames(genesprot2)=filas

rownames(smoke_alcohol)=smoke_alcohol$participant_id
rownames(family)=family$participant_id
rownames(demo)=demo$participant_id
rownames(updrs)=updrs$participant_id
rownames(mds_updrs1)=mds_updrs1$participant_id
rownames(medical_history)=medical_history$participant_id

# smoke alcohol
smoke_alcohol2=merge(x=genesprot2,y = smoke_alcohol,by = "row.names")
rownames(smoke_alcohol2)=smoke_alcohol2$Row.names
smoke_alcohol2=smoke_alcohol2[,-1]
smoketest=c()
for(i in c("tobacco_ever_used","tobacco_current_use","alcohol_ever_used","alcohol_current_use","tobacco_recent_use","tobacco_prior_use","tobacco_product_type","cigarettes_per_day","cigarettes_packs_per_day","alcohol_recent_use","alcohol_prior_use","alcohol_use_frequency","alcohol_drinks_daily_range","alcohol_six_more_drinks_frequency","alcohol_related_hospitalization")){
  smoketest=c(smoketest,chisq.test(table(smoke_alcohol2$icluster,smoke_alcohol2[,i]))$p.value)
}
library(corrplot)
corrplot(chisq.test(table(smoke_alcohol2$icluster,smoke_alcohol2[,"tobacco_prior_use"]))$residuals, is.cor = FALSE)
balloonplot(table(smoke_alcohol2$icluster,smoke_alcohol2[,"tobacco_prior_use"]), main ="housetasks", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)
# fviz_pca_ind(wpca, label="none",habillage = as.factor(smoke[-c(24,11,47,46,53,41),]$icluster))
# fviz_pca_ind(wpca, label="none",habillage = as.factor(smoke[-c(24,11,47,46,53,41),]$tobacco_ever_used))
# fviz_pca_ind(wpca, label="none",habillage = as.factor(smoke[-c(24,11,47,46,53,41),]$tobacco_current_use))
# fviz_pca_ind(wpca, label="none",habillage = as.factor(smoke[-c(24,11,47,46,53,41),]$tobacco_product_type))
# fviz_pca_ind(wpca, label="none",habillage = as.factor(smoke[-c(24,11,47,46,53,41),]$alcohol_use_frequency))

# family
family2=merge(x=tcpms2,y = tClustMx2,by = "participant_id")
family2=merge(x=family2,y = family,by = "participant_id")

wpca=PCA(family2[-c(24,11,47,46,53,41),-c(1,2828,4298:4304)])
fviz_pca_ind(wpca,label="none",habillage = as.factor(family2[-c(24,11,47,46,53,41),]$icluster),addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.95)
fviz_pca_ind(wpca,label="none",habillage = as.factor(family2[-c(24,11,47,46,53,41),]$biological_mother_with_pd),addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.95)
fviz_pca_ind(wpca,label="none",habillage = as.factor(family2[-c(24,11,47,46,53,41),]$biological_father_with_pd),addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.95)
fviz_pca_ind(wpca,label="none",habillage = as.factor(family2[-c(24,11,47,46,53,41),]$other_relative_with_pd),addEllipses = TRUE, ellipse.type = "confidence",ellipse.level = 0.95)

# demo
demo2=merge(x=tcpms2,y = tClustMx2,by = "participant_id")
demo2=merge(x=demo2,y = demo,by = "participant_id")
wpca=PCA(demo2[-c(24,11,47,46,53,41),-c(1,2828,4298:4306)])
fviz_pca_ind(wpca,label="none",habillage = as.factor(demo2[-c(24,11,47,46,53,41),]$icluster))
fviz_pca_ind(wpca,label="none",habillage = as.factor(demo2[-c(24,11,47,46,53,41),]$age_at_baseline))
fviz_pca_ind(wpca,label="none",habillage = as.factor(demo2[-c(24,11,47,46,53,41),]$sex))
fviz_pca_ind(wpca,label="none",habillage = as.factor(demo2[-c(24,11,47,46,53,41),]$ethnicity))
fviz_pca_ind(wpca,label="none",habillage = as.factor(demo2[-c(24,11,47,46,53,41),]$race))

#updrs
mds_updrs1.2=merge(x=tcpms2,y = tClustMx2,by = "participant_id")
mds_updrs1.2=merge(x=mds_updrs1.2,y = mds_updrs1,by = "participant_id")
wpca=PCA(mds_updrs1.2[-c(24,11,47,46,53,41),-c(1,2828,4297:4332)])
fviz_pca_ind(wpca,label="none",habillage = as.factor(mds_updrs1.2[-c(24,11,47,46,53,41),]$icluster))
fviz_pca_ind(wpca,label="none",habillage = as.factor(mds_updrs1.2[-c(24,11,47,46,53,41),]$upd2102_hallucinations_and_psychosis))
fviz_pca_ind(wpca,label="none",habillage = as.factor(mds_updrs1.2[-c(24,11,47,46,53,41),]$upd2103_depressed_mood))
fviz_pca_ind(wpca,label="none",habillage = as.factor(mds_updrs1.2[-c(24,11,47,46,53,41),]$upd2104_anxious_mood))

# only genes and proteins
genesprot=merge(x=tcpms,y=tClustMx,by="row.names")
genesprot$icluster=output[[5]]$fit[[5]]$clusters
# wpca cluster 1
clust1=PCA(genesprot[genesprot$icluster==1,-1])
clust2=PCA(genesprot[genesprot$icluster==2,-1])
clust3=PCA(genesprot[genesprot$icluster==3,-1])
clust4=PCA(genesprot[genesprot$icluster==4,-1])
clust5=PCA(genesprot[genesprot$icluster==5,-1])
clust6=PCA(genesprot[genesprot$icluster==6,-1])
clust7=PCA(genesprot[genesprot$icluster==7,-1])
fviz_pca_ind(clust4,label="none",habillage = as.factor(genesprot[genesprot$icluster==4,-1]$icluster))

# medical history
rownames(genesprot)=genesprot$Row.names
genesprot=genesprot[,-1]
# rename columns
filas=c()
for(i in rownames(genesprot)){
  fila=substr(i,start = 0,stop = gregexpr("-",i)[[1]][2]-1)
  filas=c(filas,fila)
}
rownames(genesprot)=filas
medical2=merge(x=genesprot,y = medical_history,by="row.names")
rownames(medical2)=medical2$Row.names
medical2=medical2[,-1]
wpca=PCA(medical2[-c(11,24,53,47,46,41,67),-c(1,4297:4318)])
medical3=medical2[-c(11,24,53,47,46,41,67),]
fviz_pca_ind(wpca,label="none",habillage = as.factor(medical3$icluster))
fviz_pca_ind(wpca,label="none",habillage = as.factor(medical3$on_levodopa))
fviz_pca_ind(wpca,label="none",habillage = as.factor(medical3$on_dopamine_agonist))


################################################################################
genesprot=merge(x=tcpms,y=tClustMx,by="row.names")
rownames(genesprot)=genesprot$Row.names
genesprot=genesprot[,-1]
filas=c()
for(i in rownames(genesprot)){
  fila=substr(i,start = 0,stop = gregexpr("-",i)[[1]][2]-1)
  filas=c(filas,fila)
}
rownames(genesprot)=filas
genesprot$icluster2=icluster2
medical2=merge(x=genesprot,y = medical_history,by="row.names")

# -c(11,24,53,47,46,41,67)
wpca=PCA(medical2[,-c(1,4297:4318)])
fviz_pca_ind(wpca,label="none",habillage = as.factor(medical2$icluster2))
fviz_pca_ind(wpca,label="none",habillage = as.factor(medical2$on_levodopa))
fviz_pca_ind(wpca,label="none",habillage = as.factor(medical2$on_dopamine_agonist))

# cafe
caffeine <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Caffeine_history.csv")
caffeine=caffeine[caffeine$visit_name=="M0",]
caffeine=caffeine[grepl("PD-",caffeine$participant_id)|grepl("PP-",caffeine$participant_id),]
rownames(caffeine)=caffeine$participant_id
cafe=merge(x = genesprot,y=caffeine,by="row.names")

# sleep
sleep <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Epworth_Sleepiness_Scale.csv")
sleep=sleep[sleep$visit_name=="M0",]
sleep=sleep[grepl("PD-",sleep$participant_id)|grepl("PP-",sleep$participant_id),]
rownames(sleep)=sleep$participant_id
sleep1=merge(x = genesprot,y=sleep,by="row.names")
