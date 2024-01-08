# clinical clustering
# read metadata
mri <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MRI.csv")
# preprocessing
mri=mri[mri$visit_name=="M0",]
mri=mri[grepl("PD-",mri$participant_id)|grepl("PP-",mri$participant_id),]
rownames(mri)=mri$participant_id
mri2=merge(x=genesprot2,y = mri,by = "row.names")
rownames(mri2)=mri2$Row.names
mri2=mri2[,-1]
mritest=c()
for(i in c("mri_results")){
  mritest=c(mritest,chisq.test(table(mri2$icluster,mri2[,i]))$p.value)
}


sleep <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_REM_Sleep_Behavior_Disorder_Questionnaire_Mayo.csv")
# preprocessing
sleep=sleep[sleep$visit_name=="M0",]
sleep=sleep[grepl("PD-",sleep$participant_id)|grepl("PP-",sleep$participant_id),]
rownames(sleep)=sleep$participant_id
sleep2=merge(x=genesprot2,y = sleep,by = "row.names")
rownames(sleep2)=sleep2$Row.names
sleep2=sleep2[,-1]
sleeptest=c()
for(i in c("msq01_act_out_dreams","msq02_legs_jerk","msq03_restless_legs","msq04b_walked_asleep","msq05_snorted_awake","msq06_stop_breathing","msq07_leg_cramps","msq08.2")){
  sleeptest=c(sleeptest,chisq.test(table(sleep2$icluster,sleep2[,i]))$p.value)
}
sleep_aov=aov(sleep2$msq08_rate_of_alertness ~ sleep2$icluster, data = sleep2)
summary(sleep_aov)
leveneTest(sleep2$msq08_rate_of_alertness ~ as.factor(sleep2$icluster), data = sleep2)
plot(sleep_aov, 2)
shapiro.test(x = sleep_aov$residuals )
kruskal.test(sleep2$msq08_rate_of_alertness ~ sleep2$icluster, data = sleep2)
##################################################################################

mds_updrs1 <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MDS_UPDRS_Part_I.csv")
mds_updrs1=mds_updrs1[grepl("PD-",mds_updrs1$participant_id)|grepl("PP-",mds_updrs1$participant_id),]
mds_updrs1=mds_updrs1[mds_updrs1$visit_name=="M0",]
rownames(mds_updrs1)=mds_updrs1$participant_id
mds_updrs12=merge(x=genesprot2,y = mds_updrs1,by = "row.names")
rownames(mds_updrs12)=mds_updrs12$Row.names
mds_updrs12=mds_updrs12[,-1]
dim(mds_updrs12)
mds_updrs1test=c()
for(i in colnames(mds_updrs12[,3939:3966])){
  mds_updrs1test=c(mds_updrs1test,chisq.test(table(mds_updrs12$icluster,mds_updrs12[,i]))$p.value)
}
summary(aov(mds_updrs12[,3968] ~ mds_updrs12$icluster, data = mds_updrs12))



mds_updrs2 <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MDS_UPDRS_Part_II.csv")
mds_updrs2=mds_updrs2[grepl("PD-",mds_updrs2$participant_id)|grepl("PP-",mds_updrs2$participant_id),]
mds_updrs2=mds_updrs2[mds_updrs2$visit_name=="M0",]
rownames(mds_updrs2)=mds_updrs2$participant_id
mds_updrs22=merge(x=genesprot2,y = mds_updrs2,by = "row.names")
rownames(mds_updrs22)=mds_updrs22$Row.names
mds_updrs22=mds_updrs22[,-1]
mds_updrs2test=c()
for(i in colnames(mds_updrs22[,3939:ncol(mds_updrs22)])){
  mds_updrs2test=c(mds_updrs2test,chisq.test(table(mds_updrs22$icluster,mds_updrs22[,i]))$p.value)
}
summary(aov(mds_updrs22[,3965] ~ mds_updrs22$icluster, data = mds_updrs22))
kruskal.test(mds_updrs22[,3965] ~ mds_updrs22$icluster, data = mds_updrs22)

mds_updrs3 <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MDS_UPDRS_Part_III.csv")
mds_updrs3=mds_updrs3[grepl("PD-",mds_updrs3$participant_id)|grepl("PP-",mds_updrs3$participant_id),]
mds_updrs3=mds_updrs3[mds_updrs3$visit_name=="M0",]
rownames(mds_updrs3)=mds_updrs3$participant_id
mds_updrs32=merge(x=genesprot2,y = mds_updrs3,by = "row.names")
rownames(mds_updrs32)=mds_updrs32$Row.names
mds_updrs32=mds_updrs32[,-1]
mds_updrs3test=c()
for(i in colnames(mds_updrs32[,3938:4009])){
  mds_updrs3test=c(mds_updrs3test,chisq.test(table(mds_updrs32$icluster,mds_updrs32[,i]))$p.value)
}
# significant pvalue<0.05
#"code_upd2314_body_bradykinesia","code_upd2317c_rest_tremor_amplitude_right_lower_extremity"
#"upd2314_body_bradykinesia","upd2317c_rest_tremor_amplitude_right_lower_extremity" 
library(corrplot)
corrplot(chisq.test(table(mds_updrs32$icluster,mds_updrs32[,"code_upd2314_body_bradykinesia"]))$residuals, is.cor = FALSE)
corrplot(chisq.test(table(mds_updrs32$icluster,mds_updrs32[,"code_upd2317c_rest_tremor_amplitude_right_lower_extremity"]))$residuals, is.cor = FALSE)
corrplot(chisq.test(table(mds_updrs32$icluster,mds_updrs32[,"upd2314_body_bradykinesia"]))$residuals, is.cor = FALSE)
corrplot(chisq.test(table(mds_updrs32$icluster,mds_updrs32[,"upd2317c_rest_tremor_amplitude_right_lower_extremity"]))$residuals, is.cor = FALSE)




mds_updrs4 <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MDS_UPDRS_Part_IV.csv")
mds_updrs4=mds_updrs4[grepl("PD-",mds_updrs4$participant_id)|grepl("PP-",mds_updrs4$participant_id),]
mds_updrs4=mds_updrs4[mds_updrs4$visit_name=="M0",]
rownames(mds_updrs4)=mds_updrs4$participant_id
mds_updrs42=merge(x=genesprot2,y = mds_updrs4,by = "row.names")
rownames(mds_updrs42)=mds_updrs42$Row.names
mds_updrs42=mds_updrs42[,-1]
mds_updrs4test=c()
for(i in colnames(mds_updrs42[,c(3938:3950)])){
  mds_updrs4test=c(mds_updrs4test,chisq.test(table(mds_updrs42$icluster,mds_updrs42[,i]))$p.value)
}
# test assumptions
upd4_aov=aov(mds_updrs42[,3950] ~ mds_updrs42$icluster, data = mds_updrs42)
summary(upd4_aov)
leveneTest(mds_updrs42[,3950] ~ as.factor(mds_updrs42$icluster), data = mds_updrs42)
plot(upd4_aov, 2)
shapiro.test(x = upd4_aov$residuals )



soma_plasma <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Biospecimen_analyses_SomaLogic_plasma.csv")
soma_plasma=soma_plasma[grepl("PD-",soma_plasma$participant_id)|grepl("PP-",soma_plasma$participant_id),]
soma_plasma=soma_plasma[soma_plasma$visit_name=="M0",]
rownames(soma_plasma)=soma_plasma$participant_id


adl <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Modified_Schwab___England_ADL.csv")
adl=adl[grepl("PD-",adl$participant_id)|grepl("PP-",adl$participant_id),]
adl=adl[adl$visit_name=="M0",]
rownames(adl)=adl$participant_id
adl2=merge(x=genesprot2,y = adl,by = "row.names")
rownames(adl2)=adl$Row.names
adl2=adl2[,-1]
# test assumptions
adl_aov=aov(adl2$mod_schwab_england_pct_adl_score ~ adl2$icluster, data = adl2)
summary(adl_aov)
leveneTest(adl2$mod_schwab_england_pct_adl_score ~ as.factor(adl2$icluster), data = adl2)
plot(adl_aov, 2)
shapiro.test(x = adl_aov$residuals )




pdq <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_PDQ_39.csv")
pdq=pdq[grepl("PD-",pdq$participant_id)|grepl("PP-",pdq$participant_id),]
pdq=pdq[pdq$visit_name=="M0",]
rownames(pdq)=pdq$participant_id
pdq2=merge(x=genesprot2,y = pdq,by = "row.names")
rownames(pdq2)=pdq$Row.names
pdq2=pdq2[,-1]
pdqtest=c()
for(i in colnames(pdq2[,3977:3984])){
  pdq_aov=aov(pdq2[,i] ~ pdq2$icluster, data = pdq2)
  pdqtest=c(pdqtest,summary(pdq_aov)[[1]]$`Pr(>F)`[1])
}

socialscore_aov=aov(pdq2$pdq39_social_score ~ pdq2$icluster, data = pdq2)
summary(socialscore_aov)
leveneTest(pdq2$pdq39_social_score ~ as.factor(pdq2$icluster), data = pdq2)
plot(socialscore_aov, 2)
shapiro.test(x = socialscore_aov$residuals )
kruskal.test(pdq2$pdq39_social_score ~ as.factor(pdq2$icluster), data = pdq2)
library("ggpubr")
ggboxplot(pdq2, x = "icluster", y = "pdq39_social_score", 
          color = "icluster",
          order = c("1", "2", "3","4"),
          ylab = "pdq39_social_score", xlab = "cluster")
pairwise.t.test(pdq2$pdq39_social_score, pdq2$icluster,
                p.adjust.method = "BH")
#"pdq39_20_angry","pdq39_36_felt_ignored"
# corrplot(chisq.test(table(pdq2$icluster,pdq2[,"pdq39_20_angry"]))$residuals, is.cor = FALSE)
# corrplot(chisq.test(table(pdq2$icluster,pdq2[,"pdq39_36_felt_ignored"]))$residuals, is.cor = FALSE)




datscan <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_DaTSCAN_visual_interpretation.csv")
datscan=datscan[grepl("PD-",datscan$participant_id)|grepl("PP-",datscan$participant_id),]
datscan=datscan[datscan$visit_name=="LOG",]
rownames(datscan)=datscan$participant_id
datscan2=merge(x=genesprot2,y = datscan,by = "row.names")
rownames(datscan2)=datscan$Row.names
datscan2=datscan2[,-1]
datscantest=c()
for(i in c("datscan_visual_interpretation")){
  datscantest=c(datscantest,chisq.test(table(datscan2$icluster,datscan2[,i]))$p.value)
}
corrplot(chisq.test(table(datscan2$icluster,datscan2[,"datscan_visual_interpretation"]))$residuals, is.cor = FALSE)


biospecimen <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Biospecimen_analyses_other.csv")
biospecimen=biospecimen[grepl("PD-",biospecimen$participant_id)|grepl("PP-",biospecimen$participant_id),]
biospecimen=biospecimen[biospecimen$visit_name=="M0",]
# H-Ferritin and serum
biospecimenH=biospecimen[biospecimen$test_name=="H-Ferritin"&biospecimen$sample_type=="Serum",]
rownames(biospecimenH)=biospecimenH$participant_id
biospecimenH2=merge(x=genesprot2,y = biospecimenH,by = "row.names")
rownames(biospecimenH2)=biospecimenH2$Row.names
biospecimenH2=biospecimenH2[,-1]
# check asumptions
library(car)
Hserum_aov=aov(biospecimenH2$test_value ~ biospecimenH2$icluster, data = biospecimenH2)
summary(Hserum_aov)
leveneTest(biospecimenH2$test_value ~ as.factor(biospecimenH2$icluster), data = biospecimenH2)
plot(Hserum_aov, 2)
shapiro.test(x = Hserum_aov$residuals )
# remove outliers
#biospecimenH2=biospecimenH2[!(rownames(biospecimenH2) %in% c("PD-PDRF073JTG","PD-PDJH549RB6","PD-PDBE349TEN","PD-PDKT482ZER","PD-PDVV728UB6")),]
################################################################################
# Actin beta
biospecimenAB=biospecimen[biospecimen$test_name=="Actin Beta"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimenAB)=biospecimenAB$participant_id
biospecimenAB2=merge(x=genesprot2,y = biospecimenAB,by = "row.names")
rownames(biospecimenAB2)=biospecimenAB2$Row.names
biospecimenAB2=biospecimenAB2[,-1]
actinbeta_aov=aov(biospecimenAB2$test_value ~ biospecimenAB2$icluster, data = biospecimenAB2)
#################################################################################
biospecimen_adhesion=biospecimen[biospecimen$test_name=="Basal Cell Adhesion Molecule"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_adhesion)=biospecimen_adhesion$participant_id
biospecimen_adhesion2=merge(x=genesprot2,y = biospecimen_adhesion,by = "row.names")
rownames(biospecimen_adhesion2)=biospecimen_adhesion2$Row.names
biospecimen_adhesion2=biospecimen_adhesion2[,-1]
#################################################################################
biospecimen_coatomer=biospecimen[biospecimen$test_name=="Coatomer Protein Complex Subunit Zeta 1"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_coatomer)=biospecimen_coatomer$participant_id
biospecimen_coatomer2=merge(x=genesprot2,y = biospecimen_coatomer,by = "row.names")
rownames(biospecimen_coatomer2)=biospecimen_coatomer2$Row.names
biospecimen_coatomer2=biospecimen_coatomer2[,-1]
################################################################################
biospecimen_cop9=biospecimen[biospecimen$test_name=="COP9 Signalosome Subunit 7A"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_cop9)=biospecimen_cop9$participant_id
biospecimen_cop92=merge(x=genesprot2,y = biospecimen_cop9,by = "row.names")
rownames(biospecimen_cop92)=biospecimen_cop92$Row.names
biospecimen_cop92=biospecimen_cop92[,-1]
################################################################################
biospecimen_fos8=biospecimen[biospecimen$test_name=="Dual Specificity Phosphatase 8"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_fos8)=biospecimen_fos8$participant_id
biospecimen_fos82=merge(x=genesprot2,y = biospecimen_fos8,by = "row.names")
rownames(biospecimen_fos82)=biospecimen_fos82$Row.names
biospecimen_fos82=biospecimen_fos82[,-1]
################################################################################
biospecimen_faty2=biospecimen[biospecimen$test_name=="Fatty Acid Hydroxylase Domain Containing 2"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_faty2)=biospecimen_faty2$participant_id
biospecimen_faty22=merge(x=genesprot2,y = biospecimen_faty2,by = "row.names")
rownames(biospecimen_faty22)=biospecimen_faty22$Row.names
biospecimen_faty22=biospecimen_faty22[,-1]
################################################################################
biospecimen_gly3=biospecimen[biospecimen$test_name=="Glyceraldehyde-3-Phosphate Dehydrogenase"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_gly3)=biospecimen_gly3$participant_id
biospecimen_gly32=merge(x=genesprot2,y = biospecimen_gly3,by = "row.names")
rownames(biospecimen_gly32)=biospecimen_gly32$Row.names
biospecimen_gly32=biospecimen_gly32[,-1]
#################################################################################
biospecimen_lst8=biospecimen[biospecimen$test_name=="MTOR Associated Protein, LST8 Homolog"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_lst8)=biospecimen_lst8$participant_id
biospecimen_lst82=merge(x=genesprot2,y = biospecimen_lst8,by = "row.names")
rownames(biospecimen_lst82)=biospecimen_lst82$Row.names
biospecimen_lst82=biospecimen_lst82[,-1]
##################################################################################
biospecimen_nico=biospecimen[biospecimen$test_name=="Nicotinamide Phosphoribosyltransferase"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_nico)=biospecimen_nico$participant_id
biospecimen_nico2=merge(x=genesprot2,y = biospecimen_nico,by = "row.names")
rownames(biospecimen_nico2)=biospecimen_nico2$Row.names
biospecimen_nico2=biospecimen_nico2[,-1]
dim(biospecimen_nico2)
##################################################################################
biospecimen_poly=biospecimen[biospecimen$test_name=="Polypyrimidine Tract Binding Protein 1"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_poly)=biospecimen_poly$participant_id
biospecimen_poly2=merge(x=genesprot2,y = biospecimen_poly,by = "row.names")
rownames(biospecimen_poly2)=biospecimen_poly2$Row.names
biospecimen_poly2=biospecimen_poly2[,-1]
dim(biospecimen_poly2)
##################################################################################
biospecimen_tyr=biospecimen[biospecimen$test_name=="Protein Tyrosine Phosphatase, Non-Receptor Type 1"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_tyr)=biospecimen_tyr$participant_id
biospecimen_tyr2=merge(x=genesprot2,y = biospecimen_tyr,by = "row.names")
rownames(biospecimen_tyr2)=biospecimen_tyr2$Row.names
biospecimen_tyr2=biospecimen_tyr2[,-1]
dim(biospecimen_tyr2)
##################################################################################
biospecimen_tyc=biospecimen[biospecimen$test_name=="Protein Tyrosine Phosphatase, Receptor Type C"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_tyc)=biospecimen_tyc$participant_id
biospecimen_tyc2=merge(x=genesprot2,y = biospecimen_tyc,by = "row.names")
rownames(biospecimen_tyc2)=biospecimen_tyc2$Row.names
biospecimen_tyc2=biospecimen_tyc2[,-1]
dim(biospecimen_tyc2)
##################################################################################
biospecimen_urea=biospecimen[biospecimen$test_name=="Urea transporter 1"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_urea)=biospecimen_urea$participant_id
biospecimen_urea2=merge(x=genesprot2,y = biospecimen_urea,by = "row.names")
rownames(biospecimen_urea2)=biospecimen_urea2$Row.names
biospecimen_urea2=biospecimen_urea2[,-1]
dim(biospecimen_urea2)
##################################################################################
biospecimen_bcl2=biospecimen[biospecimen$test_name=="BCL2, Apoptosis Regulator"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_bcl2)=biospecimen_bcl2$participant_id
biospecimen_bcl22=merge(x=genesprot2,y = biospecimen_bcl2,by = "row.names")
rownames(biospecimen_bcl22)=biospecimen_bcl22$Row.names
biospecimen_bcl22=biospecimen_bcl22[,-1]
dim(biospecimen_bcl22)
##################################################################################
biospecimen_zinc=biospecimen[biospecimen$test_name=="Zinc Finger Protein 134"&biospecimen$test_units=="Ct Avg",]
rownames(biospecimen_zinc)=biospecimen_zinc$participant_id
biospecimen_zinc2=merge(x=genesprot2,y = biospecimen_zinc,by = "row.names")
rownames(biospecimen_zinc2)=biospecimen_zinc2$Row.names
biospecimen_zinc2=biospecimen_zinc2[,-1]
dim(biospecimen_zinc2)
##################################################################################
biospecimen_glu=biospecimen[biospecimen$test_name=="Total Glucosylceramide",]
rownames(biospecimen_glu)=biospecimen_glu$participant_id
biospecimen_glu2=merge(x=genesprot2,y = biospecimen_glu,by = "row.names")
rownames(biospecimen_glu2)=biospecimen_glu2$Row.names
biospecimen_glu2=biospecimen_glu2[,-1]
dim(biospecimen_glu2)
# test assumptions
biospecimen_glu2=biospecimen_glu2[!(rownames(biospecimen_glu2) %in% c("PD-PDDF170BDZ","PD-PDVV728UB6","PD-PDXK386DCB","PD-PDDM051VVT")),]
glu_aov=aov(biospecimen_glu2$test_value ~ biospecimen_glu2$icluster, data = biospecimen_glu2)
summary(glu_aov)
leveneTest(biospecimen_glu2$test_value ~ as.factor(biospecimen_glu2$icluster), data = biospecimen_glu2)
plot(glu_aov, 2)
shapiro.test(x = glu_aov$residuals )

#################################################################################
#demographics
demo <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Demographics.csv")
demo=demo[grepl("PD-",demo$participant_id)|grepl("PP-",demo$participant_id),]
demo=demo[demo$visit_name=="M0",]
rownames(demo)=demo$participant_id
demo2=merge(x=genesprot2,y = demo,by = "row.names")
rownames(demo2)=demo2$Row.names
demo2=demo2[,-1]
dim(demo2)
# test assumptions
demo2=demo2[!(rownames(demo2) %in% c("PD-PDNA262BHX","PD-PDWK436LWD","PD-PDUA877LP2")),]
demo_aov=aov(demo2$age_at_baseline ~ demo2$icluster, data = demo2)
summary(demo_aov)
leveneTest(demo2$age_at_baseline ~ as.factor(demo2$icluster), data = demo2)
plot(demo_aov, 2)
shapiro.test(x = demo_aov$residuals )
kruskal.test(demo2$age_at_baseline ~ demo2$icluster, data = demo2)
demo_test=c()
for(i in c("sex","race","education_level_years","ethnicity")){
  demo_test=c(demo_test,chisq.test(table(demo2$icluster,demo2[,i]))$p.value)
}
################################################################################


mut <- read.csv("AMP PD Data v3 2022/clinical_metadata/tier2/releases_2022_v3release_1115_amp_pd_participant_mutations.csv")
mut=mut[grepl("PD-",mut$participant_id)|grepl("PP-",mut$participant_id),]
rownames(mut)=mut$participant_id
mut2=merge(x=genesprot2,y = mut,by = "row.names")
rownames(mut2)=mut2$Row.names
mut2=mut2[,-1]
dim(mut2)
mut_test=c()
for(i in colnames(mut2[,3935:3939])){
  mut_test=c(mut_test,chisq.test(table(mut2$icluster,mut2[,i]))$p.value)
}
corrplot(chisq.test(table(mut2$icluster,mut2[,"has_known_SNCA_mutation_in_WGS"]))$residuals, is.cor = FALSE)

#################################################################################

mutape <- read.csv("AMP PD Data v3 2022/clinical_metadata/tier2/releases_2022_v3release_1115_amp_pd_participant_apoe_mutations.csv")
mutape=mutape[grepl("PD-",mutape$participant_id)|grepl("PP-",mutape$participant_id),]
rownames(mutape)=mutape$participant_id
mutape2=merge(x=genesprot2,y = mutape,by = "row.names")
rownames(mutape2)=mutape2$Row.names
mutape2=mutape2[,-1]
dim(mutape2)
mutape_test=c()
for(i in colnames(mutape2[,3935:3936])){
  mutape_test=c(mutape_test,chisq.test(table(mutape2$icluster,mutape2[,i]))$p.value)
}

################################################################################
genstat <- read.csv("AMP PD Data v3 2022/clinical_metadata/tier2/clinical/releases_2022_v3release_1115_clinical_Clinically_Reported_Genetic_Status.csv")
genstat=genstat[grepl("PD-",genstat$participant_id)|grepl("PP-",genstat$participant_id),]
rownames(genstat)=genstat$participant_id
genstat2=merge(x=genesprot2,y = genstat,by = "row.names")
rownames(genstat2)=genstat2$Row.names
genstat2=genstat2[,-1]
dim(genstat2)
genstat_test=c()
for(i in c("genetic_status_enrollment","genetic_status_wgs")){
  genstat_test=c(genstat_test,chisq.test(table(genstat2$icluster,genstat2[,i]))$p.value)
}
#################################################################################

abeta <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Biospecimen_analyses_CSF_abeta_tau_ptau.csv")
abeta=abeta[grepl("PD-",abeta$participant_id)|grepl("PP-",abeta$participant_id),]
abeta=abeta[abeta$visit_name=="M0",]
# measure: abeta
abeta=abeta[abeta$test_name=="p-Tau",]
rownames(abeta)=abeta$participant_id
abeta2=merge(x=genesprot2,y = abeta,by = "row.names")
rownames(abeta2)=abeta2$Row.names
abeta2=abeta2[,-1]
dim(abeta2)
# test assumptions
#abeta2=abeta2[!(rownames(abeta2) %in% c("PD-PDDF170BDZ","PD-PDVV728UB6","PD-PDXK386DCB","PD-PDDM051VVT")),]
abeta_aov=aov(abeta2$test_value ~ abeta2$icluster, data = abeta2)
summary(abeta_aov)
leveneTest(abeta2$test_value ~ as.factor(abeta2$icluster), data = abeta2)
plot(abeta_aov, 2)
shapiro.test(x = abeta_aov$residuals )
kruskal.test(abeta2$test_value ~ as.factor(abeta2$icluster), data = abeta2)

################################################################################

updr <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_UPDRS.csv")
updr=updr[grepl("PD-",updr$participant_id)|grepl("PP-",updr$participant_id),]
updr=updr[updr$visit_name=="M0",]
rownames(updr)=updr$participant_id
updr2=merge(x=genesprot2,y = updr,by = "row.names")
rownames(updr2)=updr2$Row.names
updr2=updr2[,-1]
updrtest=c()
for(i in colnames(updr2[,3939:ncol(updr2)])){
  updrtest=c(updrtest,chisq.test(table(updr2$icluster,updr2[,i]))$p.value)
}

################################################################################

glucer <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Biospecimen_analyses_CSF_beta_glucocerebrosidase.csv")

################################################################################

kolster <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv")
kolster=kolster[grepl("PD-",kolster$participant_id)|grepl("PP-",kolster$participant_id),]
kolster=kolster[kolster$visit_name=="M0",]
rownames(kolster)=kolster$participant_id
kolster2=merge(x=genesprot2,y = kolster,by = "row.names")
rownames(kolster2)=kolster2$Row.names
kolster2=kolster2[,-1]
kolstest=c()
for(i in colnames(kolster2[,3939:ncol(kolster2)])){
  kolstest=c(kolstest,chisq.test(table(kolster2$icluster,kolster2[,i]))$p.value)
}

##################################################################################

epw <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Epworth_Sleepiness_Scale.csv")
epw=epw[grepl("PD-",epw$participant_id)|grepl("PP-",epw$participant_id),]
epw=epw[epw$visit_name=="M0",]
rownames(epw)=epw$participant_id
epw2=merge(x=genesprot2,y = epw,by = "row.names")
rownames(epw2)=epw2$Row.names
epw2=epw2[,-1]
dim(epw2)
epwtest=c()
for(i in colnames(epw2[,3939:3954])){
  epwtest=c(epwtest,chisq.test(table(epw2$icluster,epw2[,i]))$p.value)
}
summary(aov(epw2[,3955] ~ epw2$icluster, data = epw2))

#################################################################################

global <- read.csv("AMP PD Data v3 2022/clinical_metadata/releases_2022_v3release_1115_amp_pd_global_sample_inventory.csv")

#################################################################################

cafe <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Caffeine_history.csv")
cafe=cafe[grepl("PD-",cafe$participant_id)|grepl("PP-",cafe$participant_id),]
cafe=cafe[cafe$visit_name=="M0",]
rownames(cafe)=cafe$participant_id
cafe2=merge(x=genesprot2,y = cafe,by = "row.names")
rownames(cafe2)=cafe2$Row.names
cafe2=cafe2[,-1]
dim(cafe2)
#################################################################################

dti <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_DTI.csv")
dti=dti[grepl("PD-",dti$participant_id)|grepl("PP-",dti$participant_id),]
dti=dti[dti$visit_name=="LOG",]
# eigen1
dti=dti[dti$dti_measure=="Fractional Anisotropy",]
rownames(dti)=dti$participant_id
dti2=merge(x=genesprot2,y = dti,by = "row.names")
rownames(dti2)=dti2$Row.names
dti2=dti2[,-1]
dim(dti2)

dtitest=c()
for(i in colnames(dti2[,3940:3947])){
  dti_aov=aov(dti2[,i] ~ dti2$icluster, data = dti2)
  dtitest=c(dtitest,summary(dti_aov)[[1]]$`Pr(>F)`[1])
}
for(i in colnames(dti2[,3938:3939])){
  print(chisq.test(table(dti2$icluster,dti2[,i]))$p.value)
}
corrplot(chisq.test(table(dti2$icluster,dti2[,3938]))$residuals, is.cor = FALSE)
corrplot(chisq.test(table(dti2$icluster,dti2[,3939]))$residuals, is.cor = FALSE)

##################################################################################3

upsit <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_UPSIT.csv")
upsit=upsit[grepl("PD-",upsit$participant_id)|grepl("PP-",upsit$participant_id),]
upsit=upsit[upsit$visit_name=="M0",]
rownames(upsit)=upsit$participant_id
upsit2=merge(x=genesprot2,y = upsit,by = "row.names")
rownames(upsit2)=upsit2$Row.names
upsit2=upsit2[,-1]
dim(upsit2)
upsittest=c()
for(i in colnames(upsit2[,3938:3943])){
  upsittest=c(upsittest,print(chisq.test(table(upsit2$icluster,upsit2[,i]))$p.value))
}
##################################################################################

mmse <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MMSE.csv")
mmse=mmse[grepl("PD-",mmse$participant_id)|grepl("PP-",mmse$participant_id),]
mmse=mmse[mmse$visit_name=="M0",]
mmse=mmse[mmse$participant_id==unique(mmse$participant_id), ]
rownames(mmse)=mmse$participant_id
mmse2=merge(x=genesprot2,y = mmse,by = "row.names")
rownames(mmse2)=mmse2$Row.names
mmse2=mmse2[,-1]
dim(mmse2)
################################################################################

moca <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_MOCA.csv")
moca=moca[grepl("PD-",moca$participant_id)|grepl("PP-",moca$participant_id),]
moca=moca[moca$visit_name=="M0",]
rownames(moca)=moca$participant_id
moca2=merge(x=genesprot2,y = moca,by = "row.names")
rownames(moca2)=moca2$Row.names
moca2=moca2[,-1]
dim(moca2)
mocatest=c()
for(i in c("moca_visuospatial_executive_subscore","moca_naming_subscore","moca_attention_digits_subscore","moca11_attention_vigilance","moca12_attention_serial_7s","moca13_sentence_repetition","moca15_verbal_fluency","moca_language_subscore","moca16_abstraction","moca_abstraction_subscore","moca_delayed_recall_subscore_optnl_cat_cue","moca_delayed_recall_subscore_optnl_mult_choice","moca_delayed_recall_subscore","moca_orientation_subscore","code_education_12years_complete")){
  mocatest=c(mocatest,(chisq.test(table(moca2$icluster,moca2[,i]))$p.value))
}
corrplot(chisq.test(table(moca2$icluster,moca2[,"moca_delayed_recall_subscore"]))$residuals, is.cor = FALSE,cl.align.text="l")

##################################################################################

smoke <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Smoking_and_alcohol_history.csv")
smoke=smoke[grepl("PD-",smoke$participant_id)|grepl("PP-",smoke$participant_id),]
smoke=smoke[smoke$visit_name=="M0",]
rownames(smoke)=smoke$participant_id
smoke2=merge(x=genesprot2,y = smoke,by = "row.names")
rownames(smoke2)=smoke2$Row.names
smoke2=smoke2[,-1]
dim(smoke2)

smoketest=c()
for(i in c("tobacco_ever_used","tobacco_current_use","alcohol_ever_used","alcohol_current_use","tobacco_recent_use","tobacco_prior_use","tobacco_product_type","cigarettes_per_day","cigarettes_packs_per_day","alcohol_recent_use","alcohol_prior_use","alcohol_use_frequency","alcohol_drinks_daily_range","alcohol_six_more_drinks_frequency")){
  smoketest=c(smoketest,(chisq.test(table(smoke2$icluster,smoke2[,i]))$p.value))
}
corrplot(chisq.test(table(smoke2$icluster,smoke2[,"tobacco_recent_use"])[,-1])$residuals, is.cor = FALSE,cl.align.text="l")
corrplot(chisq.test(table(smoke2$icluster,smoke2[,"tobacco_prior_use"])[,-1])$residuals, is.cor = FALSE,cl.align.text="l")
summary(aov(smoke2[,"tobacco_start_age"] ~ smoke2$icluster, data = smoke2))

#################################################################################

family <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_Family_History_PD.csv")
family=family[grepl("PD-",family$participant_id)|grepl("PP-",family$participant_id),]
family=family[family$visit_name=="M0",]
rownames(family)=family$participant_id
family2=merge(x=genesprot2,y = family,by = "row.names")
rownames(family2)=family2$Row.names
family2=family2[,-1]
dim(family2)
famtest=c()
for(i in c("biological_mother_with_pd","biological_father_with_pd","other_relative_with_pd")){
  famtest=c(famtest,(chisq.test(table(family2$icluster,family2[,i]))$p.value))
}

#################################################################################

med <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_PD_Medical_History.csv")
med=med[grepl("PD-",med$participant_id)|grepl("PP-",med$participant_id),]
med=med[med$visit_name=="M0",]
rownames(med)=med$participant_id
med2=merge(x=genesprot2,y = med,by = "row.names")
rownames(med2)=med2$Row.names
med2=med2[,-1]
dim(med2)
medtest=c()
for(i in c("on_levodopa","on_dopamine_agonist","on_other_pd_medications")){
  medtest=c(medtest,(chisq.test(table(med2$icluster,med2[,i]))$p.value))
}

################################################################################

lbd <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_LBD_Cohort_Clinical_Data.csv")
lbd=lbd[grepl("PD-",lbd$participant_id)|grepl("PP-",lbd$participant_id),]
lbd=lbd[lbd$visit_name=="M0",]
rownames(lbd)=lbd$participant_id
lbd2=merge(x=genesprot2,y = lbd,by = "row.names")
rownames(lbd2)=lbd2$Row.names
lbd2=lbd2[,-1]
dim(lbd2)
################################################################################

lbd2 <- read.csv("AMP PD Data v3 2022/clinical_metadata/clinical/releases_2022_v3release_1115_clinical_LBD_Cohort_Path_Data.csv")
