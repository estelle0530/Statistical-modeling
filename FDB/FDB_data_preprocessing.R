setwd("/Users/estelleyao/Desktop/FDB_updated_data/")

library(data.table)
library(impute)
library(stringr)

##CRISPR DepMap Public 22Q1 Primary Files
crispr_effect<-fread("CRISPR_gene_effect.csv") 
crispr_info<-read.csv('sample_info.csv')
colnames(crispr_effect)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(crispr_effect)))
crispr_effect<-merge(crispr_info, crispr_effect, by="DepMap_ID")

rnai_effect_achillies<-fread("RNAi_(Achilles,_DEMETER2).csv")
colnames(rnai_effect_achillies)=sapply(colnames(rnai_effect_achillies), function(x) strsplit(x, "_")[[1]][1])
rnai_effect_achillies$V1<-gsub("\\s*\\([^\\)]+\\)","",as.character(rnai_effect_achillies$V1))
colnames(rnai_effect_achillies)[1]<-"DepMap_ID" 
rnai_effect_achillies = data.frame(rnai_effect_achillies)

###filter out genes with more than 50% NAs across cell lines
delete_gene_rnai<-colnames(rnai_effect_achillies)[(which(colSums(is.na(rnai_effect_achillies))/nrow(rnai_effect_achillies)>0.5))] 
print(length(delete_gene_rnai))
delete_gene_crispr<-colnames(crispr_effect)[(which(colSums(is.na(crispr_effect))/nrow(crispr_effect)>0.5))] 
print(length(delete_gene_crispr))

###Filtered RNAi file
rnai_effect_achillies = rnai_effect_achillies[,-(which(colSums(is.na(rnai_effect_achillies))/nrow(rnai_effect_achillies)>0.5))]

###CRISPR, RNAi merge 
intergene<-intersect(colnames(crispr_effect)[which(colnames(crispr_effect)=="culture_type")+1:ncol(crispr_effect)],
                     colnames(rnai_effect_achillies))
intersample<-intersect(crispr_effect$DepMap_ID,rnai_effect_achillies$DepMap_ID)

setdiff(rnai_effect_achillies$DepMap_ID, crispr_effect$DepMap_ID) 

m_gene<-na.omit(match(intergene, colnames(rnai_effect_achillies)))
m_sample<-na.omit(match(intersample,rnai_effect_achillies$DepMap_ID))
achilles_crispr<-rnai_effect_achillies[m_sample, m_gene]

m_gene<-na.omit(match(intergene,colnames(crispr_effect)))
m_sample<-na.omit(match(intersample,crispr_effect$DepMap_ID))
crispr_info = crispr_effect[m_sample, 1:which(colnames(crispr_effect)=="culture_type")]
crispr_achilles = crispr_effect[m_sample, m_gene]

ncol(achilles_crispr) == sum(colnames(crispr_achilles) == colnames(achilles_crispr))
ncol(crispr_achilles) == sum(colnames(crispr_achilles) == colnames(achilles_crispr))
nrow(achilles_crispr) == nrow(crispr_achilles)  

###Imputation
rnai_impute<-impute.knn(as.matrix(t(achilles_crispr)))
rnai_impute<-rnai_impute$data
rnai_impute<-as.data.frame(t(rnai_impute))
rownames(rnai_impute)<-crispr_info$DepMap_ID

crispr_impute<-impute.knn(as.matrix(t(crispr_achilles)))
crispr_impute<-crispr_impute$data
crispr_impute<-as.data.frame(t(crispr_impute))
rownames(crispr_impute)<-crispr_info$DepMap_ID
colnames(crispr_impute)<-colnames(crispr_achilles)

fwrite(rnai_impute, "achilles_crispr_impute.csv")
fwrite(crispr_impute, "crispr_achilles_impute.csv")

###Drug data ###AUC data

####################CTRP##########################
ctrp = fread("Drug_sensitivity_AUC_(CTD^2).csv")
ctrp_info = fread("ctrp_info.csv")
colnames(ctrp)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(ctrp)))
colnames(ctrp)[1] = "DepMap_ID"

####Filter for drug that have at least 80% data over all cell lines and drugs 
delete_drug<-(which( colSums(is.na(ctrp))/nrow(ctrp)>0.8))
delete_sample<-(which( rowSums(is.na(ctrp))/ncol(ctrp)>0.8))
ctrp = data.frame(ctrp)
ctrp = ctrp[-delete_sample, -delete_drug]
ctrp_info = ctrp_info[-(delete_drug-1),]

fwrite(ctrp_info,"ctrp_info_v2.csv")
setdiff(ctrp$DepMap_ID, crispr_effect$DepMap_ID)

####Match full CRISPR with CTRP
intersample = intersect(crispr_effect$DepMap_ID, 
                        ctrp$DepMap_ID)
m_sample<-na.omit(match(intersample,crispr_effect$DepMap_ID))
crispr_info_ctrp = crispr_effect[m_sample, 1:which(colnames(crispr_effect)=="culture_type")]
crispr_ctrp = crispr_effect[m_sample, (which(colnames(crispr_effect)=="culture_type")+1):ncol(crispr_effect)]

m_sample<-na.omit(match(intersample,ctrp$DepMap_ID))
ctrp_crispr = ctrp[m_sample, 2:ncol(ctrp)]

####Match intersected CRISPR, RNAi with CTRP
intersample = intersect(crispr_info$DepMap_ID, 
                        ctrp$DepMap_ID)
m_sample<-na.omit(match(intersample,crispr_info$DepMap_ID))
crispr_achilles_info_ctrp = crispr_info[m_sample, 1:which(colnames(crispr_info)=="culture_type")]
crispr_achilles_ctrp = crispr_achilles[m_sample, ]
achilles_crispr_ctrp = achilles_crispr[m_sample, ]

m_sample<-na.omit(match(intersample,ctrp$DepMap_ID))
ctrp_crispr_achilles = ctrp[m_sample, 2:ncol(ctrp) ]

####Imputation 
ctrp_impute<-impute.knn(as.matrix(t(ctrp_crispr)))
ctrp_impute<-ctrp_impute$data
ctrp_impute<-as.data.frame(t(ctrp_impute))
rownames(ctrp_impute)<-crispr_info_ctrp$DepMap_ID

crispr_impute<-impute.knn(as.matrix(t(crispr_ctrp)))
crispr_impute<-crispr_impute$data
crispr_impute<-as.data.frame(t(crispr_impute))
rownames(crispr_impute)<-crispr_info_ctrp$DepMap_ID

fwrite(crispr_ctrp, "crispr_ctrp.csv")
fwrite(ctrp_crispr, "ctrp_crispr.csv")
fwrite(ctrp_impute, "ctrp_crispr_impute.csv")
fwrite(crispr_impute, "crispr_ctrp_impute.csv")
fwrite(crispr_info_ctrp, "crispr_info_ctrp.csv")

ctrp_impute<-impute.knn(as.matrix(t(ctrp_crispr_achilles)))
ctrp_impute<-ctrp_impute$data
ctrp_impute<-as.data.frame(t(ctrp_impute))
rownames(ctrp_impute)<-crispr_achilles_info_ctrp$DepMap_ID

rnai_impute<-impute.knn(as.matrix(t(achilles_crispr_ctrp)))
rnai_impute<-rnai_impute$data
rnai_impute<-as.data.frame(t(rnai_impute))
rownames(rnai_impute)<-crispr_achilles_info_ctrp$DepMap_ID

crispr_impute<-impute.knn(as.matrix(t(crispr_achilles_ctrp)))
crispr_impute<-crispr_impute$data
crispr_impute<-as.data.frame(t(crispr_impute))
rownames(crispr_impute)<-crispr_achilles_info_ctrp$DepMap_ID

fwrite(crispr_achilles_ctrp, "crispr_achilles_ctrp.csv")
fwrite(ctrp_crispr_achilles, "ctrp_crispr_achilles.csv")
fwrite(achilles_crispr_ctrp, "achilles_crispr_ctrp.csv")
fwrite(crispr_achilles_info_ctrp,"crispr_achilles_info_ctrp.csv")

fwrite(ctrp_impute, "ctrp_crispr_achilles_impute.csv")
fwrite(rnai_impute, "achilles_crispr_ctrp_impute.csv")
fwrite(crispr_impute, "crispr_achilles_ctrp_impute.csv")

####################GDSC##########################
gdsc1 = fread("Drug_sensitivity_AUC_(Sanger_GDSC1).csv")
gdsc2 = fread("Drug_sensitivity_AUC_(Sanger_GDSC2).csv")
colnames(gdsc2)[1] = "DepMap_ID"
colnames(gdsc1)[1] = "DepMap_ID"

colnames(gdsc1)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(gdsc1)))
colnames(gdsc2)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(gdsc2)))
gdsc_info = fread("screened_compunds_rel_8.2.csv")

length(intersect(colnames(gdsc1), colnames(gdsc2)))
length(intersect(gdsc1$DepMap_ID,  gdsc2$DepMap_ID))

gdsc_info$lower_drug = tolower(gdsc_info$DRUG_NAME)

###combine both GDSC1, GDSC2

# s = fread("sanger-viability.csv")
# d = fread("sanger-dose-response.csv")

delete_drug<-(which( colSums(is.na(gdsc2))/nrow(gdsc2)>0.8))
delete_sample<-(which( rowSums(is.na(gdsc2))/ncol(gdsc2)>0.8))
gdsc2 = data.frame(gdsc2)
gdsc2 = gdsc2[-delete_sample, -delete_drug]
setdiff(gdsc2$DepMap_ID, crispr_effect$DepMap_ID)

gdsc_name = data.frame('drug' = tolower(colnames(gdsc2)))
fwrite(gdsc_name, "gdsc_name.csv")
d = merge(data.frame('drug' = tolower(colnames(gdsc2))), 
          gdsc_info, by.x = "drug", by.y = "lower_drug")
d= d[!duplicated(d$drug),]
fwrite(d, "gdsc_info.csv")

####Match full CRISPR with CTRP
intersample = intersect(crispr_effect$DepMap_ID, 
                        gdsc2$DepMap_ID)
m_sample<-na.omit(match(intersample,crispr_effect$DepMap_ID))
crispr_info_gdsc = crispr_effect[m_sample, 1:which(colnames(crispr_effect)=="culture_type")]
crispr_gdsc = crispr_effect[m_sample, (which(colnames(crispr_effect)=="culture_type")+1):ncol(crispr_effect)]

m_sample<-na.omit(match(intersample,gdsc2$DepMap_ID))
gdsc_crispr = gdsc2[m_sample, 2:ncol(gdsc2)]

###Blood cell lines are present in all intersections 
# View(table(crispr_info_gdsc$primary_disease))
# View(table(crispr_info_ctrp$primary_disease))
# View(table(crispr_achilles_info_ctrp$primary_disease))
# View(table(crispr_info$primary_disease))

####Match intersected CRISPR, RNAi with CTRP
intersample = intersect(crispr_info$DepMap_ID, 
                        gdsc2$DepMap_ID)
m_sample<-na.omit(match(intersample,crispr_info$DepMap_ID))
crispr_achilles_info_gdsc = crispr_info[m_sample, 1:which(colnames(crispr_info)=="culture_type")]
crispr_achilles_gdsc = crispr_achilles[m_sample, ]
achilles_crispr_gdsc = achilles_crispr[m_sample, ]

m_sample<-na.omit(match(intersample,gdsc2$DepMap_ID))
gdsc_crispr_achilles = gdsc2[m_sample, 2:ncol(gdsc2) ]

####Imputation 
gdsc_impute<-impute.knn(as.matrix(t(gdsc_crispr)))
gdsc_impute<-gdsc_impute$data
gdsc_impute<-as.data.frame(t(gdsc_impute))
rownames(gdsc_impute)<-crispr_info_gdsc$DepMap_ID

crispr_impute<-impute.knn(as.matrix(t(crispr_gdsc)))
crispr_impute<-crispr_impute$data
crispr_impute<-as.data.frame(t(crispr_impute))
rownames(crispr_impute)<-crispr_info_gdsc$DepMap_ID

fwrite(crispr_gdsc, "crispr_gdsc.csv")
fwrite(gdsc_crispr, "gdsc_crispr.csv")
fwrite(crispr_info_gdsc, "crispr_info_gdsc.csv")
fwrite(gdsc_impute, "gdsc_crispr_impute.csv")
fwrite(crispr_impute, "crispr_gdsc_impute.csv")

gdsc_impute<-impute.knn(as.matrix(t(gdsc_crispr_achilles)))
gdsc_impute<-gdsc_impute$data
gdsc_impute<-as.data.frame(t(gdsc_impute))
rownames(gdsc_impute)<-crispr_achilles_info_gdsc$DepMap_ID

rnai_impute<-impute.knn(as.matrix(t(achilles_crispr_gdsc)))
rnai_impute<-rnai_impute$data
rnai_impute<-as.data.frame(t(rnai_impute))
rownames(rnai_impute)<-crispr_achilles_info_gdsc$DepMap_ID

crispr_impute<-impute.knn(as.matrix(t(crispr_achilles_gdsc)))
crispr_impute<-crispr_impute$data
crispr_impute<-as.data.frame(t(crispr_impute))
rownames(crispr_impute)<-crispr_achilles_info_gdsc$DepMap_ID

fwrite(crispr_achilles_gdsc, "crispr_achilles_gdsc.csv")
fwrite(gdsc_crispr_achilles, "gdsc_crispr_achilles.csv")
fwrite(achilles_crispr_gdsc, "achilles_crispr_gdsc.csv")
fwrite(crispr_achilles_info_gdsc,"crispr_achilles_info_gdsc.csv")

fwrite(gdsc_impute, "gdsc_crispr_achilles_impute.csv")
fwrite(rnai_impute, "achilles_crispr_gdsdc_impute.csv")
fwrite(crispr_impute, "crispr_achilles_gdsc_impute.csv")

########################PRISM############################
prism = fread("Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_19Q4.csv")
colnames(prism)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(prism)))
colnames(prism)[1] = "DepMap_ID"

delete_drug<-(which( colSums(is.na(prism))/nrow(prism)>0.8))
delete_sample<-(which( rowSums(is.na(prism))/ncol(prism)>0.8))
prism = data.frame(prism)
# prism = prism[-delete_sample, -delete_drug]
# setdiff(prism$DepMap_ID, crispr_effect$DepMap_ID)

####Match full CRISPR with PRISM
intersample = intersect(crispr_effect$DepMap_ID, 
                        prism$DepMap_ID)
m_sample<-na.omit(match(intersample,crispr_effect$DepMap_ID))
crispr_info_prism = crispr_effect[m_sample, 1:which(colnames(crispr_effect)=="culture_type")]
crispr_prism = crispr_effect[m_sample, (which(colnames(crispr_effect)=="culture_type")+1):ncol(crispr_effect)]

m_sample<-na.omit(match(intersample,prism$DepMap_ID))
prism_crispr = prism[m_sample, 2:ncol(prism)]

####Match intersected CRISPR, RNAi with PRISM
intersample = intersect(crispr_info$DepMap_ID, 
                        prism$DepMap_ID)
m_sample<-na.omit(match(intersample,crispr_info$DepMap_ID))
crispr_achilles_info_prism = crispr_info[m_sample, 
                                           1:which(colnames(crispr_info)=="culture_type")]
crispr_achilles_prism = crispr_achilles[m_sample, ]
achilles_crispr_prism = achilles_crispr[m_sample, ]

m_sample<-na.omit(match(intersample,prism$DepMap_ID))
prism_crispr_achilles = prism[m_sample, 2:ncol(prism) ]

####Imputation 
prism_impute<-impute.knn(as.matrix(t(prism_crispr)))
prism_impute<-prism_impute$data
prism_impute<-as.data.frame(t(prism_impute))
rownames(prism_impute)<-crispr_info_prism$DepMap_ID

crispr_impute<-impute.knn(as.matrix(t(crispr_prism)))
crispr_impute<-crispr_impute$data
crispr_impute<-as.data.frame(t(crispr_impute))
rownames(crispr_impute)<-crispr_info_prism$DepMap_ID

fwrite(crispr_prism, "crispr_prism.csv")
fwrite(prism_crispr, "prism_crispr.csv")
fwrite(crispr_info_prism, "crispr_info_prism.csv")
fwrite(prism_impute, "prism_crispr_impute.csv")
fwrite(crispr_impute, "crispr_prism_impute.csv")

prism_impute<-impute.knn(as.matrix(t(prism_crispr_achilles)))
prism_impute<-prism_impute$data
prism_impute<-as.data.frame(t(prism_impute))
rownames(prism_impute)<-crispr_achilles_info_prism$DepMap_ID

rnai_impute<-impute.knn(as.matrix(t(achilles_crispr_prism)))
rnai_impute<-rnai_impute$data
rnai_impute<-as.data.frame(t(rnai_impute))
rownames(rnai_impute)<-crispr_achilles_info_prism$DepMap_ID

crispr_impute<-impute.knn(as.matrix(t(crispr_achilles_prism)))
crispr_impute<-crispr_impute$data
crispr_impute<-as.data.frame(t(crispr_impute))
rownames(crispr_impute)<-crispr_achilles_info_prism$DepMap_ID

fwrite(crispr_achilles_prism, "crispr_achilles_prism.csv")
fwrite(prism_crispr_achilles, "prism_crispr_achilles.csv")
fwrite(achilles_crispr_prism, "achilles_crispr_prism.csv")
fwrite(crispr_achilles_info_prism,"crispr_achilles_info_prism.csv")

fwrite(prism_impute, "prism_crispr_achilles_impute.csv")
fwrite(rnai_impute, "achilles_crispr_prism_impute.csv")
fwrite(crispr_impute, "crispr_achilles_prism_impute.csv")


