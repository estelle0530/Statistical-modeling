library(data.table)
library(pcaMethods)
#install.packages("glmnet")
library(glmnet)
library(tidyverse)
#install.packages("caret")
library(caret)

#### Load data + Data Cleaning ####
##Expression file 
##Downloaded from https://depmap.org/portal/download/
ccle<-fread("/Users/yaoestelle/Downloads/CCLE_expression (2).csv")
ccle_info<-fread("/Users/yaoestelle/Downloads/sample_info (3).csv")
colnames(ccle)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(ccle)))
ccle<-as.data.frame(ccle)
colnames(ccle)[1]<-"DepMap_ID"
m<-na.omit(match(ccle$DepMap_ID, ccle_info$DepMap_ID))
ccle_info<-ccle_info[m,]
m<-na.omit(match( ccle_info$DepMap_ID, ccle$DepMap_ID))
ccle<-ccle[m,]

###Copy number files
###Downloaded from https://depmap.org/portal/download/ 
copynumber<-fread("/Users/yaoestelle/Downloads/CCLE_gene_cn (1).csv")
colnames(copynumber)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(copynumber)))
colnames(copynumber)[1]<-"DepMap_ID"

##Mutation file 
##Downloaded from https://portals.broadinstitute.org/ccle/data 
##ccle_mutant<-fread("/Users/yaoestelle/Downloads/CCLE_mutations.csv")
mut_ccle_wide<-fread("wide_mutation_ccle_20Q3.csv")

###Match cell lines between all files
intersect<-intersect(ccle$DepMap_ID, copynumber$DepMap_ID)
intersect<-intersect(intersect, mut_ccle_wide$DepMap_ID)

m<-na.omit(match( intersect,mut_ccle_wide$DepMap_ID))
mut_ccle<-mut_ccle_wide[m,]

m<-na.omit(match(intersect ,copynumber$DepMap_ID ))
copy_ccle_mut<-copynumber[m,]
copy_ccle_mut<-as.data.frame(copy_ccle_mut)

m<-na.omit(match( intersect,ccle$DepMap_ID))
ccle_mut<-ccle[m,]

m<-na.omit(match(intersect ,ccle_info$DepMap_ID))
ccle_info<-ccle_info[m,]

##common gene sequence 
common_gene<-intersect(colnames(ccle_mut), colnames(copy_ccle_mut))
common_gene<-intersect(common_gene, colnames(mut_ccle_wide))

m<-na.omit(match(common_gene ,colnames(ccle_mut)))
ccle_mut<-ccle_mut[,m]

m<-na.omit(match(common_gene ,colnames(copy_ccle_mut)))
copy_ccle_mut<-copy_ccle_mut[,m]

mut_ccle<-as.data.frame(mut_ccle)
m<-na.omit(match(common_gene ,colnames(mut_ccle)))
mut_ccle_wide<-mut_ccle[,m]

####TP53 network genes
dnv1<-read.delim("/Users/yaoestelle/Downloads/geneset (6).txt",stringsAsFactors = F, header = F)
upv1<-read.delim("/Users/yaoestelle/Downloads/geneset (7).txt",stringsAsFactors = F,header = F)
dnv2<-read.delim("/Users/yaoestelle/Downloads/geneset (8).txt",stringsAsFactors = F,header = F)
upv2<-read.delim("/Users/yaoestelle/Downloads/geneset (9).txt",stringsAsFactors = F,header = F)
network<-rbind(dnv1,upv1,dnv2, upv2)
network<-unique(network$V1)


###Function for lasso regression 
eval_results <- function(true, predicted,df) {
  
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}
lasso_cv<-function( X_matrix, response, alpha)
{
  lambdas <- 10^seq(3, -2, by = -0.1);
  cvfit <- cv.glmnet(
    x=as.matrix(X_matrix),
    y=response,
    alpha = alpha,
    lambda = lambdas
  );
  #plot(cvfit);
  fit.min <- glmnet(
    x=as.matrix(X_matrix),
    y=response,
    alpha = alpha,
    lambda = cvfit$lambda.min
  );
  lasso_betas <- coef(fit.min) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column('Variables')%>%
    filter(s0 != 0 & Variables !="(Intercept)") %>%
    arrange(desc(s0));
  
  results<-eval_results(true=response,
               predicted=predict(cvfit, s = "lambda.min", newx = as.matrix(X_matrix) ) ,
               df=X_matrix)
  
  return(list(betas=lasso_betas, model_perform=results))
  
  # lasso_intercept<-coef(fit.min) %>%
  #   as.matrix() %>%
  #   as.data.frame() %>%
  #   rownames_to_column('Variables')%>%
  #   filter(Variables =="(Intercept)")
}

####1) Lasso regression Gene expression ~ Copynumber 
mod1<-lasso_cv(copy_ccle_mut[,which(colnames(copy_ccle_mut)%in%network)] , ccle_mut[,which(colnames(ccle_mut)=="TP53")], alpha=1)

####2) Lasso regression Gene expression ~ Mutation 
mod2<-lasso_cv(mut_ccle_wide[,which(colnames(mut_ccle_wide)%in%network)] , ccle_mut[,which(colnames(ccle_mut)=="TP53")], alpha=1)

####3) Lasso regression Gene expression ~ Copynumber + Mutation 
train_mut_network<-mut_ccle_wide[,which(colnames(mut_ccle_wide)%in%network)]
train_copy_network<-copy_ccle_mut[,which(colnames(copy_ccle_mut)%in%network)]
colnames(train_mut_network)<-paste0("Mutation_",colnames(train_mut_network))
colnames(train_copy_network)<-paste0("Copy_",colnames(train_copy_network))
train_mut_copy_network<-cbind(train_mut_network, train_copy_network)

mod3<-lasso_cv(train_mut_copy_network , ccle_mut[,which(colnames(ccle_mut)=="TP53")], alpha=1)

####4) Lasso regression Gene expression ~ Copynumber + Mutation + Expression 

train_exp_network<-ccle_mut[,which(colnames(ccle_mut)%in% network[-which(network=="TP53")]) ]
colnames(train_exp_network)<-paste0("Exp_",colnames(train_exp_network))
train_mut_copy_exp_network<-cbind(train_mut_network, train_copy_network,train_exp_network)

mod4<-lasso_cv(train_mut_copy_exp_network , ccle_mut[,which(colnames(ccle_mut)=="TP53")], alpha=1)

####5) Lasso regression Gene expression ~ Copynumber + Mutation + Copynumber * Mutation

###Generate pairwise interaction between genes selected in model1 
mod1_beta<-mod1$betas
inter_mut<-mut_ccle_wide[,which(colnames(mut_ccle_wide)%in%mod1_beta$Variables)]
inter_copy<-copy_ccle_mut[,which(colnames(copy_ccle_mut)%in%mod1_beta$Variables)]
colnames(inter_mut)<-paste0("Mutation_",colnames(inter_mut))
colnames(inter_copy)<-paste0("Copy_",colnames(inter_copy))
inter_mut_copy<-cbind(inter_mut, inter_copy)

mut_copy_network_inter<-model.matrix( ~.^2, data=inter_mut_copy)[,-1]
mod5<-lasso_cv(mut_copy_network_inter , ccle_mut[,which(colnames(ccle_mut)=="TP53")], alpha=1)

####6) Lasso regression Gene expression ~ Copynumber + Mutation + Expression + Copynumber * Mutation

inter_mut_copy_exp_network<-cbind(mut_copy_network_inter, train_exp_network)
mod6<-lasso_cv(inter_mut_copy_exp_network , ccle_mut[,which(colnames(ccle_mut)=="TP53")], alpha=1)

###Comparisons of model fit by rMSE, r square 
comp<-rbind(data.frame(model= "Copynumber", RMSE=mod1$model_perform[1], Rsquare=mod1$model_perform[2]),
            data.frame(model= "Mutation", RMSE=mod2$model_perform[1], Rsquare=mod2$model_perform[2]),
            data.frame(model= "Copynumber + Mutation", RMSE=mod3$model_perform[1], Rsquare=mod3$model_perform[2]),
            data.frame(model= "Copynumber + Mutation + Expression", RMSE=mod4$model_perform[1], Rsquare=mod4$model_perform[2]),
            data.frame(model= "Copynumber + Mutation + Copy*Mut", RMSE=mod5$model_perform[1], Rsquare=mod5$model_perform[2]),
            data.frame(model= "Copynumber + Mutation + Expression + Copy*Mut", RMSE=mod6$model_perform[1], Rsquare=mod6$model_perform[2]))

g1<-ggplot(comp, aes(x=model, y=RMSE))+
  geom_bar(stat="identity")+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=20),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))

g2<-ggplot(comp, aes(x=model, y=Rsquare))+
  geom_bar(stat="identity")+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=20),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))

ggsave("bst210_modelrmse.png",g1,width=20, height=10)
ggsave("bst210_modelrsquare.png",g2,width=20, height=10)

####plot coefficients for the top 2 models

plot_mod4<-mod4$betas
plot_mod4$type<-ifelse(grepl("Copy", plot_mod4$Variables), "Copy number",
                       ifelse(grepl("Mutation", plot_mod4$Variables), "Mutation", "Expression"))
plot_mod4$Variables<-factor(plot_mod4$Variables, levels=plot_mod4$Variables[order(plot_mod4$s0)])
g1<-ggplot(plot_mod4, aes(x=Variables, y=s0))+
  geom_bar(stat="identity", aes(fill=type))+
  scale_fill_manual("Omic type", values = c("Copy number"="purple","Expression"="pink", "Mutation"="navyblue" ))+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=15),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))+
  labs(y="Coefficient", title="TP53 Expression~ Copy number + Mutation + Expression")+
  coord_flip()
ggsave("cn_mut_exp_coef.png", g1, width=10, height=10)

plot_mod6<-mod6$betas
plot_mod6$type<-""
for (i in 1:nrow(plot_mod6))
{

  if(grepl(":", plot_mod6$Variables[i]))
  {
    if(lengths(regmatches(plot_mod6$Variables[i], gregexpr("Copy", plot_mod6$Variables[i])))==2 ){
      plot_mod6$type[i]="Copy number interaction"
    }else if(lengths(regmatches(plot_mod6$Variables[i], gregexpr("Mutation", plot_mod6$Variables[i])))==2 ){
        plot_mod6$type[i]="Mutation interaction"
    }else{
          plot_mod6$type[i]="Copy number and Mutation interaction"
    }
  }else{
    if(grepl("Copy",plot_mod6$Variables[i])){
      plot_mod6$type[i]="Copy number"
    }else if(grepl("Mutation",plot_mod6$Variables[i])){
      plot_mod6$type[i]="Mutation"
    }else{
      plot_mod6$type[i]="Expression"
    }
  }
}

plot_mod6$Variables<-factor(plot_mod6$Variables, levels=plot_mod6$Variables[order(plot_mod6$s0)])
g2<-ggplot(plot_mod6, aes(x=Variables, y=s0))+
  geom_bar(stat="identity", aes(fill=type))+
  scale_fill_manual("Omic type", values = c("Copy number"="purple","Expression"="pink", "Mutation"="navyblue",
                                            "Copy number interaction"="lightgreen", "Copy number and Mutation interaction"="orange","Mutation interaction"="skyblue"))+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=14),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))+
  labs(y="Coefficient", title="TP53 Expression~ Copy number + Mutation + Copy*Mutation + Expression")+
  coord_flip()
ggsave("cn_mut_exp_inter_coef.png", g2, width=12, height=16)


####KRAS network genes
dnv1<-read.delim("/Users/yaoestelle/Downloads/geneset (10).txt",stringsAsFactors = F, header = F)
upv1<-read.delim("/Users/yaoestelle/Downloads/geneset (11).txt",stringsAsFactors = F,header = F)
dnv2<-read.delim("/Users/yaoestelle/Downloads/geneset (12).txt",stringsAsFactors = F,header = F)
upv2<-read.delim("/Users/yaoestelle/Downloads/geneset (13).txt",stringsAsFactors = F,header = F)
network<-rbind(dnv1,upv1,dnv2, upv2)
network<-unique(network$V1)

####1) Lasso regression Gene expression ~ Copynumber 
mod1<-lasso_cv(copy_ccle_mut[,which(colnames(copy_ccle_mut)%in%network)] , ccle_mut[,which(colnames(ccle_mut)=="KRAS")], alpha=1)

####2) Lasso regression Gene expression ~ Mutation 
mod2<-lasso_cv(mut_ccle_wide[,which(colnames(mut_ccle_wide)%in%network)] , ccle_mut[,which(colnames(ccle_mut)=="KRAS")], alpha=1)

####3) Lasso regression Gene expression ~ Copynumber + Mutation 
train_mut_network<-mut_ccle_wide[,which(colnames(mut_ccle_wide)%in%network)]
train_copy_network<-copy_ccle_mut[,which(colnames(copy_ccle_mut)%in%network)]
colnames(train_mut_network)<-paste0("Mutation_",colnames(train_mut_network))
colnames(train_copy_network)<-paste0("Copy_",colnames(train_copy_network))
train_mut_copy_network<-cbind(train_mut_network, train_copy_network)

mod3<-lasso_cv(train_mut_copy_network , ccle_mut[,which(colnames(ccle_mut)=="KRAS")], alpha=1)

####4) Lasso regression Gene expression ~ Copynumber + Mutation + Expression 

train_exp_network<-ccle_mut[,which(colnames(ccle_mut)%in% network) ]
colnames(train_exp_network)<-paste0("Exp_",colnames(train_exp_network))
train_mut_copy_exp_network<-cbind(train_mut_network, train_copy_network,train_exp_network)

mod4<-lasso_cv(train_mut_copy_exp_network , ccle_mut[,which(colnames(ccle_mut)=="KRAS")], alpha=1)

####5) Lasso regression Gene expression ~ Copynumber + Mutation + Copynumber * Mutation

###Generate pairwise interaction between genes selected in model1 
mod1_beta<-mod1$betas
inter_mut<-mut_ccle_wide[,which(colnames(mut_ccle_wide)%in%mod1_beta$Variables)]
inter_copy<-copy_ccle_mut[,which(colnames(copy_ccle_mut)%in%mod1_beta$Variables)]
colnames(inter_mut)<-paste0("Mutation_",colnames(inter_mut))
colnames(inter_copy)<-paste0("Copy_",colnames(inter_copy))
inter_mut_copy<-cbind(inter_mut, inter_copy)

mut_copy_network_inter<-model.matrix( ~.^2, data=inter_mut_copy)[,-1]
mod5<-lasso_cv(mut_copy_network_inter , ccle_mut[,which(colnames(ccle_mut)=="KRAS")], alpha=1)

####6) Lasso regression Gene expression ~ Copynumber + Mutation + Expression + Copynumber * Mutation

inter_mut_copy_exp_network<-cbind(mut_copy_network_inter, train_exp_network)
mod6<-lasso_cv(inter_mut_copy_exp_network , ccle_mut[,which(colnames(ccle_mut)=="KRAS")], alpha=1)

###Comparisons of model fit by rMSE, r square 
comp<-rbind(data.frame(model= "Copynumber", RMSE=mod1$model_perform[1], Rsquare=mod1$model_perform[2]),
            data.frame(model= "Mutation", RMSE=mod2$model_perform[1], Rsquare=mod2$model_perform[2]),
            data.frame(model= "Copynumber + Mutation", RMSE=mod3$model_perform[1], Rsquare=mod3$model_perform[2]),
            data.frame(model= "Copynumber + Mutation + Expression", RMSE=mod4$model_perform[1], Rsquare=mod4$model_perform[2]),
            data.frame(model= "Copynumber + Mutation + Copy*Mut", RMSE=mod5$model_perform[1], Rsquare=mod5$model_perform[2]),
            data.frame(model= "Copynumber + Mutation + Expression + Copy*Mut", RMSE=mod6$model_perform[1], Rsquare=mod6$model_perform[2]))

g1<-ggplot(comp, aes(x=model, y=RMSE))+
  geom_bar(stat="identity")+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=20),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))

g2<-ggplot(comp, aes(x=model, y=Rsquare))+
  geom_bar(stat="identity")+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=20),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))

ggsave("bst210_KRAS_modelrmse.png",g1,width=20, height=10)
ggsave("bst210_KRAS_modelrsquare.png",g2,width=20, height=10)

####plot coefficients for the top 2 models

plot_mod4<-mod4$betas
plot_mod4$type<-ifelse(grepl("Copy", plot_mod4$Variables), "Copy number",
                       ifelse(grepl("Mutation", plot_mod4$Variables), "Mutation", "Expression"))
plot_mod4$Variables<-factor(plot_mod4$Variables, levels=plot_mod4$Variables[order(plot_mod4$s0)])
g1<-ggplot(plot_mod4, aes(x=Variables, y=s0))+
  geom_bar(stat="identity", aes(fill=type))+
  scale_fill_manual("Omic type", values = c("Copy number"="purple","Expression"="pink", "Mutation"="navyblue" ))+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=15),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))+
  labs(y="Coefficient", title="KRAS Expression~ Copy number + Mutation + Expression")+
  coord_flip()
ggsave("cn_mut_exp_coef_KRAS.png", g1, width=10, height=16)

plot_mod6<-mod6$betas
plot_mod6$type<-""
for (i in 1:nrow(plot_mod6))
{
  
  if(grepl(":", plot_mod6$Variables[i]))
  {
    if(lengths(regmatches(plot_mod6$Variables[i], gregexpr("Copy", plot_mod6$Variables[i])))==2 ){
      plot_mod6$type[i]="Copy number interaction"
    }else if(lengths(regmatches(plot_mod6$Variables[i], gregexpr("Mutation", plot_mod6$Variables[i])))==2 ){
      plot_mod6$type[i]="Mutation interaction"
    }else{
      plot_mod6$type[i]="Copy number and Mutation interaction"
    }
  }else{
    if(grepl("Copy",plot_mod6$Variables[i])){
      plot_mod6$type[i]="Copy number"
    }else if(grepl("Mutation",plot_mod6$Variables[i])){
      plot_mod6$type[i]="Mutation"
    }else{
      plot_mod6$type[i]="Expression"
    }
  }
}

plot_mod6$Variables<-factor(plot_mod6$Variables, levels=plot_mod6$Variables[order(plot_mod6$s0)])
g2<-ggplot(plot_mod6, aes(x=Variables, y=s0))+
  geom_bar(stat="identity", aes(fill=type))+
  scale_fill_manual("Omic type", values = c("Copy number"="purple","Expression"="pink", "Mutation"="navyblue",
                                            "Copy number interaction"="lightgreen", "Copy number and Mutation interaction"="orange","Mutation interaction"="skyblue"))+
  theme( axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         text = element_text(size=14),
         axis.title.x = element_text(size=20,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
         axis.title.y = element_text(size=20,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10))+
  labs(y="Coefficient", title="KRAS Expression~ Copy number + Mutation + Copy*Mutation + Expression")+
  coord_flip()
ggsave("cn_mut_exp_inter_coef_KRAS.png", g2, width=12, height=16)



