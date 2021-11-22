#' Run regularized canonical correlation analysis 
#' 
#' Input: commonvariate score from cca_read output
#' Input: projected loading score from matrix B cca_read output
#' 
#' Note: 
#' 1) one can also use only variate score from matrix A or matrix B  
#' 2) common variates are the average fo variate score of matrix A and matrix B rotaions 
#' 
#' Output: 
#' Bar plots of enriched types in each component based on KS test p value 
#' 
#' cca_variate_enrich
#' @param variate dataframe, from cca_read output 
#' @param category dataframe, variable cluster specification (if null, no color labels) 
#' @param name string, name of output files
#' 
#' Note: 
#' 1) variate and category must match by sample names 
#' 
#' cca_variate_scatter
#' @param variate dataframe, from cca_read output 
#' @param category dataframe, variable cluster specification (if null, no color labels) 
#' @param number_plot int, numbers of components to plot 
#' @param name string, name of output files
#' 
#' @export
#' 

###Load 2 sided KS test and name the function ks.test.2 from https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/ks.test.R

require(ggplot2)
require(ggpubr)
require(ggrepel)

###This approach can also applied to analyze enriched types by PCA scores 
cca_variate_enrich<-function( variate, category, name )
{
  ###disease contains cancer type 
  type_table<-data.frame(table(category$disease))
  ks_type<-type_table[type_table$Freq>1,]
  
  common_variate<-data.frame(variate, category)
  
  ##result 1 stores ks test p value for each cancer type for each component 
  result1<-data.frame()
  for (i in 1:nrow(ks_type))
  {
    result<-data.frame()
    test<-variate[variate$disease==ks_type$Var1[i],]
    background<-variate[variate$disease!=ks_type$Var1[i],]
    for (t in 1:5)
    {
      ks<-ks.test.2(test[,t], background[,t], alternative="two.sided")
      df<-data.frame(Type=ks_type$Var1[i], Component=t, ks.pval=ks$p.value, ESscore=ks$ES )
      result<-rbind(result,df)
    }
    result1<-rbind(result1,result)
  }
  result1$log10<-(-log10(result1$ks.pval))
  
  ##PC1-5
  g<-ggplot(data=result1, aes(x=Type, y=log10, fill=as.factor(Component))) +
    geom_bar(stat="identity", position=position_dodge())+
    scale_fill_brewer(palette="Paired")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=60, hjust=1, size=8))+
    geom_hline(yintercept=(-log10(0.05)))+
    labs(x="Cancer types", y="-log10 transformed Kolmogorov Simorov test p value")
  
  ##PC1-2
  plot<-result1[result$Component==1|result$Component==2,]
  g2<-ggplot(data=plot, aes(x=Type, y=log10, fill=as.factor(Component))) +
    geom_bar(stat="identity", position=position_dodge())+
    scale_fill_brewer(palette="Paired")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=60, hjust=1, size=8))+
    geom_hline(yintercept=(-log10(0.05)))+
    labs(x="Cancer types", y="-log10 transformed Kolmogorov Simorov test p value")
  
  ggsave(paste0(name, "_variate_enrich_first5comp.png"), g1)
  ggsave(paste0(name, "_variate_enrich_first2comp.png"), g2)
  
}


cca_variate_scatter<-function( variate, category, number_plot ,name )
{
  common_variate<-data.frame(variate, category)
  
  for( i in 1:number_plot)
  {
    
    g<-ggplot(file, aes_string(x=paste0("X",i), y= paste0("X",i+1)))+
      geom_point( size=I(2),aes(color=disease), show.legend = F) 
      theme(legend.position="none",legend.title = element_text(size=4),
            axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
      labs(title = paste0("CCA Component ",i,"-",i+1) , 
           x=paste0("Component ",i), y=paste0("Component ",i+1))+
      theme(text = element_text(size=8),
            axis.title.x = element_text(size=8,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
            axis.title.y = element_text(size=8,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      geom_text(data = g_file, mapping = aes(label = plotgene),hjust=0, vjust=0, check_overlap = T, size = 1.5)+
      geom_hline(yintercept= 0, linetype="dashed")+
      geom_vline(xintercept = 0, linetype="dashed")+
      coord_fixed()+
      scale_x_continuous(breaks = round(seq(-0.8,0.8, by = 0.1),1), limits = c(-0.8,0.8))+
      scale_y_continuous(breaks = round(seq(-0.8,0.8, by = 0.1),1), limits = c(-0.8,0.8))+
      geom_abline(intercept = 0, linetype="dashed", slope=(-tan(pi/6)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/6)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(-tan(pi/3)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/3)) )
    
    ggsave(paste0(name,i,"-",i+1 ,"_variate_score.png"), g)
   }
  
}
