#' Run regularized canonical correlation analysis 
#' 
#' Input: projected loading score from matrix A cca_read output
#' Input: projected loading score from matrix B cca_read output
#' 
#' Note: 
#' 1) projected loading avoids collinearity of loadings 
#' 
#' Output: 
#' Scatter plots from specified number of components for matrix A
#' Scatter plots from specified number of components for matrix B  
#' 
#' 
#' cca_loading_graph
#' @param projx dataframe, from cca_read output 
#' @param projy dataframe, from cca_read output
#' @param modulecolor dataframe, variable cluster specification (if null, no color labels) 
#' @param number_plot int, numbers of component to plot
#' @param name string, name of output files
#' 
#' @export
#' 

require(ggplot2)
require(ggpubr)
require(ggrepel)

cca_loading_graph<-function(projx, projy, modulecolor , number_plot, name)
{
  df1=fread(projx)
  color=fread(modulecolor)
  colnames(color)[1]="var"
  
  file<-merge(df1,color, by="var")
  
  circles <- data.frame(x0 = c(0,0,0),y0 = c(0,0,0),
                        r = c(0.2, 0.3, 0.4))
  
  tb<-as.data.frame(table(file$moduleColors))
  file$alpha<-ifelse(file$moduleColors=="grey",0.8, 1)
  
  for ( i in 1:number_plot)
  {
    g<-ggplot(file, aes_string(x=paste0("X",i), y= paste0("X",i+1)))+
      geom_point( size=I(2),aes(color=moduleColors, alpha=alpha), show.legend = F) +
      scale_color_manual(values=as.character(as.factor(tb_g$Var1)))+
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
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/3)) )+
      geom_circle(aes(x0=x0, y0=y0, r=r),inherit.aes = F, data=circles, linetype=2)
    
    ggsave(paste0("CCA Component projx_",name,i,"-",i+1,".png"), g, height = 10, width=10)
    
  }
  
  
  df2=fread(projy)
  file<-data.frame()
  file<-merge(df2, color , by="var")
  
  tb<-as.data.frame(table(file$moduleColors))
  file$alpha<-ifelse(file$moduleColors=="grey",0.8, 1)
  
  for ( i in 1:number_plot)
  {
    g<-ggplot(file, aes_string(x=paste0("X",i), y= paste0("X",i+1)))+
      geom_point( size=I(2),aes(color=moduleColors, alpha=alpha), show.legend = F) +
      scale_color_manual(values=as.character(as.factor(tb_g$Var1)))+
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
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/3)) )+
      geom_circle(aes(x0=x0, y0=y0, r=r),inherit.aes = F, data=circles, linetype=2)
    
    ggsave(paste0("CCA Component projy_",name,i,"-",i+1,".png"), g, height = 10, width=10)
    
  }
  
}

