#' Read in GSEA results and leading edge genes from terminal GSEA output
#' 
#' Returns combined list of standard GSEA outputs with member genes  
#' based on https://www.gsea-msigdb.org/gsea/index.jsp
#' 
#' @param path directory of GSEA outputs are stored 
#' @param name name of the output file
#' @param enriched boolean, True = only returned enriched genes from each gene set; False = returned all member genes from each gene set
#' 
#' @export
#' 

require(data.table)

read_gsea<-function(path, name, enriched)
{
  setwd(path)
  temp<-list.files(pattern=".xls")
  files<-lapply(temp, fread)
  names(files)<-temp
  names(files)<-sub(".xls","", names(files))
  
  if (enriched==T)
  {
    final<-data.frame()
    for(i in 1:length(temp))
    {
      r<-files[[i]]
      if(ncol(r)==9)
      {
        r<-r[,-c(1,9)]
        r$Geneset<-names(files)[i]
        r<-r[which(r$`CORE ENRICHMENT`=="Yes"),]
        final<-rbind(final, r)
      }
      else
      {next}
    }
    
    gseapos<-files[[which(grepl("gsea", names(files)))[1]]]
    gseaneg<-files[[which(grepl("gsea", names(files)))[2]]]
    
    if(nrow(gseapos)==1)
    {bind<-gseaneg}else{
      bind<-rbind(gseapos, gseaneg)
    }
    
    finalfile<-merge(bind, final, by.x="NAME", by.y="Geneset", all.x=T)
    
    fwrite(finalfile, paste0(path, name, ".csv"))
  }
  
  else
  {
    final<-data.frame()
    for(i in 1:length(temp))
    {
      r<-files[[i]]
      if(ncol(r)==9)
      {
        r<-r[,-c(1,9)]
        r$Geneset<-names(files)[i]
        final<-rbind(final, r)
      }
      else
      {next}
    }
    
    gseapos<-files[[which(grepl("gsea", names(files)))[1]]]
    gseaneg<-files[[which(grepl("gsea", names(files)))[2]]]
    
    if(nrow(gseapos)==1)
    {bind<-gseaneg}else{
      bind<-rbind(gseapos, gseaneg)
    }
    
    finalfile<-merge(bind, final, by.x="NAME", by.y="Geneset", all.x=T)
    
    fwrite(finalfile, paste0(path, name, ".csv"))
    
  }
  
}
