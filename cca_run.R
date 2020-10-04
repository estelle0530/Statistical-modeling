#' Run regularized canonical correlation analysis 
#' 
#' Input: matrix A (rows as samples, cols as variables)
#' Input: matrix B (rows as samples, cols as variables)
#' 
#' Note: 
#' 1) matrix A and B need to be matched by samples  
#' 2) numbers of variables (columns) > numbers of samples (rows)
#' 
#' Output: 
#' cca_run saves cca object from mixOmics package 
#' cca_read writes to file average variate scores and projected loadings of decomposed matrices that maximize the correlation between A and B
#' 
#' based on http://mixomics.org/methods/rcca/
#' 
#' cca_run 
#' @param file1 directory of matrix A
#' @param file2 directory of matrix B
#' @param ncomp int, numbers of components to output
#' @param name string, name of the output file
#' 
#' cca_read
#' @param output_directory directory to output projection scores and variate 
#' @param rdata_path directory of RData stored by cca_run 
#' @param demographic_path directory of sample information, must matched to matrixA, matrix B samples
#' @param name string, name of the output file
#'
#' 
#' @export
#' 

require(mixOmics)
require(data.table)

cca_run<-function(file1, file2, ncomp,name)
{
  df1=fread(file1)
  df2=fread(file2)
  
  cca<-rcc(df1, df2, method="shrinkage", ncomp=ncomp)
  save(cca, file=paste0(name, ".RData" ))
}

cca_read<-function(output_directory,rdata_path, demographic_path, name )
{ 
  setwd(output_directory)
  
  load(rdata_path)
  
  demographic<-fread(demographic_path)
  
  projx<-data.frame(var=colnames(cca$X), cor(cca$X, cca$variates$X+cca$variates$Y, use="pairwise"))
  projy<-data.frame(var=colnames(cca$Y), cor(cca$Y, cca$variates$X+cca$variates$Y, use="pairwise"))
  commonvariate<-data.frame(demographic,cca$variates$X+cca$variates$Y)
  
  write.table(projx, paste0(name,"_projx.txt"), sep="\t", row.names = F, col.names = T, quote = F)
  write.table(projy, paste0(name,"_projy.txt"), sep="\t", row.names = F, col.names = T, quote = F)
  write.table(commonvariate, paste0(name,"_commonvariate.txt"), sep="\t", row.names = F, col.names = T, quote = F)
}


