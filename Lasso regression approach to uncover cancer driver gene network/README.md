# Lasso regression approach to uncover cancer driver gene network

## Background

Correlation-based analysis has shown significant correlations between each gene’s expression level to the listed genomic features (Cai et al). Past studies have shown that linear regression models reveal strong correlation between copy number variation of a gene to their gene expression pattern (Shao et al.). For example, copy number loss of a gene is strongly correlated with its expression for tumor suppressor genes (Zhao et al.).  Multiple linear regression analysis has also been applied to study co-expression pathways controlling for their copy number variation. However, to date, correlation studies are feature specific without including the interplay of copy number variation, mutation, and expression patterns. 

## Method

To advance understanding of how different genomic features are associated with oncogenes and tumor suppressor gene expression patterns across cancer types, this study proposes six lasso regression models that uncover associations between cancer drivers’ expression with various omi-types. We further construct a multi-omit network for three genes of interest: TP53, KRAS, and MYC to reveal uncovered genomic associations in both pan-cancer and cancer type specific settings.

A regularized model is selected in light of the collinearities among genetic features. With the primary goal to point out signals for top features associated with each driver oncogene of interest, we further choose the lasso regularization approach, which uses L1 penalties  to shrink coefficients to zeros, resulting in a parsimonious model. This project proposes six variations of a lasso model for each central gene of interest, and uses the mean squared error and R squared of the model to choose the best-performed model as a base to construct and analyze the multi-omic network with the central gene. 

## Result

We showcased differential performances for lasso regression based model across driver oncogenes and cancer types
 
![alt_text](https://github.com/estelleyao0530/Statistical-modeling/blob/main/Figure/lasso_result.png)
