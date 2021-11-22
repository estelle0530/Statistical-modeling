
# **Canonical Correlation Analysis Functions**  
## Background
Multivariate analysis such as principle component analysis (PCA) and partial least square regression (PLS) have been utilized to understand genetic data. While PCA selects variables by maximizing variance within a data set and PLS selects variables by maximizing covariance between two data sets, Canonical correlation analysis (CCA) selects variables by maximizing correlation between two data sets. This is a useful multivariate model to understand high dimensional data when the underlying assumption is multivariable correlations between two data sets of interest. 

## CCA interpretation
### 1. CCA variate score  
CCA variate score represents the rotation of the sample space to achieve to maximized correlations between two data sets. One can think of CCA variate score synonymous to PCA scores.   

### 2. CCA projection loading 
The advantage of using CCA to understand multivariate variable selection is its interpretability as one can project the rotation score back to the two data sets' original space which then yield correlation coefficients that are bounded between -1 to 1. One can interpret variables by these projected scores as how much a variable contribute to maximize correlation between two data sets of interest. 

## Overview of CCA functions
### 1. CCA_run 
Input two matrices that are matched by samples into the function. This is a regularized CCA, implemented for cases where numbers of variables exceed the numbers of samples. For computation efficiency, the default parameter is shrinkage for relaxed hyperparameter tuning. User can also change it to cross-validation which will result in a more robust parameter tuning. 

### 2. CCA_loading_graph
Using the loading outputs of CCA_run, one can graph scatter plots of the projected loadings for each data set. This is a great way to visualize the variables that are contributing to maximizing correlations.   

### 3. CCA_variate_graph 
Using the common variate  output from CCA_run, one can test for type enrichment for each component. One can also visualize by scatter plot for variate score distirbution across components. 


