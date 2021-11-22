# Bayesian NMF-EM framework

## Background
[Bayesian non-negative matrix factorization model](https://arxiv.org/pdf/1111.6085.pdf) models single cell counts with exponential prior and infer genetic signatures in single cell RNA-seq transcriptomic data. However, background rate can confound true signal when conducting downstream analysis

## Method 
Derived and implemented exponential-gaussian mixture model to infer true expression activities
![alt text](https://github.com/estelleyao0530/Statistical-modeling/blob/main/Figure/EM_example.png)

## Results
Increased power in detecting differentially expresed genetic signatures across conditions
