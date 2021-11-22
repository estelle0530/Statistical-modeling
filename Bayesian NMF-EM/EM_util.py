import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import random
import string
import subprocess
import anndata as an
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy 
import warnings
from collections import Counter
from google.cloud import storage
import re
from scipy.stats import norm, expon

def EM(w, total_iter):

    ###initialize CRS +, CRS- distribution proportion by cell counts at D7
    pi_n0 = Counter(d7.CRS)[1.0]/d7.shape[0]
    pi_e0 = Counter(d7.CRS)[0.0]/d7.shape[0]
    pi_n0, pi_e0

    ###initialize normal distritbution parameter per signature by MLE fitting normal distribution on the entire dataset
    mu_0, std_0 = norm.fit(d7[f'W{w}'])
    ###initialize exponential distritbution parameter per signature by MLE fitting exp distribution on the entire dataset
    loc, lambda_0 = expon.fit(d7[f'W{w}'])
    scale_0 = 0 ###center at 0 (not tuned)

    mu = mu_0
    std = std_0
    sigma = std**2
    pi_n = pi_n0
    lambda_current = lambda_0

    data = np.asarray(d7[f'W{w}'].values)
    N = d7.shape[0]

    def calc_rnorm(pi_n, mu, std,x, evidence):

        return (pi_n * norm.pdf(x, mu, std))/evidence

    def log_likelihood(data, mu, sigma, lambda_current, pi_n):

        return np.sum( np.log(pi_n * norm.pdf(data, mu, sigma**0.5) + (1-pi_n) * expon.pdf(data, loc = loc, scale = lambda_current) ) )

    pi_tracker, mu_tracker, sigma_tracker, lambda_tracker, likelihood_tracker = [],[],[],[],[]

    # pi_tracker.append(pi_n)
    # mu_tracker.append(mu)
    # sigma_tracker.append(sigma)
    # lambda_tracker.append(lambda_current)
    # likelihood_tracker.append(log_likelihood(data, mu, sigma, lambda_current, pi_n))
    
    i = 0
    threshold = 1e-50
    
    pi_diff = 1
    mu_diff = 1
    sigma_diff = 1
    lambda_diff = 1
    
    while i < total_iter and mu_diff > threshold and pi_diff > threshold and sigma_diff > threshold and lambda_diff > threshold:

        ###E step posterior of mixture gaussian and exponential 
        evidence = pi_n * norm.pdf(data, mu, std) + (1-pi_n) * expon.pdf(data, loc =loc, scale = lambda_current)
        r_n = (pi_n * norm.pdf(data, mu, std)) / evidence
        #r_e = (1-pi_n) * (expon.pdf(x, scale = lamd)) / evidence

        ### M step - maximize ELBO for mu, sigma, lambda, pi 

        pi_new = np.sum(r_n) / N ###pi_exp = 1-pi_new

        mu_new = np.sum(data * r_n) /np.sum(r_n)
        #mu_new = np.sum([x*calc_rnorm(x) for x in data])/ np.sum(r_n)
        sigma_new = np.sum([(data-mu_new)**2 * r_n])/np.sum(r_n)
        #sigma_new = np.sum([(x-mu)**2*calc_rnorm(x) for x in data])/np.sum(r_n)

        ###this lambda = scale parameter in exponential function for scipy
        lambda_new = np.sum(data*(1-r_n))/ np.sum(1-r_n)
        #lambda_new = np.sum([x* (1-calc_rnorm(x)) for x in data])/np.sum(1-r_n)
        
        #convergence check
        mu_diff = np.linalg.norm(mu_new - mu)
        pi_diff = np.linalg.norm(pi_new - pi_n)
        sigma_diff = np.linalg.norm(sigma_new - sigma)
        lambda_diff = np.linalg.norm(lambda_new - lambda_current)

        pi_n = pi_new
        mu = mu_new
        sigma = sigma_new
        lambda_current = lambda_new

        likelihood = log_likelihood(data, mu, sigma, lambda_current, pi_n)

        pi_tracker.append(pi_n)
        mu_tracker.append(mu)
        sigma_tracker.append(sigma)
        lambda_tracker.append(lambda_current)
        likelihood_tracker.append(likelihood)

        i+=1
        
    return pi_tracker, mu_tracker,sigma_tracker, lambda_tracker, likelihood_tracker

def EM_diagnostic(w ,pi_tracker, mu_tracker,sigma_tracker, lambda_tracker, likelihood_tracker):
    
    fig,ax = plt.subplots(1,5, figsize = (35,5))

    ax[0].plot(likelihood_tracker)
    ax[0].set_xlabel("Iterations")
    ax[0].set_ylabel("Log likelihood")

    ax[1].plot(pi_tracker, label = "normal proportion")
    ax[1].plot(1- np.asarray(pi_tracker), label = "exponential proportion")

    ax[1].set_xlabel("Iterations")
    ax[1].set_ylabel("% estimate for \n normal/exponential distribution")
    ax[1].legend()

    ax[2].plot(mu_tracker, label = "mu estimate")
    ax[3].plot(sigma_tracker, label = "sigma estimate")
    ax[4].plot(lambda_tracker, label = "lambda estimate")

    ax[2].set_xlabel("Iterations")
    ax[2].set_ylabel("Estimated mu")
    ax[3].set_xlabel("Iterations")
    ax[3].set_ylabel("Estimated sigma")
    ax[4].set_xlabel("Iterations")
    ax[4].set_ylabel("Estimated lambda")
    
    fig.suptitle(f"EM for W{w}")

    plt.show()
    
def EM_plot(data, w , pi_new, mu, sigma, lambda_current):
    
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    x = np.linspace(data.min(), data.max(), 100)
    ax.hist(data, bins=50, density=True,color='gray', alpha=0.5, label='histogram of original data')
    ax.plot(x, pi_new * norm(mu, sigma**0.5).pdf(x), color='red', label='Gaussian component')
    ax.plot(x, (1-pi_new) * expon(scale = lambda_current).pdf(x), color='blue', label='Exponential component')
    
    ax.axvline(x = expon(scale = lambda_current).ppf(0.5), ls = "--", color = "green", label ="50th percentile exponential distribution")
    ax.axvline(x = norm.ppf(0.5, mu, sigma**0.5), ls = "--", color = "purple",label ="50th percentile normal distribution")

    ax.set_title(f'Gaussian-Exponential mixture EM estimate for W{w}')
    ax.legend()
    lg = ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

    plt.show()