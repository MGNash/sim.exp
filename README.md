# sim.exp
R code that simulates proteomics experiments for prospective power analysis.

I wrote this code for a paid internship with the Biostatistics and Design Program at Oregon Health and Science University. I am its sole author. I am sharing it with the Github community with their permission.

The function 'sim.exp' simulates the results of experiments measuring associations between amounts of specific proteins in biological samples and a phenotype measurement, but could be used for prospective power analysis for any type of experiment in which one is testing for pairwise associations between one measure and many others. A permutation test is used to assign p-values to sample correlations based on a null hypothesis of no association between protein and phenotype. The Benjamini-Hochberg procedure is used to select significant associations between phenotype and protein measures at a designated false discovery rate (FDR). 

The number of subjects, proteins, true positives, and times the experiment is repeated, the desired FDR, the strength of the association for those proteins associated with the phenotype, and the type of correlation being measured (Spearman or Pearson) can be set by arguments to the function 'sim.exp'. Simulated data are bivariate normal (i.e. normal distribution of protein and phenotype with linear association) with no missing values. For each simulated experiment, sim.exp returns the number of true and false positives, and it is easy to calculate from these the observed proportion of false positives.

Files contained in this repository:
simulate.PDF - contains an example showing how to use 'sim.exp' to simulate results for experiments, repeating much of the information found here.  
simulate.Rmd - contains the code used to generate 'simulate.PDF'  
simulate.R - contains definitions for two functions: 'sim.exp' and 'getcdf'. 'getcdf' is used by 'sim.exp' to assign p values to sample correlations.  

'sim.exp' takes the following arguments:

nrep: number of simulated experiments  
nsubj: number of subjects in each experiment  
nprot: number of protein measurements in each experiment  
ntrue: number of protein measurements associated with the phenotype in the population (from zero to nprot)  
FDR: the false discovery rate one wishes to use (from 0 to 1)  
spearman: correlations are measured using Spearman correlation if TRUE, Pearson otherwise  

The defaults are as follows: nrep=1,nsubj = 16,nprot=600,ntrue=30,FDR=.5,rho=.5,spearman=TRUE

'sim.exp' returns a data frame containing the number of true and false positives for each experiment, labeled as 'true_positive' and 'false_positive', respectively. To get basic summary statistics for results of simulated experiments without saving results, use the following syntax:

summary(sim.exp())

I encourage researchers to use my code freely, as long as they credit its author, Michael G. Nash.
