---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# DBHC

Package DBHC is an implementation of a sequence clustering algorithm that uses a 
mixture of discrete-output hidden Markov models (HMMs), the Discrete Bayesian 
HMM Clustering (DBHC) algorithm. The algorithm uses heuristics based on the 
Bayesian Information Criterion (BIC) to search for the optimal number of hidden 
states in each HMM and the optimal number of clusters. The packages provides 
functions for finding clusters in discrete sequence data with the DBHC algorithm 
and for plotting heatmaps of the probability matrices that are estimated in the 
cluster models. 

## Example

Below a basic example of how to use package DBHC for obtaining sequence clusters
for the Swiss Household data in package TraMineR:


```r
library(DBHC)
library(TraMineR)

## Swiss Household Data
data("biofam", package = "TraMineR")

# Clustering algorithm
new.alphabet <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
sequences <- seqdef(biofam[,10:25], alphabet = 0:7, states = new.alphabet)

# Code below takes long time to run
res <- hmm.clust(sequences)

# Heatmaps
cluster <- 1  # display heatmaps for cluster 1
transition.heatmap(res$partition[[cluster]]$transition_probs,
                   res$partition[[cluster]]$initial_probs)
emission.heatmap(res$partition[[cluster]]$emission_probs)

## A smaller example, which takes less time to run
subset <- sequences[sample(1:nrow(sequences), 20, replace = FALSE),]

# Clustering algorithm
res <- hmm.clust(subset, K.max = 3)

# Number of clusters
print(res$n.clusters)

# Table of cluster memberships
table(res$memberships[,"cluster"])

# BIC for each number of clusters
print(res$bic)

# Heatmaps
cluster <- 1  # display heatmaps for cluster 1
transition.heatmap(res$partition[[cluster]]$transition_probs,
                   res$partition[[cluster]]$initial_probs)
emission.heatmap(res$partition[[cluster]]$emission_probs)
```
