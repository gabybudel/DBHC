---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# DBHC

Package DBHC is an implementation of the DBHC algorithm, an HMM sequence 
clustering algorithm that finds a mixture of discrete-output HMMs. The algorithm 
uses heuristics based on BIC to search for the optimal number of hidden states 
ineach HMM and the optimal number of clusters. The packages provides functions 
for finding clusters in discrete sequence data with the DBHC algorithm and for 
plotting heatmaps of the probability matrices that are estimated in the cluster 
models. 

## Example

Below a basic example of how to use package DBHC for obtaining sequence clusters
for the Swiss Household data:


```r
## Swiss Household Data
data("biofam", package = "TraMineR")

# Clustering algorithm
new.alphabet <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
sequences <- seqdef(biofam[,10:25], alphabet = 0:7, states = new.alphabet)
#> Error in seqdef(biofam[, 10:25], alphabet = 0:7, states = new.alphabet): could not find function "seqdef"
res <- hmm.clust(sequences)
#> Error in "stslist" %in% class(sequences): object 'sequences' not found

# Heatmaps
cluster <- 1  # display heatmaps for cluster 1
transition.heatmap(res$partition[[cluster]]$transition_probs,
                   res$partition[[cluster]]$initial_probs)
#> Error in transition.heatmap(res$partition[[cluster]]$transition_probs, : object 'res' not found
emission.heatmap(res$partition[[cluster]]$emission_probs)
#> Error in emission.heatmap(res$partition[[cluster]]$emission_probs): object 'res' not found
```
