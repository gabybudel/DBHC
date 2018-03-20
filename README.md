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
library(TraMineR)
#> 
#> TraMineR stable version 2.0-8 (Built: 2018-03-15)
#> Website: http://traminer.unige.ch
#> Please type 'citation("TraMineR")' for citation information.

## Swiss Household Data
data("biofam", package = "TraMineR")

# Clustering algorithm
new.alphabet <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
sequences <- seqdef(biofam[,10:25], alphabet = 0:7, states = new.alphabet)
#>  [>] state coding:
#>        [alphabet]  [label]  [long label]
#>      1             0P        P
#>      2             1L        L
#>      3             2M        M
#>      4             3LM       LM
#>      5             4C        C
#>      6             5LC       LC
#>      7             6LMC      LMC
#>      8             7D        D
#>  [>] 2000 sequences in the data set
#>  [>] min/max sequence length: 16/16
res <- hmm.clust(sequences)
#> Error in hmm.clust(sequences): could not find function "hmm.clust"

# Heatmaps
cluster <- 1  # display heatmaps for cluster 1
transition.heatmap(res$partition[[cluster]]$transition_probs,
                   res$partition[[cluster]]$initial_probs)
#> Error in transition.heatmap(res$partition[[cluster]]$transition_probs, : could not find function "transition.heatmap"
emission.heatmap(res$partition[[cluster]]$emission_probs)
#> Error in emission.heatmap(res$partition[[cluster]]$emission_probs): could not find function "emission.heatmap"
```
