## Test environments
* local Windows 10 install, R 3.4.4

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking R code for possible problems ... NOTE
  emission.heatmap: no visible binding for global variable 'symbol_names'
  emission.heatmap: no visible binding for global variable 'state_names'
  emission.heatmap: no visible binding for global variable 'value'
  transition.heatmap: no visible binding for global variable 'to'
  transition.heatmap: no visible binding for global variable 'from'
  transition.heatmap: no visible binding for global variable 'value'
  Undefined global functions or variables:
  from state_names symbol_names to value

## The above NOTE arises from the fact that I built plotting functions for 
(heatmaps of) the probability matrices that come out of the seqHMM package. 
These variable names will always be exactly as specified if the probability 
matrices come from the seqHMM package, but I could not find a way to make these 
variable names binding.