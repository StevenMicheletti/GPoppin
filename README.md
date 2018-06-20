# GPoppin
GPoppin, an R package for editing, simulating, and analyzing genotypes in Genepop format

### Functions

prog2gp - Convert Progeny csv/txt exports (1 allele per columns format) into Genepop format.  
fixformat - Fix special characters and locus names in Genepop that analysis tools have issues with.  
popgen - Perform a suite of analyses with adegenet, generate summaries, QC loci.  
sortloci - Remove or keep loci based on external list.  
Buildpops -  Define populations and individuals that belong to them. Good for subsetting individuals.  
matchloci - Rearrange locus order and content based on a user-provided list.  
genosim - Simulate individual genotypes based on allele frequencies. Can incorporate Ho-He.  

## Installation

GPoppin is intended to be installed on desktop computers, but also functions on Linux-based servers.  

### Dependencies
GPoppin is not available on CRAN. Therefore, you must pre-install the following packages if you do not have devtools installed:
    
  >install.packages(c('adegenet', 'data.table', 'ggplot2', 'hierfstat', 'pegas'))
  
### Installing in R
set your working directory and :

> install.packages("GPoppin.tar", repos=NULL, type="source")   
> library("GPoppin")

### Alternative Install with Devtools

Unpack GPoppin.tar into a directory and: 

> install("GPoppin",dependencies = TRUE)  

devtools will automatically install required dependencies

## Example files 

There are multiple example files built into GPoppin which include:

ex_prog : example of progeny input from a .csv  
ex_loci: example loci input for filtering loci  
ex_match: example loci input for matching locus order  
ex_pop : example input for rearranging populations  

## Troubleshooting

GPoppin relies on multiple external dependencies. Any problems are likely to stem from these external packages.  

adegenet, in particular, is paired with many dependencies. On occasion an associated dependency will not properly install. If this occurs, installation of GPoppin will exit with an error, indicating the name of the package that failed. Usually this issue can be resolved by simply installing the package that failed. If prompted "Do you want to install from sources the package which needs compilation?," select no. 
