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
GPoppin is not available on CRAN. Therefore, you must pre-install the following packages:
    
    > install.packages('adegenet') 
    > install.packages('data.table')  
    > install.packages('ggplot2')  
    > install.packages('hierfstat')  
    > install.packages('pegas')  

### Installing in R
set your working directory and :

> install.packages("GPoppin.tar", repos=NULL, type="source")   
> library("GPoppin")

## Example files 

There are multiple example files built into GPoppin which include:

ex_prog : example of progeny input from a .csv  
ex_loci: example loci input for filtering loci  
ex_match: example loci input for matching locus order  
ex_pop : example input for rearranging populations  
