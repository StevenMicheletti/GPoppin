# GPoppin
GPoppin a R package for editing, simulation, and analysis using Genepop format

Functions

buildpops - Specify which populations individuals belong to. Also good for subsetting individuals  
fixformat - Fix special characters and locus names that analysis tools frequency have issues with  
genosim - Simulate individual genotypes based on allele frequencies. Can incorporate He-Ho   
matchloci - Rearrange loci based on a provided list  
popgen - Perform a suite of analyses with adegenet, generate summaries  
prog2gp - Convert Progeny exports to genepop format  
sortloci - Remove or keep loci based on external list.  

## Installation

### Dependencies
GPoppin is not yet available on CRAN. Therefore you must have pre-installed the following packages:
    
    > install.packages('adegenet') 
    > install.packages('data.table')  
    > install.packages('ggplot2')  
    > install.packages('hierfstat')  
    > install.packages('pegas')  

### Installing in R
set your working directory and 

> install.packages("GPoppin", repos=NULL, type="source")   
> library("GPoppin")

## Example files 

There are multiple example files built into GPoppin which include:

ex_prog : example of progent input from a .csv
ex_loci: example loci input for filtering loci
ex_match: example loci input for matching locus order
ex_pop : example input for rearranging populations
