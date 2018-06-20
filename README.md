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
It can be installed locally or using install_github from devtools 

### Local Installation
Download GPoppin.tar into a directory. Prior to installing, you will need to have the following dependencies:  
    
  >install.packages(c('adegenet', 'data.table', 'ggplot2', 'hierfstat', 'pegas'))
  
Set your working director to the location of GPoppin.tar and :

> install.packages("GPoppin.tar", repos=NULL, type="source")   
> library("GPoppin")

### devtools/Github Install

Load devtools and run install_github

> library("Devtools")
> install_github("StevenMicheletti/GPoppin/install")  

OR

Unpack GPoppin.tar as a directory

> library("Devtools")  
> install("GPoppin",dependencies = TRUE)


## Example files 

There are multiple example files built into GPoppin which include:

ex_prog : example of progeny input from a .csv  
ex_loci: example loci input for filtering loci  
ex_match: example loci input for matching locus order  
ex_pop : example input for rearranging populations  

The easiest way to utilize example files is to begin with prog2gp:

> prog2gp(examp=T, split.pop=T, rem.thres=0.5)

Then use the produced file: ex_prog_out.gen for other functions.

## Troubleshooting

GPoppin relies on multiple external dependencies. Any problems are likely to stem from these external packages.  

adegenet, in particular, is paired with many dependencies. On occasion an associated dependency will not properly install. If this occurs, installation of GPoppin will exit with an error, indicating the name of the package that failed. Usually this issue can be resolved by simply installing the package that failed. If prompted "Do you want to install from sources the package which needs compilation?," select no. 
