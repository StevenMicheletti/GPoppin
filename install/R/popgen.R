#' popgen
#'
#' Perform multiple population genetic analyses on a Genepop file and flag loci that don't meet criteria
#' @param gpfile Genepop filename. REQUIRED.
#' @param maf Create minor allele frequency table for each locus. Default=TRUE.
#' @param geno.s Determine genotype success of each locus. Default=TRUE.
#' @param fail.pops Flag loci that aren't sequenced in at least N population. Default=TRUE.
#' @param pop.thres Number of populations a locus must be present in. If fail.pops = TRUE. Default = 1.
#' @param loc.stats Calculate stats such as He, Ho, He-Ho, FIS for each locus. Default = TRUE.
#' @param maf.fz Flag loci that have a minor allele frequency less than N. Default = 0.05.
#' @param pop.table Create a summary table of population level analyses (Ho,He, Ho-He, FIS, GenoS,N). Default=TRUE.
#' @param hwe.calc Calculate Hardy-Weinberg equilibrium for each locus. Default = FALSE.
#' @param iter Number of iterations to perform in HWE analysis. If hwe.calc = TRUE. Default = 100.
#' @param p.pop Flag loci that are out of HWE in N proportion of populations. If hwe.calc = TRUE. Default = 0.33.
#' @param p.thres p-value threshold for HWE tests. If hwe.calc = TRUE. Default = 0.05.
#' @param calc.fst Calculate FST matrix (by population). Default = FALSE.
#' @param pca Perform an individual-level principal component analysis. Default = TRUE.
#' @param fz.table Create an allele frequency table for each locus. Default = TRUE.
#' @keywords popgen
#' @export
#' @examples
#' popgen(gpfile="Populations.gen", maf= TRUE, geno.s= TRUE, fail.pops= TRUE, pop.thres=1,
#' loc.stats= TRUE, maf.fz=0.05, pop.table= TRUE, hwe.calc=TRUE,
#' iter=100, p.pop=0.33, p.thres=0.05, calc.fst= TRUE, pca=TRUE, fz.table= TRUE)

popgen <- function (gpfile,
                    maf = TRUE,
                    geno.s = TRUE,
                    fail.pops = TRUE,
                    pop.thres = 1,
                    loc.stats = TRUE,
                    maf.fz = 0.05,
                    fail.loci = TRUE,
                    pop.table = TRUE,
                    hwe.calc = FALSE,
                    iter = 100,
                    p.pop = 0.33,
                    p.thres = 0.05,
                    calc.fst = FALSE,
                    pca=TRUE,
                    fz.table = TRUE
){

  
  #Install these packages if they aren't currently installed
  suppressMessages(require("adegenet"))
  suppressMessages(require("data.table"))
  suppressMessages(require("pegas"))
  suppressMessages(require("ggplot2"))
  suppressMessages(require("hierfstat"))
  
print ("Reading original genepop file")
system.time(GP <- fread(gpfile,stringsAsFactors=FALSE, sep="\t", header= F))
GP <- as.data.frame(GP[-1,])
colnames(GP)[1] <- "GP"
zis <-grep('pop', GP$GP, ignore.case =TRUE)
zis <-as.data.frame(zis)
zis <- as.data.frame(zis[1,1])
L = as.numeric((zis))
system.time(GP1 <- fread(gpfile,stringsAsFactors=FALSE, nrows= L, sep=" ", fill=TRUE, header= F, quote= "\""))
GP1 <-as.data.frame(GP1)
GP1 <- (GP1[-1, ])
GP1 <- as.data.frame(GP1)
colnames(GP1)[1] <- "V1"
keeps <- "V1"
GP1 = GP1[keeps]
print("Grabbing genotype information from genepop file")

if (ncol(GP) > 1){
  GPH <- fread(gpfile, header = F, skip = L , sep="\t", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
  GP2 <- fread(gpfile, header = F, skip = L+1 , sep= "\t", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
  GP2 <- rbind(GPH[1,],GP2, fill=TRUE)
  GP2 [is.na(GP2)] <- ""
  
}

if (ncol(GP) < 2){
  GPH <- fread(gpfile, header = F, skip = L , sep=" ", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
  GP2 <- fread(gpfile, header = F, skip = L+1 , sep= " ", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
  GP2 <- rbind(GPH[1,],GP2, fill=TRUE)
  GP2 [is.na(GP2)] <- ""
  
}

if (ncol(GP2) < 2){
  GPH <- fread(gpfile, header = F, skip = L , sep="\t", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
  GP2 <- fread(gpfile, header = F, skip = L+1 , sep= "\t", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
  GP2 <- rbind(GPH[1,],GP2, fill=TRUE)
  GP2 [is.na(GP2)] <- ""
}

getnum <- GP2$V3
gpcode = nchar(getnum[2])/2

ptm <- proc.time()
#Load file
gdata <- read.genepop(gpfile, ncode= gpcode)

#Split file into populations
pops <-seppop(gdata)
varlist <-names(pops)

#Get N population and N loci
P = as.numeric(length(pops))
L = as.numeric(length(gdata$loc.n.all))


if(P < 2) stop ("You need >1 population to run popgen")


if (fz.table == TRUE){

  gpop2 <- genind2genpop(gdata)
  temp1 <-makefreq(gpop2, missing="mean")
  temp2 <- t(temp1)
  temp2 <- cbind(Loci= rownames(temp2), temp2)

  #Save the pop frequency file
  fzout= sub('\\..*', '', gpfile)
  fzout = paste0(fzout,"_freqtable.txt")
  fzout = paste(fzout, collapse=" ")

  write.table(temp2, file= "fztemp", sep='\t',
              row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
  file.rename("fztemp", fzout)
  print ("Look for freq table")

}

if( pca == TRUE){

  #pop structure

  gpop3 <- genind2genpop(gdata)
  temp1 <-makefreq(gpop3, missing="mean")
  ind <- seq(2, ncol(temp1), by=2)
  temp2 <- temp1[ ,ind]
  temp1<-as.data.frame(temp1)
  temp3<- as.data.frame(temp2)
  corsnps<-gsub( "\\..*", "", names(temp3))
  corsnps2<- as.data.frame(corsnps)
  I =nrow(corsnps2)
  I = as.numeric(I)
  corpops<-gsub( "_.*$", "", rownames(temp3))
  corpops<-as.data.frame(corpops)
  row.names(temp3) <- corpops$corpops
  names(temp3) <-corsnps
  popsz <-(row.names(temp3))
  popsz <-as.factor(popsz)
  temp2<-as.matrix(temp3)
  pca1 <- dudi.pca(temp2, cent=TRUE, scale=FALSE, scannf=FALSE, nf=2)
  #Visualize eigenvalues for PCA axes
  #barplot(pca1$eig[1:50],main="PCA eigenvalues", xlab= "PC", col=heat.colors(50))
  #in general, eig represent the amount of genetic diversity. First two axes have a bunch of indo

  #Get PCA table
  #write.table(pca1$li, file= "pcatab.txt", sep='\t',
  #            row.names=TRUE, quote=FALSE)

  ax1 <-round((pca1$eig[1]/sum(pca1$eig)*100), digits=1)
  ax2 <-round((pca1$eig[2]/sum(pca1$eig)*100), digits=1)
  ax1 <- paste ("PC1 ","(", ax1,"%",")", sep= "")
  ax2 <- paste ("PC2 ","(", ax2,"%",")", sep= "")


  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_pop_PCA.pdf")
  gpout = paste(gpout, collapse=" ")

  #Plot PCA
  pdf(gpout, width = 8, height = 6)
      colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab= ax1, ylab= ax2, col= rainbow(2), col.lab= "black", col.axis= "Black")
      text(pca1$li, labels= rownames(pca1$li), cex= 0.5, pos=3)
      myCol <- colorplot(pca1$li, pca1$li, cex=3, transp=TRUE, add=TRUE)
  dev.off()

  
  #Get loading scores squared
  top_pop = as.data.frame(pca1$c1)
  top_pop$CS1 = (top_pop$CS1)^2
  top_pop = as.data.frame(top_pop[order(-top_pop[[1]]),])
  top_pop_out = as.data.frame(row.names(top_pop))
  top_pop_out[[2]] = top_pop[[1]]
  colnames(top_pop_out) <- c("Locus", "SNP_score2")
  
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_pop_PCA_score.txt")
  gpout = paste(gpout, collapse=" ")
  
  write.table(top_pop_out, file= gpout, sep='\t', col.names = TRUE,
              row.names=FALSE, quote=FALSE)
  
  ##individual PCA

  X <- tab(gdata, NA.method="mean")
  x <- X[complete.cases(X), ]
  X <- na.omit(X)
  
  #Remove monomorphic
  #X<- X[sapply(X, function(x) length(unique(x))>1)]

  pca1 <- dudi.pca(X, cent= FALSE, scale =FALSE, scannf= FALSE)

  #plot pca

  #s.class(pca$scores, pop(gdat)
  #abline(h=0,v=0, col="grey")

  ax1 <-round((pca1$eig[1]/sum(pca1$eig)*100), digits=1)
  ax2 <-round((pca1$eig[2]/sum(pca1$eig)*100), digits=1)
  ax1 <- paste ("PC1 ","(", ax1,"%",")", sep= "")
  ax2 <- paste ("PC2 ","(", ax2,"%",")", sep= "")

  ax3<-round((pca1$eig[3]/sum(pca1$eig)*100), digits=1)

  #Plot PCA
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_ind_PCA.pdf")
  gpout = paste(gpout, collapse=" ")

  pdf(gpout, width = 14, height = 12)
      colorplot(pca1$li, pca1$li, transp=TRUE, cex=2, xlab= ax1, ylab= ax2, col= rainbow(2), col.lab= "black", col.axis= "Black")
      text(pca1$li, labels= rownames(pca1$li), cex= 0.25, pos=3)
  dev.off()

  #Get loading scores squared
  top_ind = as.data.frame(pca1$c1)
  top_ind$CS1 = (top_ind$CS1)^2
  top_ind = as.data.frame(top_ind[order(-top_ind$CS1),])
  top_ind_out = as.data.frame(row.names(top_ind))
  top_ind_out[[2]] = top_ind[[1]]
  colnames(top_ind_out) <- c("Locus.Allele", "SNP_score2")
  
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_ind_PCA_score.txt")
  gpout = paste(gpout, collapse=" ")
  
  write.table(top_ind_out, file= gpout, sep='\t', col.names = TRUE,
              row.names=FALSE, quote=FALSE)  
  
}


# Check number of alleles per SNP (should be 2)
allele <- nAll(gdata)
max(allele)
min(allele)
###

#MAF per locus
mfreq <- minorAllele(gdata)
mfreq <-as.data.frame(mfreq)

#MAF per pop
maftable <- data.frame(matrix(ncol= P, nrow = L))
colnames(maftable) <- c(1:P)

print ('Calculating minor allele frequencies')
for (i in (1:P)) {
  mfreq2 <- minorAllele(pops[[i]])
  mfreq2 <- as.data.frame(mfreq2)
  maftable[[i]] = mfreq2$mfreq2
  popn <- names(pops)[i]
  colnames(maftable)[i] = popn
}

#add locus names
rownames(maftable) <- levels(gdata$loc.fac)
#calculate mean
maftable$mean <- mfreq$mfreq
maftable[maftable=="numeric(0)"] <-NA

if (maf == TRUE) {
  write.table(as.matrix(maftable), file= "tempAF", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_maft.txt")
  gpout = paste(gpout, collapse=" ")
  file.rename("tempAF", gpout)
  print ("Minor allele frequencies printed (_maft.txt)")
}

########

print ("Calculating summary stats")
pgen <- summary(gdata)

#Ho vs He Plots

#plot(pgen$Hobs, xlab="Loci number", ylab="Observed Heterozygosity",
#     main="Observed heterozygosity per locus")

gpout= sub('\\..*', '', gpfile)
gpout = paste0(gpout,"_he_vs_ho.pdf")
gpout = paste(gpout, collapse=" ")

pdf(gpout, width = 8, height = 6)
  plot(pgen$Hobs,pgen$Hexp, xlab="Locus Hobs", ylab=" Locus Hexp",
     main="Expected heterozygosity vs. Observed heterozygosity", col = "red")
dev.off()

#print(" Does He = Ho?")
#bartlett.test(list(pgen$Hexp, pgen$Hobs))

#Get basic popgen stats
print ("Getting popgen stats")

basicstat <- basic.stats(gdata, diploid = TRUE, digits = gpcode)
# Number of individuals genotyped for each loci at each pop
inds <- basicstat$n.ind.samp
inds[is.na(inds)] <- 0
inds<- as.data.frame(inds)
colnames(inds) <- c(1:P)

suctable <- data.frame(matrix(ncol= P, nrow = L))
colnames(suctable) <- c(1:P)

print ("Calculating genotype success %")
for (i in (1:P)) {
  exp <-nrow(pops[[i]]@tab)
  exp <-as.numeric(exp)
  suctable[[i]] <- inds[[i]]/exp
  popn <- names(pops)[i]
  colnames(suctable)[i] = popn
}


nsum <- as.data.frame (rowSums(inds,))
nsum$total <- pgen$n
nsum$success <- nsum$`rowSums(inds, )`/nsum$total
rownames(nsum) <- levels(gdata$loc.fac)
rownames(suctable) <- levels(gdata$loc.fac)
suctable2 <- suctable
suctable2$mean <- nsum$success


if (geno.s == TRUE) {
  write.table(as.matrix(suctable2), file= "tempgs", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_genos.txt")
  file.rename("tempgs", gpout)
  print ("Done with Genotype success % (_genos.txt)")
}


fail <- as.data.frame (rowSums (suctable == 0))
colnames(fail) <- ("V1")

if (sum(fail$V1 == 0)) {
  print("Congratulations, no loci are completely missing in any population")
}

#hist(fail$V1, main = "Non-genotyped pop distribution", xlab = "Number of pops", ylab = "number of Snps")

failed <- (which(fail$V1 > pop.thres))
failed <- as.data.frame(failed)
fail$order <-1:nrow(fail)
failedl <- subset(fail, fail$order %in% failed$failed)
failedn <- as.data.frame (rownames(failedl))

if(fail.pops == TRUE){
  write.table(failedn, file= "tempf", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_failed_geno.txt")
  file.rename("tempf", gpout)
  print ("Made list of SNPs that didn't genotype in  pops (_faileD_geno.txt")
}

# # Individuals

indtable <- data.frame(matrix(ncol= P, nrow = 1))
colnames(suctable) <- 1

for (i in (1:P)) {
  nind <- nrow(pops[[i]]@tab)
  nind <-as.numeric(nind)
  popn <- names(pops)[i]
  indtable[[i]] <- nind
  colnames(indtable)[i] = popn
}

indtable <- as.data.frame(t(indtable))
colnames(indtable) = "N"

#write.table(as.matrix(indtable), file= "temp", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
#gpout = c(gfile,".Nind")
#gpout = paste(gpout, collapse=" ")
#file.rename("temp", gpout)
#print ("Done with individual number")

## Get FIS , Ho, and He
fist <- as.data.frame(basicstat$Fis)
fist$mean <- rowMeans(fist, na.rm =TRUE)
rownames(fist) <- levels(gdata$loc.fac)

hot <- as.data.frame(basicstat$Ho)
hot$mean <- rowMeans(hot, na.rm =TRUE)
rownames(hot) <- levels(gdata$loc.fac)

het <- as.data.frame(basicstat$Hs)
het$mean <- rowMeans(het, na.rm =TRUE)
rownames(het) <- levels(gdata$loc.fac)


if (loc.stats==TRUE) {
  write.table(fist, file= "temp1", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_loc_FIS.txt")
  file.rename("temp1", gpout)

  write.table(hot, file= "temp2", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)

  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_loc_ho.txt")
  file.rename("temp2", gpout)

  write.table(het, file= "temp3", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)

  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_loc_he.txt")
  file.rename("temp3", gpout)

  print ("Basic diversity measures for lic calculated (_ho, _he, _FIS")

  print ("Overall stats")
}


basicstat$overall

#boot.ppfis(gdata)

#Calculate HW disequilbrium
distable <- hot - het

if (loc.stats==TRUE) {
  write.table(distable, file= "temp4", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_loc_dis.txt")
  file.rename("temp4", gpout)
  print ("Ho - He (disequilibrium) table wrote(_dis)")

}

#Fin plot

print ('Plot: MAF Vs Ho-He')

gpout= sub('\\..*', '', gpfile)
gpout = paste0(gpout,"_fin_plot.pdf")
gpout = paste(gpout, collapse=" ")

pdf(gpout, width = 8, height = 6)
    plot(maftable$mean, distable$mean, xlab = "MAF", ylab = "Hardy-Weinberg disequilibrium", 
         main ="FIN Plot", col = "blue")
dev.off()

#Fin plot - remove HIGH HWE DIS

fintable <- data.frame(matrix(ncol= 1, nrow = L))
colnames(fintable)[1] <- "V1"
fintable$V1 <- rownames(distable)
fintable$V2 <- distable$mean
fintable$V3 <- maftable$mean
# fintable$V4 <- fintable$V2 * fintable$V3

# Additional MAF filter

mafilt <- subset(fintable,fintable$V3 < maf.fz)


alert1 <- paste0(nrow(mafilt), " loci targeted after maf filter of ", maf.fz)
print(alert1)

if(fail.loci == TRUE) {
  write.table(mafilt$V1, file= "temfp", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,".MAF_filt_fail.txt")
  file.rename("temfp", gpout)



  print ("Fin table filter")

  ####  #     #     #
  #     #    # #   #
  ###   #   #   # #
  #     #  #    ##

  fintab1 <- subset(fintable,fintable$V3 < 0.1 & fintable$V2 > 0.1 | fintable$V3 < 0.1 & fintable$V2 < - 0.15  )
  fintab2 <- subset(fintable,fintable$V3 < 0.2 & fintable$V2 > 0.175 | fintable$V3 < 0.2 & fintable$V2 < - 0.225  )
  fintab3 <- subset(fintable,fintable$V3 < 0.3 & fintable$V2 > 0.225 | fintable$V3 < 0.3 & fintable$V2 < - 0.25  )
  fintab4 <- subset(fintable,fintable$V3 < 0.4 & fintable$V2 > 0.25 | fintable$V3 < 0.4 & fintable$V2 < - 0.275  )
  fintab5 <- subset(fintable,fintable$V3 < 0.5 & fintable$V2 > 0.3 | fintable$V3 < 0.5 & fintable$V2 < - 0.3  )
  fintab6 <- subset(fintable,fintable$V3 < 0.01)

  ftabout <- rbind(fintab1, fintab2, fintab3, fintab4, fintab5, fintab6)


  write.table(ftabout$V1, file= "tempfy", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_fin_fail.txt")
  file.rename("tempfy", gpout)

  alert2<- paste0(nrow(ftabout), " loci targeted for failing fin filter")
  print (alert2)

  #Get singletons
  singleton <- data.frame(matrix(ncol= 1, nrow = L))
  colnames(singleton)[1] <- "V1"
  singleton$V1 <- rownames(distable)
  singleton$V2 <- hot$mean
  singleton$V3 <- maftable$mean
  sintab1 <- subset(singleton,singleton$V2 ==0 & singleton$V3 < .02)
  sintab2 <- subset(singleton,singleton$V2 > 0.25  & singleton$V3 < .02)

  sintab <- rbind(sintab1, sintab2)

  write.table(sintab$V1, file= "temps", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)

  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_singleton.txt")
  file.rename("temps", gpout)

  alert2<- paste0(nrow(sintab), " loci targeted as singletons")
  print (alert2)

}

#Create Population table summary


poptable <- data.frame(matrix(ncol= 1, nrow = P))
colnames(poptable) <- "V1"

poptable$V1 = rownames(indtable)
hot$mean = NULL
het$mean = NULL
distable$mean = NULL
fist$mean = NULL
suctable2$mean = NULL


poptable$V2 = colMeans (hot, na.rm = TRUE, dims = 1)
poptable$V3 = colMeans (het, na.rm = TRUE, dims = 1)
poptable$V4 = colMeans (distable, na.rm = TRUE, dims = 1)
poptable$V5 = colMeans (fist, na.rm = TRUE, dims = 1)
poptable$V6 = colMeans (suctable2, na.rm = TRUE, dims = 1)
poptable$V7 = indtable$N

colnames(poptable) <- c("Pop", "Ho", "He", "Ho-He", "FIS", "GenoS", "N" )


if (pop.table == TRUE) {
  write.table(poptable, file= "temppt", sep='\t',row.names= FALSE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_pop_table.txt")
  file.rename("temppt", gpout)
  print("Population table created")
}


#  #  #     # #####
#  #  #     # #
####  #  #  # #####
#  #  #  #  # #
#  #  ####### #####

#Calulate HWE and make a list of loci that do not meet threshold

if (hwe.calc == TRUE) {

  print('HWE exact test for all populations... this will take awhile')
  htable <- data.frame(matrix(ncol= P, nrow = L))
  colnames(htable) <- c(1:P)

  for (i in (1:P)) {
    HWE <-hw.test(pops[[i]], B= iter)
    HWE <-as.data.frame(HWE)
    popn <- names(pops)[i]
    htable[[i]] = HWE$Pr.exact
    colnames(htable)[i] = popn
    progress = (i / P) * 100
    progress = round(progress, digits=2)
    print ( progress)
  }

  rownames(htable) <- levels(gdata$loc.fac)

  write.table(htable, file= "temphw", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_hwe.txt")
  file.rename("temphw", gpout)
  print("Uncorrected HWE Pvalue table printed")


  ### P-value Corrections



  btable <- data.frame(matrix(ncol= P, nrow = L))
  colnames(btable) <- c(1:P)

  for (i in (1:P)) {
    b <-p.adjust(htable[[i]], "bonferroni")
    b <-as.data.frame(b)
    popn <- names(pops)[i]
    btable[[i]] = b$b
    colnames(btable)[i] = popn
  }

  rownames(btable) <- levels(gdata$loc.fac)


  ftable <- data.frame(matrix(ncol= P, nrow = L))
  colnames(ftable) <- c(1:P)

  for (i in (1:P)) {
    b <-p.adjust(htable[[i]], "fdr")
    b <-as.data.frame(b)
    popn <- names(pops)[i]
    ftable[[i]] = b$b
    colnames(ftable)[i] = popn
  }

  rownames(ftable) <- levels(gdata$loc.fac)


  ytable <- data.frame(matrix(ncol= P, nrow = L))
  colnames(ytable) <- c(1:P)

  for (i in (1:P)) {
    b <-p.adjust(htable[[i]], "BY")
    b <-as.data.frame(b)
    popn <- names(pops)[i]
    ytable[[i]] = b$b
    colnames(ytable)[i] = popn
  }

  rownames(ytable) <- levels(gdata$loc.fac)

  print ('P-value corrections calculated')


  htable[is.na(htable)] <- " "
  ftable[is.na(ftable)] <- " "
  btable[is.na(btable)] <- " "
  ytable[is.na(ytable)] <- " "

  ###
  per <- as.data.frame (rowSums (htable < p.thres))
  htable$out <- per$`rowSums(htable < p.thres)` / P
  htable$order<- 1:nrow(htable)
  htable$name <- rownames(htable)


  per <- as.data.frame (rowSums (ftable < p.thres))
  ftable$out <- per$`rowSums(ftable < p.thres)` / P
  ftable$order<- 1:nrow(ftable)
  ftable$name <- rownames(ftable)

  per <- as.data.frame (rowSums (btable <p.thres))
  btable$out <- per$`rowSums(btable < p.thres)` / P
  btable$order<- 1:nrow(btable)
  btable$name <- rownames(btable)


  per <- as.data.frame (rowSums (ytable <p.thres))
  ytable$out <- per$`rowSums(ytable < p.thres)` / P
  ytable$order<- 1:nrow(ytable)
  ytable$name <- rownames(ytable)

  #Uncorrected loci that don't meet thredhold
  badh <- (which(htable$out > p.pop))
  badh <- as.data.frame(badh)
  badhl <- subset(htable, htable$order %in% badh$badh)

  #Bonferoni loci that don't meet threshold
  badb <- (which(btable$out > p.pop))
  badb <- as.data.frame(badb)
  badbl <- subset(btable, btable$order %in% badb$badb)

  #BY loci that don't meet thredhols
  bady <- (which(ytable$out > p.pop))
  bady <- as.data.frame(bady)
  badyl <- subset(ytable, ytable$order %in% bady$bady)

  #FDR loci that don't meet thredhols
  badf <- (which(ftable$out > p.pop))
  badf <- as.data.frame(badf)
  badfl <- subset(ftable, ftable$order %in% badf$badf)

  write.table(badhl$name, file= "temp", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  gpout= sub('\\..*', '', gpfile)
  gpout = c(gpout,"_badHWE.txt")
  gpout = paste0(gpout, collapse="")
  file.rename("temp", gpout)

  
  write.table(badbl$name, file= "temp", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  gpout= sub('\\..*', '', gpfile)
  gpout = c(gpout,"_badBON.txt")
  gpout = paste0(gpout, collapse="")
  file.rename("temp", gpout)


  write.table(badyl$name, file= "temp", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  gpout= sub('\\..*', '', gpfile)
  gpout = c(gpout,"_badBY.txt")
  gpout = paste0(gpout, collapse="")
  file.rename("temp", gpout)
  
  write.table(badfl$name, file= "temp", sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  gpout= sub('\\..*', '', gpfile)
  gpout = c(gpout,"_badFDR.txt")
  gpout = paste0(gpout, collapse="")
  file.rename("temp", gpout)

  print ("Bad HWE loci lists created")

  ctable <- data.frame(matrix(ncol= 4, nrow = L))
  colnames(ctable) <- c('HWEout', 'BONout', 'BYout', 'FDRout')

  ctable$HWEout <- htable$out
  ctable$BONout <- btable$out
  ctable$BYout <- ytable$out
  ctable$FDRout <- ftable$out
  ccomb <- as.data.frame(rowSums(ctable > p.pop))
  ctable$calls <-ccomb$`rowSums(ctable > p.pop)`
  ctable$SNP <- rownames(htable)

  write.table(ctable, file= "temp", sep='\t',row.names= FALSE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = c(gpout,"_HWEcomp.txt")
  gpout = paste0(gpout, collapse="")
  file.rename("temp", gpout)

  print ("HWE propotion out comparison table created")

}



if (calc.fst == TRUE){
  #FST
  # boostrapped FST (Takes a long time)

  print ("Calculating fst matrix")
  fst.mat <-pairwise.fst(gdata,res.type="matrix")

  write.table(fst.mat, file= "temp", sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
  gpout= sub('\\..*', '', gpfile)
  gpout = c(gpout,".fst")
  gpout = paste0(gpout, collapse="")
  file.rename("temp", gpout)

  print ("FST created")
}

proc.time() - ptm

}
