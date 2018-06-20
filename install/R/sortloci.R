#' sortloci
#'
#' Keep or remove a list of loci from a Genepop file
#' @param gpfile Name of Genepop file that will be filtered. REQUIRED.
#' @param locus.file Name of file that is a list of locus names. REQUIRED.
#' @param option Indicate if you want to keep or remove loci in the locus.file. Either "keep" or "remove". REQUIRED.
#' @keywords sortloci
#' @export
#' @examples
#' sortloci(gpfile = "Population.gen", locus.file= "significant_snps.txt", option = "keep")


sortloci <- function (gpfile,
                      locus.file,
                      option
) {
  suppressMessages(require(data.table))

  ptm <- proc.time()

  print ("Reading original genepop file")
  system.time(GP <- fread(gpfile,stringsAsFactors=FALSE, sep="\t", header= F))
  GP <- as.data.frame(GP[-1,])
  colnames(GP)[1] <- "GP"
  zis <-grep('pop', GP$GP, ignore.case =TRUE)
  zis <-as.data.frame(zis)
  zis <- as.data.frame(zis[1,1])
  L = as.numeric((zis))

  #Loading loci
  loci= fread(locus.file, header = F)
  # preparing output name
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_loci_sort.gen")


  # Read in genepop file in two chunks: 1) list of loci, 2) individuals x loci
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
  
  GP2 <-as.data.frame(GP2)

  #check format and fix
  gptest <- GP2
  gptest$V2[gptest$V2 ==""] <-0
  check <-tail(names(sort(table(gptest$V2))), 1)
  newnames <- colnames(gptest)
  newnames <- head(newnames,-1)
  check2 <- tail(names(sort(table(gptest$V3))), 1)

  if (check == ",") {
    print ("Genepop Format is Bad... Fixing")
    GP2$V1 <- paste (GP2$V1, GP2$V2, sep="")
    GP2$V2 =NULL
    colnames(GP2) <- newnames
  }

  if (check2 == "") {
    print ("Genepop Format is Bad... Fixing")
    GP2$V3 = NULL
    colnames(GP2) <- newnames
  }

  if (check == "0") {
    print ("Genepop Format is Bad... Fixing")
    GP2$V2 = NULL
    colnames(GP2) <- newnames
  }

  #Determine loci position in genepop file
  GP1$V2 <- 1:nrow(GP1)
  GP1$V2 <- 1 + GP1$V2
  GPR <- GP1 [! GP1$V1 %in% loci$V1, ]
  GPR <- subset(GP1, GP1$V1 %in% loci$V1)
  GPR$V2 <- sub("","V",GPR$V2)
  GPR <- as.character(GPR$V2)
  GPR <- c(keeps, GPR)

  #Keep loci or remove
  if (option == "keep") {
    GP1$V2 <- 1:nrow(GP1)
    GP1$V2 <- 1 + GP1$V2
    GPR <- GP1 [! GP1$V1 %in% loci$V1, ]
    GPR <- subset(GP1, GP1$V1 %in% loci$V1)
    GPR$V2 <- sub("","V",GPR$V2)
    GPR <- as.character(GPR$V2)
    GPR <- c(keeps, GPR)
    locilist <- subset(GP1, GP1$V1 %in% loci$V1)
    locilist <- as.data.frame(locilist)
    locilist2 <- GP2[,names(GP2) %in% GPR]
  }

  if (option == "remove") {
    GP1$V2 <- 1:nrow(GP1)
    GP1$V2 <- 1 + GP1$V2
    GPR <- GP1 [! GP1$V1 %in% loci$V1, ]
    GPR <- subset(GP1, GP1$V1 %in% loci$V1)
    GPR$V2 <- sub("","V",GPR$V2)
    GPR <- as.character(GPR$V2)
    locilist <- GP1 [! GP1$V1 %in% loci$V1, ]
    locilist <- as.data.frame(locilist)
    locilist2 <- GP2[,!names(GP2) %in% GPR]
  }

  locilist<- locilist[,-2]
  locilist<- as.matrix(locilist)

  if (sum (na.omit(as.numeric(locilist2[,2]))) < 1) {
    print ("ERROR, something wrong with the file, check warnings and errors above")
  }

  if (sum (na.omit(as.numeric(locilist2[,2]))) == "NA") {
    print ("ERROR, something wrong with the file, check warnings and errors above")
  }

  if (sum (na.omit(as.numeric(locilist2[,2]))) > 1) {
    locilist2<- as.matrix(locilist2)
    #Reconstruct GenePop file
    header = "Sorted by sort.loci on"
    header = paste(header, Sys.time())
    print("Creating new file")
    cat(header, '\n',  file = "temp")
    write.table(locilist, file="temp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(locilist2, file="temp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    file.rename("temp", gpout)
    print ("Done! Look for loci_sort.gen file")
  }

  if (sum (na.omit(as.numeric(locilist2[,2]))) < 1) stop ("Something is wrong with the input file... check it")

  proc.time() - ptm
}
