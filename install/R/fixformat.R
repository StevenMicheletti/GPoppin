#' fixformat
#'
#' Fixes common whitespace and invalid character errors in genepop files
#' @param gpfile Genepop file. REQUIRED.
#' @keywords fixformat
#' @export
#' @examples
#' fixformat(gpfile="Populations.gen")


fixformat <- function (gpfile
) {

  suppressMessages(require("data.table"))

  #Read in genepop file and get dimensions
  print ("Reading original genepop file")
  system.time(GP <- fread(gpfile,stringsAsFactors=FALSE, sep="\t", header= F))
  GP <- as.data.frame(GP[-1,])
  colnames(GP)[1] <- "GP"
  zis <-grep("\\pop\\b|\\Pop\\b|\\POP\\b", GP$GP, ignore.case =TRUE)
  zis <-as.data.frame(zis)
  zis <- as.data.frame(zis[1,1])
  L = as.numeric((zis))

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

  #Remove problematic characters from locus names
  GP1<-as.data.frame(gsub("Pop|POP|pop|PoP|poP", "pep", GP1$V1))
  colnames(GP1) <-"name"
  GP1<-as.data.frame(gsub("`|'|[.]|:", "$", GP1$name))
  colnames(GP1) <-"name"
  print("Problematic characters in locus names were replaced with $")
  print("Any occurance of pop in locus names were replaced with pep")

  locilist <- GP1
  locilist <- as.data.frame(locilist)
  locilist2 <- GP2

  #Reconstruct GenePop file
  header = "Fixed by fixformat on"
  header = paste(header, Sys.time())

  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_fixed.gen")

  print("Creating new file")
  cat(header, '\n',  file = "tempff")
  write.table(locilist, file="tempff", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(locilist2, file="tempff", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)

  file.rename("tempff", gpout)
  print ("Done! Look for _fixed file")

}
