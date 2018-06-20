#' matchloci 
#'
#' Reorganizes loci in Genepop file to match the order of another list of loci.
#' @param gpfile Genepop filename to be reorganized. REQUIRED.
#' @param loci.file Filename of list of loci in desired order. REQUIRED. 
#' @keywords matchloci
#' @export
#' @examples
#' matchloci(gpfile="Population.gen", loci.file="Loci_by_significance.txt")


matchloci <- function (gpfile,
                       loci.file
){

  ptm <- proc.time()
  suppressMessages (require(data.table))
  print ("Reading original genepop file")
  system.time(GP <- fread(gpfile,stringsAsFactors=FALSE, sep="\t", header= F))
  GP <- as.data.frame(GP[-1,])
  colnames(GP)[1] <- "GP"
  zis <-grep('pop', GP$GP, ignore.case =TRUE)
  zis <-as.data.frame(zis)
  zis <- as.data.frame(zis[1,1])
  L = as.numeric((zis))
  #Loading order file
  loci= read.delim(loci.file, header = F)
  #preparing output name
  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_matched.gen")
  # Read in genepop file in two chunks: 1) list of loci, 2) individuals x loci
  system.time(GP1 <- fread(gpfile,stringsAsFactors=FALSE, nrows= L, sep=" ", fill=TRUE, header= F, quote= "\""))
  GP1 <-as.data.frame(GP1)
  GP1 <- (GP1[-1, ])
  GP1 <- as.data.frame(GP1)
  colnames(GP1)[1] <- "V1"

  keeps <- "V1"
  GP1 = GP1[keeps]
  print("Grabbing genotype information from genepop file")
  #OLS# GP2 <- read.table(gpfile, header = F, skip = L , sep=" ", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)

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
  check2 <- tail(names(sort(table(gptest$V3))), 1)
  if (check == ",") {
    print ("Genepop Format is Bad... Fixing")
    GP2$V1 <- paste (GP2$V1, GP2$V2, sep="")
    GP2$V2 =NULL
  }

  if (check2 == "") {
    print ("Genepop Format is Bad... Fixing")
    GP2$V3 = NULL
  }

  if (check == "0") {
    print ("Genepop Format is Bad... Fixing")
    GP2$V2 = NULL
  }

  #Determine loci positions in genepop file
  GP1$V2 <- 1:nrow(GP1)
  GP1$V2 <- 1 + GP1$V2
  GPR <- GP1 [! GP1$V1 %in% loci$V1, ]
  GPR <- subset(GP1, GP1$V1 %in% loci$V1)
  GPR$V2 <- sub("","V",GPR$V2)
  GPX <- as.data.frame(GPR)
  GPR <- as.character(GPR$V2)
  GPR <- c(keeps, GPR)

  if (nrow(GPX) != nrow(loci)) {
    print ("WARNING: Loci list differs from Loci in Genepop File")
  }

  #if (nrow(GPX) != nrow(loci))


  #determine desired loci positions
  loci$V2 <- rownames(loci)
  loci$V2 <- as.numeric(loci$V2)
  loci$V3 <- loci$V2 +1




  mergzX <- merge(loci,GPX,by="V1",all.x=TRUE)
  mergz <- merge(loci, GPX, by="V1")


  if ( (nrow (loci) == nrow(GP1)) & (nrow (mergzX[complete.cases(mergzX), ]) != nrow(mergzX)))   {
    print("Same number of loci, but they don't all match. Removing non-matching row")
    mergT <- merge(GP1, loci, by="V1", all.x=TRUE)
    mergT2 <- mergT[rowSums(is.na(mergT)) > 0,]
    mloc<- c(mergT2[,2])
    GP2[,mloc] <-NULL
    mloc2 = mloc -1
    remrow <- function(x, rows) x[-rows,, drop = FALSE]
    GP1<-remrow(GP1,mloc2)
  }


  if (nrow (mergzX[complete.cases(mergzX), ]) != nrow(mergzX)) {
    print("Loci list has loci not present in Genepop file, these will be populated with 0's")
    extraZ <- mergzX[rowSums(is.na(mergzX)) > 0,]
    extraG <- as.character(extraZ[[1]])
    print("Additional Genes are:")
    print(extraG)

    #sort
    mergz$V1 <-as.character(mergz$V1)
    GP1$V1 <-as.character(GP1$V1)
    mergz$V1[order(match(mergzX$V1,GP1$V1))]
    mergz <- mergz[match(GP1$V1,  mergz$V1),]
    mergz <-  mergz[complete.cases(mergz), ]

    t2<- as.character(loci$V1)
    mergzX$V1 <-as.character(mergzX$V1)
    mergzX <- mergzX[match(t2, mergzX$V1),]

    neworder <- mergzX$V2.x + 1
    neworder <- as.numeric(neworder)
    neworder <- c(1, neworder)
    ncharz<- as.numeric(nchar(GP2[2,2]))
    NoID= paste0(replicate(ncharz, "0"))
    NoID= paste(NoID, sep="", collapse="")

    prepORD <- rbind(mergz,extraZ)

    corder <- c(1, prepORD$V3)

    mergT <- merge(GP1, loci, by="V1", all.x=TRUE)
    mergT2 <- mergT[rowSums(is.na(mergT)) > 0,]
    mloc<- c(mergT2[,2])
    GP2[,mloc] <-NULL
    mloc2 = mloc -1
    remrow <- function(x, rows) x[-rows,, drop = FALSE]
    GP1<-remrow(GP1,mloc2)



    for (i in (ncol(GP2)+1)  : (ncol(GP2) + length(extraG))) {
      GP2[,i]= NoID
    }


    colnames(GP2) <- corder
    X <- as.numeric(colnames(GP2))
    sprintf("%08i",X)
    colnames(GP2) <- sprintf("%08i",X)

    GP2 <- GP2[ , order(names(GP2))]

    #GPT<-GP2[!(GP2[,1]=="Pop"),]


    toMatch<-c("\\bPop\\b","\\bPOP\\b","\\bpop\\b")
    matcheZ <- unique (grep(paste(toMatch,collapse="|"),
                            GP2[,1]))


    for (i in matcheZ) {
      GP2[i,2:ncol(GP2)] =""
    }

    header = "Matched by match.loci on"
    header = paste(header, Sys.time())

    print("Creating new file")
    cat(header, '\n',  file = "Xtemp")
    write.table(loci$V1, file="Xtemp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(GP2, file="Xtemp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    file.rename("Xtemp", gpout)
    print ("Done! Look for _sort_matched file")




  }



  if (nrow (mergzX[complete.cases(mergzX), ]) == nrow(mergzX)) {
    #sort it
    mergz <- mergz[order(mergz$V3),]
    neworder <- sub("V","",mergz$V2.y)
    neworder <- as.numeric(neworder)
    neworder <- c(1, neworder)

    #reorder
    rearg <- GP2[,neworder]
    if (sum (na.omit(as.numeric(rearg[,2]))) < 1) {
      print ("ERROR, something wrong with the file")
    }

    if (sum (na.omit(as.numeric(rearg[,2]))) == "NA") {
      print ("ERROR, something wrong with the file")
    }

    if (sum (na.omit(as.numeric(rearg[,2]))) > 1) {
      #Reconstruct GenePop file
      header = "Matched by match.loci on"
      header = paste(header, Sys.time())

      print("Creating new file")
      cat(header, '\n',  file = "Ztemp")
      write.table(loci$V1, file="Ztemp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
      write.table(rearg, file="Ztemp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
      file.rename("Ztemp", gpout)
      print ("Done! Look for _sort_matched file")
    }

    proc.time() - ptm
  }
}

