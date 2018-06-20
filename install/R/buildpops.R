#' buildpops
#'
#' Sort, filter, and rearrange genepop populations by individual names 
#' @param gpfile Name of genepop file. REQUIRED.
#' @param pop.file Name of table with population designation, 1 pop per row. REQUIRED.
#' @keywords buildpops
#' @export
#' @examples
#' buildpops (gpfile="Populations.gen", pop.file="5_pops.txt")

buildpops <- function (gpfile,
                       pop.file
) {

  suppressMessages(require(data.table))
  print( "Reading  genepop file to get population and locus information")
  GP <- fread(gpfile, header = FALSE, sep="\t", fill = TRUE, quote="\"")
  GP <- GP[-1,]
  GP <-as.data.frame(GP)
  colnames(GP) <- 'GP'

  if (ncol(GP) > 1){
    zis <-grep('pop', GP$V1, ignore.case =TRUE)
    zis <-as.data.frame(zis)
    zis <- as.data.frame(zis[1,1])
    L = as.numeric((zis))
  }

  if (ncol(GP) < 2){
    zis <-grep('pop', GP$GP, ignore.case =TRUE)
    zis <-as.data.frame(zis)
    zis <- as.data.frame(zis[1,1])
    L = as.numeric((zis))
  }

  gpops= read.csv(pop.file, header = F, sep='\t', row.name = 1, na.strings=c("","NA"))
  n <- as.numeric(nrow(gpops))
  pop.list <- split(gpops, seq(nrow(gpops)))
  pop.list <- lapply(pop.list, function(x) x[!is.na(x)])

  gpout= sub('\\..*', '', gpfile)
  gpout = paste0(gpout,"_pop_sort.gen")

  GP1 <- fread(gpfile, header = F, sep=" ", nrows = L  , fill = TRUE, quote= "\"")
  GP1 <-as.data.frame(GP1)
  GP1 <- (GP1[-1, ])
  GP1 <- as.data.frame(GP1)
  colnames(GP1)[1] <- "V1"
  keeps <- "V1"
  GP1 = GP1[keeps]
  print( "Reading and splitting genepop file")
  
  
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
  
  GP3 <-GP2
  GP2$ord <- 1:nrow(GP2)
  sl <- as.numeric(ncol(GP2) -1)
  temprow <- matrix(c(rep.int("",length(GP3))),nrow=1,ncol=length(GP3))
  newrow <- data.frame(temprow, stringsAsFactors=FALSE)
  colnames(newrow) <- colnames(GP3)
  newrow[1, 1] = "Pop"
  say = "Making"
  say2 = "Population(s)"
  say3 = paste(say,n,say2)
  print(say3)

  for (i in (pop.list)) {
    rm<-as.array(i)
    rm2 <-paste0("\\b",rm,"\\b")
    rm <-paste(rm2, collapse = "|")
    rml <- grep(rm, GP2$V1)
    rml <-as.data.frame(rml)
    rmt<- unique(rml)
    rmt<- as.data.frame(rmt)
    GPR <- subset(GP2, GP2$ord %in% rmt$rml)
    GPR$ord=NULL
    GPR <- GPR[!GPR$V1 == "pop", ]
    GPR <- GPR[!GPR$V1 == "Pop", ]
    GPR <- GPR[!GPR$V1 == "POP", ]
    GPR <- rbind(newrow,GPR)
    print("Compiling populations:")
    print (i)
    GPR <-as.matrix(GPR)
    write.table(GPR, file="stemp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  print( "Loading compiled list")
  GPR <- fread(file="stemp", header = F, sep=" ", stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
  GPR <- na.omit(GPR)
  GPR <- as.matrix(GPR)
  file.remove("stemp")

  GP1 <- as.matrix(GP1)
  #Reconstruct GenePop file
  header = "Sorted by pop.sort on"
  header = paste(header,Sys.time())
  cat(header, '\n',  file = "temp")
  write.table(GP1, file="temp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(GPR, file="temp", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  file.rename("temp", gpout)
  print( "Done! look for pop_sort file")

}
