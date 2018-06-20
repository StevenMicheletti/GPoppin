#' prog2gp
#'
#' Convert a progeny export csv, 1 allele per column format, to Genepop format.
#' @param infile Input progeny export table. csv or txt. REQUIRED.
#' @param rem.indels Remove indels from dataset. Default = FALSE.
#' @param rem.sex Remove sex markers from dataset. Default = FALSE.
#' @param fix.loc Fix problematic locus names that contain invalid characters for R packages. Default = TRUE.
#' @param rem.thres Remove loci that are not genotyped in N percent of individuals. Default = 0.
#' @param split.pop Attempt to split Genepop file into pops based on individual names. Default = FALSE.
#' @param rem.dup Remove duplicated individuals by name. Default = TRUE.
#' @param first.col Genotypes start at the second column, after individual names. Default = FALSE.
#' @param examp Use example file. Default = FALSE.
#' @keywords prog2gp
#' @export
#' @examples
#' prog2gp(infile = "progeny_export.csv", rem.indels=FALSE, rem.sex=FALSE,
#' fix.loc=TRUE, rem.thres=0, split.pop=TRUE, rem.dup=TRUE, first.col = FALSE, examp=FALSE)


prog2gp <- function (infile,
                     rem.indels = FALSE,
                     rem.sex = FALSE,
                     fix.loc = TRUE,
                     rem.thres = 0,
                     split.pop = FALSE,
                     rem.dup = TRUE,
                     first.col = FALSE,
                     examp = FALSE

) {


  suppressMessages(require('data.table'))
  #Begin and read file
  print("Loading file and determining locus start column...")
  ptm <- proc.time()


  if (examp == FALSE) {
  gpops <- fread(infile,stringsAsFactors=FALSE, na.strings =c(""," ","NA"), colClasses = c("character"), header = T)
  gpops <- as.data.frame(gpops)
  }

  else {
  gpops <-as.data.frame(ex_prog)
  }
  
   #remove blank rows
  gpops <-gpops[rowSums(is.na(gpops)) != ncol(gpops),]

  #Now blank columns
  gpops <- Filter(function(x)!all(is.na(x)), gpops)

  #make all columns characters
  gpops <- data.frame(lapply(gpops, as.character), stringsAsFactors=FALSE)

  #make NAs 0
  gpops[is.na(gpops)] <- 0

  #Fix weird genotype calls to 0
  gpops <- as.data.frame(sapply(gpops, function(x) gsub("Invalid", "0", x)))
  gpops <- as.data.frame(sapply(gpops, function(x) gsub("No Call", "0", x)))
  gpops <- data.frame(lapply(gpops, as.character), stringsAsFactors=FALSE)

  
  if (first.col == TRUE) {
    print("Auto-detection off, genotypes should start on first column")
    GC=2
  } else { 
  
    #Auto Determine starting genotype column by a set of criteria
    # Has to have one character per row
    autod <- as.data.frame(sapply(gpops, function(x) sum(nchar(x))))
    autod$RN <- 1:ncol(gpops)
    colnames(autod) <- c("V1", "RN")
    autod$V2 <- (sapply(gpops, function(x) tail(names(sort(table(x))), 1)))
    colnames(autod) <- c("V1","RN","V2")

    #Feature to allow for combined allele columns in the future
    AN=1
    autopos <- subset(autod, V1 == nrow(gpops))
      #if (AN == 1) {
      #  autopos <- subset(autod, V1 == nrow(gpops))
      #}

     #if (AN == 2) {
      #  autopos <- subset(autod, V1 == 2*nrow(gpops))
      #not tested yet
      #}

    # Has to consist of an allele character
    autopos$V2 <- as.character(autopos$V2)
    autopos <-autopos[  autopos$V2 %in% c('A','T','C','G','X','Y','-','0'), ]
   # Has to be next to other genotype columns
    autopos$V3 <- as.numeric(1:nrow(autopos))
    autopos <-data.table(autopos)
    setkey(autopos,V1)
    autopos[,V4:=c(NA,diff(RN)), by=V1]
    autopos[is.na(autopos)] <- 1

    bc <- autopos[which.max(autopos$V4),]
    bc <- as.numeric(bc[1,4])
    shart <- autopos[bc:nrow(autopos),]

   GC <- as.numeric(shart[1,2])


   LCHECK=(ncol(gpops[, GC:ncol(gpops)]))/2
    if ((LCHECK%%1==0) == FALSE) {
        GC = GC - 1
    }
  }
    START <- colnames(gpops)[GC]
  
  
  alert2<- paste0("Starting locus is ", START, " at column ", GC, " of gpops")
  print(alert2)

  #Check for individual name column
  #Ensure individual name case is correct
  if ((!"Individual.name" %in% colnames(gpops)) & (!"Individual.Name" %in% colnames(gpops)) & (!"Individual.name" %in% colnames(gpops)) & (!"Individual name" %in% colnames(gpops)) & (!"Individual Name" %in% colnames(gpops)) & (!"individual.name" %in% colnames(gpops)) & (!"individual name" %in% colnames(gpops))) stop("Individual Name column not found, make sure this heading is present!")
  
  if ("Individual.name" %in% colnames(gpops)) {
    indZ <- data.frame(gpops$`Individual.name`)
  }
  
  if ("Individual name" %in% colnames(gpops)) {
    indZ <- data.frame(gpops$`Individual name`)
  }
  
  if ("Individual.Name" %in% colnames(gpops)) {
    indZ <- data.frame(gpops$`Individual.Name`)
  }
  
  if ("Individual Name" %in% colnames(gpops)) {
    indZ <- data.frame(gpops$`Individual Name`)
  }
  
  if ("individual name" %in% colnames(gpops)) {
    indZ <- data.frame(gpops$`individual name`)
  }
  
  if ("individual.name" %in% colnames(gpops)) {
    indZ <- data.frame(gpops$`individual.name`)
  }
  
  if (nrow(indZ) == 0 ) stop("Problem identifying individual column. Check headings")
  
  colnames(indZ) <- "IND"
  genoZ <- as.matrix(gpops[, GC:ncol(gpops)])


  if (nrow(genoZ) != nrow(indZ))  stop("Individual rows and number of genotypes do not match. Check input file")

  #Replace these progeny characters with these numbers
  from = c("A", "C", "G", "T", "-","0","X","Y")
  to = c("102", "104", "106", "108", "100", "000", "006", "009")

  map = setNames(to, from)
  genoZ[] = map[genoZ]

  #Prepare data
  genoZ <- as.data.frame(genoZ)
  sloci <- ncol(genoZ)/2

  alert1<-paste0('Starting loci number = ', sloci)
  print(alert1)

  
  #split genotype file into separate allele 1 and 2 data frames
  A1 <-as.matrix(genoZ[c(T,F)])
  A2 <-as.matrix(genoZ[c(F,T)])

  #Merge Allele 1 with Allele 2 to get 4 character Genepop format
  combmatrix <- data.frame(matrix(, nrow = nrow(genoZ), ncol = ncol(A1)))
  for (i in (1:ncol(A1))){
    combmatrix[,i] <- as.data.frame(paste0(A1[,i],A2[,i]))
  }
  locusNA <- colnames(A1)
  c.names <- gsub(".A1", "", locusNA)
  colnames(combmatrix) <-c.names
  combmatrix2<-combmatrix

  #If error pops up, make sure the correct start column is being called


  #Optional removal of sex markers
  if (rem.sex == TRUE) {
    cols <- colSums(mapply('==', '006009', combmatrix2))
    combmatrix2 <- combmatrix2[,which(cols == 0)]
    cols <- colSums(mapply('==', '006006', combmatrix2))
    combmatrix2 <- combmatrix2[,which(cols == 0)]
    dloci <-ncol(combmatrix2)
    dloci <- sloci-dloci
    sx<- paste( dloci, 'sex markers removed')
    print(sx)
  }

  #Optional removal of insertion markers
  if (rem.indels == TRUE) {
    cols <- colSums(mapply('==', '100100', combmatrix2))
    combmatrix2 <- combmatrix2[,which(cols == 0)]
    eloci <-ncol(combmatrix2)
    eloci <- sloci-eloci-1
    ins <-paste( eloci, 'insertion markers removed')
    print(ins)
  }

  #Number of final loci retained
  st3 <- paste(ncol(combmatrix2), "loci retained after any filters")
  print(st3)

  #Get population names (determine splitting character)
  locilist <- as.data.frame(colnames(combmatrix2))

  com.character <- gsub('[0-9]+', '', as.character(indZ$IND))
  com.character <- gsub('[a-zA-Z]+', '', com.character)
  com.character <- gsub('[,]+', '', com.character)
  fun1 <- function(InVec) {
    names(which.max(table(InVec)))
  }
  fun2 <- function(InVec) {
    m <- which.max(table(as.character(InVec)))
    as.character(InVec)[m]
  }
  com.character <- fun2(as.vector(com.character))
  com.character <- as.name(com.character)
  com.character <- paste0(com.character, ".*$")
  colnames(locilist)  <-"name"
  ##If com.character fails, put in custom
  pops <- gsub(com.character, "", indZ$IND)
  #pops <- gsub("-.*$", "", indZ$IND)

  #Add comma for Genepop format

  indZ$IND <- paste0(indZ$IND,',')
  indZ$POP <- pops

  #Add populations and get names of unique populations
  genoZ <- cbind(indZ, combmatrix2)
  pop.list<- unique(indZ$POP)

  st4<- paste(nrow(genoZ), "individuals initially")
  print(st4)

  #remove duplicate individuals
  badg <- as.data.frame(genoZ$IND)
  dupz <- as.data.frame(as.numeric(duplicated(genoZ$IND) +1))
  badg$dup <-dupz$`as.numeric(duplicated(genoZ$IND) + 1)`
  badg2<-badg[!(badg$dup < 2),]
  ind <- as.character(badg2$`genoZ$IND`)
  st6<- paste(nrow(badg2), "individual IDs are duplicates and have been resolved")
  print(st6)
  alert3 <- paste0("Duplicated individual: ", ind)
  print(alert3)

  if (rem.dup == TRUE) {
    dupz<-dupz[!(dupz$`as.numeric(duplicated(genoZ$IND))` > 0 ),]
    badg3<-badg[!(badg$dup > 1),]
    retain <- rownames(badg3)
    genoZ  <- as.data.frame(genoZ[c(retain), ])
    rownames(genoZ) <-c(1:nrow(genoZ))
  }

  #Remove low genotype individuals
  badg <- as.matrix(genoZ[, c(3:ncol(genoZ))])
  class(badg) <-"numeric"
  # If there is a warning here, check your genotypes for odd characters (i.e, 0 = 0:00)
  badg[badg > 0] <- 1
  sumz <- as.data.frame(rowSums(badg))
  sumz$success <- sumz$`rowSums(badg)`/ncol(badg)
  sumz<-sumz[!(sumz$success < rem.thres),]
  prop = 100-(rem.thres * 100)
  if (rem.thres > 0) {
    st5<- paste0(nrow(genoZ) - nrow(sumz), " individuals removed due to low genotyping (Have " , prop, "% missing data or more)" )
    print(st5)
  }
  retain <- rownames(sumz)
  genoZ  <- as.data.frame(genoZ[c(retain), ])

  #Create 'Pop' heading for each population
  temprow <- matrix(c(rep.int("",length(genoZ))),nrow=1,ncol=length(genoZ))
  newrow <- data.frame(temprow, stringsAsFactors=FALSE)
  colnames(newrow) <- colnames(genoZ)
  newrow[1, 1] = "Pop"

  if (split.pop == TRUE) {
    #Split pop genotype data frames to process separately
    splitpops <- split(genoZ, genoZ$POP)

    #Add pop heading to each list
    VaR <- list()
    for (i in (1:length(splitpops))) {
      test<- as.data.frame(splitpops[i])
      colnames(newrow) <- colnames(test)
      VaR[[i]]<- as.data.frame(rbind(newrow,test))
      colnames(VaR[[i]]) <- colnames(genoZ)
    }

    #Reassemble list into single data.frame
    reassemble <- do.call("rbind", VaR)
    reassemble$POP <- NULL
    st7 <- paste(length(VaR), "populations created")
  }

  if (split.pop == FALSE) {
    reassemble <- as.data.frame(rbind(newrow,genoZ))
    colnames(reassemble) <- colnames(genoZ)
    reassemble$POP <- NULL
    st7 <- paste("1 population created")
  }

  #Remove problematic characters from locus names
  if (fix.loc == TRUE) {
    locilist<-as.data.frame(gsub("Pop|POP|pop|PoP|poP", "pep", locilist$name))
    colnames(locilist) <-"name"
    locilist<-as.data.frame(gsub("`|'|[.]|:", "$", locilist$name))
    colnames(locilist) <-"name"
    print("Problematic characters in locus names were replaced with $")
    print("Any occurrence of pop in locus names were replaced with pep")
  }

  #Reconstruct GenePop file

  if (examp == FALSE){
   gpout= sub('\\..*', '', infile)
   gpout = paste0(gpout,".gen")
    gpout = paste(gpout, collapse=" ")
  }
  else{
    gpout="ex_prog_out.gen"
  }


    if(file.exists(gpout)) {
     print("WARNING. A previous version of the output file existed and was overwritten")
   }

  
  reassemble <- reassemble[complete.cases(reassemble), ]

  header = "A=102,C=104,G=106,T=108,-=100,X=006,Y=009; Created by prog2gp on"
  header = paste(header,Sys.time())
  cat(header, '\n',  file = "tempZ")
  write.table(locilist, file="tempZ", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(reassemble, file="tempZ", append = TRUE, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  alert7 <- paste0("Writing file: ", gpout)
  print(alert7)
  file.rename("tempZ", gpout)



  print("Finished without errors")

  proc.time() - ptm
  #rm(list=setdiff(ls(), "ind"))
}
