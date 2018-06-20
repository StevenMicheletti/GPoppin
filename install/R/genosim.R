#' genosim
#'
#' Simulate Genepop genotypes using allele frequency tables and HW disequilibrium 
#' @param infile Input allele frequency table (from GPoppin::popgen). REQUIRED.
#' @param iter Number of individuals to simulate per population. Default = 100.
#' @param hwefile Input Hardy-Weinberg disequilbrium table (from GPoppin::popgen) Default = NULL.
#' @keywords genosim
#' @export
#' @examples
#' genosim(infile="Populations_freqtable.txt", iter = 200, hwefile = "Populations_loc_dis.txt")

genosim <- function (infile, iter = 100, hwefile = NULL) 
{
  if (is.null(hwefile)) {
    print("No disequilibrium file provided, will skip HWE correction")
  }
  freqz <- read.csv(infile, sep = "\t", as.is = TRUE, header = TRUE, 
                    row.names = 1, check.names = FALSE)
  freqz <- as.data.frame(t(freqz))
  
  check <- data.frame(min=sapply(freqz,min),max=sapply(freqz,max))
  check <- data.frame(check$max - check$min)
  check$V2 <- rownames(check)
  check <- check[check[[1]] == 0,]
  
  if (nrow(check) > 0) {
      print("Monomorphic loci detected and removed")
      freqz<-freqz[,-c(as.numeric(check$V2))]
  }
    
  coL <- as.numeric(ncol(freqz))
  coL2 <- coL/2
  roW <- as.numeric(nrow(freqz))
  headz <- seq(1, ncol(freqz), by = 2)
  headz2 <- freqz[, headz]
  headz2 <- colnames(headz2)
  headz2 = sub("\\..*", "", headz2)
  start <- seq(1, by = 2, length = coL/2)

  splitter <- lapply(start, function(i, freqz) freqz[i:(i + 1)], freqz = freqz)
  
  
  if (!is.null(hwefile)) {
    hwedis <- read.csv(hwefile, sep = "\t", as.is = TRUE, 
                       header = TRUE, row.names = 1)
    hwedis$mean = NULL
    
    if (nrow(check) > 0) {
      hwedis <- hwedis[-c(as.numeric(check$V2)),]
    }
    
    hwedis <- as.data.frame(t(hwedis))
    start2 <- seq(1, by = 1, length = coL2)
    hwesplit <- lapply(start2, function(i, hwedis) hwedis[i:i], 
                       hwedis = hwedis)
  }
  P = coL2
  L = iter * roW
  if (file.exists("sim_temp1")) {
    file.remove("sim_temp1")
  }
  if (file.exists("sim_temp2")) {
    file.remove("sim_temp2")
  }
  if (file.exists("sim_temp3")) {
    file.remove("sim_temp3")
  }
  outtable <- data.frame(matrix(ncol = P, nrow = L))
  colnames(outtable) <- c(1:P)
  print("Simulation is running. Status...")
  ptm <- proc.time()
  for (f in (1:coL2)) {
    spliz = as.matrix(splitter[[f]])
    for (i in (1:roW)) {
      for (j in (1:iter)) {
        k = i + 1
        p = spliz[[i, 1]]
        q = spliz[[i, 2]]
        AA = p^2
        Aa = 2 * p * q
        aa = q^2
        n = j
        if (!is.null(hwefile) & p > 0 & p < 1 & q > 0 & 
            q < 1) {
          hx = hwesplit[[f]][[i, 1]]
          if (is.na(hx)) {
            hx=0
            } else {
            ha = hx/2
            }
          Aa = Aa + hx
          if (Aa < 0) {
            Aa <- 0
          } else {
            Aa=Aa
          }
          if (Aa > 1) {
            Aa <- 1
          } else {
            Aa=Aa
          }
          AA = AA - ha
          if (AA < 0) {
            AA <- 0
          }
          if (AA > 1) {
            AA <- 1
          }
          aa = aa - ha
          if (aa > 1) {
            aa <- 1
          }
          if (aa < 0) {
            aa <- 0
          }
        }
        rname <- rownames(splitter[[f]])
        popname <- rname[i]
        cname <- colnames(splitter[[f]])
        locname <- cname[1]
        al1 <- as.character(sub(".*\\.", "", colnames(spliz)[1]))
        al2 <- as.character(sub(".*\\.", "", colnames(spliz)[2]))
        g1 <- paste0(al1, al1)
        g2 <- paste0(al1, al2)
        g3 <- paste0(al2, al2)
        geno <- sample(c(g1, g2, g3), size = 1, replace = TRUE, 
                       prob = c(AA, Aa, aa))
        idtable <- paste(popname, j, ", ", sep = "_")
        idtable <- noquote(idtable)
        write.table(geno, file = "sim_temp1", append = TRUE, 
                    sep = " ", row.names = FALSE, col.names = FALSE)
        write.table(idtable, file = "sim_temp2", append = TRUE, 
                    sep = " ", row.names = FALSE, col.names = FALSE)
      }
    }
    tabwrite <- read.table("sim_temp1", as.is = TRUE, header = FALSE, 
                           colClasses = "character")
    outtable[[f]] <- tabwrite$V1
    file.remove("sim_temp1")
    e = (f/coL2) * 100
    print(e)
  }
  rownamez <- read.table("sim_temp2", as.is = TRUE, header = FALSE, 
                         nrows = L, colClasses = "character")
  rn2 <- gsub("_, ", ",", rownamez$V1)
  rownames(outtable) <- rn2
  file.remove("sim_temp2")
  print("Results are in, making into genepop format")
  temprow <- matrix(c(rep.int("", coL2)), nrow = 1, ncol = coL2)
  newrow <- data.frame(temprow, stringsAsFactors = FALSE)
  newrow <- newrow[rep(seq_len(roW), each = 1), ]
  seqn = iter + 1
  seqx = (seqn * roW)
  sec1 <- seq(from = 1, to = seqx, by = seqn)
  popprint <- rep("Pop", roW)
  newrow$name = popprint
  newrow$seq = seq(from = 1, to = seqx, by = seqn)
  newrow[is.na(newrow)] <- " "
  seql = nrow(outtable) + roW
  seqt <- seq(from = 1, to = seql, by = 1)
  sec2 <- seqt[-(sec1)]
  outtable$name <- rownames(outtable)
  outtable$seq <- sec2
  colnames(newrow) <- colnames(outtable)
  comtest <- rbind(outtable, newrow)
  comtest <- comtest[order(as.numeric(as.character(comtest$seq))), 
                     ]
  cat <- as.data.frame(comtest$name)
  comtestcat = cbind(cat, comtest)
  comtestcat$name = NULL
  comtestcat$seq = NULL
  loci <- as.data.frame(headz2)
  header = "Simulated individuals by gpsim: "
  header = paste0(header, iter, " individuals made on ", Sys.time())
  print("Creating new file")
  cat(header, "\n", file = "sim_temp3")
  write.table(loci, file = "sim_temp3", append = TRUE, sep = " ", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(comtestcat, file = "sim_temp3", append = TRUE, 
              sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  gpout = sub("\\..*", "", infile)
  gpout = c(gpout, "_", iter, "_sim.gen")
  gpout = paste0(gpout, collapse = "")
  file.rename("sim_temp3", gpout)
  print("Simulation complete, look for _sim.gen file")
  proc.time() - ptm
}

