
readVCF <- function(infilename, direct="./", gform="reference", sampthres = 0.01){
  outname <- VCFtoRA(infilename,direct=direct)
  return(readRA(outname[[1]], gform=gform, sampthres=sampthres))
}

VCFtoRA <- function(infilename, direct="./", makePed=F){
  
  ## Do some checks
  if(!is.character(infilename) || length(infilename) !=1)
    stop("The input file name is not a string of length 1.")
  if(!file.exists(infilename))
    stop("Input file does not exist. Check your wording or the file path.")
  if(!is.character(direct) || length(direct) != 1)
    stop("Invalid input for the path to the directory where the RA file is to be written.")
  
  cat("Processing VCF file: Converting to RA format.\n\n")
  
  outfilename <- paste0(tail(strsplit(infilename,split=.Platform$file.sep)[[1]],1),".ra.tab")
  outpath <- dts(normalizePath(direct, winslash=.Platform$file.sep, mustWork=T))
  
  headerlist = c('CHROM', 'POS')
  
  ## Read in the lines of the file
  Lines <- readLines(infilename)
  
  ## entries for empty genotypes
  empty_genotypes <- c("./.",".,.",".",".|.")
  
  start = which(unlist(lapply(Lines, function(x) substr(x,1,2) == "##")))  ## starting position
  start = max(start,0) + 1
  
  ## create the file
  outfile = file.path(outpath,outfilename)
  file.create(outfile,showWarnings = F)
  
  ## write the heading line
  line = strsplit(Lines[start], split="\t")[[1]]
  headerlist = c(headerlist, line[10:length(line)])
  ## Write the header to the RA file
  newLines <- list(paste(headerlist, collapse="\t"))
  cat("Found",length(headerlist),"samples\n")
  
  ## Now write the SNPs
  for(i in (start+1):length(Lines)){
    line = trimws(Lines[[i]])
    
    line = strsplit(line, split="\t")[[1]]
    chrom = line[1]
    pos = line[2]
    ref = line[4]
    alt = line[5]
    ## Filter out indels and multiple alternative alleles
    if (ref == "." || alt == "." || length(alt) > 1)
      next
    else{
      ## Check that there is allelic depth in the VCF file
      format = strsplit(line[9], split= ":")[[1]]
      if("AD" %in% format)
        ad_pos = which(format == "AD")
      else if(all(c("RO","AO") %in% format)){
        ro_pos = which(format == "RO")
        ao_pos = which(format == "AO")
      }
      else if("DP4" %in% format)
        dp4_pos = which(format == "DP4")
      else
        stop("We can't use this vcf file. AD (alleleic depth or RO (Reference allele observation count) and AO (Alternate allele observation count) information is needed.\n")
      ## Extract the alleles depths
      newline = c()
      #for(j in line[10:length(line)]){
      for(j in line[10:length(line)]){
        if (j %in% empty_genotypes)
          newline = c(newline,"0,0")
        else{
          j = strsplit(j, split = ":")[[1]]
          if( length(ad_pos) > 0 ){              #If ad_pos is not null, i.e., there is a value for it, append it to outlist
            if( j[ad_pos] %in% empty_genotypes ) #gatk vcf4.2 will fill out genotype fields, even for uncovered data
              newline = c(newline,"0,0")
            else
              newline = c(newline, j[ad_pos])
          }
          else if(length(ro_pos) > 0 && length(ao_pos) > 0){ #ELSE IF ro_pos and ao_pos are not null or are equal to 0
            ad = paste0(j[ro_pos], ",", j[ao_pos])
            newline = c(newline, ad)
          }
          else if( dp4_pos ){
            counts = strsplit(j[dp4_pos], split=",")
            allele1 = as.integer(counts[1]) + as.integer(counts[2])
            allele2 = as.integer(counts[3]) + as.integer(counts[4])
            ad = paste0(allele1, ",", allele2)
            newline = c(newline, ad)
          }
          else ##Should never really get here, but if AD, AO and RO are all null, it will break the script
            stop("Can't Find either Allele Depth (AD) or RO (Reference allele observation count) and AO (Alternate allele observation count) at this position.\n")
        }
      }
      newLines[[i-(start)+1]] <- paste0(c(chrom,pos,newline), collapse = "\t")
    }
  }
  ## open the connection to the file
  con <- file(outfile)
  writeLines(unlist(newLines), con=con)
  close(con)
  ## output the information
  cat(length(Lines)-start,"SNPs written\n\n")
  cat("Name of RA file:    ",outfilename,"\n")
  cat("Location of RA file: ",outpath,"/\n\n",sep="")
  ## Initialize the pedigree file
  if(makePed){
    pedfile <- paste0(strsplit(paste0(tail(strsplit(infilename,split=.Platform$file.sep)[[1]],1)),split="\\.")[[1]][1],"_ped.csv")
    pedpath <- file.path(outpath,pedfile)
    if(!file.exists(pedpath)){
      cat("A pedigree file has been initialized.\n")
      cat("Name of pedigree file:     ",pedfile,"\n")
      cat("Location of pedigree file: ",outpath,"/\n\n",sep="")
      nSamp <- length(headerlist) - 2
      write.csv(cbind("SampleID"=c(headerlist[-c(1,2)]),"IndividualID"=rep("",nSamp),"Mother"=rep("",nSamp),
                      "Father"=rep("",nSamp), "Family"=rep("",nSamp)), file = pedpath, row.names=F)
    }
  }  
  return(invisible(outfile))
}

#### Some functions from the kutils package for removing trailing spaces for filenames.
dts <- function (name) 
  gsub("/$", "", dms(name))

dms <- function(name)
  gsub("(/)\\1+", "/", name)

#### Function for reading in RA data and converting to genon and depth matrices.
readRA <- function(genofile, gform, sampthres = 0.01, excsamp = NULL){
  
  if(!is.character(genofile) || length(genofile) != 1)
    stop("File name of RA data set is not a string of length one")
  if(!is.character(gform) || length(gform) != 1 || !(gform %in% c("reference","uneak")))
    stop("gform argument must be either 'reference' ot 'uneak'")
  
  ## separate character between reference and alternate allele count
  gsep <- switch(gform, uneak = "|", reference = ",")
  ## Process the individuals info
  ghead <- scan(genofile, what = "", nlines = 1, sep = "\t")
  
  ## Read in the data
  # If reference based
  if (gform == "reference"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = "", coord = 0), rep(list(""), length(ghead) - 2)))
    chrom <- genosin[[1]]
    pos <- genosin[[2]]
    SNP_Names <- paste(genosin[[1]],genosin[[2]],sep="_")
    indID <- ghead[3:length(ghead)]
    AFrq <- NULL
  }
  else if (gform == "uneak"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = ""), rep(list(""), length(ghead) - 6), list(hetc1 = 0, hetc2 = 0, acount1 = 0, acount2 = 0, p = 0)))
    SNP_Names <- genosin[[1]]
    indID <- ghead[2:(length(ghead)-5)]
    AFrq <- genosin[[length(genosin)]]
    chrom <- pos <- NULL
  }
  
  ## compute dimensions
  nSnps <- length(SNP_Names)
  nInd <- length(ghead) - switch(gform, reference=2, uneak=6)
  
  ## generate the genon and depth matrices
  depth_Ref <- depth_Alt <- matrix(0, nrow = nInd, ncol = nSnps)
  start.ind <- switch(gform, uneak=1, reference=2)
  for (i in 1:nInd){ 
    depths <- strsplit(genosin[[start.ind+i]], split = gsep, fixed = TRUE)
    depth_Ref[i, ] <- as.numeric(unlist(lapply(depths,function(z) z[1])))
    depth_Alt[i, ] <- as.numeric(unlist(lapply(depths,function(z) z[2])))
  }
  genon <- (depth_Ref > 0) + (depth_Alt == 0)
  genon[which(depth_Ref == 0 & depth_Alt == 0)] <- NA
  
  ## Check that the samples meet the minimum sample treshold
  sampDepth <- rowMeans(depth_Ref + depth_Alt)
  badSamp <- which(sampDepth < sampthres)
  if(length(badSamp) > 0){
    cat("Samples removed due to having a minimum sample threshold below ",sampthres,":\n",sep="")
    cat(paste0(indID[badSamp],collapse = "\n"),"\n\n")
    excsamp <- unique(c(excsamp,indID[badSamp]))
  }
  ## Remove any sample which we don't want
  if(!is.null(excsamp)){
    toRemove <- which(indID %in% excsamp)
    if(length(excsamp) > 0){
      depth_Ref <- depth_Ref[-toRemove,]
      depth_Alt <- depth_Alt[-toRemove,]
      genon <- genon[-toRemove,]
      indID <- indID[-toRemove]
      nInd <- length(indID)
    }
  }
  
  ## Create the R6 object
  obj <- list(genon = genon, depth_Ref = depth_Ref, depth_Alt = depth_Alt, chrom = chrom, pos = pos,
         SNP_Names = SNP_Names, indID = indID, nSnps = nSnps, nInd = nInd, gform = gform, AFrq = AFrq)
  
  return(obj)
}




