
vcf <- deerData()
rafile <- VCFtoRA(vcf)
radata <- readRA(rafile,gform = "reference")

deer <- makeUR(radata, ploid = 2, filter = list(MAF=0.01, MISS=0.9))

## test out GUSLD

LDmat <- GUSLD(deer,nClust = 3, file="test",sf=6)
LDmat <- GUSLD(deer,SNPpairs = matrix(c(1:5,2:6), ncol=2),
               nClust = 3, file="test_sub", sf=6, LDmeasure = c("r2","LD_r"))
