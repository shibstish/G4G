#' @param pheno file with numeric phenotypic values
#' @param geno data.frame with genotype calls coded as 0,1,2.
#' @param covar numeric data.frame with covariates values
#' @param map genetic map of data with chr and position of each SNP
#' @param PCA.m number of principal components to use default is 3

#Step 1. Start of G4G function
G4G <- function(pheno = NULL,geno = NULL, covar, PCA.m = 3, map = NULL){
  start_time <- Sys.time() #start recording the amount of time it takes to run this function
  m <- ncol(geno)
  n <- nrow(geno)
  y <- pheno
  P <- array(dim = m)

  #1.1. PC analysis
  PCA <- prcomp(geno)
  PCs <- PCA$x
  PCs <- as.matrix(PCs[,1:3])
  for(i in 1:m){
    x = geno[,i]
  if(max(x)==min(x)){p=1
 Pvalues <- p
  }else{
    #1.2. Check PC dependencies

    # if covariates are not provided by the user:
    if(missing(covar)==TRUE){
    X <- cbind(1, x, PCA$x[,1:PCA.m])}else{

    #if covariates are provided:
    #check linear dependence between the PCAs and the given covariates. Rank doesn't change if you drop linearly dependent data. Don't include the PCs that are dependent
       COVPC <- as.matrix(cbind(covar, PCs))
      rankifremoved <- sapply(3:ncol(COVPC), function (x) qr(COVPC[,-x])$rank) #3 refers to the 3rd column which is the 1st PC
      depend =which(rankifremoved == max(rankifremoved)) #tells which data to remove
      PC_nodep = PCs[,-depend]

      X <- as.matrix(cbind(1,x, covar, PC_nodep))
  }}
    #1.3. LINEAR REGRESSION MATRIX
    LHS <- t(X)%*%X
    C <- solve(LHS, tol = 1e-19)
    RHS <- t(X)%*%y
    b <- C%*%RHS
    yb <- X%*%b
    e <- y-yb
    n <- length(y)
    ve <- sum(e^2)/(n-1)
    vt <- C*ve
    t <- b/sqrt(diag(vt)) #t-statistic
    Pvalues <- 2*(1-pt(abs(t),n-2))#here is multiplied by 2 because of symmetry of the Gaussian dist.
    P[i] <-  Pvalues[2]

    #done with p-values
  }
  #1.4. Add manhattan plot
  color.vector <- rep(c("cyan2","darkgoldenrod1","chartreuse3","deeppink"),10)
  maplength=nrow(map)
  plot(-log10(P)~seq(1:maplength),col=color.vector[map[,2]], main = "manhattan plot")

  #1.5. Add QQ plot
  length2=length(P)
  p.uni=runif(length2,0,1)
  order.obs=order(P)
  order.uni=order(p.uni)

  plot(-log10(p.uni[order.uni]),-log10(P[order.obs]), main = "QQ plot", col = "deepskyblue4")
  abline(a = 0, b = 1, col = "red")

  #1.6. Processing time and returning p-values
  end_time <- Sys.time()
  finalTime <- end_time - start_time

  P.time <- list(P, finalTime)
  return(P.time)
  }






