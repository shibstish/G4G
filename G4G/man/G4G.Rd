\name{G4G}
\alias{G4G}

\title{
Calculate GWAS using the GLM model
}
\description{
G4G calculates p-values for SNPs associated with phenotypes by including cofactors and principle components. The result is
}
\usage{
G4G(pheno, geno, covar, PCA.m, map)
}
\arguments{
  \item{pheno}{
  numeric phenotype data
  *REQUIRED*
  }
  \item{geno}{
  data.frame with genotype calls coded as 0,1,2.
  *REQUIRED
  }
  \item{covar}{
  numeric data.frame with covariates values
  }
  \item{PCA.m}{
  default is 3
  }
  \item{map}{
  genetic map of data with chr and position of each SNP. used to generate the manhattan plot
  *REQUIRED*
  }
}
\value{
  \item{comp1 }{[[1]] list of SNP p-values
  [[2]] the amount of time the function took}
  \item{comp2 }{manhattan plot}
  \item{comp3 }{QQ plot}

}
\author{
Cristina Ocana Gallegos and Shelby Tisinai
}

