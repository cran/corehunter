#' Core Hunter 3
#'
#' Core Hunter is a tool to sample diverse, representative subsets from large germplasm
#' collections, with minimum redundancy. Such so-called core collections have applications in plant
#' breeding and genetic resource management in general. Core Hunter can construct cores based on
#' genetic marker data, phenotypic traits or precomputed distance matrices, optimizing one of many
#' provided evaluation measures depending on the precise purpose of the core (e.g. maximum diversity,
#' representativeness, or allelic richness). In addition, multiple measures can be simultaneously
#' optimized as part of a weighted index to bring the different perspectives closer together.
#' The Core Hunter library is implemented in Java 8 as an open source project
#' (see \url{http://www.corehunter.org}).
#'
#' @examples
#' \dontrun{
#' # sample core based on genetic marker data
#' geno.file <- system.file("extdata", "genotypes.csv", package = "corehunter")
#' geno <- genotypes(file = geno.file)
#' sampleCore(geno)
#'
#' # sample core based on phenotypic traits
#' pheno.file <- system.file("extdata", "phenotypes.csv", package = "corehunter")
#' pheno <- phenotypes(file = pheno.file)
#' sampleCore(pheno)
#'
#' # sample core based on precomputed distance matrix
#' dist.file <- system.file("extdata", "distances.csv", package = "corehunter")
#' dist <- distances(file = dist.file)
#' sampleCore(dist)
#'
#' # sample core from genotypes with custom objective (allelic richness)
#' sampleCore(geno, obj = objective("HE"))
#'
#' # sample core from genotypes with custom size and objective (representativeness)
#' sampleCore(geno, obj = objective("AN", "MR"), size = 0.1)
#'
#' # sample core from genotypes with custom size and stop condition
#' sampleCore(geno, size = 0.1, impr.time = 2)
#'
#' # sample core based on both genotypes and phenotypes
#' geno.pheno <- coreHunterData(geno, pheno)
#' sampleCore(geno.pheno)
#' }
#'
#' @seealso \code{\link{coreHunterData}}, \code{\link{genotypes}},
#'  \code{\link{phenotypes}}, \code{\link{distances}},
#'  \code{\link{sampleCore}}, \code{\link{evaluateCore}},
#'  \code{\link{objective}}
#'
#' @docType package
#' @name corehunter
NULL
