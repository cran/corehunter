# ------------ #
# EXAMPLE DATA #
# ------------ #

#' Small example dataset with 218 individuals.
#'
#' Data was genotyped using 190 SNP markers and 4 quantitative traits were recorded.
#' Includes a precomputed distance matrix read from \code{"extdata/distances.csv"},
#' genotypes read from \code{"extdata/genotypes-biparental.csv"} and phenotypes read
#' from \code{"extdata/phenotypes.csv"}.
#' The distance matrix is computed from the genotypes (Modified Rogers' distance).
#'
#' Data was taken from the CIMMYT Research Data Repository (Study Global ID
#' hdl:11529/10199; real data set 5, cycle 0).
#'
#' @source Cerón-Rojas, J. Jesús ; Crossa, José; Arief, Vivi N.; Kaye Basford;
#'         Rutkoski, Jessica; Jarquín, Diego ; Alvarado, Gregorio; Beyene, Yoseph;
#'         Semagn, Kassa ; DeLacy, Ian, 2015-06-04, "Application of a Genomics
#'         Selection Index to Real and Simulated Data",
#'         \url{http://hdl.handle.net/11529/10199} V10
#'
#' @examples
#' exampleData()
#'
#' @return Core Hunter data of class \code{chdata}
#' @export
exampleData <- function(){
  getFile <- function(file){
    system.file("extdata", file, package = "corehunter")
  }
  coreHunterData(
    distances = distances(file = getFile("distances.csv")),
    genotypes = genotypes(file = getFile("genotypes-biparental.csv"), format = "biparental"),
    phenotypes = phenotypes(file = getFile("phenotypes.csv"))
  )
}

# ---------------- #
# CORE HUNTER DATA #
# ---------------- #

#' Initialize Core Hunter data.
#'
#' The data may contain genotypes, phenotypes and/or a precomputed distance matrix.
#' All provided data should describe the same individuals which is verified by comparing
#' the item ids and names.
#'
#' @param genotypes Genetic marker data (\code{chgeno}).
#' @param phenotypes Phenotypic trait data (\code{chpheno}).
#' @param distances Precomputed distance matrix (\code{chdist}).
#'
#' @return Core Hunter data (\code{chdata}) with elements
#' \describe{
#'  \item{\code{geno}}{Genotype data of class \code{chgeno} if included.}
#'  \item{\code{pheno}}{Phenotype data of class \code{chpheno} if included.}
#'  \item{\code{dist}}{Distance data of class \code{chdist} if included.}
#'  \item{\code{size}}{Number of individuals in the dataset.}
#'  \item{\code{ids}}{Unique item identifiers.}
#'  \item{\code{names}}{Item names. Names of individuals to which no explicit name
#'    has been assigned are equal to the unique \code{ids}.}
#'  \item{\code{java}}{Java version of the data object.}
#' }
#'
#' @return Core Hunter data of class \code{chdata}.
#'
#' @examples
#' geno.file <- system.file("extdata", "genotypes.csv", package = "corehunter")
#' pheno.file <- system.file("extdata", "phenotypes.csv", package = "corehunter")
#' dist.file <- system.file("extdata", "distances.csv", package = "corehunter")
#'
#' my.data <- coreHunterData(
#'   genotypes(file = geno.file, format = "default"),
#'   phenotypes(file = pheno.file),
#'   distances(file = dist.file)
#' )
#'
#' @seealso \code{\link{genotypes}}, \code{\link{phenotypes}}, \code{\link{distances}}
#'
#' @import rJava
#' @importFrom methods is
#' @export
coreHunterData <- function(genotypes, phenotypes, distances){

  # check arguments
  if(!missing(genotypes) && !is(genotypes, "chgeno")){
    stop("Argument 'genotypes' should contain Core Hunter genotype data of class 'chgeno'.")
  }
  if(!missing(phenotypes) && !is(phenotypes, "chpheno")){
    stop("Argument 'phenotypes' should contain Core Hunter phenotype data of class 'chpheno'.")
  }
  if(!missing(distances) && !is(distances, "chdist")){
    stop("Argument 'distances' should contain Core Hunter distance matrix data of class 'chdist'.")
  }
  if(missing(genotypes) && missing(phenotypes) && missing(distances)){
    stop("Please specify at least one type of data (genotypes, phenotypes and/or distances).")
  }

  # create data
  j.geno <- .jnull(ch.genotypes())
  j.pheno <- .jnull(ch.phenotypes())
  j.dist <- .jnull(ch.distances())
  if(!missing(genotypes)){
    j.geno <- genotypes$java
  }
  if(!missing(phenotypes)){
    j.pheno <- phenotypes$java
  }
  if(!missing(distances)){
    j.dist <- distances$java
  }
  j.data <- new(ch.data(), j.geno, j.pheno, j.dist)

  # create R object
  data <- list()
  size.ids.names <- c("size", "ids", "names")
  if(!missing(genotypes)){
    data$geno <- genotypes
    tmp <- genotypes
  }
  if(!missing(phenotypes)){
    data$pheno <- phenotypes
    tmp <- phenotypes
  }
  if(!missing(distances)){
    data$dist <- distances
    tmp <- distances
  }
  data[size.ids.names] <- tmp[size.ids.names]
  data$java <- j.data
  # set class
  class(data) <- c("chdata", class(data))

  return(data)

}

#' @export
print.chdata <- function(x, ...){
  cat("Core Hunter data\n")
  cat("----------------\n\n")

  cat("Number of accessions =", x$size, "\n")
  cat("Ids:")
  str(x$ids)
  if(!isTRUE(all.equal(x$ids, x$names))){
    cat("Names:")
    str(x$names)
  }

  for(data.type in c("geno", "pheno", "dist")){
    if(!is.null(x[[data.type]])){
      cat("\n")
      print(x[[data.type]], include.size = FALSE)
    }
  }
}

# -------------------- #
# DISTANCE MATRIX DATA #
# -------------------- #

#' Create Core Hunter distance data from matrix or file.
#'
#' Specify either a symmetric distance matrix or the file from which to read the matrix.
#' See \url{www.corehunter.org} for documentation and examples of the distance matrix
#' file format used by Core Hunter.
#'
#' @param data Symmetric distance matrix. Unique row and column headers are required,
#'  should be the same and are used as item ids. Can be a \code{numeric} matrix or a data frame.
#'  The data frame may optionally include a first column \code{NAME} used to assign names to some
#'  or all individuals. The remaining columns should be \code{numeric}.
#' @param file File from which to read the distance matrix.
#'
#' @return Distance matrix data of class \code{chdist} with elements
#' \describe{
#'  \item{\code{data}}{Distance matrix (\code{numeric} matrix).}
#'  \item{\code{size}}{Number of individuals in the dataset.}
#'  \item{\code{ids}}{Unique item identifiers.}
#'  \item{\code{names}}{Item names. Names of individuals to which no explicit name
#'    has been assigned are equal to the unique \code{ids}.}
#'  \item{\code{java}}{Java version of the data object.}
#'  \item{\code{file}}{Normalized path of file from which data was read (if applicable).}
#' }
#'
#' @examples
#' # create from distance matrix
#' m <- matrix(runif(100), nrow = 10, ncol = 10)
#' diag(m) <- 0
#' # make symmetric
#' m[lower.tri(m)] <- t(m)[lower.tri(m)]
#' # set headers
#' rownames(m) <- colnames(m) <- paste("i", 1:10, sep = "-")
#'
#' dist <- distances(m)
#'
#' # read from file
#' dist.file <- system.file("extdata", "distances.csv", package = "corehunter")
#' dist <- distances(file = dist.file)
#'
#' @import rJava
#' @export
distances <- function(data, file){

  # check input
  if(missing(data) && missing(file)){
    stop("Please specify matrix or file.")
  }
  if(!missing(data) && !missing(file)){
    stop("Please specify either matrix or file, not both.")
  }

  api <- ch.api()

  check.matrix <- function(matrix){
    if(!is.numeric(matrix)){
      stop("Distance matrix should be numeric")
    }
    if(is.null(rownames(matrix)) || is.null(colnames(matrix))){
      stop("Row and column names are required.")
    }
    if(!isSymmetric(matrix)){
      stop("Distance matrix should be symmetric.")
    }
  }

  if(!missing(file)){

    # read from file

    # check file path
    if(!is.character(file)){
      stop("Argument 'file' should be a file path (character).")
    }
    if(!file.exists(file)){
      stop("File 'file' does not exist.")
    }
    file <- normalizePath(file)

    # create Java object from file
    j.data <- api$readDistanceMatrixData(file)
    # read matrix as data frame
    data <- read.autodelim(file)
    # extract distance matrix
    matrix <- extract.matrix(data)

  } else {

    # in memory

    # check type
    if(!is.data.frame(data) && !is.matrix(data)){
      stop("Argument 'matrix' should be a matrix or a data frame.")
    }

    # extract matrix
    if(is.matrix(data)){
      names <- as.character(rep(NA, nrow(data)))
      matrix <- data
    } else {
      # extract names and convert to matrix
      names <- extract.names(data)
      matrix <- extract.matrix(data)
    }
    check.matrix(matrix)

    j.matrix <- .jarray(matrix, dispatch = TRUE)
    j.ids <- .jarray(rownames(data))
    j.names <- .jarray(names)
    j.data <- api$createDistanceMatrixData(j.matrix, j.ids, j.names)

  }

  # obtain ids and names from Java object
  ids <- api$getIds(j.data)
  names <- api$getNames(j.data)

  # create R object
  dist <- list(
    data = matrix,
    size = j.data$getSize(),
    ids = ids,
    names = names,
    java = j.data
  )
  if(!missing(file)){
    dist$file = file
  }
  class(dist) <- c("chdist", "chdata", class(dist))

  return(dist)

}

#' @export
print.chdist <- function(x, include.size = TRUE, ...){
  n <- x$size
  cat("# Precomputed distance matrix\n")
  if(include.size){
    cat(sprintf("\nNumber of accessions = %d\n", n))
    cat("Ids:")
    str(x$ids)
    if(!isTRUE(all.equal(x$ids, x$names))){
      cat("Names:")
      str(x$names)
    }
  }
  if(!is.null(x$file)){
    cat(sprintf("\nRead from file: \"%s\"\n", x$file))
  }
}

# ------------- #
# GENOTYPE DATA #
# ------------- #

#' Create Core Hunter genotype data from data frame, matrix or file.
#'
#' Specify either a data frame or matrix, or a file from which to read the genotypes.
#' See \url{www.corehunter.org} for documentation and examples of the genotype data
#' file format used by Core Hunter.
#'
#' @param data Data frame or matrix containing the genotypes (individuals x markers)
#'  depending on the chosen format:
#'  \describe{
#'    \item{\code{default}}{
#'      Data frame. One row per individual and one or more columns per marker.
#'      Columns contain the names, numbers, references, ... of observed alleles.
#'      Unique row names (item ids) are required and columns should be named
#'      after the marker to which they belong, optionally extended with an
#'      arbitrary suffix starting with a dot (\code{.}), dash (\code{-}) or
#'      underscore (\code{_}) character.
#'    }
#'    \item{\code{biparental}}{
#'      Numeric matrix or data frame. One row per individual and one column per marker.
#'      Data consists of 0, 1 and 2 coding for homozygous (AA), heterozygous (AB) and
#'      homozygous (BB), respectively. Unique row names (item ids) are required and
#'      optionally column (marker) names may be included as well.
#'    }
#'    \item{\code{frequency}}{
#'      Numeric matrix or data frame. One row per individual (or bulk sample) and multiple
#'      columns per marker. Data consists of allele frequencies, grouped per marker in
#'      consecutive columns named after the corresponding marker, optionally extended
#'      with an arbitrary suffix starting with a dot (\code{.}), dash (\code{-}) or
#'      underscore (\code{_}) character.. The allele frequencies of each marker should
#'      sum to one in each sample. Unique row names (item ids) are required.
#'    }
#'    In case a data frame is provided, an optional first column \code{NAME}
#'    may be included to specify item names. The remaining columns should follow
#'    the format as described above.
#'    See \url{www.corehunter.org} for more details about the supported genotype formats.
#'    Note that both the \code{frequency} and \code{biparental} format syntactically also
#'    comply with the \code{default} format but with different semantics, meaning that it
#'    is very important to specify the correct format. Some checks have been built in that
#'    raise warnings in case it seems that the wrong format might have been specified based
#'    on an inspection of the data. If you are sure that you have selected the correct format
#'    these warnings, if any, can be safely ignored.
#'  }
#' @param alleles Allele names per marker (\code{character} vector).
#'  Ignored except when creating \code{frequency} data from a matrix or data frame.
#'  Allele names should be ordered in correspondence with the data columns.
#' @param file File containing the genotype data.
#' @param format Genotype data format, one of \code{default}, \code{biparental} or \code{frequency}.
#'
#' @return Genotype data of class \code{chgeno} with elements
#' \describe{
#'  \item{\code{data}}{Genotypes. Data frame for default format, \code{numeric} matrix for other formats.}
#'  \item{\code{size}}{Number of individuals in the dataset.}
#'  \item{\code{ids}}{Unique item identifiers (\code{character}).}
#'  \item{\code{names}}{Item names (\code{character}). Names of individuals to which no explicit name
#'    has been assigned are equal to the unique \code{ids}.}
#'  \item{\code{markers}}{Marker names (\code{character}).
#'    May contain \code{NA} values in case only some or no marker names were specified.
#'    Marker names are always included for the \code{default} and \code{frequency} format
#'    but are optional for the \code{biparental} format.}
#'  \item{\code{alleles}}{List of character vectors with allele names per marker.
#'    Vectors may contain \code{NA} values in case only some or no allele names were
#'    specified. For \code{biparental} data the two alleles are name \code{"0"} and
#'    \code{"1"}, respectively, for all markers. For the \code{default} format allele
#'    names are inferred from the provided data. Finally, for \code{frequency} data
#'    allele names are optional and may be specified either in the file or through
#'    the \code{alleles} argument when creating this type of data from a matrix or
#'    data frame.}
#'  \item{\code{java}}{Java version of the data object.}
#'  \item{\code{format}}{Genotype data format used.}
#'  \item{\code{file}}{Normalized path of file from which data was read (if applicable).}
#' }
#'
#' @examples
#' # create from data frame or matrix
#'
#' # default format
#' geno.data <- data.frame(
#'  NAME = c("Alice", "Bob", "Carol", "Dave", "Eve"),
#'  M1.1 = c(1,2,1,2,1),
#'  M1.2 = c(3,2,2,3,1),
#'  M2.1 = c("B","C","D","B",NA),
#'  M2.2 = c("B","A","D","B",NA),
#'  M3.1 = c("a1","a1","a2","a2","a1"),
#'  M3.2 = c("a1","a2","a2","a1","a1"),
#'  M4.1 = c(NA,"+","+","+","-"),
#'  M4.2 = c(NA,"-","+","-","-"),
#'  row.names = paste("g", 1:5, sep = "-")
#' )
#' geno <- genotypes(geno.data, format = "default")
#'
#' # biparental (e.g. SNP)
#' geno.data <- matrix(
#'  sample(c(0,1,2), replace = TRUE, size = 1000),
#'  nrow = 10, ncol = 100
#' )
#' rownames(geno.data) <- paste("g", 1:10, sep = "-")
#' colnames(geno.data) <- paste("m", 1:100, sep = "-")
#' geno <- genotypes(geno.data, format = "biparental")
#'
#' # frequencies
#' geno.data <- matrix(
#'  c(0.0, 0.3, 0.7, 0.5, 0.5, 0.0, 1.0,
#'    0.4, 0.0, 0.6, 0.1, 0.9, 0.0, 1.0,
#'    0.3, 0.3, 0.4, 1.0, 0.0, 0.6, 0.4),
#'  byrow = TRUE, nrow = 3, ncol = 7
#' )
#' rownames(geno.data) <- paste("g", 1:3, sep = "-")
#' colnames(geno.data) <- c("M1", "M1", "M1", "M2", "M2", "M3", "M3")
#' alleles <- c("M1-a", "M1-b", "M1-c", "M2-a", "M2-b", "M3-a", "M3-b")
#' geno <- genotypes(geno.data, alleles, format = "frequency")
#'
#' # read from file
#'
#' # default format
#' geno.file <- system.file("extdata", "genotypes.csv", package = "corehunter")
#' geno <- genotypes(file = geno.file, format = "default")
#'
#' # biparental (e.g. SNP)
#' geno.file <- system.file("extdata", "genotypes-biparental.csv", package = "corehunter")
#' geno <- genotypes(file = geno.file, format = "biparental")
#'
#' # frequencies
#' geno.file <- system.file("extdata", "genotypes-frequency.csv", package = "corehunter")
#' geno <- genotypes(file = geno.file, format = "frequency")
#'
#' @import rJava
#' @export
genotypes <- function(data, alleles, file, format){

  # check input
  if(missing(data) && missing(file)){
    stop("Please specify data or file.")
  }
  if(!missing(data) && !missing(file)){
    stop("Please specify either data or file, not both.")
  }
  if(missing(format) || is.null(format)){
    stop("Please specify data format.")
  }
  format <- match.arg(format, c("default", "biparental", "frequency"))

  api <- ch.api()

  if(missing(file)){

    # create from data
    if(format == "default"){

      #---------#
      # default #
      #---------#

      # check data
      if(!is.data.frame(data)){
        stop("Default genotype format: data should be a data frame.")
      } else {
        names <- extract.names(data)
        # erase names
        data$NAME <- NULL
        # convert to character matrix
        matrix <- data
        for(c in 1:ncol(matrix)){
          matrix[[c]] <- as.character(matrix[[c]])
        }
        matrix <- extract.matrix(matrix)
      }
      # check matrix
      if(is.null(rownames(matrix)) || is.null(colnames(matrix))){
        stop("Unique row names (item ids) and column names (marker names) are required.")
      }

      # create data
      j.matrix <- .jarray(matrix, dispatch = TRUE)
      j.ids <- .jarray(rownames(matrix))
      j.names <- .jarray(names)
      j.column.names <- .jarray(colnames(matrix))
      j.data <- api$createDefaultGenotypeData(j.matrix, j.ids, j.names, j.column.names)

    } else if(format == "biparental"){

      #------------#
      # biparental #
      #------------#

      # check data
      if(is.matrix(data)){
        names <- as.character(rep(NA, nrow(data)))
      } else if(is.data.frame(data)){
        # extract names and matrix
        names <- extract.names(data)
        data <- extract.matrix(data)
      } else {
        stop("Biparental genotype format: data should be either a matrix or data frame.")
      }
      # check matrix
      values <- unique(as.vector(data))
      values <- values[!is.na(values)]
      expected <- c(0,1,2)
      if(!is.numeric(data) || length(setdiff(values, expected)) > 0){
        stop("Marker matrix should be numeric (0, 1, 2).")
      }
      # check row names
      if(is.null(rownames(data))){
        stop("Unique row names are required (item ids).")
      }

      # create data
      missing.allele.score <- ch.constants()$MISSING_ALLELE_SCORE
      j.matrix <- .jbyte(data)
      missing <- is.na(j.matrix)
      if(any(missing)){
        j.matrix[missing] <- missing.allele.score
      }
      j.matrix <- .jarray(j.matrix, dispatch = TRUE)
      j.ids <- .jarray(rownames(data))
      j.names <- .jarray(names)
      marker.names <- colnames(data)
      if(is.null(marker.names)){
        marker.names <- as.character(rep(NA, ncol(data)))
      }
      j.marker.names <- .jarray(marker.names)
      j.data <- api$createBiparentalGenotypeData(j.matrix, j.ids, j.names, j.marker.names)

    } else {

      #-----------#
      # frequency #
      #-----------#

      # check data
      if(is.matrix(data)){
        names <- as.character(rep(NA, nrow(data)))
      } else if(is.data.frame(data)){
        # extract names and matrix
        names <- extract.names(data)
        data <- extract.matrix(data)
      } else {
        stop("Frequency genotype format: data should be either a matrix or data frame.")
      }
      # check matrix
      if(!is.numeric(data) || any(data < 0.0, na.rm = TRUE) || any(data > 1.0, na.rm = TRUE)){
        stop("Frequencies should be numeric values between 0.0 and 1.0.")
      }
      # check row and column names
      if(is.null(rownames(data)) || is.null(colnames(data))){
        stop("Unique row names (item ids) and column names (marker names) are required.")
      }
      # check allele names
      if(missing(alleles)){
        alleles <- rep(NA, ncol(data))
      }
      if(!is.vector(alleles) || length(alleles) != ncol(data)){
        stop("Alleles should be a vector of length equal to the number of data columns.")
      }

      # create data
      j.freqs <- .jarray(data, dispatch = TRUE)
      j.ids <- .jarray(rownames(data))
      j.names <- .jarray(names)
      j.column.names <- .jarray(colnames(data))
      j.alleles <- .jarray(as.character(alleles))
      j.data <- api$createFrequencyGenotypeData(j.freqs, j.ids, j.names, j.column.names, j.alleles)

    }

  } else {

    # read from file

    # check file path
    if(!is.character(file)){
      stop("Argument 'file' should be a file path (character).")
    }
    if(!file.exists(file)){
      stop("File '", file, "' does not exist.")
    }
    file <- normalizePath(file)

    # read from file
    j.data <- api$readGenotypeData(file, format)
    # read raw data
    data <- read.autodelim(file)
    # drop names
    data$NAME <- NULL
    # clean frequency format
    if(format == "frequency"){
      # drop allele name row
      data <- data[rownames(data) != "ALLELE", ]
      # convert to numeric
      for(c in 1:ncol(data)){
        data[[c]] <- as.numeric(data[[c]])
      }
    }
    # convert to matrix for frequency and biparental format
    if(format == "frequency" || format == "biparental"){
      data <- as.matrix(data)
    }

  }

  # obtain ids, names and marker names from Java object
  ids <- api$getIds(j.data)
  names <- api$getNames(j.data)
  markers <- api$getMarkerNames(j.data)
  # obtain allele names per marker from Java object
  alleles <- .jevalArray(api$getAlleles(j.data), simplify = TRUE)
  # convert to list of vectors
  if(is.list(alleles)){
    # different number of alleles per marker
    alleles <- lapply(alleles, as.vector)
  } else {
    # same number of alleles for each marker
    alleles <- split(t(alleles), rep(1:nrow(alleles), each = ncol(alleles)))
  }
  # assign marker names to alleles (if specified)
  names(alleles) <- NULL
  if(!all(is.na(markers))){
    names(alleles) <- markers
  }

  # create R object
  geno <- list(
    data = data,
    size = j.data$getSize(),
    ids = ids,
    names = names,
    markers = markers,
    alleles = alleles,
    java = j.data,
    format = format
  )
  if(!missing(file)){
    geno$file <- file
  }
  class(geno) <- c("chgeno", "chdata", class(geno))

  # some checks to detect likely misspecified default format
  if(format == "default"){
    if(length(setdiff(unique(unlist(alleles)), c(0,1,2))) == 0){
      warning("Genotypes seem to be in 'biparental' format, not 'default'. Please check.")
    } else {
      # function to check whether all allele names are in fact frequencies (suspicious!)
      alleles.are.frequencies <- function(a){
        a <- suppressWarnings(as.numeric(unlist(a)))
        if(any(is.na(a))){
          # not all allele names could be converted to numeric
          return(FALSE)
        }
        # all allele names successfully converted to numeric: check range
        return(all(a >= 0.0 & a <= 1.0))
      }
      if(ids[1] == "ALLELE" || alleles.are.frequencies(alleles)){
        warning("Genotypes seem to be in 'frequency' format, not 'default'. Please check.")
      }
    }
  }

  return(geno)

}

#' Get Allele frequency matrix.
#'
#' @param data Core Hunter data containing genotypes
#'
#' @return allele frequency matrix
#'
#' @import rJava
#' @export
getAlleleFrequencies <- function(data){
  if(!is(data, "chdata") && !is(data, "chgeno")){
    stop("Data should be of class 'chdata' or 'chgeno'.")
  }
  if(is(data, "chdata") && !is(data, "chgeno")){
    data <- data$geno
  }
  if(is.null(data)){
    stop("No genotypes available in given data.")
  }
  api <- ch.api()
  freqs <- .jevalArray(api$getAlleleFrequencies(data$java), simplify = TRUE)
  rownames(freqs) <- data$ids
  return(freqs)
}

#' @export
print.chgeno <- function(x, include.size = TRUE, ...){
  n <- x$size
  m <- length(x$markers)
  a <- range(sapply(x$alleles, length))
  a <- ifelse(a[1] == a[2],
              as.character(a[1]),
              sprintf("%d-%d", a[1], a[2])
  )

  cat("# Genotypes\n")

  cat(sprintf("\nFormat = %s\n", x$format))

  if(include.size){
    cat(sprintf("\nNumber of accessions = %d\n", n))
    cat("Ids:")
    str(x$ids)
    if(!isTRUE(all.equal(x$ids, x$names))){
      cat("Names:")
      str(x$names)
    }
  }

  cat(sprintf("\nNumber of markers = %d\n", m))
  if(!all(is.na(x$markers))){
    cat("Marker names:")
    str(x$markers)
  }

  cat(sprintf("Number of alleles per marker = %s\n", a))

  if(!is.null(x$file)){
    cat(sprintf("\nRead from file: \"%s\"\n", x$file))
  }
}

# -------------- #
# PHENOTYPE DATA #
# -------------- #

#' Create Core Hunter phenotype data from data frame or file.
#'
#' Specify either a data frame containing the phenotypic trait observations
#' or a file from which to read the data. See \url{www.corehunter.org} for
#' documentation and examples of the phenotype data format used by Core Hunter.
#'
#' @param data Data frame containing one row per individual and one column per trait.
#'   Unique row and column names are required and used as item and trait ids, respectively.
#'   The data frame may optionally include a first column \code{NAME} used to assign names
#'   to some or all individuals.
#'
#' @param types Variable types (optional).
#'   Vector of characters, each of length one or two.
#'   Ignored when reading from file.
#'
#'   The first letter indicates the scale type and should be one of \code{N} (nominal),
#'   \code{O} (ordinal), \code{I} (interval) or \code{R} (ratio).
#'
#'   The second letter optionally indicates the variable encoding (in Java) and should
#'   be one of \code{B} (boolean), \code{T} (short), \code{I} (integer), \code{L} (long),
#'   \code{R} (big integer), \code{F} (float), \code{D} (double), \code{M} (big decimal),
#'   \code{A} (date) or \code{S} (string). The default encoding is \code{S} (string)
#'   for nominal variables, \code{I} (integer) for ordinal and interval variables
#'   and \code{D} (double) for ratio variables. Interval and ratio variables are
#'   limited to numeric encodings.
#'
#'   If no explicit variable types are specified these are automatically inferred from
#'   the data frame column types and classes, whenever possible. Columns of type
#'   \code{character} are treated as nominal string encoded variables (\code{N}).
#'   Unordered \code{factor} columns are converted to \code{character} and also
#'   treated as string encoded nominals. Ordered factors are converted to
#'   integer encoded interval variables (\code{I}) as described below.
#'   Columns of type \code{logical} are taken to be asymmetric binary variables (\code{NB}).
#'   Finally, \code{integer} and more broadly \code{numeric} columns are treated as integer
#'   encoded interval variables (\code{I}) and double encoded ratio variables (\code{R}),
#'   respectively.
#'
#'   Boolean encoded nominals (\code{NB}) are treated as asymmetric binary variables.
#'   For symmetric binary variables just use the default string encoding (\code{N}
#'   or \code{NS}). Other nominal variables are converted to factors.
#'
#'   Ordinal variables of class \code{ordered} are converted to integers respecting
#'   the order and range of the factor levels and subsequently treated as integer
#'   encoded interval variables (\code{I}). This conversion allows to model the
#'   full range of factor levels also when some might not occur in the data. For other
#'   ordinal variables it is assumed that each value occurs at least once and that
#'   values follow the natural ordering of the chosen data type (in Java).
#'
#'   If explicit types are given for some variables others can still be automatically inferred
#'   by setting their type to \code{NA}.
#'
#' @param min Minimum values of interval or ratio variables (optional).
#'   Numeric vector. Ignored when reading from file.
#'   If undefined for some variables the respective minimum is inferred from the data.
#'   If the data exceeds the minimum it is also updated accordingly.
#'   For nominal and ordinal variables just put \code{NA}.
#' @param max Maximum values of interval or ratio variables (optional).
#'   Numeric vector. Ignored when reading from file.
#'   If undefined for some variables the respective maximum is inferred from the data.
#'   If the data exceeds the maximum it is also updated accordingly.
#'   For nominal and ordinal variables just put \code{NA}.
#'
#' @param file File containing the phenotype data.
#'
#' @return Phenotype data of class \code{chpheno} with elements
#' \describe{
#'  \item{\code{data}}{Phenotypes (data frame).}
#'  \item{\code{size}}{Number of individuals in the dataset.}
#'  \item{\code{ids}}{Unique item identifiers.}
#'  \item{\code{names}}{Item names. Names of individuals to which no explicit name
#'    has been assigned are equal to the unique \code{ids}.}
#'  \item{\code{types}}{Variable types and encodings.}
#'  \item{\code{ranges}}{Variable ranges, when applicable (\code{NA} elsewhere).}
#'  \item{\code{java}}{Java version of the data object.}
#'  \item{\code{file}}{Normalized path of file from which the data was read (if applicable).}
#' }
#'
#' @examples
#' # create from data frame
#' pheno.data <- data.frame(
#'  season = c("winter", "summer", "summer", "winter", "summer"),
#'  yield = c(34.5, 32.6, 22.1, 54.12, 43.33),
#'  size = ordered(c("l", "s", "s", "m", "l"), levels = c("s", "m", "l")),
#'  resistant = c(FALSE, TRUE, TRUE, FALSE, TRUE)
#' )
#' pheno <- phenotypes(pheno.data)
#'
#' # explicit types
#' pheno <- phenotypes(pheno.data, types = c("N", "R", "O", "NB"))
#' # treat last column as symmetric binary, auto infer others
#' pheno <- phenotypes(pheno.data, types = c(NA, NA, NA, "NS"))
#'
#' # explicit ranges
#' pheno <- phenotypes(pheno.data, min = c(NA, 20.0, NA, NA), max = c(NA, 60.0, NA, NA))
#'
#' # read from file
#' pheno.file <- system.file("extdata", "phenotypes.csv", package = "corehunter")
#' pheno <- phenotypes(file = pheno.file)
#'
#' @import rJava
#' @importFrom utils write.csv
#' @export
phenotypes <- function(data, types, min, max, file){

  # check input
  if(missing(data) && missing(file)){
    stop("Data frame or file path is required.")
  }
  if(!missing(data) && !missing(file)){
    stop("Please specify either data frame or file, not both.")
  }

  # utility function
  read.raw.data <- function(file){
    # read raw data
    data <- read.autodelim(file)
    # drop names
    data$NAME <- NULL
    # extract variable types and encodings
    if("TYPE" %in% rownames(data)){
      types <- as.character(data["TYPE",])
    } else {
      stop("Variable types are required.")
    }
    # drop type, min, max
    data <- data[!(rownames(data) %in% c("TYPE", "MIN", "MAX")), , drop = FALSE]
    # convert columns accordingly
    for(c in 1:ncol(data)){
      data[[c]] <- convert.column(data[[c]], types[c])
    }
    # prepare and return result
    result <- list(data = data, types = types)
    return(result)
  }

  api <- ch.api()

  if(missing(file)){

    # write data to file to be read in Java

    # check data type
    if(!is.data.frame(data)){
      stop("Argument 'data' should be a data frame.")
    }

    # check row and column names
    if(is.null(rownames(data)) || is.null(colnames(data))){
      stop("Unique row and column names are required.")
    }

    # extract item names
    names <- extract.names(data)
    data$NAME <- NULL

    # check and/or infer types and ignore bounds for non-numeric variables
    if(missing(types)){
      types <- rep(NA, ncol(data))
    } else if(length(types) != ncol(data)){
      stop("Number of variable types does not correspond to number of data columns.")
    }
    if(missing(min)){
      min = as.numeric(rep(NA, length(types)))
    }
    if(missing(max)){
      max = as.numeric(rep(NA, length(types)))
    }
    for(t in 1:length(types)){
      if(!is.na(types[t])){
        if(!is.character(types[t]) || (nchar(types[t]) != 1 && nchar(types[t]) != 2)){
          stop("Types should consist of one or two characters.")
        }
        if(!is.numeric(data[[t]])){
          min[t] <- max[t] <- NA
        }
      } else {
        # infer type
        col <- data[[t]]
        if(is.character(col)){
          # character treated as nominal string
          types[t] <- "N"
          min[t] <- max[t] <- NA
        } else if(is.factor(col)){
          if(!is.ordered(col)){
            # unordered factor treated as nominal string
            types[t] <- "N"
            min[t] <- max[t] <- NA
          } else {
            # ordered factor: will be converted to integer (see below)
            types[t] <- "O"
          }
        } else if(is.logical(col)){
          types[t] <- "NB"
          min[t] <- max[t] <- NA
        } else if(is.integer(col)){
          types[t] <- "I"
        } else if(is.numeric(col)){
          types[t] <- "R"
        } else {
          stop("Could not automatically infer variable type.")
        }
      }
      # convert ordinal columns of class ordered to integer (interval)
      if(substr(types[t], 1, 1) == "O" && is.ordered(data[[t]])){
        orig.col <- data[[t]]
        data[[t]] <- as.integer(orig.col)
        types[t] <- "I"
        min[t] <- 1
        max[t] <- length(levels(orig.col))
      }
    }

    # convert all columns to characters
    ids <- rownames(data)
    data <- data.frame(lapply(data, as.character), stringsAsFactors = FALSE, check.names = FALSE)
    rownames(data) <- ids

    # add min and max rows (bottom to top)
    if(!is.numeric(max)){
      stop("Maximum values should be numeric.")
    }
    if(length(max) != ncol(data)){
      stop("Number of maximum values does not correspond to number of data columns.")
    }
    data <- rbind(MAX = max, data)
    if(!is.numeric(min)){
      stop("Maximum values should be numeric.")
    }
    if(length(min) != ncol(data)){
      stop("Number of minimum values does not correspond to number of data columns.")
    }
    data <- rbind(MIN = min, data)
    # add type row
    data <- rbind(TYPE = types, data)

    # reinsert names
    names <- c(rep("", sum(rownames(data) %in% c("TYPE", "MIN", "MAX"))), names)
    data <- cbind(NAME = names, data)
    # make row headers first column (ID)
    data <- cbind(ID = rownames(data), data)

    # write temporary file
    tmp <- tempfile(fileext = ".csv")
    write.csv(data, file = tmp, quote = F, row.names = F, na = "")

    # read data into Core Hunter from file
    j.data <- api$readPhenotypeData(tmp)
    # read raw data
    raw.data <- read.raw.data(tmp)

    # remove temporary file
    file.remove(tmp)

  } else {

    # read from file

    # check file path
    if(!is.character(file)){
      stop("Argument 'file' should be a file path (character).")
    }
    if(!file.exists(file)){
      stop("File 'file' does not exist.")
    }
    file <- normalizePath(file)

    # read from file
    j.data <- api$readPhenotypeData(file)
    # read raw data
    raw.data <- read.raw.data(file)

  }

  # obtain ids, names and variable ranges from Java object
  ids <- api$getIds(j.data)
  names <- api$getNames(j.data)
  ranges <- .jevalArray(api$getRanges(j.data), simplify = TRUE)

  # create R object
  pheno <- list(
    data = raw.data$data,
    size = j.data$getSize(),
    ids = ids,
    names = names,
    types = raw.data$types,
    ranges = ranges,
    java = j.data
  )
  if(!missing(file)){
    pheno$file <- file
  }
  class(pheno) <- c("chpheno", "chdata", class(pheno))

  return(pheno)

}

convert.column <- function(col, type){
  enc <- substr(type, 2, 2)
  type <- substr(type, 1, 1)
  # default encoding if not specified
  if(enc == ""){
    enc <- switch(type,
                  "N" = "S", "O" = "I", "I" = "I", "R" = "D",
                  stop(sprintf("Unsupported variable type '%s'.", type))
           )
  }
  # convert to proper encoding
  if(enc == "B"){
    # cfr. Java: boolean
    col <- as.logical(col)
  } else if(enc %in% c("T", "I", "L", "R")){
    # cfr. Java: short, integer, long, big integer
    col <- as.integer(col)
  } else if(enc %in% c("F", "D", "M")){
    # cfr. Java: float, double, big decimal
    col <- as.numeric(col)
  } else if(enc == "S"){
    # cfr. Java: string
    col <- as.character(col)
  } else if(enc == "A"){
    # cfr. Java: date
    col <- as.Date(col, format = "%Y%m%d%H%M%S%z")
  } else {
    stop(sprintf("Unsupported variable encoding '%s'.", enc))
  }
  # convert to proper type
  if(type == "N" && enc != "B"){
    col <- as.factor(col)
  } else if(type == "O"){
    col <- as.ordered(col)
  }
  return(col)
}

#' @export
print.chpheno <- function(x, include.size = TRUE, ...){
  n <- x$size
  traits <- colnames(x$data)
  types <- x$types
  m <- length(traits)
  scales <- sapply(types, substr, 1, 1)
  # find qualitative traits (nominal)
  qual <- traits[which(scales == "N")]
  # find quantitative traits (ordinal, interval, ratio)
  quan <- traits[which(scales %in% c("O", "I", "R"))]

  cat("# Phenotypes\n")

  if(include.size){
    cat(sprintf("\nNumber of accessions = %d\n", n))
    cat("Ids:")
    str(x$ids)
    if(!isTRUE(all.equal(x$ids, x$names))){
      cat("Names:")
      str(x$names)
    }
  }

  format.traits <- function(traits){
    if(length(traits) > 0){
      # wrap in double quotes
      traits <- sprintf('"%s"', traits)
      # collapse
      traits <- paste(traits, collapse = " ")
    } else {
      traits <- "n/a"
    }
    return(traits)
  }

  cat(sprintf("\nNumber of traits = %d\n", m))
  cat(sprintf("Traits: %s\n", format.traits(traits)))

  cat(sprintf("Quantitative traits: %s\n", format.traits(quan)))
  cat(sprintf("Qualitative traits: %s\n", format.traits(qual)))

  if(!is.null(x$file)){
    cat(sprintf("\nRead from file: \"%s\"\n", x$file))
  }
}

# --- #
# I/O #
# --- #

#' Read delimited file.
#'
#' Delegates to \code{\link{read.delim}} where the separator is inferred from the file extension (CSV or TXT).
#' For CSV files the delimiter is set to \code{","} while for TXT file \code{"\t"} is used. Also sets
#' some default argument values as used by Core Hunter.
#'
#' @param file File path.
#' @param ... Further arguments to be passed to  \code{\link{read.delim}}.
#' @inheritParams utils::read.table
#'
#' @return Data frame.
#'
#' @importFrom utils read.delim
#' @export
read.autodelim <- function(file, quote = "'\"",
                           row.names = 1,
                           na.strings = "",
                           check.names = FALSE,
                           strip.white = TRUE,
                           stringsAsFactors = FALSE,
                           ...){
  sep <- switch(tolower(tools::file_ext(file)),
                "csv" = ",",
                "txt" = "\t")
  read.delim(file, sep = sep, quote = quote,
             row.names = row.names, na.string = na.strings, check.names = check.names,
             strip.white = strip.white, stringsAsFactors = stringsAsFactors,
             ...)
}

# ----------------- #
# PRIVATE UTILITIES #
# ----------------- #

#' Wrap distances, genotypes or phenotypes in Core Hunter data.
#'
#' If the given data does not match any of these three classes
#' it is returned unchanged.
#'
#' @param data of class \code{chgeno}, \code{chpheno} or \code{chdist}
#' @return Core Hunter data of class \code{chdata}
#'
#' @importFrom methods is
wrapData <- function(data){
  if(is(data, "chdist")){
    data <- coreHunterData(distances = data)
  }
  if(is(data, "chgeno")){
    data <- coreHunterData(genotypes = data)
  }
  if(is(data, "chpheno")){
    data <- coreHunterData(phenotypes = data)
  }
  if(!is(data, "chdata")){
    stop("Argument 'data' should be of class 'chdata' (see function 'coreHunterData').")
  }
  return(data)
}

# Extract NAME column from data frame. If no NAME column is included
# the row names (unique ids) are used as names.
extract.names <- function(data){
  names <- data$NAME
  if(is.null(names)){
    names <- rep(NA, nrow(data))
  }
  names <- as.character(names)
  # replace blank names with NAs
  names[names == ""] <- NA
  return(names)
}

# extract matrix from data frame with initial NAME column followed by
# the columns of the matrix
extract.matrix <- function(data){
  # store IDs
  ids <- rownames(data)
  # discard names (if set)
  data$NAME <- NULL
  # extract matrix
  matrix <- as.matrix(data)
  if(!is.null(colnames(matrix)) && all(is.na(colnames(matrix)))){
    colnames(matrix) <- NULL
  }
  rownames(matrix) <- ids
  return(matrix)
}



