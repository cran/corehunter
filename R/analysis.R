# -------- #
# ANALYSIS #
# -------- #

#' Evaluate a core collection using the specified objective.
#'
#' @param core A core collection of class \code{chcore}, or a
#'   numeric or character vector indicating the indices or ids,
#'   respectively, of the individuals in the evaluated core.
#' @param data Core Hunter data (\code{chdata}) containing genotypes,
#'   phenotypes and/or a precomputed distance matrix. Can also be an
#'   object of class \code{chdist}, \code{chgeno} or \code{chpheno}
#'   if only one type of data is provided.
#' @param objective Objective function (\code{chobj}) used to evaluate the core.
#'
#' @return Value of the core when evaluated with the chosen objective (numeric).
#'
#' @examples
#' \donttest{
#' data <- exampleData()
#' core <- sampleCore(data, objective("EN", "PD"))
#' evaluateCore(core, data, objective("EN", "PD"))
#' evaluateCore(core, data, objective("AN", "MR"))
#' evaluateCore(core, data, objective("EE", "GD"))
#' evaluateCore(core, data, objective("CV"))
#' evaluateCore(core, data, objective("HE"))
#' }
#'
#' @seealso \code{\link{coreHunterData}}, \code{\link{objective}}
#'
#' @import rJava
#' @export
evaluateCore <- function(core, data, objective){
  UseMethod("evaluateCore")
}

#' @export
evaluateCore.chcore <- function(core, data, objective){
  evaluateCore(core$sel, data, objective)
}

#' @export
evaluateCore.character <- function(core, data, objective){
  # convert ids to indices
  api <- ch.api()
  data <- wrapData(data)
  core <- toRIndices(api$getIndicesFromIds(data$java, .jarray(core)))
  evaluateCore(core, data, objective)
}

#' @export
evaluateCore.numeric <- function(core, data, objective){

  data <- wrapData(data)
  # check arguments
  if(!is(data, 'chdata')){
    stop("Argument 'data' should be of class 'chdata' (see function 'coreHunterData').")
  }
  if(!is(objective, 'chobj')){
    stop("Argument 'objective' should be of class 'chobj' (see function 'objective').")
  }

  # convert objective to Java
  j.objective <- ch.objectives(list(objective))[[1]]

  # evaluate core
  api <- ch.api()
  j.core <- .jarray(toJavaIndices(core))
  j.data <- data$java
  value <- api$evaluateCore(j.core, j.data, j.objective)

  return(value)

}
