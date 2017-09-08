# ------------- #
# CORE SAMPLING #
# ------------- #

#' Determine normalization ranges of all objectives in a multi-objective configuration.
#'
#' Executes an independent stochastic hill-climbing search (random descent) per objective
#' to approximate the optimal solution for each objective, from which a suitable normalization
#' range is inferred based on the Pareto minima/maxima. These normalization searches are
#' executed in parallel.
#'
#' For an objective that is being maximized, the upper bound is set to the value of the best
#' solution for that objective, while the lower bound is set to the Pareto minimum, i.e. the
#' minimum value obtained when evaluating all optimal solutions (for each single objective)
#' with the considered objective. For an objective that is being minimized, the roles of
#' upper and lower bound are interchanged, and the Pareto maximum is used instead.
#'
#' Because Core Hunter uses stochastic algorithms, repeated runs may produce different
#' results. To eliminate randomness, you may set a random number generation seed using
#' \code{\link{set.seed}} prior to executing Core Hunter. In addition, when reproducible
#' results are desired, it is advised to use step-based stop conditions instead of the
#' (default) time-based criteria, because runtimes may be affected by external factors,
#' and, therefore, a different number of steps may have been performed in repeated runs
#' when using time-based stop conditions.
#'
#' @param data Core Hunter data (\code{chdata}) containing genotypes,
#'   phenotypes and/or a precomputed distance matrix. Can also be an
#'   object of class \code{chdist}, \code{chgeno} or \code{chpheno}
#'   if only one type of data is provided.
#' @param obj List of objectives (\code{chobj}).
#'   If no objectives are specified Core Hunter maximizes a weighted
#'   index including the default entry-to-nearest-entry distance
#'   (\code{EN}) for each available data type.
#'   For genotypes, the Modified Roger's distance (\code{MR}) is
#'   used. For phenotypes, Gower's distance (\code{GD}) is applied.
#' @param size Desired core subset size (numeric). If larger than one the value
#'   is used as the absolute core size after rounding. Else it is used as the
#'   sampling rate and multiplied with the dataset size to determine the size of
#'   the core. The default sampling rate is 0.2.
#' @param always.selected vector with indices (integer) or ids (character) of
#'   items that should always be selected in the core collection
#' @param never.selected vector with indices (integer) or ids (character) of
#'   items that should never be selected in the core collection
#' @param mode Execution mode (\code{default} or \code{fast}). In default mode,
#'   the normalization searches terminate when no improvement is found for ten
#'   seconds. In fast mode, searches terminate as soon as no improvement is
#'   made for two seconds. These stop conditions can be overridden using arguments
#'   \code{time}, \code{impr.time}, \code{steps} and/or \code{impr.steps}. In
#'   \code{default} mode, the value of the latter two, step-based conditions is
#'   multiplied with 500, in line with the behaviour of \code{\link{sampleCore}}
#'   when executed in \code{default} mode.
#' @param time Absolute runtime limit in seconds. Not used by default (\code{NA}).
#'   If used, it should be a strictly positive value, which is rounded to the
#'   nearest integer.
#' @param impr.time Maximum time without improvement in seconds. If no explicit
#'   stop conditions are specified, the maximum time without improvement defaults
#'   to ten or two seconds, when executing Core Hunter in \code{default} or
#'   \code{fast} mode, respectively. If a custom improvement time is specified,
#'   it should be strictly positive and is rounded to the nearest integer.
#' @param steps Maximum number of search steps. Not used by default (\code{NA}).
#'              If used, it should be a strictly positive value, which is rounded
#'              to the nearest integer. In \code{default} mode, the value is
#'              multiplied with 500, in line with the behaviour of
#'              \code{\link{sampleCore}} when executed in \code{default} mode.
#' @param impr.steps Maximum number of steps without improvement. Not used by
#'                   default (\code{NA}). If used, it should be a strictly
#'                   positive value, which is rounded to the nearest integer.
#'                   In \code{default} mode, the value is multiplied with 500,
#'                   in line with the behaviour of \code{\link{sampleCore}}
#'                   when executed in \code{default} mode.
#'
#' @return Numeric matrix with one row per objective and two columns:
#' \describe{
#'  \item{\code{lower}}{Lower bound of normalization range.}
#'  \item{\code{upper}}{Upper bound of normalization range.}
#' }
#'
#' @examples
#' \donttest{
#' data <- exampleData()
#'
#' # maximize entry-to-nearest-entry distance between genotypes and phenotypes (equal weight)
#' objectives <- list(objective("EN", "MR"), objective("EN", "GD"))
#' # get normalization ranges for default size (20%)
#' ranges <- getNormalizationRanges(data, obj = objectives, mode = "fast")
#'
#' # set normalization ranges and sample core
#' objectives <- lapply(1:2, function(o){setRange(objectives[[o]], ranges[o,])})
#' core <- sampleCore(data, obj = objectives)
#' }
#'
#' @seealso \code{\link{coreHunterData}}, \code{\link{objective}}
#'
#' @import rJava
#' @importFrom methods is
#' @export
getNormalizationRanges <- function(data, obj, size = 0.2,
                                   always.selected = integer(0),
                                   never.selected = integer(0),
                                   mode = c("default", "fast"),
                                   time = NA, impr.time = NA,
                                   steps = NA, impr.steps = NA){

  # check mode
  mode <- match.arg(mode)
  # check and process stop conditions
  time <- checkTimeOrSteps(time, "Time limit")
  impr.time <- checkTimeOrSteps(impr.time, "Maximum time without improvement")
  steps <- checkTimeOrSteps(steps, "Maximum number of search steps")
  impr.steps <- checkTimeOrSteps(impr.steps, "Maximum number of steps without improvement")

  # check objectives or set default
  obj <- defaultObjectives(data, obj)
  # create arguments
  j.args <- createArguments(data, obj, size, always.selected, never.selected, normalize = TRUE)

  # run Core Hunter normalization
  api <- ch.api()
  ranges <- .jevalArray(
    api$getNormalizationRanges(j.args, mode,
                               time, impr.time, .jlong(steps), .jlong(impr.steps),
                               genSeed()),
    simplify = TRUE
  )
  obj.ids <- sapply(obj, function(o){
    id <- o$type
    if(!is.null(o$meas)){
      id <- paste(id, o$meas, sep = "-")
    }
    return(id)
  })
  rownames(ranges) <- obj.ids
  colnames(ranges) <- c("lower", "upper")

  # return result
  return(ranges)

}

#' Sample a core collection.
#'
#' Sample a core collection from the given data.
#'
#' Because Core Hunter uses stochastic algorithms, repeated runs may produce different
#' results. To eliminate randomness, you may set a random number generation seed using
#' \code{\link{set.seed}} prior to executing Core Hunter. In addition, when reproducible
#' results are desired, it is advised to use step-based stop conditions instead of the
#' (default) time-based criteria, because runtimes may be affected by external factors,
#' and, therefore, a different number of steps may have been performed in repeated runs
#' when using time-based stop conditions.
#'
#' @param data Core Hunter data (\code{chdata}) containing genotypes,
#'   phenotypes and/or a precomputed distance matrix. Typically the
#'   data is obtained with \code{\link{coreHunterData}}. Can also be
#'   an object of class \code{chdist}, \code{chgeno} or \code{chpheno}
#'   if only one type of data is provided.
#' @param obj Objective or list of objectives (\code{chobj}).
#'   If no objectives are specified Core Hunter maximizes a weighted
#'   index including the default entry-to-nearest-entry distance
#'   (\code{EN}) for each available data type, with equal weight.
#'   For genotypes, the Modified Roger's distance (\code{MR}) is
#'   used. For phenotypes, Gower's distance (\code{GD}) is applied.
#' @param size Desired core subset size (numeric). If larger than one the value
#'   is used as the absolute core size after rounding. Else it is used as the
#'   sampling rate and multiplied with the dataset size to determine the size of
#'   the core. The default sampling rate is 0.2.
#' @param always.selected vector with indices (integer) or ids (character) of
#'   items that should always be selected in the core collection
#' @param never.selected vector with indices (integer) or ids (character) of
#'   items that should never be selected in the core collection
#' @param mode Execution mode (\code{default} or \code{fast}). In default mode,
#'   Core Hunter uses an advanced parallel tempering search algorithm and terminates
#'   when no improvement is found for ten seconds. In fast mode, a simple stochastic
#'   hill-climbing algorithm is applied and Core Hunter terminates as soon as no
#'   improvement is made for two seconds. Stop conditions can be overridden with
#'   arguments \code{time} and \code{impr.time}.
#' @param normalize If \code{TRUE} (default), the applied objectives in a multi-objective
#'   configuration (two or more objectives) are automatically normalized prior to execution.
#'   For single-objective configurations, this argument is ignored.
#'
#'   Normalization requires an independent preliminary search per objective (fast stochastic
#'   hill-climber, executed in parallel for all objectives). The same stop conditions, as
#'   specified for the main search, are also applied to each normalization search. In
#'   \code{default} execution mode, however, any step-based stop conditions are multiplied
#'   by 500 for the normalization searches, because in that case the main search (parallel
#'   tempering) executes 500 stochastic hill-climbing steps per replica, in a single step
#'   of the main search.
#'
#'   Normalization ranges can also be precomputed (see \code{\link{getNormalizationRanges}})
#'   or manually specified in the objectives to save computation time when sampling core
#'   collections. This is especially useful when multiple cores are sampled for the same
#'   objectives, with possibly varying weights.
#' @param time Absolute runtime limit in seconds. Not used by default (\code{NA}).
#'   If used, it should be a strictly positive value, which is rounded to the
#'   nearest integer.
#' @param impr.time Maximum time without improvement in seconds. If no explicit
#'   stop conditions are specified, the maximum time without improvement defaults
#'   to ten or two seconds, when executing Core Hunter in \code{default} or
#'   \code{fast} mode, respectively. If a custom improvement time is specified,
#'   it should be strictly positive and is rounded to the nearest integer.
#' @param steps Maximum number of search steps. Not used by default (\code{NA}).
#'              If used, it should be a strictly positive value, which is rounded
#'              to the nearest integer. The number of steps applies to the main
#'              search. Details of how this stop condition is transferred to
#'              normalization searches, in a multi-objective configuration, are
#'              provided in the description of the argument \code{normalize}.
#' @param impr.steps Maximum number of steps without improvement. Not used by
#'                   default (\code{NA}). If used, it should be a strictly
#'                   positive value, which is rounded to the nearest integer.
#'                   The maximum number of steps without improvement applies
#'                   to the main search. Details of how this stop condition is
#'                   transferred to normalization searches, in a multi-objective
#'                   configuration, are provided in the description of the argument
#'                   \code{normalize}.
#' @param indices If \code{TRUE}, the result contains the indices instead of ids
#'   (default) of the selected individuals.
#' @param verbose If \code{TRUE}, search progress messages are printed to the console.
#'   Defaults to \code{FALSE}.
#'
#' @return Core subset (\code{chcore}). It has an element \code{sel}
#'  which is a character or numeric vector containing the sorted ids or indices,
#'  respectively, of the selected individuals (see argument \code{indices}).
#'  In addition the result has one or more elements that indicate the value
#'  of each objective function that was included in the optimization.
#'
#' @examples
#' \donttest{
#' data <- exampleData()
#'
#' # default size, maximize entry-to-nearest-entry Modified Rogers distance
#' obj <- objective("EN", "MR")
#' core <- sampleCore(data, obj)
#'
#' # fast mode
#' core <- sampleCore(data, obj, mode = "f")
#' # absolute size
#' core <- sampleCore(data, obj, size = 25)
#' # relative size
#' core <- sampleCore(data, obj, size = 0.1)
#'
#' # other objective: minimize accession-to-nearest-entry precomputed distance
#' core <- sampleCore(data, obj = objective(type = "AN", measure = "PD"))
#' # multiple objectives (equal weight)
#' core <- sampleCore(data, obj = list(
#'  objective("EN", "PD"),
#'  objective("AN", "GD")
#' ))
#' # multiple objectives (custom weight)
#' core <- sampleCore(data, obj = list(
#'  objective("EN", "PD", weight = 0.3),
#'  objective("AN", "GD", weight = 0.7)
#' ))
#'
#' # custom stop conditions
#' core <- sampleCore(data, obj, time = 5, impr.time = 2)
#' core <- sampleCore(data, obj, steps = 300)
#'
#' # print progress messages
#' core <- sampleCore(data, obj, verbose = TRUE)
#' }
#'
#' @seealso \code{\link{coreHunterData}}, \code{\link{objective}}, \code{\link{getNormalizationRanges}}
#'
#' @import rJava naturalsort
#' @importFrom methods is
#' @export
sampleCore <- function(data, obj, size = 0.2,
                       always.selected = integer(0), never.selected = integer(0),
                       mode = c("default", "fast"), normalize = TRUE,
                       time = NA, impr.time = NA,
                       steps = NA, impr.steps = NA,
                       indices = FALSE, verbose = FALSE){
  # check mode
  mode <- match.arg(mode)
  # check and process stop conditions
  time <- checkTimeOrSteps(time, "Time limit")
  impr.time <- checkTimeOrSteps(impr.time, "Maximum time without improvement")
  steps <- checkTimeOrSteps(steps, "Maximum number of search steps")
  impr.steps <- checkTimeOrSteps(impr.steps, "Maximum number of steps without improvement")

  # check logicals
  if(!is.logical(normalize)){
    stop("Argument 'normalize' should be a logical.")
  }
  if(!is.logical(indices)){
    stop("Argument 'indices' should be a logical.")
  }
  if(!is.logical(verbose)){
    stop("Argument 'verbose' should be a logical.")
  }

  # check objectives or set default
  obj <- defaultObjectives(data, obj)
  # create arguments
  j.args <- createArguments(data, obj, size, always.selected, never.selected, normalize)

  # run Core Hunter
  api <- ch.api()
  sel <- api$sampleCore(j.args, mode,
                        time, impr.time, .jlong(steps), .jlong(impr.steps),
                        genSeed(), !verbose)
  if(indices){
    sel <- toRIndices(sel)
  } else {
    # convert indices to ids
    sel <- api$getIdsFromIndices(j.args$getData(), .jarray(sel))
  }
  # sort selection
  sel <- naturalsort(sel)

  # wrap result
  core <-list(
    sel = sel
  )
  # add objective function values
  for(o in obj){
    value <- evaluateCore(sel, data, o)
    if(is.null(o$meas)){
      core[[o$type]] <- value
    } else {
      if(is.null(core[[o$type]])){
        core[[o$type]] <- list()
      }
      core[[o$type]][[o$meas]] <- value
    }
  }
  # set class and return
  class(core) <- c("chcore", class(core))
  return(core)

}

checkTimeOrSteps <- function(value, description){
  if(!is.na(value)){
    if(!is.numeric(value)){
      stop(sprintf("%s should be numeric.", description))
    }
    value <- as.integer(round(value))
    if(value <= 0){
      stop(sprintf("%s should be a positive number.", description))
    }
  } else {
    value <- as.integer(-1)
  }
  return(value)
}

defaultObjectives <- function(data, obj){
  # wrap and check data class
  data <- wrapData(data)
  # set default objectives or check given objectives
  j.data <- data$java
  # create CH API
  api <- ch.api()
  # set default objectives
  if(missing(obj)){
    obj <- api$createDefaultObjectives(j.data)
    obj <- lapply(obj, function(o){
      objective(
        type = o$getObjectiveType()$getAbbreviation(),
        measure = o$getMeasure()$getAbbreviation(),
        weight = o$getWeight()
      )
    })
  }
  # wrap single objective in list
  if(is(obj, 'chobj')){
    obj <- list(obj)
  }
  # check objectives
  if(!all(sapply(obj, is, 'chobj'))){
    stop("Objectives should be of class 'chobj'.")
  }
  if(length(unique(obj)) != length(obj)){
    stop("Duplicate objectives.")
  }
  return(obj)
}

createArguments <- function(data, obj, size,
                            always.selected,
                            never.selected,
                            normalize){

  # wrap and check data class
  data <- wrapData(data)
  # set and check size
  if(!is.numeric(size)){
    stop("Core 'size' should be numeric.")
  }
  n <- data$size
  if(size > 0 && size < 1){
    size <- size * n
  }
  size <- round(size)
  if(size < 2 || size >= n){
    stop(sprintf("Core 'size' should be >= 2 and < %d (dataset size). Got: %d.", n, size))
  }

  # init API
  api <- ch.api()

  # check always/never selected and convert to Java indices
  convert <- function(subset){
    if(is.character(subset)){
      ind <- api$getIndicesFromIds(data$java, .jarray(subset))
    } else if(is.integer(subset) || (is.numeric(subset) && all(subset == as.integer(subset)))){
      ind <- toJavaIndices(subset)
    } else {
      stop("Arguments 'always.selected' and 'never.selected' ",
           "should be of type integer (indices) or character (ids).")
    }
    return(ind)
  }
  always.selected <- convert(always.selected)
  never.selected <- convert(never.selected)

  # convert objectives to Java objects
  j.obj <- ch.objectives(obj)

  # create Core Hunter arguments
  j.data <- data$java
  j.size <- as.integer(size)
  j.obj.array <- .jarray(j.obj, contents.class = ch.obj()@name)
  j.always.selected <- .jarray(always.selected)
  j.never.selected <- .jarray(never.selected)
  j.args <- api$createArguments(
    j.data, j.size, j.obj.array,
    j.always.selected, j.never.selected,
    normalize
  )

  return(j.args)

}

#' @importFrom stats runif
genSeed <- function(){
  .jlong(ceiling(runif(1, 0, 2^31-1)))
}

# ---------- #
# OBJECTIVES #
# ---------- #

#' Create Core Hunter objective.
#'
#' The following optimization objectives are supported by Core Hunter:
#' \describe{
#'  \item{\code{EN}}{
#'    Average entry-to-nearest-entry distance (default). Maximizes the average distance
#'    between each selected individual and the closest other selected item
#'    in the core. Favors diverse cores in which each individual is sufficiently
#'    different from the most similar other selected item (low redundancy).
#'    Multiple distance measures are provided to be used with this objective (see below).
#'  }
#'  \item{\code{AN}}{
#'    Average accession-to-nearest-entry distance. Minimizes the average distance
#'    between each individual (from the full dataset) and the closest selected item
#'    in the core (which can be the individual itself). Favors representative cores
#'    in which all items from the original dataset are represented by similar individuals
#'    in the selected subset. Multiple distance measures are provided to be used with this
#'    objective (see below).
#'  }
#'  \item{\code{EE}}{
#'    Average entry-to-entry distance. Maximizes the average distance between
#'    each pair of selected individuals in the core. This objective is related to
#'    the entry-to-nearest-entry (EN) distance but less effectively avoids redundant,
#'    similar individuals in the core. In general, use of \code{EN} is preferred.
#'    Multiple distance measures are provided to be used with this objective (see below).
#'  }
#'  \item{\code{SH}}{
#'    Shannon's allelic diversity index. Maximizes the entropy, as used in information
#'    theory, of the selected core. Independently takes into account all allele frequencies,
#'    regardless of the locus (marker) where to which the allele belongs. Requires genotypes.
#'  }
#'  \item{\code{HE}}{
#'    Expected proportion of heterozygous loci. Maximizes the expected proportion of heterozygous
#'    loci in offspring produced from random crossings within the selected core. In contrast to
#'    Shannon's index (\code{SH}) this objective treats each marker (locus) with equal importance,
#'    regardless of the number of possible alleles for that marker. Requires genotypes.
#'  }
#'  \item{\code{CV}}{
#'    Allele coverage. Maximizes the proportion of alleles observed in the full dataset that are
#'    retained in the selected core. Requires genotypes.
#'  }
#' }
#' The first three objective types (\code{EN}, \code{AN} and \code{EE}) aggregate pairwise distances
#' between individuals. These distances can be computed using various measures:
#' \describe{
#' \item{\code{MR}}{
#'    Modified Rogers distance (default). Requires genotypes.
#'  }
#'  \item{\code{CE}}{
#'    Cavalli-Sforza and Edwards distance. Requires genotypes.
#'  }
#'  \item{\code{GD}}{
#'    Gower distance. Requires phenotypes.
#'  }
#'  \item{\code{PD}}{
#'    Precomputed distances. Uses the precomputed distance matrix of the dataset.
#'  }
#' }
#'
#' @param type Objective type, one of \code{EN} (default), \code{AN}, \code{EE},
#'   \code{SH}, \code{HE} or \code{CV} (see description). The former three
#'   objectives are distance based and require to choose a distance
#'   \code{measure}. By default, Modified Roger's distance is used,
#'   computed from the genotypes.
#' @param measure Distance measure used to compute the distance between two
#'   individuals, one of \code{MR} (default), \code{CE}, \code{GD} or \code{PD}
#'   (see description). Ignored when \code{type} is \code{SH}, \code{HE} or
#'   \code{CV}.
#' @param weight Weight assigned to the objective when maximizing a weighted
#'   index. Defaults to 1.0.
#' @param range Normalization range [l,u] of the objective when maximizing a weighted
#'   index. By default the range is not set (\code{NULL}) and will be determined
#'   automatically prior to execution, if normalization is enabled (default).
#'   Values are rescaled to [0,1] with the linear formula
#'   \eqn{
#'    v' = (v - l)/(u - l)
#'   }.
#'   When an explicit normalization range is set, it overrides the automatically inferred
#'   range. Also, setting the range for all included objectives reduces the computation time
#'   when sampling a multi-objective core collection. In case of repeated sampling from the
#'   same dataset with the same objectives and size, it is therefore advised to determine the
#'   normalization ranges only once using \code{\link{getNormalizationRanges}} so that
#'   they can be reused for all executions.
#'
#' @return Core Hunter objective of class \code{chobj} with elements
#' \describe{
#'  \item{\code{type}}{Objective type.}
#'  \item{\code{meas}}{Distance measure (if applicable).}
#'  \item{\code{weight}}{Assigned weight.}
#'  \item{\code{range}}{Normalization range (if specified).}
#' }
#'
#' @examples
#' objective()
#' objective(meas = "PD")
#' objective("EE", "GD")
#' objective("HE")
#' objective("EN", "MR", range = c(0.150, 0.300))
#' objective("AN", "MR", weight = 0.5, range = c(0.150, 0.300))
#'
#' @seealso \code{\link{getNormalizationRanges}}, \code{\link{setRange}}
#'
#' @export
objective <- function(type = c("EN", "AN", "EE", "SH", "HE", "CV"),
                      measure = c("MR", "CE", "GD", "PD"),
                      weight = 1.0, range = NULL){
  # check arguments
  type <- match.arg(type)
  measure <- match.arg(measure)
  if(!is.numeric(weight) || weight < 0.0){
    stop("Objective 'weight' should be a positive number.")
  }
  # create objective
  obj <- list(
    type = type
  )
  if(type %in% c("EE", "EN", "AN")){
    obj$meas <- measure
  }
  obj$weight <- weight
  class(obj) <- c("chobj", class(obj))
  # set range if specified
  if(!is.null(range)){
    obj <- setRange(obj, range)
  }

  return(obj)
}

#' Set the normalization range of the given objective.
#'
#' See argument \code{range} of \code{\link{objective}} for details.
#'
#' @param obj Core Hunter objective of class \code{chobj}.
#' @param range Normalization range [l,u].
#'              See argument \code{range} of \code{\link{objective}} for details.
#'
#' @return Objective including normalization range.
#'
#' @importFrom methods is
#' @seealso \code{\link{objective}}
#' @export
setRange <- function(obj, range){
  if(!is(obj, "chobj")){
    stop("Objective 'obj' should be of class 'chobj'.")
  }
  if(!is.numeric(range) || length(range) != 2){
    stop("Normalization range should be a numeric vector of length two.")
  }
  if(any(is.na(range))){
    stop("Normalization range not fully defined (contains NA). Set to NULL to omit range.")
  }
  if(range[1] > range[2]){
    stop("Lower bound of normalization range exceeds upper bound.")
  }
  obj$range <- range
  return(obj)
}

#' @export
print.chobj <- function(x, ...){
  prefix <- "Core Hunter objective"
  range <- "N/A"
  if(!is.null(x$range)){
    range <- sprintf("[%f, %f]", x$range[1], x$range[2])
  }
  if(!is.null(x$meas)){
    cat(sprintf("%s: %s (measure = %s, weight = %.2f, range = %s)", prefix, x$type, x$meas, x$weight, range))
  } else {
    cat(sprintf("%s: %s (weight = %.2f, range = %s)", prefix, x$type, x$weight, range))
  }
}







