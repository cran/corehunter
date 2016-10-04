# --- #
# API #
# --- #

ch.api <- function(){
  J("org.corehunter.API")
}

# ---- #
# DATA #
# ---- #

ch.data <- function(){
  J("org.corehunter.data.CoreHunterData")
}
ch.distances <- function(){
  J("org.corehunter.data.DistanceMatrixData")
}
ch.genotypes <- function(){
  J("org.corehunter.data.GenotypeData")
}
ch.phenotypes <- function(){
  J("uno.informatics.data.dataset.FeatureData")
}

# ---------- #
# OBJECTIVES #
# ---------- #

ch.obj <- function(){
  J("org.corehunter.CoreHunterObjective")
}

ch.objectives <- function(objectives){
  api <- ch.api()
  j.objectives <- lapply(objectives, function(obj){
    if(is.null(obj$range)){
      # without normalization range
      api$createObjective(obj$type, obj$meas, obj$weight)
    } else {
      # with normalization range
      api$createObjective(obj$type, obj$meas, obj$weight, obj$range[1], obj$range[2])
    }
  })
  return(j.objectives)
}

# ---------------- #
# INDEX CONVERSION #
# ---------------- #

toJavaIndices <- function(indices){
  as.integer(indices-1)
}

toRIndices <- function(indices){
  indices + 1
}

# --------------------------- #
# COMPLEX ARRAYS AND MATRICES #
# --------------------------- #

# create Integer matrix
# .jIntegerMatrix <- function(matrix){
#   .jarray(
#     apply(matrix, 1, function(row){
#       .jarray(lapply(row, function(value){
#           .jnew("java.lang.Integer", as.integer(value))
#         }), "java.lang.Integer"
#       )
#     }),"[Ljava.lang.Integer;"
#   )
# }





