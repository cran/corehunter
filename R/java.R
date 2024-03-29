# --- #
# API #
# --- #

ch.api <- function(){
  J("org.corehunter.API")
}

# --------- #
# CONSTANTS #
# --------- #

ch.constants <- function(){
  J("org.corehunter.util.CoreHunterConstants")
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
  J("org.corehunter.data.FrequencyGenotypeData")
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

# ------------------ #
# GENERAL JAVA STUFF #
# ------------------ #

java.version <- function(){
  version.string <- java.version.string()
  if(startsWith(version.string, "1.")){
    # JDK <= 8
    version <- as.integer(strsplit(version.string, ".", fixed = TRUE)[[1]][2])
  } else {
    # JDK 9+
    short.version.string <- strsplit(version.string, "[-+]")[[1]][1]
    version <- as.integer(strsplit(short.version.string, ".", fixed = TRUE)[[1]][1])
  }
  return(version)
}

java.version.string <- function(){
  version.string <- J("java.lang.System")$getProperty("java.version")
  return(version.string)
}





