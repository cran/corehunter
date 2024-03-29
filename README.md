# R Package for Core Hunter 3

### Latest release

[![CRAN version](http://www.r-pkg.org/badges/version/corehunter)](https://CRAN.R-project.org/package=corehunter)
![Build Status](https://github.com/corehunter/corehunter3-r/actions/workflows/check-standard.yaml/badge.svg)


### Development snapshot

![Build Status](https://github.com/corehunter/corehunter3-r/actions/workflows/check-standard.yaml/badge.svg?branch=develop)

Core Hunter is a tool to sample diverse, representative subsets from large germplasm collections, with minimum redundancy. Such so-called core collections have applications in plant breeding and genetic resource management in general. Core Hunter can construct cores based on genetic marker data, phenotypic traits or precomputed distance matrices, optimizing one of many provided evaluation measures depending on the precise purpose of the core (e.g. high diversity, representativeness, or allelic richness). In addition, multiple measures can be simultaneously optimized as part of a weighted index to bring the different perspectives closer together. The Core Hunter library is implemented in Java 8 as an open source project (see 
<http://www.corehunter.org>).

Version 3 has been recoded from scratch using the [JAMES framework](http://www.jamesframework.org) which provides the applied optimization algorithms.
Requirements
------------

A [Java Runtime Environment] (http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) (JRE) version 8 or later is required to run Core Hunter.

Getting started
---------------

The package `corehunter` can be installed from CRAN with

```R
> install.packages("corehunter")
```

All provided functions are documented in the package, including many examples, for example try

```R
> ?corehunter
> ?sampleCore
> ?genotypes
> ?phenotypes
> ?distances
```

For more information please visit <http://www.corehunter.org>.

Supported data types
--------------------

Core Hunter 3 supports multiple types of genetic marker data, phenotypic traits and precomputed distance matrices. See <http://www.corehunter.org/data> for more details. Data can be loaded from files, data frames and matrices.

Evaluation measures
-------------------

One of the main strengths of Core Hunter is that it can directly optimize a number of different evaluation measures. If desired, multiple measures can be simultaneously optimized as part of a weighted index. The measures included in Core Hunter 3 are listed below.

#### Distance based measures

- Average entry-to-nearest-entry distance (diversity)
- Average accession-to-nearest-entry distance (representativeness)
- Average entry-to-entry distance (provided for historical reasons, not preferred)

Gower's distance is used to compute distances from phenotypic traits, and both the Modified Roger's as well as Cavalli-Sforza & Edwards distances are supported for genetic marker data. Alternatively, a precomputed distance matrix can be used.

#### Allelic richness

- Shannon's index
- Expected heterozygosity
- Allele coverage

Available for genetic marker data only.
