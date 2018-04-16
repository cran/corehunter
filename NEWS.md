Core Hunter 3
=============

Version 3.2.1 (16/04/2018)
--------------------------

 - Added support for Java 9.

Version 3.2.0 (08/09/2017)
--------------------------

 - New option to specify a set of always and/or never selected accessions.
 - Added support for step-based in addition to time-based stop conditions.
 - Using `set.seed` prior to executing Core Hunter, in combination with step-based stop conditions, now yields reproducible results.
 - Largely reduced memory footprint of genotype and distance data.
 - Default maximum time without improvement is now only applied in case no explicit stop conditions have been specified.
 - Fixed issue when loading phenotype data with a single trait.

Version 3.1.0 (27/01/2017)
--------------------------

 - More informative output when printing Core Hunter data objects (#19, #23).
 - It is now required to explicitly specify the format when reading genotype data. Moreover, a warning is raised in case it seems that the wrong format may have been selected, through an inspection of the data that was read (#20).
 - Fixed issues when reading files with single or double quoted values.
 - Improved test coverage.

Version 3.0.1 (04/10/2016)
--------------------------

 - The available Java version is now checked when installing or loading the package. Core Hunter requires Java 8 or higher.  
 - Package title and description have been updated.

Version 3.0.0 (03/10/2016)
--------------------------

 - Initial release of Core Hunter 3 package
