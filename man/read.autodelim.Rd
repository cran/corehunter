% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{read.autodelim}
\alias{read.autodelim}
\title{Read delimited file.}
\usage{
read.autodelim(
  file,
  quote = "'\\"",
  row.names = 1,
  na.strings = "",
  check.names = FALSE,
  strip.white = TRUE,
  stringsAsFactors = FALSE,
  ...
)
}
\arguments{
\item{file}{File path.}

\item{quote}{the set of quoting characters. To disable quoting
    altogether, use \code{quote = ""}.  See \code{\link{scan}} for the
    behaviour on quotes embedded in quotes.  Quoting is only considered
    for columns read as character, which is all of them unless
    \code{colClasses} is specified.}

\item{row.names}{a vector of row names.  This can be a vector giving
    the actual row names, or a single number giving the column of the
    table which contains the row names, or character string giving the
    name of the table column containing the row names.

    If there is a header and the first row contains one fewer field than
    the number of columns, the first column in the input is used for the
    row names.  Otherwise if \code{row.names} is missing, the rows are
    numbered.

    Using \code{row.names = NULL} forces row numbering. Missing or
    \code{NULL} \code{row.names} generate row names that are considered
    to be \sQuote{automatic} (and not preserved by \code{\link{as.matrix}}).
  }

\item{na.strings}{a character vector of strings which are to be
    interpreted as \code{\link{NA}} values.  Blank fields are also
    considered to be missing values in logical, integer, numeric and
    complex fields.  Note that the test happens \emph{after}
    white space is stripped from the input, so \code{na.strings}
    values may need their own white space stripped in advance.}

\item{check.names}{logical.  If \code{TRUE} then the names of the
    variables in the data frame are checked to ensure that they are
    syntactically valid variable names.  If necessary they are adjusted
    (by \code{\link{make.names}}) so that they are, and also to ensure
    that there are no duplicates.}

\item{strip.white}{logical. Used only when \code{sep} has
    been specified, and allows the stripping of leading and trailing
    white space from unquoted \code{character} fields (\code{numeric} fields
    are always stripped).  See \code{\link{scan}} for further details
    (including the exact meaning of \sQuote{white space}),
    remembering that the columns may include the row names.}

\item{stringsAsFactors}{logical: should character vectors be converted
    to factors?  Note that this is overridden by \code{as.is} and
    \code{colClasses}, both of which allow finer control.}

\item{...}{Further arguments to be passed to  \code{\link{read.delim}}.}
}
\value{
Data frame.
}
\description{
Delegates to \code{\link{read.delim}} where the separator is inferred from the file extension (CSV or TXT).
For CSV files the delimiter is set to \code{","} while for TXT file \code{"\t"} is used. Also sets
some default argument values as used by Core Hunter.
}
