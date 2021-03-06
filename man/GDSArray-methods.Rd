% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GDSArray-methods.R
\name{dim,GDSArraySeed-method}
\alias{dim,GDSArraySeed-method}
\alias{dimnames,GDSArraySeed-method}
\alias{seed,GDSArray-method}
\alias{seed<-,GDSArray-method}
\alias{gdsfile}
\alias{gdsfile,GDSArraySeed-method}
\alias{gdsfile,GDSArray-method}
\alias{gdsfile,DelayedArray-method}
\alias{gdsfile<-}
\alias{gdsfile<-,GDSArraySeed-method}
\alias{gdsfile<-,GDSArray-method}
\title{GDSArraySeed or GDSArray related methods, slot getters and setters.}
\usage{
\S4method{dim}{GDSArraySeed}(x)

\S4method{dimnames}{GDSArraySeed}(x)

\S4method{seed}{GDSArray}(x)

\S4method{seed}{GDSArray}(x) <- value

gdsfile(object)

\S4method{gdsfile}{GDSArraySeed}(object)

\S4method{gdsfile}{GDSArray}(object)

\S4method{gdsfile}{DelayedArray}(object)

gdsfile(object) <- value

\S4method{gdsfile}{GDSArraySeed}(object) <- value

\S4method{gdsfile}{GDSArray}(object) <- value
}
\arguments{
\item{x}{the \code{GDSArray} and \code{GDSArraySeed} objects.}

\item{value}{the new \code{GDSArraySeed} for the \code{GDSArray}
object.}

\item{object}{GDSArray, GDSMatrix, GDSArraySeed, GDSFile or
SummarizedExperiment object.}
}
\value{
\code{dim}: the integer vector of dimensions for
    \code{GDSArray} or \code{GDSArraySeed} objects.

\code{dimnames}: the unnamed list of dimension names for
    \code{GDSArray} and \code{GDSArraySeed} objects.

\code{seed}: the \code{GDSArraySeed} of \code{GDSArray}
    object.

\code{gdsfile}: the character string for the gds file path.
}
\description{
\code{dim}, \code{dimnames}: dimension and dimnames of
    object contained in the GDS file.

\code{seed}: the \code{GDSArraySeed} getter for
    \code{GDSArray} object.

\code{seed<-}: the \code{GDSArraySeed} setter for
    \code{GDSArray} object.

\code{gdsfile}: on-disk location of GDS file
    represented by this object.

\code{gdsfile<-}: the setter of the gds file path for
    `GDSArraySeed` and `GDSArray`.
}
\examples{
file <- SNPRelate::snpgdsExampleFileName()
ga <- GDSArray(file, "sample.annot/pop.group")
dim(ga)
dimnames(ga)
type(ga)
seed(ga)
dim(seed(ga))
gdsfile(ga)
}
