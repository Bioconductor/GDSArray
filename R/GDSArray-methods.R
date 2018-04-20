
#' GDSArraySeed or GDSArray related methods, slot getters and setters.
#' 
#' @description \code{dim}, \code{dimnames}: dimension and dimnames of
#'     object contained in the GDS file.
#' @param x the \code{GDSArray} and \code{GDSArraySeed} objects.
#' @return \code{dim}: the integer vector of dimensions for
#'     \code{GDSArray} or \code{GDSArraySeed} objects.
#' @rdname GDSArray-methods
#' @exportMethod dim
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' ga <- GDSArray(file, "sample.annot/pop.group")
#' dim(ga)
#' dimnames(ga)
#' type(ga)
#' seed(ga)
#' dim(seed(ga))
#' gdsfile(ga)
setMethod("dim", "GDSArraySeed", function(x) x@dim)
setMethod("dim", "GDSArray", function(x) dim(seed(x)))

#' @rdname GDSArray-methods
#' @exportMethod dimnames
#' @return \code{dimnames}: the unnamed list of dimension names for
#'     \code{GDSArray} and \code{GDSArraySeed} objects.
setMethod("dimnames", "GDSArraySeed", function(x) x@dimnames)
setMethod("dimnames", "GDSArray", function(x) dimnames(seed(x)))

## setGeneric("seed", function(x) standardGeneric("seed"))
#' @rdname GDSArray-methods
#' @exportMethod seed
#' @description \code{seed}: the \code{GDSArraySeed} getter for
#'     \code{GDSArray} object.
#' @return \code{seed}: the \code{GDSArraySeed} of \code{GDSArray}
#'     object.
setMethod("seed", "GDSArray", function(x) x@seed)

#' @rdname GDSArray-methods
#' @exportMethod type
#' @description \code{type}: the data type of \code{GDSArray}.
#' @return \code{type}: the data type of \code{GDSArray} or
#'     \code{GDSArraySeed}.
setMethod("type", "GDSArray", function(x) type(seed(x)))

## setGeneric( "seed<-", function(x, value) standardGeneric("seed<-"),
##     signature="x" )
#' @rdname GDSArray-methods
#' @exportMethod "seed<-"
#' @description \code{seed<-}: the \code{GDSArraySeed} setter for
#'     \code{GDSArray} object.
#' @param value the new \code{GDSArraySeed} for the \code{GDSArray}
#'     object.
setReplaceMethod("seed", "GDSArray", function(x, value) {
    x@seed <- BiocGenerics:::replaceSlots(x, seed=value, check=FALSE)
})

#' @description \code{gdsfile}: on-disk location of GDS file
#'     represented by this object.
#' @param x GDSArray, GDSMatrix, GDSArraySeed, GDSFile or
#'     SummarizedExperiment object.
#' @return \code{gdsfile}: the character string for the gds file path.
#' @rdname GDSArray-methods
setGeneric("gdsfile", function(x) standardGeneric("gdsfile"), signature="x")

#' @rdname GDSArray-methods
#' @exportMethod gdsfile
setMethod("gdsfile", "GDSArraySeed", function(x) x@file)

#' @rdname GDSArray-methods
setMethod("gdsfile", "GDSArray", function(x) gdsfile(seed(x)))

#' @rdname GDSArray-methods
setMethod("gdsfile", "DelayedArray", function(x) gdsfile(seed(x)))

## #' @rdname GDSArray-methods
## #' @aliases GDSFile-method
## setMethod("gdsfile", "GDSFile", function(x) x@file)

#' "gdsfile<-"
#' @description \code{gdsfile<-}: the setter of the gds file path for
#'     `GDSArraySeed` and `GDSArray`.
#' @rdname GDSArray-methods
setGeneric(
    "gdsfile<-",
    function(x, value) standardGeneric("gdsfile<-"),
    signature="x")

#' @rdname GDSArray-methods
#' @exportMethod "gdsfile<-"
setReplaceMethod( "gdsfile", "GDSArraySeed", function(x, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    ## Set new path.
    BiocGenerics:::replaceSlots(x, file=value, check=FALSE)
})

#' @rdname GDSArray-methods
#' @exportMethod "gdsfile<-"
setReplaceMethod("gdsfile", "GDSArray", function(x, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    x@seed <- BiocGenerics:::replaceSlots(seed(x), file=value, check=FALSE)
    x
})
