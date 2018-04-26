
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

#' @rdname GDSArray-methods
#' @exportMethod dimnames
#' @return \code{dimnames}: the unnamed list of dimension names for
#'     \code{GDSArray} and \code{GDSArraySeed} objects.
setMethod("dimnames", "GDSArraySeed", function(x) x@dimnames)

## setGeneric("seed", function(x) standardGeneric("seed"))
#' @rdname GDSArray-methods
#' @exportMethod seed
#' @description \code{seed}: the \code{GDSArraySeed} getter for
#'     \code{GDSArray} object.
#' @return \code{seed}: the \code{GDSArraySeed} of \code{GDSArray}
#'     object.
setMethod("seed", "GDSArray", function(x) x@seed)

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
#' @param object GDSArray, GDSMatrix, GDSArraySeed, GDSFile or
#'     SummarizedExperiment object.
#' @return \code{gdsfile}: the character string for the gds file path.
#' @rdname GDSArray-methods
setGeneric("gdsfile", function(object) standardGeneric("gdsfile"),
           signature="object")

#' @rdname GDSArray-methods
#' @exportMethod gdsfile
setMethod("gdsfile", "GDSArraySeed", function(object) object@file)

#' @rdname GDSArray-methods
setMethod("gdsfile", "GDSArray", function(object) gdsfile(seed(object)))

#' @rdname GDSArray-methods
setMethod("gdsfile", "DelayedArray", function(object) gdsfile(seed(object)))

#' "gdsfile<-"
#' @description \code{gdsfile<-}: the setter of the gds file path for
#'     `GDSArraySeed` and `GDSArray`.
#' @rdname GDSArray-methods
setGeneric(
    "gdsfile<-",
    function(object, value) standardGeneric("gdsfile<-"),
    signature="object")

#' @rdname GDSArray-methods
#' @exportMethod "gdsfile<-"
setReplaceMethod( "gdsfile", "GDSArraySeed", function(object, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    ## Set new path.
    BiocGenerics:::replaceSlots(object, file=value, check=FALSE)
})

#' @rdname GDSArray-methods
#' @exportMethod "gdsfile<-"
setReplaceMethod("gdsfile", "GDSArray", function(object, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    object@seed <- BiocGenerics:::replaceSlots(seed(object),
                                               file=value,
                                               check=FALSE)
    object
})
