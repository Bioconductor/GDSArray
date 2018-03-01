
#' GDSArraySeed or GDSArray related methods, slot getters and setters.
#' 
#' @description \code{dim}, \code{dimnames}: dimension and dimnames of
#'     object contained in the GDS file.
#' @rdname GDSArray-methods
#' @exportMethod dim
#' 
setMethod("dim", "GDSArraySeed", function(x) x@dim)

#' @rdname GDSArray-methods
#' @exportMethod dimnames
setMethod("dimnames", "GDSArraySeed", function(x) x@dimnames)

#' @description \code{gdsfile}: on-disk location of GDS file
#'     represented by this object.
#' @param x GDSArray, GDSMatrix, GDSArraySeed or SummarizedExperiment
#'     object.
#' @rdname GDSArray-methods
setGeneric("gdsfile", function(x) standardGeneric("gdsfile"))

#' @rdname GDSArray-methods
#' @exportMethod gdsfile
setMethod("gdsfile", "GDSArraySeed", function(x) x@file)

#' @rdname GDSArray-methods
setMethod("gdsfile", "GDSArray", function(x) gdsfile(seed(x)))

#' @rdname GDSArray-methods
setMethod("gdsfile", "DelayedArray", function(x) gdsfile(seed(x)))

#' "gdsfile<-"
#' @description \code{gdsfile<-}: the setter of the gds file path for `GDSArraySeed` and `GDSArray`.
#' @param value gds file path.
#' @rdname GDSArray-methods
setGeneric("gdsfile<-",
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
    BiocGenerics:::replaceSlots(seed(x), file=value, check=FALSE)
})
