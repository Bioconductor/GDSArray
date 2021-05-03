###-----------------------
### GDSArraySeed class 
###-----------------------

setClass("GDSArraySeed",
         contains="Array", # from DelayedArray: A virtual class with no slots
                           # to be extended by concrete subclasses with
                           # an array-like semantic.
         slots = c(
             gds = "gds.class",
             filename = "character",  # Absolute path to the gds file so the object
                                      # doesn't break when the user changes the working
                                      # directory (e.g. with setwd()).
             varname = "character",
             dim  = "integer",
             dimnames = "list"
         )
         )
## NOTE: gds(added), file->filename, name->varname, permute(removed), first_val(removed).

## if the file is open, no action internally
.reopen <- function(x) gdsfmt:::.reopen(x)

## slot accessors
.gds      <- function(x) x@gds
.filename <- function(x) x@filename
.varname  <- function(x) x@varname
## .dim      <- function(x) x@dim  ## REMOVE: dim() works directly if is one slot.

###--------------------------------------
### show method for GDSArraySeed object
###--------------------------------------

setMethod(
    "show", "GDSArraySeed",
    function(object) {
        cat("GDSArraySeed\n",
            "File: ", .filename(object), "\n",
            "Array node: ", .varname(object), "\n",
            "Dim: ", paste(dim(object), collapse=" x "), "\n",
            sep="")
    }
)

###-----------------
### extract_array() 
###-----------------
.extract_array_from_GDSArraySeed <- function(x, index)
{
    ## check
    stopifnot(is.list(index), !anyNA(index))
    ## NOTE: must support empty (NULL) or missing (integer(0)) index. -- YES!
    ##       must support duplicate indices. -- works automatically!
    ##       cannot contain NAs or non-positive values.
    ## https://bioconductor.org/packages/release/bioc/vignettes/DelayedArray/inst/doc/02-Implementing_a_backend.html#extract_array
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    # reopen the file if needed
    .reopen(.gds(x))
    on.exit(closefn.gds(.gds(x)))
    # read
    if (any(ans_dim == 0L))
    {     
        tp <- objdesp.gdsn(index.gdsn(.gds(x), .varname(x)))$type
        ans <- switch(as.character(tp),
            Raw=raw(), Integer=integer(), Logical=logical(),
            Real=double(), String=character(),
            stop("Unsupported data type: ", tp))
        dim(ans) <- ans_dim
    } else {
        nd <- index.gdsn(.gds(x), .varname(x))
        ## ans <- readex.gdsn(nd, index, .sparse=FALSE)  ## kept from SCArray
        ans <- readex.gdsn(nd, index, simplify = "none")
        if (!is.array(ans))  # ans must be an array 
            dim(ans) <- ans_dim
    }
    ans
}

#' GDSArray constructor and coercion methods.
#'
#' @rdname GDSArray-classes
#' @name extract_array
#' @exportMethod extract_array
#' @description \code{extract_array}: the function to extract data from
#'     a \code{GDS} file, by taking \code{GDSArraySeed} as input. This
#'     function is required by the \code{DelayedArray} for the seed
#'     contract.
#' @param x the GDSArraySeed object
#' @param index An unnamed list of subscripts as positive integer
#'     vectors, one vector per dimension in \code{x}. Empty and
#'     missing subscripts (represented by \code{integer(0)} and
#'     \code{NULL} list elements, respectively) are allowed. The
#'     subscripts can contain duplicated indices. They cannot contain
#'     NAs or non-positive values.
#' @aliases extract_array,GDSArraySeed-method

setMethod("extract_array", "GDSArraySeed", .extract_array_from_GDSArraySeed)

###---------------------------
### GDSArraySeed constructor
###---------------------------

#' @import methods
#' @import gdsfmt
#' @importFrom tools file_path_as_absolute

GDSArraySeed <- function(gds, varname)
{
    ## check gds
    stopifnot(inherits(gds, "gds.class"))
    ## check varname
    stopifnot(is.character(varname), length(varname)==1L, !is.na(varname))
    nd <- index.gdsn(gds, varname, silent=TRUE)
    if (is.null(nd))
        stop("No '", varname, "'")
    # check dimension
    dp <- objdesp.gdsn(nd)
    if (!dp$is.array)
        stop("'", varname, "' is not an array.")
    dm <- dp$dim
    # output
    new2("GDSArraySeed", gds=gds, filename=gds$filename, varname=varname,
        dim=dm, dimnames=vector("list", length(dm)))
}
## FIXME: may need to remove ".get_gdsnode_dim", and
## ".get_gdsnode_dimnames", ".get_gdsnode_type",
## ".get_gdsnode_first_val", ".read_gdsnode_sampleInCol" (check for
## SNP_ARRAY and SEQ_ARRAY and maybe add the @permute slot because in
## VariantExperiment it needs to be consistent with the first 2
## dimensions).

###--------------------------------
### GDSArray and GDSMatrix objects
###--------------------------------
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate GDSArray and GDSMatrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

#' @import S4Vectors
#' @import BiocGenerics
#' @import DelayedArray
#'
#' @rdname GDSArray-classes
#' @exportClass GDSArray
#' @aliases GDSArray-class matrixClass,GDSArray-method
setClass("GDSArray", contains="DelayedArray")

#' @rdname GDSArray-classes
#' @name GDSMatrix
#' @exportClass GDSMatrix
#' @aliases GDSMatrix-class 
setClass("GDSMatrix", contains=c("DelayedMatrix", "GDSArray"))

### For internal use only.
setMethod("matrixClass", "GDSArray", function(x) "GDSMatrix")

### Automatic coercion method from GDSArray to GDSMatrix (muted for
### higher dimensions) this function works only when GDSArray is
### 2-dimensional, otherwise, it fails.

#' @rdname GDSArray-classes
#' @name coerce
#' @exportMethod coerce
#' @aliases coerce,GDSArray,GDSMatrix-method
#'     coerce,GDSMatrix,GDSArray-method coerce,ANY,GDSMatrix-method

setAs("GDSArray", "GDSMatrix", function(from) new("GDSMatrix", from))    
setAs("GDSMatrix", "GDSArray", function(from) from)
setAs(
    "ANY", "GDSMatrix",
    function(from) as(as(from, "GDSArray"), "GDSMatrix"))

.validate_GDSArray <- function(x)
{
    if (!is(x@seed, "GDSArraySeed"))
        return(wmsg("'x@seed' must be a GDSArraySeed object"))
    TRUE
}

setValidity2("GDSArray", .validate_GDSArray)

###-----------------------
### GDSArray Constructor
###-----------------------

setMethod(
    "DelayedArray", "GDSArraySeed",
    function(seed) new_DelayedArray(seed, Class="GDSArray")
)

#' @description \code{GDSArray}: The function to convert a gds file
#'     into the GDSArray data structure.
#' @rdname GDSArray-classes
#' @aliases GDSArray-method
#' @param gdsfile Can be a GDSArraySeed, a character string of gds
#'     file name, or an "gds.class" R object.
#' @param varname A character string specifying the gds array node to
#'     be read into GDSArray.
#' @return \code{GDSArray} class object.
#' @export
#' @examples
#' fn <- gdsExampleFileName("snpgds") 
#' allnodes <- gdsnodes(fn)  ## print all available gds nodes in fn.
#' allnodes
#' GDSArray(fn, "genotype")
#' GDSArray(fn, "sample.annot/pop.group")
#'
#' fn1 <- gdsExampleFileName("seqgds")
#' allnodes1 <- gdsnodes(fn1)  ## print all available gds nodes in fn1. 
#' allnodes1
#' ## GDSArray(fn1, "genotype/data")
#' GDSArray(fn1, "variant.id")
#' GDSArray(fn1, "sample.annotation/family")
#' GDSArray(fn1, "annotation/format/DP/data")
#' GDSArray(fn1, "annotation/info/DP")

GDSArray <- function(gdsfile, varname)
{
    if (is(gdsfile, "GDSArraySeed")) {
        if (!missing(varname))
            stop(wmsg(
                "GDSArray() must be called with a single argument ",
                "when passed an GDSArraySeed object"))
        seed <- gdsfile
    } else {
        if (is.character(gdsfile)) {
            ## gdsfile <- openfn.gds(gdsfile, readonly = TRUE, allow.duplicate = FALSE)
            gdsfile <- openfn.gds(gdsfile)
            on.exit(closefn.gds(gdsfile))
        }
        seed <- GDSArraySeed(gdsfile, varname)
    }
    DelayedArray(seed)   ## does the automatic coercion to GDSMatrix if 2-dim.
}

###------------------------
### GDSArray example data
###------------------------

#' @rdname GDSArray-classes
#' @description \code{GDSArray} example data
#' @aliases GDSArray-data
#' @param type the type of gds file, available are "seqgds" for
#'     \code{SeqVarGDSClass} and "snpgds" for \code{SNPGDSFileClass}.
#' @export
#' @examples
#' gdsExampleFileName("snpgds")
#' gdsExampleFileName("seqgds")

gdsExampleFileName <- function(type = c("seqgds", "snpgds"))
{
    type <- match.arg(type)
    fn <- switch(type, snpgds = SNPRelate::snpgdsExampleFileName(),
                 seqgds = SeqArray::seqExampleFileName("gds"))
    f <- openfn.gds(fn)
    on.exit(closefn.gds(f))
    fmt <- get.attr.gdsn(f$root)$FileFormat
    msg <- switch(fmt, SNP_ARRAY = "This is a SNP GDS file",
                  SEQ_ARRAY = "This is a SeqArray GDS file")
    cat(msg, "\n", sep="")
    fn
}
