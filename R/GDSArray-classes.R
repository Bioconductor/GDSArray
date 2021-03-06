### -----------------------------------------------------
### GDSArraySeed class 
###
setClass(
    "GDSArraySeed",
    contains = "Array", ## from DelayedArray: A virtual class with no slots
                        ## to be extended by concrete subclasses with
                        ## an array-like semantic.
    slots = c(
        file="character",   # Absolute path to the gds file so the object
                            # doesn't break when the user changes the working
                            # directory (e.g. with setwd()).
        name="character",   # Name of the dataset in the gds file.
        dim = "integer",
        dimnames = "list",
        permute = "logical",
        first_val = "ANY"      ## remove this slot. 
    )
)

###
## show method for GDSArraySeed object
###
setMethod(
    "show", "GDSArraySeed",
    function(object) {
        cat("GDSArraySeed\n",
            "gds file: ", object@file, "\n",
            "array data: ", object@name, "\n",
            ## "dim: ", nrow(object), " x ", ncol(object), "\n",
            "dim: ", paste(dim(object), collapse=" x "), "\n",
            sep="")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array() ## must return array. 
###
#' @importFrom SeqArray seqSetFilter seqResetFilter seqGetData
#'
.extract_array_from_GDSArraySeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)){
        tpe <- class(x@first_val)
        ans <- get(tpe)(0)
        ## ans <- x@first_val[0]  ## integer(0) / character(0)
        dim(ans) <- ans_dim
    } else {
        ff <- .get_gds_fileFormat(x@file)
        if (ff == "SEQ_ARRAY") {
            f <- seqOpen(x@file)
            on.exit(seqClose(f))
            nodegroup <- names(dimnames(x))
            variant.sel <- sample.sel <- NULL
            if (any(grepl("variant", nodegroup))) {
                variant.sel <- index[[grep("variant", nodegroup)]]
            }
            if (any(grepl("sample", nodegroup))) {
                sample.sel <- index[[grep("sample", nodegroup)]]
            }
            seqSetFilter(f, variant.sel = variant.sel,
                         sample.sel = sample.sel, verbose = FALSE)
            ans <- seqGetData(f, x@name)
            ## if (grepl("format/DP", x@name))  ## for "annotation/format/DP" only
            ##     ans <- ans$data
            ## The above commented coded was due to the SeqArray(1.27.13) commit on 3/31/2020
            if (is.factor(ans)) ans <- unfactor(ans)
            seqResetFilter(f, verbose = FALSE)
            if (x@permute)
                ans <- aperm(ans)
        } else {
            f <- openfn.gds(x@file)
            on.exit(closefn.gds(f))
            if (x@permute)
                index <- rev(index)
            ans <- readex.gdsn(index.gdsn(f, x@name), index)
            if (x@permute)
                ans <- aperm(ans)
        }
        if (!is.array(ans))  ## 'ans' must be an array
            dim(ans) <- ans_dim
    }
    ans
}

#' GDSArray constructor and coercion methods.
#'
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
#' @rdname GDSArray-classes
setMethod("extract_array", "GDSArraySeed", .extract_array_from_GDSArraySeed)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArraySeed constructor
###

#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom SeqArray seqOpen seqClose seqSummary
#' @import methods
#' @import gdsfmt
#' @importFrom tools file_path_as_absolute

GDSArraySeed <- function(file, name=NA)
{
    if (!isSingleString(file))
        stop(wmsg(
            "'file' must be a single string specifying the path to ",
            "the gds file where the dataset is located."))
    if (!isSingleStringOrNA(name))
        stop("'name' must be a single string or NA")
    file <- file_path_as_absolute(file)

    dims <- .get_gdsnode_dim(file, node = name)
    dimnames <- .get_gdsnode_dimnames(file, node = name)

    if (is.null(dims)) {
        type <- .get_gdsnode_type(file, name)
        stop(wmsg("The gds node \"", name, "\" is type: ", type,
                  ", which is not valid for constructing GDSArray"))
    }
    if(any (dims == 0L))
        stop(wmsg("The dimension of gds node \"", name, "\" is: ",
                  paste(dims, collapse=" x "),
                  ",", "\n", "which could not construct GDSArray"))
    if (!identical(lengths(dimnames, use.names=FALSE), dims)) {
        stop(wmsg(
            "the lengths of dimnames ",
            "is not consistent with data dimensions."))
    }

    first_val <- .get_gdsnode_first_val(file, node = name)

    if (length(dims) == 1)
        permute = FALSE
    else 
        permute = !.read_gdsnode_sampleInCol(file, node = name)

    if (permute) {
        dims <- rev(dims)
        dimnames <- dimnames[rev(seq_len(length(dimnames)))]
    }
    new2(
        "GDSArraySeed", file=file,
        name=name,
        dim=dims,
        dimnames = dimnames,
        permute = permute,
        first_val = first_val
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArray and GDSMatrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate GDSArray and GDSMatrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

## #' @rawNamespace import(IRanges, except="concatenateObjects")
#' @import S4Vectors
#' @import BiocGenerics
#' @import DelayedArray
#'
#' @exportClass GDSArray
#' @aliases GDSArray-class matrixClass,GDSArray-method
#' @param file the gds file name.
#' @param name the gds array node to be read into GDSArraySeed / GDSArray. For
#'     \code{GDSArray}, the default value for \code{name} is the
#'     genotype data.
#' @return \code{GDSArray} class object.
#' @rdname GDSArray-classes
setClass("GDSArray", contains="DelayedArray")

#' @name GDSMatrix
#' @exportClass GDSMatrix
#' @aliases GDSMatrix-class 
#' @rdname GDSArray-classes
setClass("GDSMatrix", contains=c("DelayedMatrix", "GDSArray"))

### For internal use only.
setMethod("matrixClass", "GDSArray", function(x) "GDSMatrix")

### Automatic coercion method from GDSArray to GDSMatrix (muted for
### higher dimensions) this function works only when GDSArray is
### 2-dimensional, otherwise, it fails.

#' @name coerce
#' @exportMethod coerce
#' @aliases coerce,GDSArray,GDSMatrix-method
#'     coerce,GDSMatrix,GDSArray-method coerce,ANY,GDSMatrix-method
#' @rdname GDSArray-classes

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod(
    "DelayedArray", "GDSArraySeed",
    function(seed) new_DelayedArray(seed, Class="GDSArray")
)

#' @description \code{GDSArray}: The function to convert a gds file
#'     into the GDSArray data structure.
#' @export
#' @aliases GDSArray-method
#' @rdname GDSArray-classes
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' allnodes <- gdsnodes(file)  ## print all available gds nodes in file.
#' allnodes
#' ## GDSArray(file) #> deactivate temporarily 3/4/20
#' GDSArray(file, "sample.annot/pop.group")
#'
#' file1 <- SeqArray::seqExampleFileName("gds")
#' allnodes1 <- gdsnodes(file1)  ## print all available gds nodes in file1. 
#' allnodes1
#' ## GDSArray(file1) #> deactivate temporarily 3/4/20
#' GDSArray(file1, "variant.id")
#' GDSArray(file1, "sample.annotation/family")
#' GDSArray(file1, "annotation/format/DP")
#' GDSArray(file1, "annotation/info/DP")

GDSArray <- function(file, name=NA)
{
    if (is(file, "GDSArraySeed")) {
        if (!missing(name))
            stop(wmsg(
                "GDSArray() must be called with a single argument ",
                "when passed an GDSArraySeed object"))
        seed <- file
    } else {
        ## ff <- .get_gds_fileFormat(file)
        ## if (ff == "SNP_ARRAY") {
        ##     if (is.na(name)) name <- "genotype"
        ## } else if (ff == "SEQ_ARRAY") {
        ##     if (is.na(name)) name <- "genotype/data"
        ## }
        if (is.na(name))
            name <- "genotype"
        seed <- GDSArraySeed(file, name)
    }
    DelayedArray(seed)   ## does the automatic coercion to GDSMatrix if 2-dim.
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArray example data
###

#' @description \code{GDSArray} example data
#' @export
#' @aliases GDSArray-data
#' @rdname GDSArray-classes
#' @param pkg the package name, which is "GDSArray" by default. 
#' @examples
#' example("GDSArray")
#' 
example <- function(pkg="GDSArray")
{
    file <- SeqArray::seqExampleFileName("gds")
    gds <- GDSArray(file, "annotation/format/DP")
    gds
}
