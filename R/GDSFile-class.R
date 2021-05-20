###----------------
### GDSFile class
###----------------

#' @rdname GDSFile-class
#' @exportClass GDSFile
#' @description \code{GDSFile}: \code{GDSFile} is a light-weight class
#'     to represent a GDS file. It has the `$` completion method to
#'     complete any possible gds nodes. If the slot of `current_path`
#'     in `GDSFile` object represent a valid gds node, it will return
#'     the `GDSArray` of that node directly. Otherwise, it will return
#'     the `GDSFile` object with an updated `current_path`.
#' @aliases GDSFile-class

setClass(
    "GDSFile",
    slots=c(
        file = "character",
        current_path = "character"
    )
)

setMethod("show", "GDSFile", function(object) {
    nodes <- gdsnodes(object)
    nodes <- nodes[startsWith(nodes, object@current_path)]
    cat(
        "class: ", class(object), "\n",
        "file: ", object@file, "\n",
        "current node: ", object@current_path, "\n",
        "subnodes:\n  ", paste(nodes, collapse="\n  "), "\n",
        sep = ""
    )
})

###----------------------
### GDSFile constructor
###----------------------

#' GDSFile constructor and methods. 
#' 
#' @name GDSFile
#' @rdname GDSFile-class
#' @aliases GDSFile-constructor
#' @description \code{GDSFile}: the \code{GDSFile} class constructor.
#' @param file the GDS file path.
#' @param current_path the current path to the closest gds node.
#' @export

GDSFile <- function(file, current_path="")
{
    new("GDSFile", file = file, current_path = current_path)
}

###------------
### accessors
###------------

#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @description \code{gdsfile}: \code{file} slot getter for
#'     \code{GDSFile} object.
#' @param object \code{GDSFile} object.
#' @return \code{gdsfile}: the file path of corresponding
#'     \code{GDSfile} object.
#' @exportMethod gdsfile
#' @examples
#' fn <- gdsExampleFileName("seqgds")
#' gf <- GDSFile(fn)
#' gdsfile(gf)

setMethod("gdsfile", "GDSFile", function(object) object@file)

#' @rdname GDSFile-class
#' @aliases GDSFile-method GDSFile,gdsfile-method
#' @description \code{gdsfile<-}: \code{file} slot setter for
#'     \code{GDSFile} object.
#' @param value the new gds file path
#' @exportMethod "gdsfile<-"

setReplaceMethod("gdsfile", "GDSFile", function(object, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    BiocGenerics:::replaceSlots(object, file=value, check=FALSE)
})

###--------------------
### dollar completion
###--------------------
#' @importFrom utils .DollarNames
#' @export
.DollarNames.GDSFile <- function(x, pattern = "") {
    nodes <- gdsnodes(x)
    nodes <- nodes[startsWith(nodes, x@current_path)]
    completions <- sub(sprintf("^%s/", x@current_path), "", nodes)
    sub("/.*", "", completions)
}

#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @param name the name of gds node
#' @return \code{$}: a \code{GDSFile} with updated \code{@current_path}, or
#'     \code{GDSArray} object if the \code{current_path} is a valid
#'     gds node.
#' @exportMethod $

setMethod("$", "GDSFile", function(x, name)
{
    if (nzchar(x@current_path)) {
        name <- paste(x@current_path, name, sep="/")
    }
    ## check if exist
    nodes <- gdsnodes(x)
    nodes <- nodes[startsWith(nodes, name)]
    pattern <- sprintf("(%s*)/.*$", name)

    if (name %in% sub(pattern, "\\1", nodes)) {
        x@current_path <- name
    } else {
        stop(wmsg("the gds path of '", name, "' does not exist"))
    }
    if (x@current_path %in% gdsnodes(x@file))
        GDSArray(gdsfile(x), x@current_path)
    else
        x
})

###------------
### methods
###------------

#' @exportMethod gdsnodes

setGeneric("gdsnodes", function(x, node) standardGeneric("gdsnodes"), signature="x")

#' @name gdsnodes
#' @rdname GDSFile-class
#' @aliases GDSFile-method gdsnodes,ANY-method gdsnodes,GDSFile-method
#' @description \code{gdsnodes}: to get the available gds nodes from a
#'     gds file name or a \code{GDSFile} object. 
#' @param x a character string for the GDS file name or a \code{GDSFile} object.
#' @param node the node name of a gds file or \code{GDSFile} object. 
#' @return \code{gdsnodes}: a character vector of all available gds
#'     nodes within the related GDS file and the specified node.
#' @examples
#' fn <- gdsExampleFileName("seqgds")
#' gdsnodes(fn)
#' gdsnodes(fn, "annotation/info")
#' fn1 <- gdsExampleFileName("snpgds")
#' gdsnodes(fn1)
#' gdsnodes(fn1, "sample.annot")
#' gf <- GDSFile(fn)
#' gdsnodes(gf)
#' gdsnodes(gf, "genotype")
#' gdsfile(gf)

setMethod("gdsnodes", "ANY", function(x, node)
{
    f <- openfn.gds(x)
    on.exit(closefn.gds(f))
    if (missing(node))
        node <- ls.gdsn(f)
    repeat {
        a <- lapply(node, function(x) ls.gdsn(index.gdsn(f, x)))
        if (all(lengths(a)==0L)) {
            break
        } else {
            a[lengths(a)==0] <- ""
            ns <- rep(node, lengths(a))
            all.gdsn <- paste(ns, unlist(a), sep="/")
            all.gdsn <- sub("/$", "", all.gdsn)
            node <- all.gdsn
        }
    }
    node
})

#' @exportMethod gdsnodes
setMethod("gdsnodes", "GDSFile", function(x, node)
{
    gdsnodes(gdsfile(x), node)
})
