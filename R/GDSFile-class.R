###----------------
### GDSFile class
###----------------
#' @exportClass GDSFile
#' @rdname GDSFile-class
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
    nodes <- gdsNodes(object@file)
    nodes <- nodes[startsWith(nodes, object@current_path)]
    cat(
        "class: ", class(object), "\n",
        "file: ", object@file, "\n",
        "current node: ", object@current_path, "\n",
        "subnodes:\n  ", paste(nodes, collapse="\n  "), "\n",
        sep = ""
    )
})

###--------------
### constructor
###--------------

#' GDSFile constructor and methods. 
#' 
#' @name GDSFile
#' @description the \code{GDSFile} constructor.
#' @param file the GDS file path.
#' @param current_path the current path to the closest gds node.
#' @export
#' @rdname GDSFile-class
#' @aliases GDSFile-constructor
GDSFile <- function(file, current_path="")
{
    new("GDSFile", file = file, current_path = current_path)
}

###------------
### accessors
###------------

#' @exportMethod gdsfile
#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @description \code{gdsfile}: \code{file} slot getter and setter for
#'     \code{GDSFile} object.
#' @return \code{gdsfile}: the file path of corresponding GDS file.

setMethod("gdsfile", "GDSFile", function(x) x@file)

#' @exportMethod "gdsfile<-"
#' @rdname GDSFile-class
#' @aliases GDSFile-method GDSFile,gdsfile-method
#' @param value the new gds file path
setReplaceMethod("gdsfile", "GDSFile", function(x, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    BiocGenerics:::replaceSlots(x, file=value, check=FALSE)
})

###--------------------
### dollar completion
###--------------------

.DollarNames.GDSFile <- function(x, pattern = "") {
    nodes <- gdsNodes(x@file)
    nodes <- nodes[startsWith(nodes, x@current_path)]
    completions <- sub(sprintf("^%s/", x@current_path), "", nodes)
    sub("/.*", "", completions)
}

#' @exportMethod $
#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @param name the name of gds node
#' @return \code{$}: a \code{GDSFile} with updated \code{@current_path}, or
#'     \code{GDSArray} object if the \code{current_path} is a valid
#'     gds node.
#' 
setMethod("$", "GDSFile", function(x, name)
{
    if (nzchar(x@current_path)) {
        name <- paste(x@current_path, name, sep="/")
    }
    ## check if exist
    nodes <- gdsNodes(x)
    nodes <- nodes[startsWith(nodes, name)]
    pattern <- sprintf("(%s*)/.*$", name)

    if (name %in% sub(pattern, "\\1", nodes)) {
        x@current_path <- name
    } else {
        stop(wmsg("the gds path of '", name, "' does not exist"))
    }
    if (x@current_path %in% gdsNodes(x@file))
        GDSArray(gdsfile(x), x@current_path)
    else
        x
})

###------------
### methods
###------------

#' @exportMethod gdsNodes

setGeneric("gdsNodes",
           function(x) standGeneric(x),
           signature="x")

#' @name gdsNodes
#' @rdname GDSFile-class
#' @aliases GDSFile-method gdsNodes,ANY-method gdsNodes,GDSFile-method
#' @description \code{gdsNodes}: to get the available gds nodes from
#'     the \code{GDSFile} object or the file path with extension of
#'     ".gds".
#' @param x a \code{GDSFile} object or GDS file path.
#' @return \code{gdsNodes}: a character vector for the available gds
#'     nodes.
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' gdsNodes(file)
#' file1 <- SeqArray::seqExampleFileName("gds")
#' gdsNodes(file1)
#' gf <- GDSFile(file)
#' gdsNodes(gf)
#' gdsfile(gf)
setMethod("gdsNodes", "ANY", function(x)
{
    f <- openfn.gds(x)
    on.exit(closefn.gds(f))
    
    names.gdsn <- ls.gdsn(f)
    repeat {
        a <- lapply(names.gdsn, function(x) ls.gdsn(index.gdsn(f, x)))
        if (all(lengths(a)==0L)) {
            break
        } else {
            a[lengths(a)==0] <- ""
            n <- rep(names.gdsn, lengths(a))
            all.gdsn <- paste(n, unlist(a), sep="/")
            all.gdsn <- sub("/$", "", all.gdsn)
            names.gdsn <- all.gdsn
        }
    }
    names.gdsn
})

#' @exportMethod gdsNodes
setMethod("gdsNodes", "GDSFile", function(x) gdsNodes(gdsfile(x)) )
