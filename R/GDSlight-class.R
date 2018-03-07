
setClass(
    "GDSlight",
    slots=c(
        file = "character",
        current_path = "character"
    )
)

#' GDSlight constructor and methods. 
#' 
#' @name GDSlight
#' @description GDSlight is a light-weight class to represent a GDS file. It has the \code{$} method defined as a convenient input for \code{GDSArray} constructor.
#' @export
#' @rdname GDSlight-class
#' @aliases GDSlight-constructor
GDSlight <- function(file, current_path="")
{
    new("GDSlight", file = file, current_path = current_path)
}

.DollarNames.GDSlight <- function(x, pattern = "") {
    if (nzchar(x@current_path)) 
        pattern <- paste(x@current_path, pattern, sep="/")
    pattern <- sprintf("^(%s[^/]+).*", pattern)
    g <- sub(pattern, "\\1", grep(pattern, gdsNodes(x@file), value=TRUE))
    sub(".*/", "", g)
}

#' @exportMethod $
#' @rdname GDSlight-class
#' @aliases GDSlight-method
#' @details \code{$} takes \code{GDSlight} object as input, and output
#'     by updating \code{@current_path}.
setMethod("$", "GDSlight", function(x, name)
{
    if (nzchar(x@current_path)) 
        name <- paste(x@current_path, name, sep="/")
    x@current_path = name
    x
})
