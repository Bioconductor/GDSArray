#' gdsNodes
#' @rdname gdsNode
#' @param file the gds file name.
#' @description the function to return all available gds nodes inside
#'     the input gds file.
#' @return a vector of characters.
#' @export
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' gdsNodes(file)
#' file1 <- SeqArray::seqExampleFileName("gds")
#' gdsNodes(file1)
#' 
gdsNodes <- function(file)
{
    f <- openfn.gds(file)
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
}

