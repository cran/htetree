#
#' honest re-estimation and change the frame of object using estimation sample
#'
#' @param fit an \code{causalTree} object
#' @param x input training data
#' @param wt optional weights
#' @param treatment treatment variable
#' @param y outcome variable
#' @returns An object of class \code{rpart}.  See \code{\link{rpart.object}}.
#'
honest.est.causalTree <- function(fit, x, wt, treatment, y)
{
    frame <- fit$frame

    nc <- frame[, c("ncompete", "nsurrogate")]
    frame$index <- 1L + c(0L, cumsum((frame$var != "<leaf>") +
                                         nc[[1L]] + nc[[2L]]))[-(nrow(frame) + 1L)]
    frame$index[frame$var == "<leaf>"] <- 0L
    vnum <- match(rownames(fit$split), colnames(x))
    if (any(is.na(vnum)))
        stop("Tree has variables not found in new data")
    temp <- .Call("honest_estimate_causalTree",
                  as.integer(dim(x)),
                  as.integer(dim(frame)[1L]),
                  as.integer(dim(fit$splits)),
                  as.integer(if (is.null(fit$csplit)) rep(0L, 2L) else dim(fit$csplit)),
                  as.integer(row.names(frame)),
                  as.integer(unlist(frame[, c("ncompete", "nsurrogate", "index")])),
                  as.integer(frame[, "n"]),
                  as.double(frame[, "wt"]),
                  as.double(frame[, "dev"]),
                  as.double(frame[, "yval"]),
                  as.integer(vnum),
                  as.double(fit$splits),
                  as.integer(fit$csplit - 2L),
                  as.integer((fit$control)$usesurrogate),
                  as.double(x),
                  as.double(wt),
                  as.double(treatment),
                  as.double(y),
                  as.integer(is.na(x)))
    return (fit)
}
