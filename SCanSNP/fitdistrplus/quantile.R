
#quantile function for fitdist objects
quantile.fitdist <- function(x, probs = seq(0.1, 0.9, by=0.1), ...)
{
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    myquantiles.fitdist(f = x, probs = probs, cens = FALSE)
}

#internal quantile function for fitdist
myquantiles.fitdist <- function(f, probs, cens)
{
    qdistname<-paste("q", f$distname, sep="")
    if (!exists(qdistname, mode="function"))
        stop(paste("The ", qdistname, " function must be defined"))

    # computation and print of quantiles using estimations of parameters
    para=c(as.list(f$estimate), as.list(f$fix.arg))
    quantiles <- do.call(qdistname, c(list(probs), as.list(para)))
    if (length(probs)>1)
        quantiles <- as.data.frame(t(quantiles))
    else
        quantiles <- as.data.frame(quantiles)
    colnames(quantiles) <- paste("p=", probs, sep="")
    rownames(quantiles) <- "estimate"

    reslist <- list(quantiles = quantiles, probs = probs)
    if(!cens)
        class(reslist) <- "quantile.fitdist"
    else
        class(reslist) <- "quantile.fitdistcens"

    reslist
}


#quantile function for fitdistcens objects
quantile.fitdistcens <- function(x, probs = seq(0.1, 0.9, by=0.1), ...)
{
    if (!inherits(x, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    myquantiles.fitdist(f = x, probs = probs, cens = TRUE)
}
