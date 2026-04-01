"lines.gcmrec"<- function (x, type.plot="surv", ...) 
{
    dostep <- function(x, y) {
        n <- length(x)
        if (n > 2) {
            dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
            n2 <- sum(dupy)
            xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
            yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
            list(x = xrep, y = yrep)
        }
        else if (n == 1) 
            list(x = x, y = y)
        else list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
    }
    
    plot.type <- charmatch(type.plot, c("survival","hazard"), 
                           nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be hazard or survival")
    }
    
    if(plot.type==1)
    {
        y.lab <- c("Baseline survivor function")
        y <- x$Surv
        lines(dostep(x$diseff, y), ...)        
    } 
    
    if(plot.type==2)
    {
        y.lab <- c("Baseline hazard function")
        y <- x$Lam
        lines(dostep(x$diseff, y), ...)        
    }
    return(invisible())
}