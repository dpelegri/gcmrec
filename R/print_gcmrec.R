"print.gcmrec" <- function (x, digits = max(options()$digits - 4, 3), ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" gmcrec failed.", x$fail, "\n")
        return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    
    label.Rho<-ifelse(x$rho.type==1,"Identity", "Alpha to k")
    
    if(x$rho.type==1)
    {
        x$coef<-c(1,x$coef)
    }
    coef <- x$coef[-1]
    
    se <- sqrt(diag(x$var))
    if(x$rho.type==1)
    {
        se<-c(NA,se)
    }
    
    if (is.null(x$Xi)) 
        se <- se[-1]
    if (!is.null(x$Xi)) 
    {
        n.Xi<-nrow(x$var)
        se <- se[-c(1,n.Xi)]
    }
    if (is.null(coef) | is.null(se)) 
        stop("Input is not valid")
    tmp <- cbind(coef, exp(coef), se, coef/se, signif(1 - pchisq((coef/se)^2, 
                                                                 1), digits - 1))
    
    if (x$se.type==1)
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
                                             "se(coef)", "z", "p"))
    if (x$se.type==2)
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
                                             "se(coef) Jacknife", "z", "p"))
    
    cat("\n")
    prmatrix(tmp)
    if (is.null(x$df)) 
        df <- sum(!is.na(coef))
    else df <- round(sum(x$df), 2)
    cat("\n")
    cat("  General class model parameter estimates", "\n")
    cat("    rho function: ", label.Rho, "\n")
    alpha <- x$coef[1]
    
    if (x$rho.type==2)
        se <- sqrt(x$var[1,1])
    if (x$rho.type==1)
        se <- "--"
    
    if (x$se.type==1)
        cat("      alpha (s.e.): ", alpha, " (", se, ")", "\n", sep = "")
    if (x$se.type==2)
        cat("      alpha (s.e. Jacknife): ", alpha, " (", se, ")", "\n", sep = "")
    
    if (!is.null(x$Xi))
    {
        se <- sqrt(x$var[n.Xi,n.Xi])  
        cat("    Frailty parameter, Xi (s.e. Jacknife): ", x$Xi," (", se, ")", "\n")
    }
    cat(" \n")
    
    if (is.null(x$Xi))
    {
        cat(paste("  log-likelihood=", round(x$loglik, 2)), "\n")
        cat("  n=", x$n, "\n")
        cat("  n times=", x$nk, "\n")
        cat("  number of iterations: ", x$kiter, " Newton-Raphson \n")
        cat("\n")
    }
    
    if (!is.null(x$Xi))
    {
        cat(paste("  Marginal log-likelihood=", round(x$loglik, 2)), "\n")
        cat("  n=", x$n, "\n")
        cat("  n times=", x$nk, "\n")
        cat("  number of iterations: ", x$kiter, " EM steps \n")
        cat("\n")
    }
    invisible()
}