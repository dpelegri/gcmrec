#' Fits General Class of Models for Recurrent Event Data
#' 
#' @param formula A formula object...
#' @param data A dataframe...
#' @return An object of class \code{gcmrec}
#' @export
#' @useDynLib gcmrec, .registration = TRUE
#' @importFrom Rcpp evalCpp
gcmrec <- function (formula, data, effageData = NULL, s, Frailty = FALSE, 
              alphaSeed, betaSeed, xiSeed, tol = 10^(-6), maxit = 100, 
              rhoFunc = "alpha to k", typeEffage = "perfect", maxXi = "Newton-Raphson", 
              se = "Information matrix", cancer = NULL) 
{
    rho.type <- charmatch(rhoFunc, c("Identity","alpha to k"), 
                          nomatch = 0)
    if (rho.type == 0) {
        stop("estimator must be 'alpha to k' or 'Identity' ")
    }
    effage.type <- charmatch(typeEffage, c("perfect", "minimal"))
    if (effage.type == 0) {
        stop("typeEffage must be perfect or minimal")
    }
    if (effage.type == 2) {
        effage.type <- 0
    }
    maxXi.type <- charmatch(maxXi, c("Newton-Raphson", "Brent"))
    if (maxXi.type == 0) {
        stop("maxXi must be Newton-Raphson or Brent")
    }
    se.type <- charmatch(se, c("Information matrix", "Jacknife"))
    if (se.type == 0) {
        stop("se must be Information matrix or Jacknife")
    }
    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) || 
        inherits(formula, "Survr")) {
        stop("formula.default(object): invalid formula")
    }
    m <- match.call(expand.dots = FALSE)
    
    m$s <- m$alphaSeed <- m$betaSeed <- m$xiSeed <- m$tol <- m$maxit <- m$rhoFunc <- m$Frailty <- m$effageData <- m$typeEffage <- m$maxXi <- m$se <- m$cancer <- m$... <- NULL
    Terms <- terms(formula, "strata")
    ord <- attr(Terms, "order")
    if (length(ord) & any(ord != 1)) 
        stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    n <- nrow(m)
    Y <- model.extract(m, "response")
    if (!is.Survr(Y)) 
        stop("Response must be a survival recurrent object")
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0) 
        rep(0, nrow(Y))
    else if (tt == 1) 
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    ll <- attr(Terms, "term.labels")
    mt <- attr(m, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, m, contrasts)
    if (ncol(X) == 1) 
        stop("model need some covariates")
    cov <- X[, -1]
    if (is.null(effageData)) {
        dataOK <- formatData(Y[, 1], Y[, 2], Y[, 3], cov, effage.type, 
                             cancer)
    }
    else {
        dataOK <- formatData.effage(Y[, 1], Y[, 2], Y[, 3], cov, 
                                    effageData)
    }
    if (max(dataOK$k) > 201) 
        stop("Procedures are implemented for number of \n recurrences less than 200")
    nvar <- ifelse(!is.null(ncol(cov)), ncol(cov), 1)
    if (missing(alphaSeed)) 
        alphaSeed <- 1
    if (missing(betaSeed)) 
        betaSeed <- rep(0, nvar)
    if (missing(xiSeed)) 
        xiSeed <- 1
    if (!Frailty) {
        if (rho.type == 2) {
            # ans <- .Fortran("newtraph", 
            #                 as.double(s), 
            #                 as.integer(dataOK$n), 
            #                 as.integer(nvar), 
            #                 as.integer(dataOK$k), 
            #                 as.integer(sum(dataOK$k)), 
            #                 as.double(dataOK$tau), 
            #                 as.double(dataOK$caltimes), 
            #                 as.double(dataOK$gaptimes), 
            #                 as.double(dataOK$censored), 
            #                 as.double(dataOK$intercepts), 
            #                 as.double(dataOK$slopes), 
            #                 as.double(dataOK$lastperrep), 
            #                 as.double(dataOK$perrepind), 
            #                 as.double(dataOK$effagebegin), 
            #                 as.double(dataOK$effage), 
            #                 as.double(dataOK$cov), 
            #                 as.double(alphaSeed), 
            #                 as.double(betaSeed), 
            #                 as.double(rep(1, dataOK$n)), 
            #                 as.double(offset), 
            #                 as.integer(rho.type), 
            #                 as.integer(rep(0, dataOK$n)), 
            #                 as.double(tol), 
            #                 as.integer(maxit), 
            #                 loglik = as.double(0), 
            #                 estim = as.double(rep(0, nvar + 1)), 
            #                 info = as.double(matrix(0, nvar + 1, nvar + 1)), 
            #                 search = as.integer(0), 
            #                 kiter = as.integer(0), 
            #                 PACKAGE = "gcmrec")
        
            
            ans <- newtraph(
                s            = as.double(s),
                n            = as.integer(dataOK$n),
                nvar         = as.integer(nvar),
                k            = as.integer(dataOK$k),
                nk           = as.integer(sum(dataOK$k)),
                tau          = as.double(dataOK$tau),
                caltimes     = as.double(dataOK$caltimes),
                gaptimes     = as.double(dataOK$gaptimes),
                censored     = as.double(dataOK$censored),
                intercepts   = as.double(dataOK$intercepts),
                slopes       = as.double(dataOK$slopes),
                lastperrep   = as.double(dataOK$lastperrep),
                perrepind    = as.double(dataOK$perrepind),
                effagebegin  = as.double(dataOK$effagebegin),
                effage       = as.double(dataOK$effage),
                cov          = as.matrix(dataOK$cov),         # Rcpp lo convierte a arma::mat
                alphaSeed    = as.double(alphaSeed),
                betaSeed     = as.double(betaSeed),
                Z            = as.double(rep(1, dataOK$n)),   # Pasamos el vector Z directamente
                offset       = as.double(offset),
                rhoFunc      = as.integer(rho.type),
                ns           = as.integer(rep(0, dataOK$n)),  # Se llenará dentro de C++
                tol          = as.double(tol),
                maxiter      = as.integer(maxit),
                loglik       = 0,                             # Pasamos un escalar (referencia en C++)
                estiEnd      = rep(0, nvar + 1),              # Vector resultado
                info         = matrix(0, nvar + 1, nvar + 1), # Matriz resultado
                searchBoth       = as.integer(1),                 # 1 para optimizar ambos (antes era 0 en tu ej)
                kkiter        = as.integer(0)                  # Contador de iteraciones
            )
        }
        if (rho.type == 1) {
            # ans <- .Fortran("nrbeta", as.double(s), as.integer(dataOK$n), 
            #                 as.integer(nvar), as.integer(dataOK$k), as.integer(sum(dataOK$k)), 
            #                 as.double(dataOK$tau), as.double(dataOK$caltimes), 
            #                 as.double(dataOK$gaptimes), as.double(dataOK$censored), 
            #                 as.double(dataOK$intercepts), as.double(dataOK$slopes), 
            #                 as.double(dataOK$lastperrep), as.double(dataOK$perrepind), 
            #                 as.double(dataOK$effagebegin), as.double(dataOK$effage), 
            #                 as.double(dataOK$cov), as.double(alphaSeed), 
            #                 as.double(betaSeed), as.double(rep(1, dataOK$n)), 
            #                 as.double(offset), as.integer(rho.type), as.integer(rep(0, 
            #                                                                         dataOK$n)), as.double(tol), as.integer(maxit), 
            #                 loglik = as.double(0), estim = as.double(rep(0, 
            #                                                              nvar)), info = as.double(matrix(0, nvar, 
            #                                                                                              nvar)), search = as.integer(0), kiter = as.integer(0), 
            #                 PACKAGE = "gcmrec")
            
            ans <- nrBeta(
                s            = as.double(s), 
                n            = as.integer(dataOK$n), 
                nvar         = as.integer(nvar), 
                k            = as.integer(dataOK$k), 
                nk           = as.integer(sum(dataOK$k)), 
                tau          = as.double(dataOK$tau), 
                caltimes     = as.double(dataOK$caltimes), 
                gaptimes     = as.double(dataOK$gaptimes), 
                censored     = as.double(dataOK$censored), 
                intercepts   = as.double(dataOK$intercepts), 
                slopes       = as.double(dataOK$slopes), 
                lastperrep   = as.double(dataOK$lastperrep), 
                perrepind    = as.double(dataOK$perrepind), 
                effagebegin  = as.double(dataOK$effagebegin), 
                effage       = as.double(dataOK$effage), 
                cov          = as.matrix(dataOK$cov), 
                alphaSeed    = as.double(alphaSeed), 
                betaSeed     = as.double(betaSeed), 
                Z            = as.double(rep(1, dataOK$n)), 
                offset       = as.double(offset), 
                rhoFunc     = as.integer(rho.type), 
                ns       = as.integer(rep(0, dataOK$n)), 
                tol          = as.double(tol), 
                maxiter        = as.integer(maxit),
                loglik      = 0,
                estiNew     = rep(0, nvar), 
                infobeta    = matrix(0, nvar, nvar),
                search      = as.integer(0),
                kiter       = as.integer(0)
            )
            
        }
        if (se.type == 2) {
            # Sustituimos el .Fortran por la función Rcpp
            seJack <- Jacknife(
                s           = as.double(s), 
                n           = as.integer(dataOK$n), 
                nvar        = as.integer(nvar), 
                k           = as.integer(dataOK$k), 
                nk          = as.integer(sum(dataOK$k)), 
                tau         = dataOK$tau, 
                caltimes    = dataOK$caltimes, 
                gaptimes    = dataOK$gaptimes, 
                censored    = dataOK$censored, 
                intercepts  = dataOK$intercepts, 
                slopes      = dataOK$slopes, 
                lastperrep  = dataOK$lastperrep, 
                perrepind   = dataOK$perrepind, 
                effagebegin = dataOK$effagebegin, 
                effage      = dataOK$effage, 
                cov         = as.matrix(dataOK$cov), 
                alphaSeed   = as.double(ans$estim[1]), 
                betaSeed    = as.double(ans$estim[-1]), 
                Z           = as.double(rep(1, dataOK$n)), 
                offset      = as.double(offset), 
                rhoFunc     = as.integer(rho.type), 
                ns          = as.integer(rep(0, dataOK$n)), 
                tol         = as.double(tol), 
                maxiter     = as.integer(maxit)
            )
            # Tu lógica original se mantiene intacta:
            JackEst <- apply(seJack, 2, mean)
            diff <- seJack - rep(1, dataOK$n) %*% matrix(JackEst, 1, nvar + 1)
            JackEstCov <- (dataOK$n - 1)/dataOK$n * t(diff) %*% diff
        }
        
        # if (se.type == 2) {
        #     jack <- .Fortran("jacknife", as.double(s), as.integer(dataOK$n), 
        #                      as.integer(nvar), as.integer(dataOK$k), as.integer(sum(dataOK$k)), 
        #                      as.double(dataOK$tau), as.double(dataOK$caltimes), 
        #                      as.double(dataOK$gaptimes), as.double(dataOK$censored), 
        #                      as.double(dataOK$intercepts), as.double(dataOK$slopes), 
        #                      as.double(dataOK$lastperrep), as.double(dataOK$perrepind), 
        #                      as.double(dataOK$effagebegin), as.double(dataOK$effage), 
        #                      as.double(dataOK$cov), as.double(ans$estim[1]), 
        #                      as.double(ans$estim[-1]), as.double(rep(1, dataOK$n)), 
        #                      as.double(offset), as.integer(rho.type), as.integer(rep(0, 
        #                                                                              dataOK$n)), as.double(tol), as.integer(maxit), 
        #                      estiEndJack = as.double(matrix(0, dataOK$n, nvar + 
        #                                                         1)), PACKAGE = "gcmrec")
        #     seJack <- matrix(jack$estiEndJack, dataOK$n, nvar + 
        #                          1)
        #     JackEst <- apply(seJack, 2, mean)
        #     diff <- seJack - rep(1, dataOK$n) %*% matrix(JackEst, 
        #                                                  1, nvar + 1)
        #     JackEstCov <- (dataOK$n - 1)/dataOK$n * t(diff) %*% 
        #         diff
        # }
        if (ans$search != 1) 
            warning("Algorithm did not converge in: ", maxit, 
                    " iterations")
        fit <- list(loglik = ans$loglik)
        fit$coef <- ans$estim
        if (se.type == 1 & rho.type == 2) 
            fit$var <- matrix(ans$info, nvar + 1, nvar + 1)
        if (se.type == 1 & rho.type == 1) 
            fit$var <- matrix(ans$info, nvar, nvar)
        if (se.type == 2) 
            fit$var <- JackEstCov
        fit$n <- ans[[2]]
        fit$nk <- ans[[5]]
        fit$kiter <- ans$kiter
        fit$Xi <- NULL
        fit$frailties <- rep(1, dataOK$n)
        fit$search <- ans$search
    }
    diseff <- sort(unique(dataOK$effage))
    ndiseff <- length(diseff)
    if (Frailty) {
        # ans <- .Fortran("estimwithfrailty", s = as.double(s),
        #                 n = as.integer(dataOK$n), nvar = as.integer(nvar),
        #                 k = as.integer(dataOK$k), nk = as.integer(sum(dataOK$k)),
        #                 as.double(dataOK$tau), as.double(dataOK$caltimes),
        #                 as.double(dataOK$gaptimes), as.double(dataOK$censored),
        #                 as.double(dataOK$intercepts), as.double(dataOK$slopes),
        #                 as.double(dataOK$lastperrep), as.double(dataOK$perrepind),
        #                 as.double(dataOK$effagebegin), as.double(dataOK$effage),
        #                 as.integer(ndiseff), as.double(diseff), as.double(dataOK$cov),
        #                 as.double(alphaSeed), as.double(betaSeed), as.double(xiSeed),
        #                 as.double(rep(1, dataOK$n)), as.double(offset), as.integer(rho.type),
        #                 as.double(tol), as.integer(maxit), as.integer(maxXi.type),
        #                 estim = as.double(rep(0, nvar + 2 + dataOK$n)),
        #                 control = as.integer(c(0, 0)),
        #                 loglik = as.double(0), PACKAGE = "gcmrec")
        
        # La función se llama directamente por su nombre en C++
        ans_cpp <- EstimWithFrailty(
            s           = s,
            n           = as.integer(dataOK$n),
            nvar        = as.integer(nvar),
            k           = as.integer(dataOK$k),
            nk          = as.integer(sum(dataOK$k)),
            tau         = dataOK$tau,
            caltimes    = dataOK$caltimes,
            gaptimes    = dataOK$gaptimes,
            censored    = dataOK$censored,
            intercepts  = dataOK$intercepts,
            slopes      = dataOK$slopes,
            lastperrep  = dataOK$lastperrep,
            perrepind   = dataOK$perrepind,
            effagebegin = dataOK$effagebegin,
            effage      = dataOK$effage,
            ndiseff     = as.integer(ndiseff),
            diseff      = diseff,
            cov         = as.matrix(dataOK$cov),
            alphaSeed   = alphaSeed,
            betaSeed    = betaSeed,
            xi          = xiSeed,
            Z           = rep(1, dataOK$n),
            offset      = offset,
            rhoFunc     = as.integer(rho.type),
            tol         = tol,
            maxiter     = as.integer(maxit),
            maxXi       = as.integer(maxXi.type)
        )
        
        # Para que el resto de tu código R siga funcionando sin cambios:
        # Rcpp devuelve una lista si la definimos así en el paso anterior.
        ans <- ans_cpp
        
        # ans <- EstimWithFrailty(
        #     alphaSeed = as.double(alphaSeed),
        #     betaSeed  = as.double(betaSeed),
        #     xiSeed    = as.double(xiSeed),
        #     Z         = as.double(rep(1, dataOK$n)), # Initial Fragility
        #     KK        = as.integer(dataOK$k),        # Events for subjet
        #     AA        = as.double(dataOK$tau),
        #     X         = as.matrix(dataOK$cov),       # Covariate matrix
        #     offset    = as.double(offset),
        #     rho_type  = as.integer(rho.type),
        #     tol       = as.double(tol),
        #     maxiter   = as.integer(maxit)
        # )
        
        
        if (se.type == 2) {
            seJack <- Jacknife2(
                s           = as.double(s), 
                n           = as.integer(dataOK$n), 
                nvar        = as.integer(nvar), 
                k           = as.integer(dataOK$k), 
                nk          = as.integer(sum(dataOK$k)), 
                tau         = dataOK$tau, 
                caltimes    = dataOK$caltimes, 
                gaptimes    = dataOK$gaptimes, 
                censored    = dataOK$censored, 
                intercepts  = dataOK$intercepts, 
                slopes      = dataOK$slopes, 
                lastperrep  = dataOK$lastperrep, 
                perrepind   = dataOK$perrepind, 
                effagebegin = dataOK$effagebegin, 
                effage      = dataOK$effage, 
                ndiseff     = as.integer(ndiseff), 
                diseff      = as.double(diseff), 
                cov         = as.matrix(dataOK$cov), 
                alphaSeed   = as.double(ans$estim[1]), 
                betaSeed    = as.double(ans$estim[2:(nvar + 1)]), 
                xi          = as.double(ans$estim[nvar + 2]), 
                Z           = as.double(rep(1, dataOK$n)), 
                offset      = as.double(offset), 
                rhoFunc     = as.integer(rho.type), 
                tol         = as.double(tol), 
                maxiter     = as.integer(maxit), 
                maxXi       = as.integer(maxXi.type)
            )
            # Tu lógica original se mantiene intacta:
            JackEst <- apply(seJack, 2, mean)
            diff <- seJack - rep(1, dataOK$n) %*% matrix(JackEst, 1, nvar + 2)
            JackEstCov <- (dataOK$n - 1)/dataOK$n * t(diff) %*% diff
        }
        
        # if (se.type == 2) {
        #     jack <- .Fortran("jacknife2", s = as.double(s), n = as.integer(dataOK$n), 
        #                      nvar = as.integer(nvar), k = as.integer(dataOK$k), 
        #                      nk = as.integer(sum(dataOK$k)), as.double(dataOK$tau), 
        #                      as.double(dataOK$caltimes), as.double(dataOK$gaptimes), 
        #                      as.double(dataOK$censored), as.double(dataOK$intercepts), 
        #                      as.double(dataOK$slopes), as.double(dataOK$lastperrep), 
        #                      as.double(dataOK$perrepind), as.double(dataOK$effagebegin), 
        #                      as.double(dataOK$effage), as.integer(ndiseff), 
        #                      as.double(diseff), as.double(dataOK$cov), as.double(ans$estim[1]), 
        #                      as.double(ans$estim[2:(nvar + 1)]), as.double(ans$estim[nvar + 
        #                                                                                  2]), as.double(rep(1, dataOK$n)), as.double(offset), 
        #                      as.integer(rho.type), as.double(tol), as.integer(maxit), 
        #                      as.integer(maxXi.type), estiEndJack = as.double(matrix(0, 
        #                                                                             dataOK$n, nvar + 2)), PACKAGE = "gcmrec")
        #     seJack <- matrix(jack$estiEndJack, dataOK$n, nvar + 
        #                          2)
        #     JackEst <- apply(seJack, 2, mean)
        #     diff <- seJack - rep(1, dataOK$n) %*% matrix(JackEst, 
        #                                                  1, nvar + 2)
        #     JackEstCov <- (dataOK$n - 1)/dataOK$n * t(diff) %*% 
        #         diff
        # }
        fit <- list(loglik = ans$loglik)
        fit$search <- ans$control[1]
        if (fit$search != 1) 
            warning("Algorithm did not converge after: ", maxit, 
                    " iterations")
        
        if (rho.type==1)  
            fit$coef <- ans$estim[2:(nvar + 1)]
        if (rho.type==2)  
            fit$coef <- ans$estim[1:(nvar + 1)]
        
        if (se.type == 1 & rho.type==1) 
            fit$var <- matrix(NA, nvar + 1, nvar + 1)
        if (se.type == 1 & rho.type==2) 
            fit$var <- matrix(NA, nvar + 2, nvar + 2)
        if (se.type == 2) 
            fit$var <- JackEstCov
        
        fit$n <- ans[[2]]
        fit$nk <- ans[[5]]
        fit$kiter <- ans$control[2]
        fit$Xi <- ans$estim[nvar + 2]
        fit$frailties <- ans$estim[(nvar + 3):(fit$n + nvar + 
                                                   2)]
    }
    fit$method <- "Newton-Raphson"
    fit$rho.type <- rho.type
    if (fit$search == 1) {
        
        EstLambSurv <- EstLambSurv(
            s = s, n = dataOK$n, nvar = nvar, k = dataOK$k, nk = sum(dataOK$k),
            tau = dataOK$tau, caltimes = dataOK$caltimes, gaptimes = dataOK$gaptimes,
            censored = dataOK$censored, intercepts = dataOK$intercepts, slopes = dataOK$slopes,
            lastperrep = dataOK$lastperrep, perrepind = dataOK$perrepind,
            effagebegin = dataOK$effagebegin, effage = dataOK$effage,
            ndiseff = ndiseff, diseff = diseff, cov = as.matrix(dataOK$cov),
            alpha = fit$coef[1], beta = fit$coef[-1], Z = fit$frailties,
            offset = offset, rhoFunc = rho.type
        )
        
        # EstLambSurv <- .Fortran("estlambsurv", as.double(s), 
        #                         as.integer(dataOK$n), as.integer(nvar), as.integer(dataOK$k), 
        #                         as.integer(sum(dataOK$k)), as.double(dataOK$tau), 
        #                         as.double(dataOK$caltimes), as.double(dataOK$gaptimes), 
        #                         as.double(dataOK$censored), as.double(dataOK$intercepts), 
        #                         as.double(dataOK$slopes), as.double(dataOK$lastperrep), 
        #                         as.double(dataOK$perrepind), as.double(dataOK$effagebegin), 
        #                         as.double(dataOK$effage), as.integer(ndiseff), as.double(diseff), 
        #                         as.double(dataOK$cov), as.double(fit$coef[1]), as.double(fit$coef[-1]), 
        #                         as.double(fit$frailties), as.double(offset), as.integer(rho.type), 
        #                         Lamb = as.double(rep(0, ndiseff)), DeltaLamb = as.double(rep(0, 
        #                                                                                      ndiseff)), Surv = as.double(rep(0, ndiseff)), 
        #                         PACKAGE = "gcmrec")
        fit$Lambda <- EstLambSurv$Lamb
        fit$Survival <- EstLambSurv$Surv
    }
    fit$se.type <- se.type
    fit$diseff <- diseff
    fit$terms <- Terms
    fit$call <- call
    if (length(fit$coef) > 2) {
        mod.X <- X[, -1]
    }
    else {
        mod.X <- matrix(X[, -1])
    }
    if (rho.type==2)
        names(fit$coef) <- c("alpha", colnames(X)[-1])
    if (rho.type==1)
        names(fit$coef) <- c(colnames(X)[-1])
    
    class(fit) <- "gcmrec"
    fit
}