"formatData" <- function (id, time, event, covariates, parameffage,cancer) 
{
    covariates <- data.frame(covariates)
    n <- length(unique(id))
    id.distinct <- unique(id)
    temp <- formatData.i(id[id == id.distinct[1]], time[id == 
                                                            id.distinct[1]], event[id == id.distinct[1]], covariates[id == 
                                                                                                                         id.distinct[1], ], parameffage,cancer[id == id.distinct[1]])
    
    
    k <- temp$k
    tau <- temp$tau
    caltimes <- temp$caltimes
    gaptimes <- temp$gaptimes
    censored <- temp$censored
    intercepts <- temp$intercepts
    slopes <- temp$slopes
    lastperrep <- temp$lastperrep
    perrepind <- temp$perrepind
    effagebegin <- temp$effagebegin
    effage <- temp$effage
    covariate <- temp$covariate
    for (i in 2:n) {
        temp <- formatData.i(id[id == id.distinct[i]], time[id == 
                                                                id.distinct[i]], event[id == id.distinct[i]], covariates[id == 
                                                                                                                             id.distinct[i], ], parameffage, cancer[id == id.distinct[i]])
        
        k <- c(k, temp$k)
        tau <- c(tau, temp$tau)
        caltimes <- c(caltimes, temp$caltimes)
        gaptimes <- c(gaptimes, temp$gaptimes)
        censored <- c(censored, temp$censored)
        intercepts <- c(intercepts, temp$intercepts)
        slopes <- c(slopes, temp$slopes)
        lastperrep <- c(lastperrep, temp$lastperrep)
        perrepind <- c(perrepind, temp$perrepind)
        effagebegin <- c(effagebegin, temp$effagebegin)
        effage <- c(effage, temp$effage)
        covariate <- cbind(covariate, temp$covariate)
    }
    ans <- list(n = n, k = k, tau = tau, caltimes = caltimes, 
                gaptimes = gaptimes, censored = censored, intercepts = intercepts, 
                slopes = slopes, lastperrep = lastperrep, perrepind = perrepind, 
                effagebegin = effagebegin, effage = effage, covariate = covariate)
    return(ans)
}