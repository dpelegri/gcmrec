"formatData.i" <- function (id, time, event, covariates, parameffage, cancer = NULL) 
{
    covariates <- data.frame(covariates)
    k <- length(id)
    tau <- sum(time)
    caltimes <- cumsum(c(0, time[event == 1]))
    gaptimes <- c(0, time[event == 1])
    censored <- time[event == 0]
    lastperrep <- c(1)
    eff <- generlmi(parameffage)
    if (is.null(cancer)) {
        intercepts <- eff$intercept
        slopes <- eff$slope
        effagebegin <- eff$intercept
        effageend <- eff$intercept
        perrepind <- c(1)
    }
    else {
        if (cancer[1] == "CR") {
            intercepts <- 0
            slopes <- 1
            effagebegin <- intercepts
            effageend <- intercepts
            perrepind <- c(1)
        }
        
        #
        # This is in the case of not CR after first treatment. We must change A_0=0
        #
        if (cancer[1] != "CR") {
            intercepts <- 0
            slopes <- 1
            effagebegin <- intercepts
            effageend <- intercepts
            perrepind <- c(1)
        }
    }
    if (k >= 2) {
        for (i in 2:k) {
            eff <- generlmi(parameffage)
            intercepts <- c(intercepts, eff$intercept)
            slopes <- c(slopes, eff$slope)
            if (is.null(cancer)) {
                perrepind <- c(perrepind, eff$perrepind)
                if (eff$perrepind == 1) {
                    lastperrep <- c(lastperrep, i)
                }
                if (eff$perrepind == 0) {
                    lastperrep <- c(lastperrep, lastperrep[i - 
                                                               1])
                }
                effagebegin <- c(effagebegin, intercepts[i - 
                                                             1] + slopes[i] * (caltimes[i] - caltimes[lastperrep[i]]))
                effageend <- c(effageend, intercepts[i - 1] + 
                                   slopes[i] * (caltimes[i] - caltimes[lastperrep[i - 
                                                                                      1]]))
            }
            else {
                if (cancer[i] == "CR") {
                    lastperrep <- c(lastperrep, i)
                    perrepind <- c(perrepind, 1)
                }
                if (cancer[i] != "CR") {
                    lastperrep <- c(lastperrep, lastperrep[i - 
                                                               1])
                    perrepind <- c(perrepind, 0)
                }
                
                if(cancer[i] == "CR")
                {
                    effagebegin <- c(effagebegin, 0)
                    effageend <- c(effageend,gaptimes[i])
                }
                
                if (cancer[i] == "SD")
                { 
                    control<-effageend[i-1]
                    effagebegin <- c(effagebegin, control+gaptimes[i])
                    effageend <- c(effageend, control+gaptimes[i])
                }
                
                if (cancer[i] == "PR")
                { 
                    control<-effageend[i-1]
                    effagebegin <- c(effagebegin, control+0.5*gaptimes[i])
                    effageend <- c(effageend, control+0.5*gaptimes[i])
                }
                
                
                
                
            }
        }
    }
    covariate <- t(covariates)
    ans <- list(k = k, tau = tau, caltimes = caltimes, gaptimes = gaptimes, 
                censored = censored, intercepts = intercepts, slopes = slopes, 
                lastperrep = lastperrep, perrepind = perrepind, effagebegin = effagebegin, 
                effage = effageend, covariate = matrix(covariate, ncol(covariates), 
                                                       k))
    return(ans)
}
