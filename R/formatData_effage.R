"formatData.effage"<- function (id, time, status, covariates,effageData) 
{
    
    # effageData: 
    #    list containing effective age information:
    #      effageData$n effageData$subject
    #
    #    effageData$subject:  
    #       list containing:
    #               intercepts, slopes, lastperrep, 
    #               perrepind, effagebegin, effage
    
    
    
    
    # some controls
    if (effageData$n!=length(unique(id))) stop('data does not match')
    
    # ...........
    
    
    id.unique<-unique(id)
    n<-length(id.unique)
    
    
    covariates<-data.frame(covariates)
    
    kk<-table(id)
    
    
    tau<-sum(time[id==id.unique[1]])
    caltimes<-cumsum(c(0,time[id==id.unique[1] & status==1]))
    gaptimes<-c(0,time[id==id.unique[1] & status==1])
    censored<-time[id==id.unique[1] & status==0]
    covariate<-t(covariates[id==id.unique[1],])
    
    
    intercepts<-effageData$subject[[1]]$intercepts
    slopes<-effageData$subject[[1]]$slopes
    lastperrep<-effageData$subject[[1]]$lastperrep
    perrepind<-effageData$subject[[1]]$perrepind
    effagebegin<-effageData$subject[[1]]$effagebegin
    effage<-effageData$subject[[1]]$effage
    
    
    for (i in 2:n)
    {
        intercepts<-c(intercepts,effageData$subject[[i]]$intercepts)
        slopes<-c(slopes,effageData$subject[[i]]$slopes)
        lastperrep<-c(lastperrep,effageData$subject[[i]]$lastperrep)
        perrepind<-c(perrepind,effageData$subject[[i]]$perrepind)
        effagebegin<-c(effagebegin,effageData$subject[[i]]$effagebegin)
        effage<-c(effage,effageData$subject[[i]]$effage)
        
        tau<-c(tau,sum(time[id==id.unique[i]]))
        caltimes<-c(caltimes,cumsum(c(0,time[id==id.unique[i] & status==1])))
        gaptimes<-c(gaptimes,c(0,time[id==id.unique[i] & status==1]))
        censored<-c(censored,time[id==id.unique[i] & status==0])
        covariate<-cbind(covariate,t(covariates[id==id.unique[i],]))
    }
    
    ans<-list(n=n, k=kk, tau=tau, caltimes=caltimes, gaptimes=gaptimes, censored=censored, 
              intercepts=intercepts, slopes=slopes, lastperrep=lastperrep, perrepind = perrepind, 
              effagebegin=effagebegin, effage=effage, covariate=covariate)
    
    return(ans)
    
}
