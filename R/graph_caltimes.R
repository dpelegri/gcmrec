"graph.caltimes"<- function(data, var=NULL, effageData=NULL, width = 2, lines = TRUE, sortevents = TRUE, ...)
{
    
    if (!is.data.frame(data))
    {
        data<-List.to.Dataframe(data)
    }
    
    id.unique<-unique(data$id)
    nsubjs<-length(id.unique)
    
    
    if (!is.null(effageData))
    {
        perfect <- lapply(effageData, function(x) x[["perrepind"]])
    }
    else
    {
        perfect <- rep(1,nrow(data))
    }
    
    tauval<-NULL
    for(i in 1:nsubjs)
    {
        tauval[i]<-sum(data$time[data$id==id.unique[i]])
    }
    
    oldcex <- par("cex")
    par(cex = 1.25)
    
    
    plot(c(0,max(tauval)), c(1,nsubjs), type = "n", xlab = "Calendar Time", ylab = "Subject",...)
    for (i in 1:nsubjs)
    {
        thiscaltimes <- cumsum(data$time[data$id==id.unique[i] & data$event==1])
        thisperfect <- perfect[[i]]
        thistau <- tauval[i]
        
        if ((!is.null(var)) && (length(unique(var))<=5))
        {  
            thiscovvals <- var[data$id==id.unique[i]]
            colors <- thiscovvals[1] 
        }
        else 
        {
            colors<-1
        } 
        
        pchs <- ifelse(thisperfect == 0, 2, 1)
        points(thiscaltimes, rep(i, length(thiscaltimes)), col = colors, pch = pchs, lwd = width) 
        points(thistau, i, col="darkgreen",pch=4, lwd = width)
        if (lines) abline(h = i, lty = 3, lwd = .5)
    }
    par(cex = oldcex)
}