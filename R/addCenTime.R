"addCenTime"<- function(datin, id=1, time=2, event=3)
{
    #
    # datin should have id,time,event variables
    # 
    
    ids <- unique(datin[,id])
    datout <- NULL
    for (i in ids)
    {
        thissubj <- datin[datin[,id] == i, , drop=FALSE]
        nrecs <- nrow(thissubj)
        if (thissubj[nrecs, event] == 1)
        {
            newrec <- thissubj[nrecs,]
            newrec[, time] <- 0
            newrec[, event] <- 0
            thissubj <- rbind(thissubj, newrec)
        }
        datout <- rbind(datout, thissubj)
    }
    datout
}