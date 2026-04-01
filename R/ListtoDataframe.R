"List.to.Dataframe"<- function (data) 
{
    time.cen <- data$subject[[1]]$tau - sum(data$subject[[1]]$gaptime)
    time <- c(data$subject[[1]]$gaptime[-1], time.cen)
    event <- c(rep(1, data$subject[[1]]$k - 1), 0)
    id <- rep(1, data$subject[[1]]$k)
    dataEnd <- cbind(id = id, time = time, event = event, t(data$subject[[1]]$cov))
    i <- 2
    while (i <= data$n) {
        time.cen <- data$subject[[i]]$tau - sum(data$subject[[i]]$gaptime)
        time <- c(data$subject[[i]]$gaptime[-1], time.cen)
        event <- c(rep(1, data$subject[[i]]$k - 1), 0)
        id <- rep(i, data$subject[[i]]$k)
        temp <- cbind(id = id, time = time, event = event, 
                      t(data$subject[[i]]$cov))
        dataEnd <- rbind(dataEnd, temp)
        i <- i + 1
    }
    dataEnd <- data.frame(dataEnd)
    nvar <- nrow(data$subject[[1]]$cov)
    names.Cov <- paste("covar.", 1:nvar, sep = "")
    names(dataEnd) <- c(names(dataEnd)[1:3], names.Cov)
    dataEnd
}