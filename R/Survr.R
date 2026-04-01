Survr <- function (id,time,event) 
{
    
    if(length(unique(id))!=length(event[event==0]))
    {
        stop("Data doesn't match. Every subject must have a censored time")
    }
    
    if(length(unique(event))>2 | max(event)!=1 | min(event)!=0)
    {
        stop("event must be 0-1")
    }
    
    ans<-cbind(id,time,event)
    
    oldClass(ans) <- "Survr"
    invisible(ans)
    
}


#' @export
"is.Survr" <- function(x)
{
    inherits(x, "Survr")
}
