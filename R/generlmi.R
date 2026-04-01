"generlmi"<- function(perrep)
{
    ll <- 0
    mm <- 1
    ii <- rbinom(1, 1, perrep)
    return(list(intercept = ll, slope = mm, perrepind = ii))
}
