sae.ul.f <-
function(...)
{
    dp <- NULL; rm(dp); #to satisfy r cmd check
    dp <- 
        ul.data.prep.f(...)
    #
    var.pars <- NULL; rm(var.pars); #to satisfy r cmd check
    var.pars <-
        ul.reml.f(
            samp.data=dp$samp.data,
            formula=dp$formula,
            samp.agg.X.pop=dp$samp.agg.X.pop,
            y.name=dp$y.name,
            X.names=dp$X.names
            )
    #
    est <-
        ul.est.f(
        samp.data = dp$samp.data
        ,samp.agg.X.pop=dp$samp.agg.X.pop
        ,X.names=dp$X.names
        ,y.name=dp$y.name
        ,beta.hat=var.pars$beta.hat
        ,cov.beta.hat=var.pars$cov.beta.hat
        ,sig.sq.e=var.pars$sig.sq.e
        ,sig.sq.v=var.pars$sig.sq.v
        ,V.bar.ee=var.pars$V.bar.ee
        ,V.bar.vv=var.pars$V.bar.vv
        ,V.bar.ve=var.pars$V.bar.ve
        ,neg.sfrac=dp$neg.sfrac
        ,resid=T
        )
    return(list(est=est, var.pars=var.pars, data=dp))
}
