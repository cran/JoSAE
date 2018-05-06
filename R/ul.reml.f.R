ul.reml.f <-
function(
        samp.data
        , formula
        , samp.agg.X.pop
        , y.name
        , X.names
        , ...
        )
    {

        fit <- lme(formula
                   , data=samp.data, random=(~1|domain.id)
                   , weights=varFixed(~k.ij^2))####SQUARE K.IJ!


    #get approximate variances of variance parameters
var <- fit$apVar
par <- attr(var, "Pars")
vc <- exp(par)^2#ranef and resid var
sig.sq.v <- vc[1]
sig.sq.e <- vc[2]


    #Rao & M 7.2.8 - 7.2.11, Approx variances
a.i.dot <- unique(samp.data[c("domain.id","a.i.dot","n.i")])
a.i.dot$alpha.i <- sig.sq.e + a.i.dot$a.i.dot * sig.sq.v

I.vv <- 0.5 * sum(a.i.dot$a.i.dot^2 * a.i.dot$alpha.i^-2)
I.ee <- 0.5 * sum((a.i.dot$n.i-1) * (1/sig.sq.e^2) + a.i.dot$alpha.i^-2)
I.ve <- 0.5 * sum(a.i.dot$a.i.dot * a.i.dot$alpha.i^-2)

I.mat <- matrix(nrow=2,ncol=2)
I.mat[1,1] <- I.vv
I.mat[2,2] <- I.ee
I.mat[2,1] <- I.mat[1,2] <- I.ve

I.mat.inv <- solve(I.mat)

V.bar.vv <- I.mat.inv[1,1]
V.bar.ee <- I.mat.inv[2,2]
V.bar.ve <- I.mat.inv[1,2]


    #calculation of beta.hat and its covariance
        #gamma.i
    samp.agg.X.pop$gamma.i <- sig.sq.v/(sig.sq.v + sig.sq.e/samp.agg.X.pop$a.i.dot)


    list(  sig.sq.e=as.vector(sig.sq.e), V.bar.ee=V.bar.ee, sqrt.V.bar.ee=sqrt(V.bar.ee)
         , sig.sq.v=as.vector(sig.sq.v), V.bar.vv=V.bar.vv, sqrt.V.bar.vv=sqrt(V.bar.vv)
         , V.bar.ve=V.bar.ve
         , beta.hat=fixed.effects(fit), cov.beta.hat=vcov(fit)
         , sd.beta.hat=coef(summary(fit))[, "Std.Error"]
         , lme.fit=fit
         )
}
