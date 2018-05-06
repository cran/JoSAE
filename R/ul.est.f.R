ul.est.f <-
function(
    samp.data
    ,
    samp.agg.X.pop
    ,
    X.names
    ,
    y.name
    ,
    beta.hat
    ,
    cov.beta.hat
    ,
    sig.sq.e
    ,
    sig.sq.v
    ,
    V.bar.ee
    ,
    V.bar.vv
    ,
    V.bar.ve
    ,
    neg.sfrac
    ,
    resid=F
    , ...
    )
{
        #gamma.i
    samp.agg.X.pop$gamma.i <- sig.sq.v/(sig.sq.v + sig.sq.e/samp.agg.X.pop$a.i.dot)
        #tau.i and tau for clacluation of transformed residuals to check the model quality. R+M p. 187
    samp.agg.X.pop$tau.i <- 1 - sqrt(1 - samp.agg.X.pop$gamma.i)
        #gamma and tau for the case with k!=1. R+M p. 193
    gamma <- sig.sq.v/(sig.sq.v+sig.sq.e)
    tau <- 1 - sqrt(1 - gamma) 
    
        #eblup-synthetic
    samp.agg.X.pop$eblup.synth <- as.matrix(samp.agg.X.pop[,paste0(X.names,".X.pop")]) %*% beta.hat
     
        #survey regression estimator
    samp.agg.X.pop$sr.e.bar <-
        as.matrix(samp.agg.X.pop[,paste0(X.names,".X.pop")] - samp.agg.X.pop[,paste0(X.names, ".xbar.i")]) %*%
        beta.hat

    samp.agg.X.pop$svreg <- samp.agg.X.pop[,paste0(y.name, ".ybar.i")] + samp.agg.X.pop$sr.e.bar
     
        #prediction of random effects. R+M 7.1.6
    samp.agg.X.pop$ran.effs.v <- samp.agg.X.pop$gamma.i *
        (samp.agg.X.pop[,paste0(y.name,".ybar.i.a")] -
         as.matrix(samp.agg.X.pop[,paste0(X.names, ".xbar.i.a")]) %*% beta.hat)
     
        #eblup estimate
    samp.agg.X.pop$eblup <- samp.agg.X.pop$eblup.synth + samp.agg.X.pop$ran.effs.v

        #residuals for model diagnostics and checking of model assumptions
            #transformed residuals in the case of k=1
    tmp <- merge(samp.data, samp.agg.X.pop[,c("domain.id","tau.i","ran.effs.v")])
    trans.resid.k1 <-
        (tmp[,y.name]           - tmp$tau.i * tmp[,paste0(y.name,".ybar.i")]) -
        as.matrix(tmp[,X.names] - tmp$tau.i * tmp[,paste0(X.names,".xbar.i")]) %*% beta.hat
            #transformed residuals. R+M p. 192. These were specific for the variance equation in Militino 2006!
    trans.resid.u.ij <-
        tmp$k.ij^-1 * 
        (
            (tmp[,y.name]           - tau * tmp[,paste0(y.name,".ybar.i")]) -
            as.matrix(tmp[,X.names] - tau * tmp[,paste0(X.names,".xbar.i")]) %*% beta.hat
            )
    trans.resid.e.ij <-
        tmp$k.ij^-1 * sqrt(sig.sq.e)^-1 *
        (tmp[,y.name]  - as.matrix(tmp[,X.names]) %*% beta.hat - tmp$ran.effs.v)
            #raw eblup residuals
    prd.eblup <- as.matrix(tmp[,X.names]) %*% beta.hat + tmp$ran.effs.v
    eblup.resid <-
        (tmp[,y.name]  - prd.eblup)
            #raw eblup residuals
    eblup.synth.resid <-
        (tmp[,y.name]  - as.matrix(tmp[,X.names]) %*% beta.hat)
    #
    resids <- data.frame(tmp[c("sample.id", y.name)], prd.eblup, eblup.resid, eblup.synth.resid, trans.resid.k1
                         , trans.resid.u.ij, trans.resid.e.ij)
        #R.sq
    ssq.y <- sum((tmp[,y.name]-mean(tmp[,y.name]))^2)
    ssq.e <- sum(eblup.resid^2)
    ssq.e.synth <- sum(eblup.synth.resid^2)
    pseudo.R.sq.eblup <- 1-ssq.e/ssq.y
    R.sq.eblup.synth <- 1-ssq.e.synth/ssq.y

    rm(tmp)
        #if residuals should not be returned                 
    if(resid==F){resids <- NULL}
     
        #mse
            #Militino 2007
    samp.agg.X.pop$g1.i.eblup <- (1- samp.agg.X.pop[,"gamma.i"]) * sig.sq.v
    samp.agg.X.pop$g1.i.synth <- sig.sq.v
    samp.agg.X.pop$g1.i.svreg <- sig.sq.e/samp.agg.X.pop[,"a.i.dot"]
     
        #g2 - uncert of model parameters A.inv is cov(beta.hat)
    X.pop.minus.gamma.X.samp.mean <- as.matrix(samp.agg.X.pop[,paste0(X.names,".X.pop")] - samp.agg.X.pop$gamma.i * samp.agg.X.pop[,paste0(X.names, ".xbar.i.a")])
    samp.agg.X.pop$g2.i.eblup <- diag(X.pop.minus.gamma.X.samp.mean %*% cov.beta.hat %*% t(X.pop.minus.gamma.X.samp.mean))
    samp.agg.X.pop$g2.i.synth <- diag(as.matrix(samp.agg.X.pop[,paste0(X.names,".X.pop")]) %*%
                                      cov.beta.hat %*% t(as.matrix(samp.agg.X.pop[,paste0(X.names,".X.pop")])))
    samp.agg.X.pop$g2.i.svreg <- diag(as.matrix(samp.agg.X.pop[,paste0(X.names,".X.pop")] - samp.agg.X.pop[,paste0(X.names, ".xbar.i")]) %*% cov.beta.hat %*% t(as.matrix(samp.agg.X.pop[,paste0(X.names,".X.pop")] - samp.agg.X.pop[,paste0(X.names, ".xbar.i")])))
     
        #calc mse based on g1 and g2
    samp.agg.X.pop$se.eblup.naive <- sqrt(samp.agg.X.pop$g1.i.eblup + samp.agg.X.pop$g2.i.eblup)
    samp.agg.X.pop$se.eblup.synth <- sqrt(samp.agg.X.pop$g1.i.synth + samp.agg.X.pop$g2.i.synth)
    samp.agg.X.pop$se.svreg <- sqrt(samp.agg.X.pop$g1.i.svreg + samp.agg.X.pop$g2.i.svreg)#similar bt not the same as in Battese

     
        #g3 and g3star - only for eblup. Has 2 components.
            #g3.c1
    samp.agg.X.pop$g3.c1 <- samp.agg.X.pop[,"a.i.dot"]^-2 * (sig.sq.v + sig.sq.e/samp.agg.X.pop[,"a.i.dot"])^-3
     
            #h.sig.sq.v.e - note that v is multiplied with ee and e is multiplied with vv!
    h.sig.sq.v.e <-
        sig.sq.v^2 * V.bar.ee + sig.sq.e^2 * V.bar.vv - 2 * sig.sq.v * sig.sq.e * V.bar.ve
     
    samp.agg.X.pop$g3.i <- samp.agg.X.pop$g3.c1 * c(h.sig.sq.v.e)
     
     
         #g3star has 3 components. Makes mse area-specific.
            #g3star.c1 - note the minor difference to g3.c1!
    samp.agg.X.pop$g3star.c1 <- samp.agg.X.pop[,"a.i.dot"]^-2 * (sig.sq.v + sig.sq.e/samp.agg.X.pop[,"a.i.dot"])^-4
     
            #g3star.c3. These are residuals based on the sample x means
    samp.agg.X.pop$g3star.c3 <- (samp.agg.X.pop[,paste0(y.name, ".ybar.i.a")] -
                  as.matrix(samp.agg.X.pop[,paste0(X.names, ".xbar.i.a")]) %*% beta.hat)^2
     
    samp.agg.X.pop$g3star.i <- samp.agg.X.pop$g3star.c1 * c(h.sig.sq.v.e) * samp.agg.X.pop$g3star.c3
     
        #global se of eblup
    samp.agg.X.pop$se.eblup.global <-
        sqrt(samp.agg.X.pop$g1.i.eblup + samp.agg.X.pop$g2.i.eblup + 2 * samp.agg.X.pop$g3.i)
        #2 area specific versions of eblup
    samp.agg.X.pop$se.eblup.area1 <-
        sqrt(samp.agg.X.pop$g1.i.eblup + samp.agg.X.pop$g2.i.eblup + 2 * samp.agg.X.pop$g3star.i)
    samp.agg.X.pop$se.eblup.area2 <-
        sqrt(samp.agg.X.pop$g1.i.eblup + samp.agg.X.pop$g2.i.eblup + samp.agg.X.pop$g3.i +  samp.agg.X.pop$g3star.i)
     
     
        #direct estimators without fpc
    srs1 <- aggregate(cbind(srs=samp.data[,y.name]), by=list(domain.id=samp.data$domain.id),mean)
    srs2 <- aggregate(cbind(se.srs=samp.data[,y.name]), by=list(domain.id=samp.data$domain.id),
                      function(x){sqrt(var(x)/length(x))})
    srs <- merge(srs1,srs2)
    samp.agg.X.pop2 <- merge(samp.agg.X.pop, srs)
     
        #GREG - needs to be calculated on unit level. Point est of greg is same as svreg
    samp.data$gr.e <- samp.data[,y.name] - as.matrix(samp.data[,X.names]) %*% beta.hat
    greg1 <- aggregate(gr.e~domain.id, samp.data, mean)
    names(greg1)[-1] <- "greg.e.bar"
        #sum of squares. variance is then given by division by n.i-p
    greg2 <- aggregate(gr.e^2~domain.id, samp.data, sum)
    names(greg2)[-1] <- "greg.sse"
    greg <- merge(greg1, greg2)
    samp.agg.X.pop3 <- merge(samp.agg.X.pop2, greg)
    samp.agg.X.pop3$greg <- samp.agg.X.pop3$eblup.synth + samp.agg.X.pop3$greg.e.bar
            #greg se can be na or inf if too few sample plots in relation to number of parameters
            #therefore, switch off warnings temporarily
    options(warn=-1)
    samp.agg.X.pop3$se.greg <- sqrt(samp.agg.X.pop3$greg.sse/(samp.agg.X.pop3$n.i-length(beta.hat)))
    options(warn=0)

    
         #df of estimates
    est <- data.frame(domain.id=samp.agg.X.pop3$domain.id
                      , n.i=samp.agg.X.pop3$n.i
                      , srs=samp.agg.X.pop3$srs
                      #, greg=samp.agg.X.pop3$greg
                      , svreg=samp.agg.X.pop3$svreg
                      , eblup.synth=samp.agg.X.pop3$eblup.synth
                      , eblup=samp.agg.X.pop3$eblup
                      , gamma.i=samp.agg.X.pop3$gamma.i
                      , ran.effs.v=samp.agg.X.pop3$ran.effs.v)
     
        #df of SEs
    se <- data.frame(domain.id=samp.agg.X.pop3$domain.id
                     , se.srs=samp.agg.X.pop3$se.srs
                     #, se.greg=samp.agg.X.pop3$se.greg
                     , se.svreg=samp.agg.X.pop3$se.svreg
                     , se.eblup.synth=samp.agg.X.pop3$se.eblup.synth, se.eblup.naive=samp.agg.X.pop3$se.eblup.naive
                     , se.eblup.global=samp.agg.X.pop3$se.eblup.global, se.eblup.area1=samp.agg.X.pop3$se.eblup.area1
                     , se.eblup.area2=samp.agg.X.pop3$se.eblup.area2)
     
        #non neggligible sampling fractions
    est.se.fpc <- NULL
    if(!neg.sfrac){     
            #add X population values with sample data removed
        samp.agg.X.pop.fpc <- samp.agg.X.pop3
     
            #eblup R+M 7.1.24
        samp.agg.X.pop.fpc$resid.eblup.synth <-
            samp.agg.X.pop.fpc[,paste0(y.name, ".ybar.i")] - samp.agg.X.pop.fpc$eblup.synth
        samp.agg.X.pop.fpc$eblup.fpc <-
            samp.agg.X.pop.fpc$eblup.synth +
            (samp.agg.X.pop.fpc$f.i + (1 - samp.agg.X.pop.fpc$f.i) * samp.agg.X.pop.fpc$gamma.i) *
            samp.agg.X.pop.fpc$resid.eblup.synth
            
            #mse of eblup.nn.sfrac: g2tilde based on non-sampled population X values
        X.pop.r.minus.gamma.X.samp.mean <-
            as.matrix(samp.agg.X.pop.fpc[,paste0(X.names, ".X.pop.r")] -
                      samp.agg.X.pop.fpc$gamma.i * samp.agg.X.pop.fpc[,paste0(X.names, ".xbar.i.a")])

        samp.agg.X.pop.fpc$g2tilde.i.eblup <-
            diag(X.pop.r.minus.gamma.X.samp.mean %*% cov.beta.hat %*% t(X.pop.r.minus.gamma.X.samp.mean))
         
         
            #mse based on non-sampled population X values
                #global 
        samp.agg.X.pop.fpc$mse.y.ir.global <- samp.agg.X.pop.fpc$g1.i.eblup + samp.agg.X.pop.fpc$g2tilde.i.eblup + 2 * samp.agg.X.pop.fpc$g3.i
                #2 area specific versions of eblup                            
        samp.agg.X.pop.fpc$mse.y.ir.area1 <- (samp.agg.X.pop.fpc$g1.i.eblup + samp.agg.X.pop.fpc$g2tilde.i.eblup + 2 * samp.agg.X.pop.fpc$g3star.i)
        samp.agg.X.pop.fpc$mse.y.ir.area2 <- (samp.agg.X.pop.fpc$g1.i.eblup + samp.agg.X.pop.fpc$g2tilde.i.eblup + samp.agg.X.pop.fpc$g3.i +  samp.agg.X.pop.fpc$g3star.i)
       
                #proportion of residual variance
        samp.agg.X.pop.fpc$res.var.prop <- samp.agg.X.pop.fpc$sum.i.k.ij.sq.r * samp.agg.X.pop.fpc$N.i^-2
         
            #SEs
        samp.agg.X.pop.fpc$se.eblup.global.fpc <-
            sqrt((1-samp.agg.X.pop.fpc$f.i)^2 * samp.agg.X.pop.fpc$mse.y.ir.global +
                 samp.agg.X.pop.fpc$res.var.prop * sig.sq.e)
        samp.agg.X.pop.fpc$se.eblup.area1.fpc <-  sqrt((1-samp.agg.X.pop.fpc$f.i)^2 * samp.agg.X.pop.fpc$mse.y.ir.area1 + samp.agg.X.pop.fpc$res.var.prop * sig.sq.e)
        samp.agg.X.pop.fpc$se.eblup.area2.fpc <-  sqrt((1-samp.agg.X.pop.fpc$f.i)^2 * samp.agg.X.pop.fpc$mse.y.ir.area2 + samp.agg.X.pop.fpc$res.var.prop * sig.sq.e)
     
            #srs with fpc
        samp.agg.X.pop.fpc$fpc <-
            (samp.agg.X.pop.fpc$N.i - samp.agg.X.pop.fpc$n.i)/samp.agg.X.pop.fpc$N.i
        samp.agg.X.pop.fpc$se.srs.fpc <- sqrt(samp.agg.X.pop.fpc$fpc * samp.agg.X.pop.fpc$se.srs^2)
        samp.agg.X.pop.fpc$se.greg.fpc <- sqrt(samp.agg.X.pop.fpc$fpc * samp.agg.X.pop.fpc$se.greg^2)
        samp.agg.X.pop.fpc$se.svreg.fpc <- sqrt(samp.agg.X.pop.fpc$fpc * samp.agg.X.pop.fpc$se.svreg^2)
            #the synth est gets also a part of the resiaual error
        samp.agg.X.pop.fpc$se.eblup.synth.fpc <- sqrt(samp.agg.X.pop.fpc$se.eblup.synth^2 +
                                                samp.agg.X.pop.fpc$res.var.prop * sig.sq.e)
     
     
        #df of estimates
        est.se.fpc <- data.frame(domain.id=samp.agg.X.pop.fpc$domain.id
                                 , n.i=samp.agg.X.pop.fpc$n.i
                                 , N.i=samp.agg.X.pop.fpc$N.i
                                 , eblup=samp.agg.X.pop.fpc$eblup.fpc
                                 , res.var.prop=samp.agg.X.pop.fpc$res.var.prop
                                 , fpc=samp.agg.X.pop.fpc$fpc
                                 , se.srs=samp.agg.X.pop.fpc$se.srs.fpc
                                 #, se.greg=samp.agg.X.pop.fpc$se.greg.fpc
                                 , se.svreg=samp.agg.X.pop.fpc$se.svreg.fpc
                                 , se.eblup.synth.fpc=samp.agg.X.pop.fpc$se.eblup.synth.fpc
                                 , se.eblup.global=samp.agg.X.pop.fpc$se.eblup.global.fpc, se.eblup.area1=samp.agg.X.pop.fpc$se.eblup.area1.fpc
                                 , se.eblup.area2=samp.agg.X.pop.fpc$se.eblup.area2.fpc)
     
        try(est.se.fpc <- est.se.fpc[order(est.se.fpc$domain.id),])
     
    }#END: if(!neg.sfrac)
     
        #output 
    list(est=est, se=se, resids=resids 
         , est.se.fpc=est.se.fpc
         , gamma=gamma, sig.sq.v=sig.sq.v, sig.sq.e=sig.sq.e
         , pseudo.R.sq.eblup=pseudo.R.sq.eblup, R.sq.eblup.synth=R.sq.eblup.synth)
}
