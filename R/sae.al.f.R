sae.al.f <-
function(
             domain.id #vector of length m with domain IDs = stand_level.mean.xy$stand.ID #Sample districts
             , n.i # domain sample size
             , psi.i # Variance of mean estimate (s^2/n)
             , formula
             , data
             , b.i  # constant for heteroskedasticity
             , type="RE" # estimation procedure for estimating among-domain variance only REML (RE) is implemented here
             , verbose=F
             , iter = 100#max # of iterations
             , a.conv=0.001#convergency rel difference between old and new
        , ...
             )
    {
        
            #calculation of some derived parameters
    #number of sample domains
m.sam <- length(domain.id) 
    #use of formula
m <- model.frame(formula, data)
        #response cannot be empty for model.matrix
y.sam <- m[,1]
        # design matrix        
X.pop <- model.matrix(formula, m)
    # n fixed parameters -- p in Rao and Goerndt
k <- c(dim(X.pop))[2]

te <- NULL# initialize try error 


## REML estimator
if(type=="RE"){
sig.sq.v <- 0
    #est variance par
te <- try(#the estimates may not always converge
for(j in 1:iter){
    if(verbose) print(sig.sq.v)
        #estimate beta~ p. 116
            #with loop
    b.list.1 <- b.list.2 <- list()
    for(i in 1:m.sam){
        b.list.1[[i]] <- X.pop[i,] %*% t(X.pop[i,]) / c(psi.i[i] + sig.sq.v * b.i[i]^2)
        b.list.2[[i]] <- X.pop[i,] * y.sam[i] / c(psi.i[i] + sig.sq.v * b.i[i]^2)
    }
    #
    b.1 <- Reduce("+", b.list.1)
    b.2 <- Reduce("+", b.list.2)
    var.beta <- solve(b.1)
    beta.tilde <- var.beta %*% b.2 
    #
    if(verbose) print(beta.tilde)
        #Rao p. 119, R+M p. 128
    B <- diag(b.i^2)
        #p. 116
    V <- diag(psi.i + sig.sq.v * b.i^2)
        #p. 101
    P <- solve(V) - solve(V) %*% X.pop %*% solve(t(X.pop) %*% solve(V) %*% X.pop) %*% t(X.pop) %*% solve(V)
        #7.1.15
    I.R <- 0.5 * sum(diag(P %*% B %*% P %*% B))
    s.R <- - 0.5 * sum(diag(P %*% B)) + 0.5 * t(y.sam) %*% P %*% B %*% P %*% y.sam
        #iteratively change sig of ranef
    sig.sq.v.new <- c(sig.sq.v + solve(I.R) * s.R)
        #stopping iterations
    if(verbose) print(sig.sq.v.new)
    if(abs(sig.sq.v-sig.sq.v.new)/abs(sig.sq.v.new+0.001)<a.conv){if(verbose){cat("Converged at iteration",j)}
                                                                  break("break")}
    sig.sq.v <- sig.sq.v.new
    if(j == iter){warning("Algorithm did not converge.\n")}
    }#end: for 1:iter
    )#end: try
    #all estimators: if between area variance is negative, it is 0
if (sig.sq.v<0) sig.sq.v=0 #print(sig.sq.v)
beta <- beta.tilde
var.beta <- var.beta
}#end: if RE
##END: REML estimator


if(sig.sq.v==0){warning("Random domain effect not significant. MSE calculation without g3.i")}

if(!is.null(te)){warning(cat("Algorithm (",type,") did not converge. Using parameter and random variance estimate of last iteration step."))}


    # Compute shrinkage estimate
gamma.i=sig.sq.v * b.i^2 /(psi.i + sig.sq.v * b.i^2)
    #EBLUP (4)
FH <- gamma.i * y.sam + (1-gamma.i) * (X.pop %*% beta)

    # Random area effects
Areaeffcets <- round(gamma.i*(y.sam - X.pop %*% beta),2)

    #eblup synthetic estimate 
eblup.synth <- X.pop %*% beta


#transformed residuals similar to R+M p. 192.
eblup.resids <- y.sam - FH#these are actually the random effects
eblup.synth.resids <- y.sam - eblup.synth
trans.resid.v.i <-
    b.i^-1 * sqrt(sig.sq.v)^-1 * eblup.resids


ssq.y <- sum((y.sam - mean(y.sam))^2)
pseudo.R.sq.eblup <- 1-sum(eblup.resids^2)/ssq.y
R.sq.eblup.synth <- 1-sum(eblup.synth.resids^2)/ssq.y



# MSE estimate (7)
# CALCULATE g1 COMPONENT OF MSE
g1.i <- gamma.i*psi.i

# CALCULATE g2 COMPONENT OF MSE
    #7.1.8
g2.i <- (1-gamma.i)^2 * diag(X.pop %*% var.beta %*% t(X.pop))#b needs to be considered in varbeta


# Calculation of assymptotic var of ran eff
if(type %in% c("ML","RE")){
    #7.1.16 
    var.sig.sq.v <- 2 / sum(b.i^4/(sig.sq.v * b.i^2 + psi.i)^2)
}


# CALCULATE g3 COMPONENT OF MSE
    #7.1.21 or R+M 6.2.2
g3.i <- (psi.i^2 * b.i^4) * (psi.i + sig.sq.v * b.i^2)^-3 * var.sig.sq.v


    #7.1.23
g.star.3.i <- (b.i^4 * psi.i^2 / (psi.i + sig.sq.v * b.i^2)^4) *
    (y.sam - eblup.synth)^2 * var.sig.sq.v


if(type %in% c("RE","MS")){
        #no bias
    b.sig.sq.v <- b.sig.sq.v.detla.g1.i <- 0
    
}


    #MSE calculation
        #7.1.22
MSE.estimator.FH <- g1.i - b.sig.sq.v.detla.g1.i + g2.i + 2 * g3.i
        #7.1.24 - area-specific versions
MSE.area.spec.FH.1 <- g1.i - b.sig.sq.v.detla.g1.i + g2.i + 2*g.star.3.i
        #7.1.25 - area-specific versions
MSE.area.spec.FH.2 <- g1.i - b.sig.sq.v.detla.g1.i + g2.i + g3.i + g.star.3.i

        #conditional MSE - 6.2.35 in SAE v2 - can be negative
MSEp.conditional <- (gamma.i * psi.i) + (1 - gamma.i)^2 * ((y.sam - eblup.synth)^2 - (psi.i + sig.sq.v))
MSEp.conditional.gt0 <- MSEp.conditional
MSEp.conditional.gt0[MSEp.conditional.gt0<0] <- Inf


                #see lines above 6.4.10 and Isabell's email
rv.comp.phi.i <- g4.i <- NA
if(!is.na(n.i[1])){
    rv.comp.phi.i <- 2 * (sig.sq.v + psi.i)^-3 * sig.sq.v^2 * (2*psi.i^2/(n.i-1))
        #Hidiroglou and You 2016: comparison of unit and area level... always very close to 0!
    g4.i <- (4/(n.i-1)) * ((sig.sq.v^2 + psi.i^2)/(sig.sq.v + psi.i)^3)
}


    #mse of eblup synthetic - R+M 6.2.14
mse.eblup.synth.parameter.uncert <- diag(X.pop %*% var.beta %*% t(X.pop))
mse.eblup.synth <- sig.sq.v + b.sig.sq.v + mse.eblup.synth.parameter.uncert



        #if raneff is not significant, g1 is zero and g3, describing the variance of the random
        #effect is not needed
if(sig.sq.v==0){MSE.estimator.FH <- g2.i}

    # return results
Result.sam<- data.frame(domain.id=domain.id, y.sam, direct.SE=sqrt(psi.i)
                        , FH
                        , FH.SE=sqrt(MSE.estimator.FH)
                        , eblup.synth, eblup.synth.se=sqrt(mse.eblup.synth), eblup.synth.parameter.uncert.se=sqrt(mse.eblup.synth.parameter.uncert)
                        , bias.sig.sq.v=b.sig.sq.v
                        , RV.SE.phi.i=sqrt(MSE.estimator.FH + rv.comp.phi.i)
                        , FH.SE.area.spec.1=sqrt(MSE.area.spec.FH.1)
                        , FH.SE.area.spec.2=sqrt(MSE.area.spec.FH.2)
                        , FH.MSE.conditional=MSEp.conditional
                        , FH.SE.conditional=sqrt(MSEp.conditional.gt0)
                        , prd.ran.effects=Areaeffcets
                        , n.i, gamma.i, g1.i, g2.i
                        , g3.i, g.star.3.i, g1.bias=b.sig.sq.v.detla.g1.i
                        , eblup.resids, trans.resid.v.i, eblup.synth.resids
                        , b.i
                        )
result.area.level <- list(results=Result.sam, type=type
                          , sig.sq.v=sig.sq.v, beta.hat=beta, var.beta.hat=var.beta,
                          try.error.converge=te#is NULL if converged
                          , des.mat=X.pop, pseudo.R.sq.eblup=pseudo.R.sq.eblup
                          , R.sq.eblup.synth=R.sq.eblup.synth
                          )
}
