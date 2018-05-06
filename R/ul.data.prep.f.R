ul.data.prep.f <-
function(
            #input parameters
                #variance parameter 1 for equal variance
        k.ij
                #sample data
        , samp.data
                #population means of explanatory variables
        , population.data
                #fixed effects formula
        , formula
                #column name of domain ID
        , domain.col
                #sample ID column name (create using rownames if not existing)
        , sample.id.col
                #are sampling fractions negligible?
        , neg.sfrac
                #population means of explanatory variables for non-sampled units
                #only needed for non-negligible sampling fractions
        , pop.r
                #sum_i of k.ij^2 for non-sampled units
                #only needed for non-negligible sampling fractions
        , sum.i.k.ij.sq.r
                #number of units per area
                #only needed for non-negligible sampling fractions
        , N.i
        , ...
        )
    {
        samp <- samp.data
        pop <- population.data
        # formula
    mod <- model.frame(formula, data=samp)
            #response
    y.sam <- data.frame(mod[,1,F])

    y.name <- make.names(colnames(y.sam))
            # design matrix        
    X.sam <- data.frame(model.matrix(formula, samp))
                #y is needed to construct design matrix
    pop.for.desmat <- pop
    pop.for.desmat[,y.name] <- 1
                #transformed x-variables must be calc'ed before!
    X.pop <- data.frame(model.matrix(formula, pop.for.desmat))
    names(X.pop) <- paste0(names(X.pop), ".X.pop")
    X.pop$domain.id  <- pop[,domain.col]
    
        #replace parentheses around Intercept
    X.names <- make.names(colnames(X.sam))
     
    df <- data.frame(domain.id=samp[,domain.col], sample.id=samp[,sample.id.col]
                     , y.sam, X.sam
                     , k.ij=k.ij, a.ij=k.ij^-2)
     
    a.i.dot <- aggregate(cbind(a.i.dot=a.ij)~domain.id, df, sum)
         
    df2 <- merge(df, a.i.dot)

    ybar.i <- aggregate(df2[,y.name,F], by=list(domain.id=df2$domain.id), mean)
    names(ybar.i)[-1] <- paste0(names(ybar.i)[-1], ".ybar.i")

    xbar.i <- aggregate(df2[,X.names], by=list(domain.id=df2$domain.id), mean)
    names(xbar.i)[-1] <- paste0(names(xbar.i)[-1], ".xbar.i")
    
    ybar.i.a <- aggregate(df2$a.ij * df2[,y.name,F] / df2$a.i.dot,
                          by=list(domain.id=df2$domain.id), sum)

    names(ybar.i.a)[-1] <- paste0(names(ybar.i.a)[-1], ".ybar.i.a")

     
    xbar.i.a <- aggregate(df2$a.ij * df2[,X.names] / df2$a.i.dot
                          , by=list(domain.id=df2$domain.id), sum)
    names(xbar.i.a)[-1] <- paste0(names(xbar.i.a)[-1], ".xbar.i.a")
     
        #number of samples
    n.i <- aggregate(cbind(n.i=k.ij)~domain.id, df2, length)
     
        #merge aggregated values with sample data
    samp.data.1 <- merge(df2, n.i)
     
    samp.data.2 <- merge(samp.data.1, ybar.i.a)

    samp.data.3 <- merge(samp.data.2, ybar.i)
     
    samp.data.4 <- merge(samp.data.3, xbar.i.a)

    samp.data <- merge(samp.data.4, xbar.i)

    
        #get mean x.ia per area
    samp.data.agg   <- unique(samp.data[   c("domain.id"
                                             , paste0(y.name, ".ybar.i.a")
                                             , paste0(y.name, ".ybar.i")
                                             , paste0(X.names, ".xbar.i.a")
                                             , paste0(X.names, ".xbar.i")
                                             ,"a.i.dot", "n.i")])

        #combine aggregated sample data and population X values
    samp.agg.X.pop <- merge(samp.data.agg, X.pop)

        #only needed for non-negligible sampling fractions
    if(!neg.sfrac){
                #y is needed to construct design matrix
        pop.for.desmat <- pop.r
        pop.for.desmat[,y.name] <- 1
        X.pop.r <- data.frame(model.matrix(formula, pop.for.desmat))
        names(X.pop.r) <- paste0(names(X.pop.r), ".X.pop.r")
        X.pop.r$domain.id <- pop.r[,domain.col]
        X.pop.r$sum.i.k.ij.sq.r <- sum.i.k.ij.sq.r
        X.pop.r$N.i <- N.i
        samp.agg.X.pop <- merge(samp.agg.X.pop, X.pop.r)
        samp.agg.X.pop$f.i <- samp.agg.X.pop$n.i/samp.agg.X.pop$N.i
    }

        #number of domains and number of observations
    n <- nrow(samp.data)
    m <- length(unique(samp.data$domain.id))

    list(samp.data=samp.data
         , samp.agg.X.pop=samp.agg.X.pop
         , m=m, n=n
         , X.names=X.names, y.name=y.name, formula=formula
#         , X.names.noIntercept = X.names.noIntercept
         , neg.sfrac=neg.sfrac
         )
}
