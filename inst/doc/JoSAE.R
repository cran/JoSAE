### R code from vignette source 'JoSAE.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: JoSAE.Rnw:67-69 (eval = FALSE)
###################################################
## install.packages("JoSAE")
## 


###################################################
### code chunk number 2: JoSAE.Rnw:75-76
###################################################
library(JoSAE)


###################################################
### code chunk number 3: JoSAE.Rnw:80-81 (eval = FALSE)
###################################################
## ?JoSAE


###################################################
### code chunk number 4: JoSAE.Rnw:100-110
###################################################
	#mean auxiliary variables for the populations in the domains
data(JoSAE.domain.data)
	#data for the sampled elements
data(JoSAE.sample.data)

#plot(biomass.ha~mean.canopy.ht,JoSAE.sample.data)

library(lattice)

print(xyplot(biomass.ha ~ mean.canopy.ht | domain.ID, data = JoSAE.sample.data))


###################################################
### code chunk number 5: JoSAE.Rnw:131-134
###################################################
    #lme model
summary(fit.lme <- lme(biomass.ha ~ mean.canopy.ht, data=JoSAE.sample.data
                       , random=~1|domain.ID))


###################################################
### code chunk number 6: JoSAE.Rnw:143-146
###################################################
    #domain data need to have the same column names as sample data or vice versa
d.data <- JoSAE.domain.data
names(d.data)[3] <- "mean.canopy.ht"


###################################################
### code chunk number 7: JoSAE.Rnw:154-155
###################################################
result <- eblup.mse.f.wrap(domain.data = d.data, lme.obj = fit.lme)


###################################################
### code chunk number 8: JoSAE.Rnw:165-174
###################################################
library(xtable)
#tmp <- result[,c(1:2,6,4,9:11)]
tmp <- result[,c("domain.ID","N.i.domain","n.i.sample", "biomass.ha.sample.mean", "GREG"
                 , "EBLUP", "Synth")]
names(tmp)[4] <- c("sample.mean")
print(xtable(tmp, label = "tab:one", caption="Number of population and sampled elements as well as simple random sample, synthetic, GREG and EBLUP estimates of the mean above-ground forest biomass within 14 Norwegian municipalities."), include.rownames=FALSE)

tmp <- result[,c("domain.ID","n.i.sample","sample.se", "GREG.se", "EBLUP.se.1", "EBLUP.se.2")]#result[,c(1:2,6,23:26)]
print(xtable(tmp, label = "tab:two", caption="Number of population and sampled elements as well as standard errors of the simple random sample, GREG and EBLUP estimates of the mean above-ground forest biomass within 14 Norwegian municipalities."), include.rownames=FALSE)


###################################################
### code chunk number 9: JoSAE.Rnw:192-218
###################################################

tmp <- result[,c("biomass.ha.sample.mean", "Synth", "GREG", "EBLUP")]
    #actual plot
tmp1 <- barplot(t(as.matrix(tmp)), beside=T
            , names.arg=result$domain.ID
            , xlab="Municipalities"
                , ylab=expression(paste("Estimated biomass (Mg ", ha^{-1}, ")" ))
                , ylim=c(0,200))
    #print n.sample plots
text(tmp1[2,]+.5, y = 50, labels = result$n.i.sample,cex=1.5)
    #error bars
tmp2<- result[,c("sample.se", "sample.se", "GREG.se", "EBLUP.se.2")]#sample.se twice to fill the column, only used once.
tmp2[is.na(tmp2)] <- 0
        #plot error bars
            #sample mean
arrows(x0=tmp1[1,], y0=tmp[,1]+tmp2[,1], x1=tmp1[1,], y1 = tmp[,1]-tmp2[,1]
           , length = 0.01, angle = 90, code = 3)
            #GREG
arrows(x0=tmp1[3,], y0=tmp[,3]+tmp2[,3], x1=tmp1[3,], y1 = tmp[,3]-tmp2[,3]
           , length = 0.01, angle = 90, code = 3)
            #EBLUP
arrows(x0=tmp1[4,], y0=tmp[,4]+tmp2[,4], x1=tmp1[4,], y1 = tmp[,4]-tmp2[,4]
           , length = 0.01, angle = 90, code = 3)
    #legend
legend(13,200, fill=grey(c(.3, .6, .8, .9)), legend=c("SRS", "Synth", "GREG", "EBLUP"), bty="n")



###################################################
### code chunk number 10: JoSAE.Rnw:236-252
###################################################


data(landsat)

    #prepare the domain data - exclude "outlying" domain
landsat.domains <- unique(landsat[-33,c(1, 7:8,10)])
        #add a numeric domain ID
landsat.domains$domain.ID <- 1:nrow(landsat.domains)
        #change names to the names in the sample data
names(landsat.domains)[2:3] <- c("PixelsCorn", "PixelsSoybeans")

    #prepare the unit-level sample data
tmp <- landsat[-33,c(2:6, 10)]
        #add numeric domain ID
landsat.sample <- merge(landsat.domains[4:5], tmp, by="CountyName")



###################################################
### code chunk number 11: JoSAE.Rnw:258-265
###################################################
summary(landsat.lme <- lme(HACorn ~ PixelsCorn + PixelsSoybeans
                           , data=landsat.sample
                           , random=~1|domain.ID))
    #obtain EBLUP estimates and MSE
result <- eblup.mse.f.wrap(domain.data = landsat.domains
             , lme.obj = landsat.lme)



###################################################
### code chunk number 12: JoSAE.Rnw:268-278
###################################################

tmp <- result[,c("CountyName.domain","n.i.sample","EBLUP", "EBLUP.se.1", "EBLUP.se.2", "GREG.se")]
#tmp <- result[,c(5,9,13,28,29,27)]
names(tmp)[1:2] <- c("County.name","n_i")#, "EBLUP.mean", "EBLUP.se.1", "EBLUP.se.2", "GREG.se")
#names(tmp) <- c("County.name","n_i", "EBLUP.mean", "EBLUP.se.1", "EBLUP.se.2", "GREG.se")
print(xtable(tmp, label = "tab:landsat"
             , caption="EBLUP estimates of county means of hectares under corn and estimated standard errors of the EBLUP and GREG estimates.")
      , include.rownames=FALSE)




