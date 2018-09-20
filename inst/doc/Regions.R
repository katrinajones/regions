## ---- warning=FALSE, message=FALSE---------------------------------------
library(regions)
data("alligator")

## ---- results = "asis",, echo=FALSE--------------------------------------
pander::pandoc.table(head(alligator[,1:5]))

## ------------------------------------------------------------------------
Xvar<-alligator[,1]
nvert<-length(Xvar)
Xvar[1:5]

## ------------------------------------------------------------------------
data<-alligator[,2:ncol(alligator)] #rest are dependent variables
data<-Missingval(data)#fill missing data

## ------------------------------------------------------------------------
data<-scale(data)

## ------------------------------------------------------------------------
pco.gower<-svdPCO(data, "gower")#PCO using svd
PCOscores<-pco.gower$scores[,1:ncol(pco.gower$scores)]

## ---- fig.cap="Variation along the column", fig.height=5, fig.width=5----
axesplot(PCOscores, 1, 2, Xvar)

## ------------------------------------------------------------------------
noregions<-3 #Set the maximum number of regions which will be calculated
regiondata<-compileregions(Xvar,PCOscores[,1:10],noregions)
pander::pandoc.table(regiondata[1:5,1:5])

## ---- fig.cap="Relationship of PCOs to regionalization", fig.height=5, fig.width=5----
plotpcoreg(eigenvals=pco.gower$eigen.val,nvert, namelabel="Alligator", regiondata, noregions)

## ------------------------------------------------------------------------
nopcos<-5
nopcos

## ---- fig.height=4, fig.width=6------------------------------------------

#bootstrapped with 100 itterations
pco.boot<-PCOcutoff(data,100, "gower")
nopcos<-pco.boot$sigpco#Select significant axes
nopcos


#Plot the eigenvalues
eigenplot(pco.boot$eigen.true, pco.boot$eigen.boot)

## ------------------------------------------------------------------------
nopcos<-length(which(pco.gower$eigen.val/sum(pco.gower$eigen.val)>0.05))#more than 5% of variance
nopcos

## ------------------------------------------------------------------------
nopcos<-PCOmax(regiondata, noregions, nvert)$pco.max
nopcos

## ------------------------------------------------------------------------
models<-modelselect(regiondata,noregions,nopcos)
pander::pandoc.table(models)

## ------------------------------------------------------------------------
support<-model_support(models,nvert, nopcos)
pander::pandoc.table(support$Model_support)

## ------------------------------------------------------------------------
rsq<-multvarrsq(Xvar,as.matrix(PCOscores[,1:nopcos]), support$Model_support)

## ---- fig.cap="Segmented regression model", fig.height=5, fig.width=6----
plotsegreg(Xvar,pcono=1, data=pco.gower$scores, modelsupport=support$Model_support)

## ---- fig.cap="Best regionalization model", fig.height=2, fig.width=6----
plot<-regionmodel(name="Alligator example", Xvar=Xvar, regiondata=support$Model_support)
  print(plot)

