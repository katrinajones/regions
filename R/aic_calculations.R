#' Select best region models from segmented regressions
#'
#' Select the best model for each number of regions by minimizing the residual sums of squares.
#'
#' @param regiondata Matrix created by \code{compileregions} containing region breaks and summed residuals
#' @param noregions Number of regions
#'@param nopcos Number of PCOs to be analyzed
#'@param startpco The pco to begin selection with
#' @return Data frame of best models for given regions
#' @export
#'
#'

modelselect<-function(regiondata, noregions, nopcos, startpco=1){

  if(any(colnames(regiondata)=="var 1")){
    colnames(regiondata)[which(colnames(regiondata)=="var 1")]<-"var.1"
  }
  pco.begin<-which(colnames(regiondata)=="var.1")
  pco.begin<-(pco.begin-1)+startpco
  if((nopcos-startpco)>0){
  sumRSS <- rowSums(regiondata[,pco.begin:(pco.begin+(nopcos-1))])
  }  else{sumRSS<-regiondata[,pco.begin]
  }
  regiondata2<-cbind(regiondata, sumRSS)
  regiondata2<-as.data.frame(regiondata2)
  models<-numeric()
  for (i in 1:noregions){
    allmodels<-subset(regiondata2, regions==i)#select only models with correct region no
    best<-allmodels[which(allmodels$sumRSS==min(allmodels$sumRSS)),]#select the lowest RSS
    best<-best[1,c(1:noregions,ncol(best))]#select columns
    models<-rbind(models, best)#fill in model table

    }
return(models)
}

#' Calculate AICc parameter for models
#'
#' Calculates corrected AIC values for least squares regression using residuals
#'
#' K is 2 (slope, intercept) multiplied by number of variables and number of regions
#'
#' No vert must be > no regions*3
#'
#'@param RSS total residual sums of squares, summed accross regions and variables
#'@param nPC Number of PCs analyzed
#'@param nvert Number of vertebrae
#'@param noregions Number of regions
#'
#' @return AICc
#' @export
#'
#'

AICcalc<-function(RSS, nPC, nvert, noregions){
  n=nPC*nvert #No of variables used
  var=RSS/n #Variance calculated ML way
  k=(2*noregions*nPC)+(noregions-1) #Based on slope, int and var estimates for each regression. Breakpoints not included as they are given.
  if(n<(k+2)) stop('ratio of variables to parameters too small. Reduce number of regions or increase variables')
  AIC=n*log(var)+(2*k) #+ n + n * log(2 * pi) #some cacls eg. r's AIC/AICc use this extra bit, shouldnt make a difference
  corr=(2*k*(k+1))/(n-k-1) #Correct for number of parameters and small sample
  AICc=AIC+corr #calculate AICc
  return(AICc)
}

#' Calculate relative support for each model
#'
#' Calculates AIC ratios and Akaike weights for each model to determine its strength of support.
#'
#' Also calculates a weighted average region score, which reflects the probability of each region model
#'
#' @param nvert No of vertebrae
#' @param nPC No of PCs used
#' @param models Data frame containing models and summed RSS
#'
#' @return Model_support Data frame containing each model and support values inlcuding AICc and Akaike weights
#' @return Region_score Mean weighted region score
#' @export
#'
#'

model_support<-function(models, nvert, nPC){
  AICc=numeric()
for (i in 1:nrow(models)){
  AICval<-AICcalc(models$sumRSS[i],nPC, nvert, models$regions[i]) #Calculate AIC score
  AICc<-rbind(AICc,AICval)
}
AICmin<-min(AICc)#AIC of best model
deltaAIC<-sapply(AICc, function(x) x-AICmin)#Calculate AIC difference
model_lik<-sapply(deltaAIC, function(x) exp(-0.5*x))#Likelihood of the model
tot_lik<-sum(model_lik)
Ak_weight<-sapply(model_lik, function(x) x/tot_lik)#Akaike Weights
AIC_models<-cbind(models, AICc, deltaAIC, model_lik, Ak_weight)
AIC_models<-AIC_models[order(AIC_models$AICc),]#Sort so best at top

weight_region<-AIC_models$regions*AIC_models$Ak_weight
Regions_score<-sum(weight_region)

return(list(Model_support=AIC_models,Region_score=Regions_score))
}


#' Run regionalization analysis
#'
#' Runs whole regionalization analysis, inlcuding running the models,
#' selecting the models and calculating the model support
#'
#' @param Xvar Vector, x variable
#' @param scores Matrix, PCO scores
#' @param noregions Number of region
#' @param nopcos Number of PCOs
#'
#' @return Model_support support for each best model
#' @return Region_scores Weighted region scores
#' @export
#'
#'
region_anal<-function(Xvar, scores, noregions, nopcos){

  nvert<-length(Xvar)
  regiondata<-compileregions(Xvar,scores,noregions)#set no regions and calculate RSS
  models<-modelselect(regiondata,noregions,nopcos)
  support<-model_support(models,nvert, nopcos)

  return(list(model_support=support$Model_support, score=support$Region_score))
}
