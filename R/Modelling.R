#' Calculate pairwise loglikelihood of models
#'
#' Test evolutionary hypotheses by comparing the distribution of their log-liklihood ratios
#'
#' Based on Boettiger et al. (2011). First likelihood ratios are generated based on data simulated under model 1, comparing
#' liklihood of model1 and model2. Then data are simulated under model2 and the liklihood ratios
#' are computed again. By comparing the distribution of these ratios with the 'true' ratio
#' the significance of the liklihood can be tested.
#'
#' Formatted to interact primarily with "Surface" models using 'surfaceSimulate' and 'startingModel'
#'
#' @param tree tree
#' @param sim1 function for simulating data under model1
#' @param sim2 function for simulating data under model2
#' @param lik1 function for fitting model1
#' @param lik2 function for fitting model2
#' @param reps no of simulations for each model
#' @param dat the 'real' data for calculating real log lik
#'
#' @return a plot comparing distributions to real ration
#' @export
#'
#'
pairwiseLoglik<-function(tree, sim1,sim2,lik1,lik2, reps, dat){
  ###Likelihood ratio test

  ##Calculate actual log lik
  olist<-convertTreeData(tree,dat)
  otree<-olist[[1]]; odata<-olist[[2]]
  real1<-lik1(otree, odata)
  real2<-lik2(otree, odata)
  actualloglik<--2*(real1[[1]]$fit$dat@loglik-real2[[1]]$fit$dat@loglik)

  ###Compare my model to BM
  likrat<-NULL
  for(i in 1:reps){
    newsim<-suppressWarnings(sim1()) #Generate null data
    olist<-convertTreeData(tree,newsim$data)
    otree<-olist[[1]]; odata<-olist[[2]]
    newbm<-lik1(otree, odata) #fit null model
    m1<-lik2(otree, odata) #fit test model
    #Likelihood ration
    likratk<--2*(newbm[[1]]$fit[[1]]@loglik-m1[[1]]$fit[[1]]@loglik)
    likrat<-c(likrat,likratk)

  }

  ###Simulate the final model
  likratmod1<-NULL
  for(i in 1:reps){
    mod1sim<-suppressWarnings(sim2()) #generate test simulation
    olist<-convertTreeData(tree,mod1sim$data)
    otree<-olist[[1]]; odata<-olist[[2]]
    mod1bm<-lik1(otree, odata)
    mod1m1<-lik2(otree, odata)
    #Likelihood ration
    likratmod1k<--2*(mod1bm[[1]]$fit[[1]]@loglik-mod1m1[[1]]$fit[[1]]@loglik)
    likratmod1<-c(likratmod1,likratmod1k)

  }

  #Plot
  d1<-density(likrat)
  d2<-density(likratmod1)
  rangedat<-range(c(d1$x,d2$x))
  plot(d1,xlim=rangedat, col="red", main="Distribution of likelihood ratios")
  polygon(d1, col=rgb(1,0,0,alpha=0.3))
  lines(d2, col="blue")
  polygon(d2, col=rgb(0,0,1,alpha=0.3))
  abline(v=actualloglik)

  ##Calculate p-value
  pval<-length(which(likrat>=actualloglik))/length(likrat)
  cutoff<-quantile(likrat,probs=c(0.95))
  power<-length(which(likratmod1>=cutoff))/length(likratmod1)

  return(list(loglik=actualloglik, pval=pval, cutoff=cutoff, power=power))
}

#' Calculate confidence interval on model based on simulation
#'
#' Estimate confidence intervals by simulating data based on true OU model. Based on Boettiger et al. (2011)
#'
#' @param tree tree
#' @param testsim simulation function
#' @param testmod  model for testing
#' @param reps no of simulations
#'@param dat no of simulations
#'
#' @return median and confidence interval on simulations
#' @export
#'
#'
confintModel<-function(tree,testsim, testmod,reps, dat){

  #Real values
  olist<-convertTreeData(tree,dat)
  otree<-olist[[1]]; odata<-olist[[2]]
  real<-testmod(otree, odata)
  real<-surfaceSummary(real)
  realdat<-c(real$n_regimes, real$alpha, real$phylhalflife, real$sigma_squared, c(real$theta))

  ###Simulate based on best model
  nregbest<-NULL
  param<-NULL
  theta<-NULL
  for(i in 1:reps){
    newsim<-suppressWarnings(testsim()) #Generate null data
    olist<-convertTreeData(tree,newsim$data)
    otree<-olist[[1]]; odata<-olist[[2]]
    bestnew<-testmod(otree, odata)
    bestnew<-surfaceSummary(bestnew)
    nregbest<-rbind(nregbest,bestnew$n_regimes)
    parami<-c(bestnew$alpha, bestnew$phylhalflife, bestnew$sigma_squared)
    names(parami)<-c("alpha", "halflife", "sigma-sq")
    param<-rbind(param, parami)
    theta<-rbind(theta, t(bestnew$theta))
   }

  nregbest<-cbind(nregbest,param, theta)
  names(realdat)<-colnames(nregbest)
  #Variation in models
  modvar<-apply(nregbest,2,median)
  modsd<-apply(nregbest,2,sd)
  modint<-apply(nregbest,2,function(x) quantile(x,probs=c(0.025,0.975)))

  return(list(true=realdat,median=modvar,confint=modint,sd=modsd))

}

#' Check for islands in surface analysis
#'
#' Rerun surface analysis with sampling of shifts. Select only credible models based on AIC threshold.
#'
#' @param tree tree
#' @param reps no of repeats
#' @param dat data
#' @param credthresh AIC cutoff for credible models. Default is 5.
#' @param samplethresh AIC cutoff for sampling shifts. Default is 5.
#'
#' @return "Real" model plus stats on variation within resampled models
#' @export
#'
#'
islandscheck<-function(tree, reps, dat, credthresh=5,samplethresh=5){

  ##Actual surface model
  real<-runSurface(tree,dat,verbose=F, only_best = T)
  realsum<-surfaceSummary(real$bwd)

  #Allow resampling of shifts to assess sensisitivity
  AICbest<-NULL
  nregbest<-NULL
  shiftsbest<-list()
  propmatch<-NULL
  for(i in 1:reps){
    best<-runSurface(tree,dat,verbose=F, sample_shifts = T, sample_threshold = samplethresh, only_best = T)
    bestnew<-surfaceSummary(best$bwd)
    AICbest<-c(AICbest, bestnew$aics[length(bestnew$aics)])
    nregbest<-rbind(nregbest,bestnew$n_regimes)
    shiftsbest[i]<-list(names(bestnew$shifts))
    propmatchk<-propRegMatch(real$bwd[[1]]$fit,best$bwd[[1]]$fit, internal=T)
    propmatch<-c(propmatch,propmatchk)
  }
  #Find credible models
  thresh<-min(AICbest)+credthresh
  credmod<-nregbest[which(AICbest<=thresh),]
  credshifts<-shiftsbest[which(AICbest<=thresh)]
  credmatch<-propmatch[which(AICbest<=thresh)]

  #Test appearance of shifts
  pshifts<-NULL
  for(a in 2:length(realsum$shifts)){
    loc<-names(realsum$shifts)[a]
    matchshift<-length(which(lapply(credshifts,function(x) any(x==loc))==TRUE))/nrow(credmod)
    pshifts<-c(pshifts,matchshift)
  }
  #Variation in models
  modvar<-apply(credmod,2,median)
  modrange<-apply(credmod,2,range)
  modAICmin<-min(AICbest)
  lowestmod<-shiftsbest[which(AICbest==modAICmin)]

  return(list(real=realsum$n_regimes,ncred=nrow(credmod),propmatch=credmatch, pshifts=pshifts, median=modvar, range=modrange, AICmin=modAICmin, bestshifts=lowestmod))

}
