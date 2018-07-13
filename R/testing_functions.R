#' Calculate single region model
#'
#' Code for examining the residuals for a single region model with specified parameters
#'
#' Useful for examining the fit of a specified region model
#'
#' @param Xvar Independent variable
#' @param data Matrix of dependent variables
#' @param noregions Number of regions
#' @param breaks Vector of breakpoints
#'
#' @importFrom stats lm
#'
#' @return RSS.results residuals sums of squares for each PC and each region
#' @export
#'
single.model<-function(Xvar,data,noregions,breaks){

  RSS.results<-matrix(0,nrow=noregions,ncol=ncol(data))

  for (i in 1: ncol(data)){
  if (noregions==1){
    fit<-lm(data[,i]~Xvar)
    RSS<-sum(fit$residuals^2)
    RSS.results[,i]<-RSS
  }
  if (noregions==2){
    fit1<-lm(data[1:which(Xvar == breaks[1]),i]~Xvar[1:which(Xvar == breaks[1])])
    RSS1<-sum(fit1$residuals^2)
    fit2<-lm(data[(which(Xvar == breaks[1])+1):length(Xvar),i]~Xvar[(which(Xvar == breaks[1])+1):length(Xvar)])
    RSS2<-sum(fit2$residuals^2)
    RSS.results[1,i]<-RSS1
    RSS.results[2,i]<-RSS2

  }
  if (noregions==3){
    fit1<-lm(data[1:which(Xvar == breaks[1]),i]~Xvar[1:which(Xvar == breaks[1])])
    RSS1<-sum(fit1$residuals^2)
    fit2<-lm(data[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2]),i]~Xvar[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2])])
    RSS2<-sum(fit2$residuals^2)
    fit3<-lm(data[(which(Xvar == breaks[2])+1):length(Xvar),i]~Xvar[(which(Xvar == breaks[2])+1):length(Xvar)])
    RSS3<-sum(fit3$residuals^2)
    RSS.results[1,i]<-RSS1
    RSS.results[2,i]<-RSS2
    RSS.results[3,i]<-RSS3
  }
  if (noregions==4){
    fit1<-lm(data[1:which(Xvar == breaks[1]),i]~Xvar[1:which(Xvar == breaks[1])])
    RSS1<-sum(fit1$residuals^2)
    fit2<-lm(data[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2]),i]~Xvar[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2])])
    RSS2<-sum(fit2$residuals^2)
    fit3<-lm(data[(which(Xvar == breaks[2])+1):which(Xvar == breaks[3]),i]~Xvar[(which(Xvar == breaks[2])+1):which(Xvar == breaks[3])])
    RSS3<-sum(fit3$residuals^2)
    fit4<-lm(data[(which(Xvar == breaks[3])+1):length(Xvar),i]~Xvar[(which(Xvar == breaks[3])+1):length(Xvar)])
    RSS4<-sum(fit4$residuals^2)
    RSS.results[1,i]<-RSS1
    RSS.results[2,i]<-RSS2
    RSS.results[3,i]<-RSS3
    RSS.results[4,i]<-RSS4
  }
  if (noregions==5){
    fit1<-lm(data[1:which(Xvar == breaks[1]),i]~Xvar[1:which(Xvar == breaks[1])])
    RSS1<-sum(fit1$residuals^2)
    fit2<-lm(data[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2]),i]~Xvar[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2])])
    RSS2<-sum(fit2$residuals^2)
    fit3<-lm(data[(which(Xvar == breaks[2])+1):which(Xvar == breaks[3]),i]~Xvar[(which(Xvar == breaks[2])+1):which(Xvar == breaks[3])])
    RSS3<-sum(fit3$residuals^2)
    fit4<-lm(data[(which(Xvar == breaks[3])+1):which(Xvar == breaks[4]),i]~Xvar[(which(Xvar == breaks[3])+1):which(Xvar == breaks[4])])
    RSS4<-sum(fit4$residuals^2)
    fit5<-lm(data[(which(Xvar == breaks[4])+1):length(Xvar),i]~Xvar[(which(Xvar == breaks[4])+1):length(Xvar)])
    RSS5<-sum(fit5$residuals^2)
    RSS.results[1,i]<-RSS1
    RSS.results[2,i]<-RSS2
    RSS.results[3,i]<-RSS3
    RSS.results[4,i]<-RSS4
    RSS.results[5,i]<-RSS5
  }
  if (noregions==6){
    fit1<-lm(data[1:which(Xvar == breaks[1]),i]~Xvar[1:which(Xvar == breaks[1])])
    RSS1<-sum(fit1$residuals^2)
    fit2<-lm(data[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2]),i]~Xvar[(which(Xvar == breaks[1])+1):which(Xvar == breaks[2])])
    RSS2<-sum(fit2$residuals^2)
    fit3<-lm(data[(which(Xvar == breaks[2])+1):which(Xvar == breaks[3]),i]~Xvar[(which(Xvar == breaks[2])+1):which(Xvar == breaks[3])])
    RSS3<-sum(fit3$residuals^2)
    fit4<-lm(data[(which(Xvar == breaks[3])+1):which(Xvar == breaks[4]),i]~Xvar[(which(Xvar == breaks[3])+1):which(Xvar == breaks[4])])
    RSS4<-sum(fit4$residuals^2)
    fit5<-lm(data[(which(Xvar == breaks[4])+1):which(Xvar == breaks[5]),i]~Xvar[(which(Xvar == breaks[4])+1):which(Xvar == breaks[5])])
    RSS5<-sum(fit5$residuals^2)
    fit6<-lm(data[(which(Xvar == breaks[5])+1):length(Xvar),i]~Xvar[(which(Xvar == breaks[5])+1):length(Xvar)])
    RSS6<-sum(fit6$residuals^2)
    RSS.results[1,i]<-RSS1
    RSS.results[2,i]<-RSS2
    RSS.results[3,i]<-RSS3
    RSS.results[4,i]<-RSS4
    RSS.results[5,i]<-RSS5
    RSS.results[6,i]<-RSS6
  }
  }
  sumRSS<-apply(RSS.results,2, sum)
  RSS.results<-rbind(RSS.results, sumRSS)
  return(RSS.results)
}
