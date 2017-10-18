#' Calculate disparity between regions
#'
#' Calculates variance of region means based on pre-defined region model
#'
#' @param data data frame with positional information in column "vert"
#' @param bestmodel vector in format (no regions, break 1, break 2...)
#'
#' @return var Variance of means
#' @return mean Region means
#' @export

disparityregions=function(data, bestmodel){

  bestmodel<-bestmodel[which(bestmodel>0)]
  for (a in 2:(length(bestmodel))){ #take into account missing verts
    bestmodel[a]<-match(bestmodel[a],data$vert,nomatch=match((bestmodel[a]+1),data$vert))
  }
  bestmodel<-as.numeric(bestmodel)
    #Calculate means
    if (bestmodel[1]==2){
      mean1=apply(data[1:bestmodel[2],3:ncol(data)],2,mean)
      mean2=apply(data[(bestmodel[2]+1):nrow(data),3:ncol(data)],2,mean)
      mean=rbind(mean1, mean2)
    }

    if (bestmodel[1]==3){
      mean1=apply(data[1:bestmodel[2],3:ncol(data)],2,mean)
      mean2=apply(data[(bestmodel[2]+1):bestmodel[3],3:ncol(data)],2,mean)
      mean3=apply(data[(bestmodel[3]+1):nrow(data),3:ncol(data)],2,mean)
      mean=rbind(mean1, mean2, mean3)
    }

    if (bestmodel[1]==4){
      mean1=apply(data[1:bestmodel[2],3:ncol(data)],2,mean)
      mean2=apply(data[(bestmodel[2]+1):bestmodel[3],3:ncol(data)],2,mean)
      mean3=apply(data[(bestmodel[3]+1):bestmodel[4],3:ncol(data)],2,mean)
      mean4=apply(data[(bestmodel[4]+1):nrow(data),3:ncol(data)],2,mean)
      mean=rbind(mean1, mean2, mean3, mean4)
    }
  if (bestmodel[1]==5){
    mean1=apply(data[1:bestmodel[2],3:ncol(data)],2,mean)
    mean2=apply(data[(bestmodel[2]+1):bestmodel[3],3:ncol(data)],2,mean)
    mean3=apply(data[(bestmodel[3]+1):bestmodel[4],3:ncol(data)],2,mean)
    mean4=apply(data[(bestmodel[4]+1):bestmodel[5],3:ncol(data)],2,mean)
    mean5=apply(data[(bestmodel[5]+1):nrow(data),3:ncol(data)],2,mean)
    mean=rbind(mean1, mean2, mean3, mean4, mean5)
  }

  if (bestmodel[1]==6){
    mean1=apply(data[1:bestmodel[2],3:ncol(data)],2,mean)
    mean2=apply(data[(bestmodel[2]+1):bestmodel[3],3:ncol(data)],2,mean)
    mean3=apply(data[(bestmodel[3]+1):bestmodel[4],3:ncol(data)],2,mean)
    mean4=apply(data[(bestmodel[4]+1):bestmodel[5],3:ncol(data)],2,mean)
    mean5=apply(data[(bestmodel[5]+1):bestmodel[6],3:ncol(data)],2,mean)
    mean6=apply(data[(bestmodel[6]+1):nrow(data),3:ncol(data)],2,mean)
    mean=rbind(mean1, mean2, mean3, mean4, mean5, mean6)
  }

    var=sum(apply(mean,2, var))
    return(list(var=var,mean=mean))

  }
