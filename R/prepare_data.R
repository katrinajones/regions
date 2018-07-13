#' Estimate and fill missing values
#'
#' Estimates missing values by averaging values from nearest adjacent vertebrae.
#'
#' In the case of a first or last vertebra, the nearest value is taken. Will only fill up to two missing data points.
#'
#' @param data data with missing elements
#'
#' @return data Data with missing values filled
#' @export
#'
#
Missingval<-function(data){

  if(any(is.na(data)==TRUE)){

    ###find strings of NAs with 2 or less missing
    for(i in 1:ncol(data)){ #for each variable
      dat<-data[,i]
      if(!any(is.na(dat)==TRUE)){ #Is there any missing data?
        next} else{
      miss.par<-which(is.na(dat))#find which ones are missing
      seqs<-split(miss.par, cumsum(c(1, diff(miss.par) != 1)))#split them into sequences
      if(length(seqs)==1){ #if theres one string
        if(length(unlist(seqs))>2){ #if the string is longer than 2 skip
          next}else{
            miss<-seqs#otherwise add to miss to fill
          }
      }
      if(length(seqs)>1){ #if theres more than one string
        l.seq<-which(sapply(seqs, length)<3)#which strings are two or less
        if(length(l.seq)==0){ next} #if no short strings skip
        seqs<-seqs[l.seq]
        miss<-seqs
      }

      for(a in 1:length(miss)){ #Fill each string
        fill<-unlist(miss[a])
        before<-min(fill)-1
        after<-max(fill)+1
        if(before<1){before<-after}#if at the beginning, use the end points
        if(after>length(dat)){after<-before}#if at the end, use beginning points
        val<-mean(c(dat[before],dat[after])) #calculate missing as mean of adjacent
        dat[fill]<-val #fill in the missing
      }

      data[,i]<-dat

      }
    }
    return(data)

  } else{return(data)}
}
