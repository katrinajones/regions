#' Cos transforms angles
#'
#' Converts angles to radians then cos transforms them.
#'
#' @param data Raw data including angles in degrees.
#' @param angles Vector of strings indicating columns which contain angles.
#' @export

angtransform <- function(data, angles) {

    # transform angles - Do I need to transform angles with Gower?? does it make a difference?
    angles <- na.omit(angles[match(colnames(data), angles)])  #find only angles that are in dataset
    if (length(angles) > 0) {
        for (m in 1:length(angles)) {
            deg2rad <- function(deg) {
                (deg * pi)/(180)
            }  ####Must convert angles to radians before cos transf
            data[, angles[m]] <- deg2rad(data[, angles[m]])
            data[, angles[m]] <- cos(data[, angles[m]])
        }
    }
    return(data)
}

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

 miss.par<-which(is.na(data), arr.ind = T)#find which ones are missing
for (i in 1:nrow(miss.par)){
  miss.point<-miss.par[i,]
  avail.points<-which(!is.na(data[,miss.point[2]]))#make list of available values for that variable
  location<-avail.points-miss.point[1]#how close are they?

  #Ensure there arent too many missing - max=2
  dist.beg<-miss.point[1]-1
dist.end<-abs(miss.point[1]-nrow(data))
  if(!any(location<0)&dist.beg>2){#break if its at the beginning and too far from first pos
    next()
  }
    if(!any(location>0)&dist.end>2){#break if its at the end and too far from last pos
      next()
    }


  if (any(location<0)){#if there are points before, use those, else use the point after
    before.val<-max(location[location<0])
    before.point<-c(avail.points[which(location==before.val)],miss.point[2])
  } else{
    before.val<-min(location[location>0])
    before.point<-c(avail.points[which(location==before.val)],miss.point[2])
  }
  if (any(location>0)){#if there are points after, use those, else use the point before
    after.val<-min(location[location>0])
    after.point<-c(avail.points[which(location==after.val)],miss.point[2])
  }      else{
    after.val<-max(location[location<0])
    after.point<-c(avail.points[which(location==after.val)],miss.point[2])
  }

  if (any(c(abs(before.val), abs(after.val))>2)){ #exclude if more than two in a row
    next()
  }
  before.num<-data[before.point[1], before.point[2]]
  after.num<-data[after.point[1], after.point[2]]
  est.point<-mean(c(before.num, after.num))
  data[miss.point[1], miss.point[2]]<-est.point[1]

}
return(data)
  }else{
  return(data)}

}
