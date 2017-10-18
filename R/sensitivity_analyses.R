#' Simulate missing data
#'
#' Randomly replace data elements with NAs, or remove whole vertebrae or variables
#'
#' If vertebrae are removed Xvar is corrected
#'
#' @param data data matrix
#' @param missing proportion missing data, from zero to one
#' @param type "elem", "vert", "var" randomly replace elements, whole rows, or whole columns
#' @param Xvar Independent variable
#'
#' @return data
#' @return xvar
#' @export
#'
#'
sim_missing<-function(data,missing,type, Xvar){

  if(type=="elem"){

    no_rm_el<-round((nrow(data)*ncol(data))*missing) #number of elements to remove
    ind<-round(runif(no_rm_el,1,(nrow(data)*ncol(data))))#generate random positions
    data[ind]<-NA # replace with NA

  }

  if(type=="vert"){

    no_rm_vert<-round(nrow(data)*missing)
    ind<-sample(c(1:nrow(data)), no_rm_vert)
    data<-data[-ind,]#drop missing vertebrae
    Xvar<-Xvar[-ind]
  }
  if(type=="var"){

    no_rm_var<-round(ncol(data)*missing)
    ind<-sample(c(1:ncol(data)), no_rm_var)
    data<-data[,-ind]

  }

  return(list(data=data,Xvar=Xvar))

}


#' Run Sensistivity analysis
#'
#' Reruns regionalization analysis on bootstrapped data
#'
#' Data will be bootstrapped by element (cells in matrix), by vertebra (rows), or by variable (cols)
#'
#' @param Xvar Positional variable
#' @param data data
#' @param noregions Number of regions
#' @param nopcos Number of PCOs
#' @param missing proportion missing, between zero and one
#' @param type "elem", "vert", or "var"
#' @param reps Number of repeats
#' @param fillmiss Fill missing TRUE or FALSE
#'
#' @return regionscore
#' @return AIC models
#' @export
#'
#'
sens_anal<-function(Xvar, data,noregions, nopcos, missing, type, reps, fillmiss=TRUE){
  AIC_models<-NULL
  regionscore<-NULL

  i=1
  while (i<=reps){

    sim_data<-regions::sim_missing(data, missing, type,Xvar)#create simulated dataset
    simdata<-sim_data$data

    if(fillmiss==TRUE){
     #fill missing based on adjacent
    simdata<-Missingval(simdata)
    }

    #test sim_data
    test<-length(which(is.na(dist(simdata))==TRUE))#check if theres any conflicting data
    if(test>0){
      next#skip if it conflicts - avoids error mg
       }

    simXvar<-sim_data$Xvar
    nvert<-length(simXvar)
    #run PCO and regionalization analysis
    simdata<-scale(simdata)#add a normalization step for bootstrapping so the variances are the same
    pco.gower<-svdPCO(simdata, "gower")

    if(nopcos=="boot"){
    pco.boot<-PCOcutoff(simdata,1000, "gower")
    nopcos2<-pco.boot$sigpco
    }  else if(nopcos=="varcut"){
    nopcos2<-length(which(pco.gower$eigen.val/sum(pco.gower$eigen.val)>0.05))#more than 5% of variance
    } else if(nopcos=="max"){
      PCOscores<-pco.gower$scores[,1:ncol(pco.gower$scores)]#calculate region models to estimate pco.max
      regiondata<-compileregions(simXvar,PCOscores,noregions)
      nopcos2<-PCOmax(regiondata, noregions, nvert)$pco.max
    }   else if (nopcos=="all"){
      nopcos2<-ncol(pco.gower$scores)
    }
    else {    nopcos2<-nopcos
    }
    if (nopcos2==0){#if none significant default to 1
      nopcos2<-1
    }
    #print(c(ncol(pco.gower$scores),nopcos2))
    PCOscores<-pco.gower$scores[,1:nopcos2]


    #Calculate sums of squares, and AIC values for all combinations of regions for all variables
    regions<-region_anal(simXvar, PCOscores,noregions, nopcos2)
    score<-cbind(i, regions$score, nopcos2)
    regionscore<-rbind(regionscore, score)
    AIC_sum<-cbind(rep(i,nrow(regions$model_support)),regions$model_support)
    AIC_models<-rbind(AIC_models, AIC_sum)

    print(i)
    i=i+1 #Move to next case
     }

  return(list(regionscore=regionscore,AIC_models=AIC_models))
}



