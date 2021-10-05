#' Plot PCO axes
#'
#' Compare any PCO axes along the series.
#'
#' @param data PCO scores or any multivariate data matrix.
#' @param x Variable for x axis.
#' @param y Variable for y axis.
#' @param Xvar Variable for labelling points, usually positional variable such as vertebral count.
#' @export

axesplot <- function(data, x, y, Xvar) {

  plot(data[, x], data[, y], col = "white", xlab = "", ylab = "")
  graphics::text(data[, x], data[, y], labels = Xvar)
  graphics::title(main = "Principal coordinates analysis", xlab = paste("PCO", x), ylab = paste("PCO",
                                                                                                y))

}

#' Graphs top region model
#'
#' Produces an illustration of serially-homologous structure divided by region, based on
#' best model from \code{model_support}.
#'
#' @param name Specimen name for graph
#' @param Xvar Positional variable such as vertebral count.
#' @param regiondata Model_support object from \code{model_support}.
#' @param plot.reg plot the data?
#' @param bestmodel option to input just the best model
#' @param first option to input the position of the first vertebral sampled (default is 3)
#'
#' @import ggplot2
#' @export

regionmodel <- function(name, Xvar, regiondata=NULL, plot.reg=TRUE, bestmodel=NULL, first=3) {

  if(is.null(bestmodel)){
    bestmodel <- as.matrix(regiondata[1, 1:7])  #bestmodel is first row from analysis
  }else{
    bestmodel<-bestmodel
  }
  firstvert <- Xvar[1]
  lastvert <- Xvar[length(Xvar)]
  break1 <- which(Xvar == bestmodel[2])  #set the breakpoints
  break2 <- which(Xvar == bestmodel[3])
  break3 <- which(Xvar == bestmodel[4])
  break4 <- which(Xvar == bestmodel[5])
  break5 <- which(Xvar == bestmodel[6])
  break6 <- which(Xvar == bestmodel[7])
  breakpoints <- rep(0, lastvert)

  # convert breakpoints into regions
  if (bestmodel[1] == 1) {
    breakpoints[Xvar]<-1
    colorlab <- c("white", "red", "black")
  }

  if (bestmodel[1] == 2) {
    breakpoints[Xvar[1:break1]] <- 1
    breakpoints[Xvar[break1 + 1:length(Xvar)]] <- 2
    colorlab <- c("white", "red", "yellow", "black")
  }
  if (bestmodel[1] == 3) {
    breakpoints[Xvar[1:break1]] <- 1
    breakpoints[Xvar[break1 + 1:break2]] <- 2
    breakpoints[Xvar[break2 + 1:length(Xvar)]] <- 3
    colorlab <- c("white", "red", "yellow", "blue", "black")
  }
  if (bestmodel[1] == 4) {
    breakpoints[Xvar[1:break1]] <- 1
    breakpoints[Xvar[break1 + 1:break2]] <- 2
    breakpoints[Xvar[break2 + 1:break3]] <- 3
    breakpoints[Xvar[break3 + 1:length(Xvar)]] <- 4
    colorlab <- c("white", "red", "darkorange", "yellow", "blue", "black")
  }
  if (bestmodel[1] == 5) {
    breakpoints[Xvar[1:break1]] <- 1
    breakpoints[Xvar[break1 + 1:break2]] <- 2
    breakpoints[Xvar[break2 + 1:break3]] <- 3
    breakpoints[Xvar[break3 + 1:break4]] <- 4
    breakpoints[Xvar[break4 + 1:length(Xvar)]] <- 5
    colorlab <- c("white", "red", "darkorange", "yellow", "green", "blue", "black")
  }
  if (bestmodel[1] == 6) {
    breakpoints[Xvar[1:break1]] <- 1
    breakpoints[Xvar[break1 + 1:break2]] <- 2
    breakpoints[Xvar[break2 + 1:break3]] <- 3
    breakpoints[Xvar[break3 + 1:break4]] <- 4
    breakpoints[Xvar[break4 + 1:break5]] <- 5
    breakpoints[Xvar[break5 + 1:length(Xvar)]] <- 6
    colorlab <- c("white", "red", "darkorange", "yellow", "green","aquamarine", "blue", "black")
  }
  if (bestmodel[1] == 7) {
    breakpoints[Xvar[1:break1]] <- 1
    breakpoints[Xvar[break1 + 1:break2]] <- 2
    breakpoints[Xvar[break2 + 1:break3]] <- 3
    breakpoints[Xvar[break3 + 1:break4]] <- 4
    breakpoints[Xvar[break4 + 1:break5]] <- 5
    breakpoints[Xvar[break5 + 1:break6]] <- 6
    breakpoints[Xvar[break6 + 1:length(Xvar)]] <- 7
    colorlab <- c("white", "red", "darkorange", "yellow", "green","aquamarine", "lavender", "blue", "black")
  }

  breakpoints[first:lastvert] <- replace(breakpoints[first:lastvert], which(breakpoints[first:lastvert] ==
                                                                              0), 8)  #make option for missing data
  breakpoints <- as.factor(breakpoints)
  x <- factor(1)
  plot.data <- data.frame(x, breakpoints)

  # Vertebra region plot
  xmin <- c(0:(lastvert - 1))
  xmax <- c(1:lastvert)
  ymin <- rep(0, lastvert)
  ymax <- rep(0.2, lastvert)
  plot.data <- cbind(plot.data, xmin, xmax, ymin, ymax)
  title <- name

  if(isTRUE(plot.reg)){

    regionplot<-ggplot2::ggplot(plot.data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + geom_rect(colour = "black",
                                                                                                                alpha = 0.5, aes(fill = breakpoints)) + scale_fill_manual(values = colorlab, guide = FALSE) +
      geom_text(aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2, label = xmax)) + theme(axis.text.y = element_blank(),
                                                                                     axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                                                                                     panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
                                                                                     panel.grid.minor.y = element_blank(), panel.background = element_rect(fill = "white")) +
      ylab("") + xlab("")+ ggtitle(title) +
      coord_fixed(ratio = 5)

    return(regionplot)

  }
  if(!isTRUE(plot.reg)){
    return(plot.data$breakpoints)
  }

}

#' Plot best segmented regression model
#'
#' Plots the segmented regression for the top model on PCO1.
#'
#' @param Xvar Positional vector such as vertebral count.
#' @param data PCO scores.
#' @param pcono Which PC to plot
#' @param modelsupport Model_support object from \code{model_support}.
#'
#' @importFrom graphics abline
#' @importFrom plotrix ablineclip
#' @export

plotsegreg <- function(Xvar, pcono, data, modelsupport) {

  bestmodel <- as.matrix(modelsupport[1, 1:7])  #bestmodel is first row from analysis
  firstvert <- Xvar[1]
  lastvert <- Xvar[length(Xvar)]
  break1 <- which(Xvar == bestmodel[2])  #set the breakpoints
  break2 <- which(Xvar == bestmodel[3])
  break3 <- which(Xvar == bestmodel[4])
  break4 <- which(Xvar == bestmodel[5])
  break5 <- which(Xvar == bestmodel[6])
  break6 <- which(Xvar == bestmodel[7])

  label.pco<-paste0("PCO",pcono)
  plot(Xvar, data[, pcono], xlab = "Vertebral position", ylab = label.pco)
  title(main = label.pco)

  if (bestmodel[1] == 1) {
    abline(lm(data[, pcono] ~ Xvar))
  }
  if (bestmodel[1] == 2) {
    ablineclip(lm(data[1:break1, pcono] ~ Xvar[1:break1]), x1 = Xvar[1], x2 = bestmodel[2])
    ablineclip(lm(data[(break1 + 1):nrow(data), pcono] ~ Xvar[(break1 + 1):nrow(data)]), x1 = bestmodel[2] +
                 1, x2 = Xvar[length(Xvar)])
    abline(v = bestmodel[2] + 0.5, col = "red", lty = 2)
  }
  if (bestmodel[1] == 3) {
    ablineclip(lm(data[1:break1, pcono] ~ Xvar[1:break1]), x1 = Xvar[1], x2 = bestmodel[2])
    ablineclip(lm(data[(break1 + 1):break2, pcono] ~ Xvar[(break1 + 1):break2]), x1 = bestmodel[2] +
                 1, x2 = bestmodel[3])
    ablineclip(lm(data[(break2 + 1):nrow(data), pcono] ~ Xvar[(break2 + 1):nrow(data)]), x1 = bestmodel[3] +
                 1, x2 = Xvar[length(Xvar)])
    abline(v = bestmodel[2] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[3] + 0.5, col = "red", lty = 2)
  }
  if (bestmodel[1] == 4) {
    ablineclip(lm(data[1:break1, pcono] ~ Xvar[1:break1]), x1 = Xvar[1], x2 = bestmodel[2])
    ablineclip(lm(data[(break1 + 1):break2, pcono] ~ Xvar[(break1 + 1):break2]), x1 = bestmodel[2] +
                 1, x2 = bestmodel[3])
    ablineclip(lm(data[(break2 + 1):break3,pcono] ~ Xvar[(break2 + 1):break3]), x1 = bestmodel[3] +
                 1, x2 = bestmodel[4])
    ablineclip(lm(data[(break3 + 1):nrow(data), pcono] ~ Xvar[(break3 + 1):nrow(data)]), x1 = bestmodel[4] +
                 1, x2 = Xvar[length(Xvar)])
    abline(v = bestmodel[2] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[3] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[4] + 0.5, col = "red", lty = 2)
  }
  if (bestmodel[1] == 5) {
    ablineclip(lm(data[1:break1, pcono] ~ Xvar[1:break1]), x1 = Xvar[1], x2 = bestmodel[2])
    ablineclip(lm(data[(break1 + 1):break2, pcono] ~ Xvar[(break1 + 1):break2]), x1 = bestmodel[2] +
                 1, x2 = bestmodel[3])
    ablineclip(lm(data[(break2 + 1):break3, pcono] ~ Xvar[(break2 + 1):break3]), x1 = bestmodel[3] +
                 1, x2 = bestmodel[4])
    ablineclip(lm(data[(break3 + 1):break4, pcono] ~ Xvar[(break3 + 1):break4]), x1 = bestmodel[4] +
                 1, x2 = bestmodel[5])
    ablineclip(lm(data[(break4 + 1):nrow(data), pcono] ~ Xvar[(break4 + 1):nrow(data)]), x1 = bestmodel[5] +
                 1, x2 = Xvar[length(Xvar)])
    abline(v = bestmodel[2] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[3] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[4] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[5] + 0.5, col = "red", lty = 2)
  }
  if (bestmodel[1] == 6) {
    ablineclip(lm(data[1:break1, pcono] ~ Xvar[1:break1]), x1 = Xvar[1], x2 = bestmodel[2])
    ablineclip(lm(data[(break1 + 1):break2, pcono] ~ Xvar[(break1 + 1):break2]), x1 = bestmodel[2] +
                 1, x2 = bestmodel[3])
    ablineclip(lm(data[(break2 + 1):break3, pcono] ~ Xvar[(break2 + 1):break3]), x1 = bestmodel[3] +
                 1, x2 = bestmodel[4])
    ablineclip(lm(data[(break3 + 1):break4, pcono] ~ Xvar[(break3 + 1):break4]), x1 = bestmodel[4] +
                 1, x2 = bestmodel[5])
    ablineclip(lm(data[(break4 + 1):break5, pcono] ~ Xvar[(break4 + 1):break5]), x1 = bestmodel[5] +
                 1, x2 = bestmodel[6])
    ablineclip(lm(data[(break5 + 1):nrow(data),pcono] ~ Xvar[(break5 + 1):nrow(data)]), x1 = bestmodel[6] +
                 1, x2 = Xvar[length(Xvar)])
    abline(v = bestmodel[2] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[3] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[4] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[5] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[6] + 0.5, col = "red", lty = 2)
  }
  if (bestmodel[1] == 7) {
    ablineclip(lm(data[1:break1, pcono] ~ Xvar[1:break1]), x1 = Xvar[1], x2 = bestmodel[2])
    ablineclip(lm(data[(break1 + 1):break2, pcono] ~ Xvar[(break1 + 1):break2]), x1 = bestmodel[2] +
                 1, x2 = bestmodel[3])
    ablineclip(lm(data[(break2 + 1):break3, pcono] ~ Xvar[(break2 + 1):break3]), x1 = bestmodel[3] +
                 1, x2 = bestmodel[4])
    ablineclip(lm(data[(break3 + 1):break4, pcono] ~ Xvar[(break3 + 1):break4]), x1 = bestmodel[4] +
                 1, x2 = bestmodel[5])
    ablineclip(lm(data[(break4 + 1):break5, pcono] ~ Xvar[(break4 + 1):break5]), x1 = bestmodel[5] +
                 1, x2 = bestmodel[6])
    ablineclip(lm(data[(break5 + 1):break6, pcono] ~ Xvar[(break5 + 1):break6]), x1 = bestmodel[6] +
                 1, x2 = bestmodel[7])
    ablineclip(lm(data[(break6 + 1):nrow(data),pcono] ~ Xvar[(break6 + 1):nrow(data)]), x1 = bestmodel[7] +
                 1, x2 = Xvar[length(Xvar)])
    abline(v = bestmodel[2] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[3] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[4] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[5] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[6] + 0.5, col = "red", lty = 2)
    abline(v = bestmodel[7] + 0.5, col = "red", lty = 2)
  }
}


#' Make a scree plot of variables against species
#'
#' @param yvar yvar
#' @param labels labels
#' @param header header
#'
#' @importFrom graphics axis
#' @export
#'
#'
plotvar<-function(yvar, labels, header){

  plot.data<-data.frame(yvar,labels)
  plot.data<-plot.data[order(plot.data$yvar),]
  plot(plot.data$yvar,xaxt="n", xlab="", ylab=header)
  title(main=header)
  axis(1,c(1:nrow(plot.data)),plot.data$labels, cex.axis=0.9, las=2)

}

#' Plot the scree distribution of regionscores across pcos
#'
#' Plot indicating how regionscore is influenced by including more PCOs
#'
#' Plots cumulative variance (line), regionscore based on each separate PCO (black) and regionscore based
#' on cumulative PCOs (red)
#'
#' @param eigenvals Eigenvalues from PCO
#' @param nvert Number of vertebrae
#' @param namelabel Label for plot (specimen name)
#' @param regiondata Segmented regression models for specimen from \code{compileregions}
#' @param noregions Maximum number of regions
#'
#' @importFrom graphics points lines plot title
#' @export
#'
#'
#'
plotpcoreg<-function(eigenvals,nvert, namelabel, regiondata, noregions){

  nvar<-ncol(regiondata[,which(colnames(regiondata)=="var 1"|colnames(regiondata)=="var.1"):ncol(regiondata)])
  #Make pco distribution plots
  var.exp<-cumsum(sapply(eigenvals[1:nvar], function(x) (x/sum(eigenvals))*100))

  pco.no.test<-data.frame(matrix(data=NA,nrow=nvar,ncol=5, dimnames=list(c(1:nvar),c("species","pco", "ind", "cum", "var.exp"))))

  for (a in 1:nvar){

    #Run for individual PCs
    models.ind<-modelselect(regiondata,noregions,a, startpco = a)
    support.ind<-model_support(models.ind,nvert, 1)

    #Run for cumulative PCs
    models.cum<-modelselect(regiondata,noregions,a, startpco = 1)
    support.cum<-model_support(models.cum,nvert, a)

    pco.no.test[a,1]<-namelabel
    pco.no.test[a,2]<-a
    pco.no.test[a,3]<-support.ind$Region_score
    pco.no.test[a,4]<-support.cum$Region_score
  }

  pco.no.test[,5]<-var.exp/(100/noregions)

  plot(pco.no.test$pco, pco.no.test$ind, xlab="PCO", ylab="Regionscore", ylim=c(1,noregions))
  points(pco.no.test$pco, pco.no.test$cum, col="red")
  lines(pco.no.test$pco, pco.no.test$var.exp)
  title(main=namelabel)

}


#' Estimate loadings on PCOS
#'
#' Estimates loadings based on correlation coefficients of PCOs against raw data
#'
#' @param data Raw data
#' @param PCOscore PCO scores
#'
#' @return loadings correlation coefficients and direction of correlation
#'
#' @export pco.load
#'
#'
pco.load<-function(data, PCOscore){

  ###Figure out loadings using correlation
  vert.size<-apply(data,1,mean)
  data<-cbind(data, vert.size)
  load.pco<-stats::cor(data[,3:ncol(data)],PCOscore, use="pairwise.complete.obs")
  abs<-abs(load.pco)
  sign<-sign(load.pco[,1])
  load.pco<-cbind(load.pco,abs, sign)

  load.pco<-load.pco[order(abs, decreasing = TRUE),]


  return(loadings=load.pco)

}


#' Plots heatmap of region boundaries
#'
#' Relative probability of break and each position based on RSS for a given number of regions
#'
#' Also calculates the standard deviation of each breakpoint in the top 5% of models based on RSS
#' If there are fewer then five models in the top 5%, the top five are used.
#'
#' @param regiondata region data
#' @param noregions number of regions for which to generate heat score
#' @param nopcos No of PCOs
#' @param Xvar Positional variable
#'
#' @return heatscore- probability of breaks at a given poistion
#' @return stdev - standard deviation of each breakpoint
#' @importFrom stats sd
#'
regionheat<-function(regiondata, noregions,nopcos, Xvar){

  #Calculate RSS
  subdata<-regiondata[which(regiondata[,1]==noregions),]#select models with correct no of regions
  pco.begin<-which(colnames(subdata)=="var.1"|colnames(subdata)=="var 1")#find the RSS values
  if(nopcos==1){sumRSS=subdata[,pco.begin:(pco.begin+(nopcos-1))]
  }else{
    sumRSS <- rowSums(subdata[,pco.begin:(pco.begin+(nopcos-1))]) #calculate total RSS
  }
  subdata<-cbind(subdata, sumRSS)
  subdata<-as.data.frame(subdata)

  #How much worse than best model
  RSSmin<-min(subdata$sumRSS)#RSS of best model
  deltaRSS<-sapply(subdata$sumRSS, function(x) x-RSSmin)#Calculate RSS difference
  subdata<-cbind(subdata, deltaRSS)

  ##Sum deltaRSS for each position
  mods<-subdata[,2:(pco.begin-1)]
  positions<-matrix(data=NA, nrow=(length(Xvar)-3), ncol=2)
  for(i in 2:(length(Xvar)-2)){
    vert<-Xvar[i]
    rowvert<-which(mods==vert, arr.ind = T)
    rowvert<-c(rowvert[,1])
    vertscore<-sum(subdata$deltaRSS[rowvert])
    positions[i-1,1]<-vert
    positions[i-1,2]<-vertscore
  }

  ###Calculate error associated with each break
  subdata<-subdata[order(subdata$deltaRSS),]#sort by RSS
  cutoff<-0.05*(range(subdata$deltaRSS)[2])#based on RSS, doesnt work for few regions
  best<-subdata[which(subdata$deltaRSS<cutoff),]
  if(nrow(best<5)){
    best<-subdata[1:5,]#minimum possible is five
  }

  sd1<-sd(best$breakpoint1)#Standard deviation of breakpoints
  sd2<-sd(best$breakpoint2)
  sd3<-sd(best$breakpoint3)
  sd4<-sd(best$breakpoint4)
  sd5<-sd(best$breakpoint5)
  sd6<-sd(best$breakpoint6)

  stdev<-list(b1=sd1, b2=sd2, b3=sd3, b4=sd4, b5=sd5, b6 = sd6)
  stdev<-stdev[1:(noregions-1)]#remove empty

  return(list(heatscore=positions,stdev=stdev))

}

