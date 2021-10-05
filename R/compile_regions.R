#' Calculates region model fit, edited for 7 region capability. SMSmith 10.05.2021
#'
#' Uses \code{calcregions} to calculate region models for each variable separately, then
#' exports the summed residual sums of squares.
#'
#' This function feeds data into \code{calcregions}, compiles the results and calculates model fitting
#' metrics.
#'
#' @param Xvar Positional variable e.g., vertebral count.
#' @param data Matrix of dependent variables.
#' @param noregions Maximum number of regions to be computed, up to 6.
#'
#' @return Returns a data matrix with one row for each model. Columns include model parameters
#'  (breakpoints) and residual sums of squares for each dependent variable
#'
#' @export
#' @examples
#' data("alligator")
#' Xvar=alligator[,1]
#' data=alligator[,2:ncol(alligator)]
#' regiondata=compileregions(Xvar, data, 2)

compileregions <- function(Xvar, data, noregions=2) {

    data <- as.matrix(data)

    if (ncol(data) == 1) {
        # for one variable only
        first <- calcregions(Xvar, data, noregions)
        first <- first[which(first[, 1] <= noregions), ]  #remove unused rows
        restab <- first[, 1:8]  #cols inlcuding RSS
        rsq <- first[, 9]  #RSQ col

    } else {
        # make original results table
        first <- calcregions(Xvar, data[, 1], noregions)
        first <- first[which(first[, 1] <= noregions), ]  #remove unused rows
        restab <- first[, 1:8]  #cols inlcuding RSS
        rsq <- first[, 9]  #RSQ col

        # run for the other PCs
        for (i in 2:ncol(data)) {
            more <- calcregions(Xvar, data[, i], noregions)
            more <- more[which(more[, 1] <= noregions), ]
            restab <- cbind(restab, more[, 8])
            rsq <- cbind(rsq, more[, 9])
        }
        colnames(rsq) <- paste("rsq", c(1:ncol(data)), sep = "")
    }
    # rename the columns
    names <- colnames(first)
    PCOvars <- c(1:ncol(data))
    varnames <- sapply(PCOvars, function(x) paste("var", x))
    names <- append(names[1:7], varnames)
    colnames(restab) <- names


    return(restab)
}
