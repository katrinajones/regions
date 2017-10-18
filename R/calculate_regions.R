#' Calculate segmented regression models
#'
#'\code{calcregions} fits segmented regression models to all possible regions, up to noregions number of regions.
#'
#' @param Xvar Positional variable e.g., vertebral count.
#' @param Yvar Dependent variable.
#' @param noregions Maximum number of regions to be computed, up to 6.
#'
#' @return Returns a matrix containing residual sums of squares and r-squared values for each model.



calcregions <- function(Xvar, Yvar, noregions) {
    # Calculates sums of squares for each possible combination of regions test Estimate number
    # of models
    n <- length(Yvar)
    mod <- 1 + (n - 3) + ((n - 5) * (n - 4)/2) + ((n - 7) * (n - 6) * (n - 5)/6) + ((n - 9) *
        (n - 8) * (n - 7) * (n - 6)/24) + ((n - 11) * (n - 10) * (n - 9) * (n - 8) * (n - 7)/120)  #3 region= triangular number, 4= tetrahedral number

    # make an empty row for the seed of the results table
    regions <- matrix(data = NA, nrow = mod+1, ncol = 8)
    colhead <- c("regions", "breakpoint1", "breakpoint2", "breakpoint3", "breakpoint4", "breakpoint5",
        "RSS", "r squared")
    colnames(regions) <- colhead
    noverts <- length(Yvar)
    modno <- 1

    # one region
    line1 <- lm(Yvar ~ Xvar)
    RSSline1 <- sum(line1$residuals^2)
    rsq <- summary(line1)$r.squared
    regions[modno, ] <- c(1, 0, 0, 0, 0, 0, RSSline1, rsq)
    modno <- modno + 1

    if (noregions == 1) {
      return(regions)
    }

    # two regions
    for (x in 2:(noverts - 2)) {
        line1 <- lm(Yvar[1:x] ~ Xvar[1:x])
        RSSline1 <- sum(line1$residuals^2)
        line2 <- lm(Yvar[(x + 1):noverts] ~ Xvar[(x + 1):noverts])
        RSSline2 <- sum(line2$residuals^2)
        totRSS <- RSSline1 + RSSline2
        rsq <- mean(c(summary(line1)$r.squared, summary(line2)$r.squared))

        regions[modno, ] <- c(2, Xvar[x], 0, 0, 0, 0, totRSS, rsq)
        modno <- modno + 1

    }
    if (noregions == 2) {
        return(regions)
    }

    # three regions
    for (x in 2:(noverts - 4)) {
        for (y in (x + 2):(noverts - 2)) {
            line1 <- lm(Yvar[1:x] ~ Xvar[1:x])
            RSSline1 <- sum(line1$residuals^2)
            line2 <- lm(Yvar[(x + 1):y] ~ Xvar[(x + 1):y])
            RSSline2 <- sum(line2$residuals^2)
            line3 <- lm(Yvar[(y + 1):noverts] ~ Xvar[(y + 1):noverts])
            RSSline3 <- sum(line3$residuals^2)
            totRSS <- RSSline1 + RSSline2 + RSSline3
            rsq <- mean(c(summary(line1)$r.squared, summary(line2)$r.squared, summary(line3)$r.squared))

            regions[modno, ] <- c(3, Xvar[x], Xvar[y], 0, 0, 0, totRSS, rsq)
            modno <- modno + 1
        }
    }

    if (noregions == 3) {
        return(regions)
    }

    # four regions
    for (x in 2:(noverts - 6)) {
        for (y in (x + 2):(noverts - 4)) {
            for (z in (y + 2):(noverts - 2)) {
                line1 <- lm(Yvar[1:x] ~ Xvar[1:x])
                RSSline1 <- sum(line1$residuals^2)
                line2 <- lm(Yvar[(x + 1):y] ~ Xvar[(x + 1):y])
                RSSline2 <- sum(line2$residuals^2)
                line3 <- lm(Yvar[(y + 1):z] ~ Xvar[(y + 1):z])
                RSSline3 <- sum(line3$residuals^2)
                line4 <- lm(Yvar[(z + 1):noverts] ~ Xvar[(z + 1):noverts])
                RSSline4 <- sum(line4$residuals^2)
                totRSS <- RSSline1 + RSSline2 + RSSline3 + RSSline4
                rsq <- mean(c(summary(line1)$r.squared, summary(line2)$r.squared, summary(line3)$r.squared,
                  summary(line4)$r.squared))

                regions[modno, ] <- c(4, Xvar[x], Xvar[y], Xvar[z], 0, 0, totRSS, rsq)
                modno <- modno + 1
            }
        }
    }

    if (noregions == 4) {
        return(regions)
    }

    # five regions
    for (x in 2:(noverts - 8)) {
        for (y in (x + 2):(noverts - 6)) {
            for (z in (y + 2):(noverts - 4)) {
                for (a in (z + 2):(noverts - 2)) {

                  line1 <- lm(Yvar[1:x] ~ Xvar[1:x])
                  RSSline1 <- sum(line1$residuals^2)
                  line2 <- lm(Yvar[(x + 1):y] ~ Xvar[(x + 1):y])
                  RSSline2 <- sum(line2$residuals^2)
                  line3 <- lm(Yvar[(y + 1):z] ~ Xvar[(y + 1):z])
                  RSSline3 <- sum(line3$residuals^2)
                  line4 <- lm(Yvar[(z + 1):a] ~ Xvar[(z + 1):a])
                  RSSline4 <- sum(line4$residuals^2)
                  line5 <- lm(Yvar[(a + 1):noverts] ~ Xvar[(a + 1):noverts])
                  RSSline5 <- sum(line5$residuals^2)
                  totRSS <- RSSline1 + RSSline2 + RSSline3 + RSSline4 + RSSline5
                  rsq <- mean(c(summary(line1)$r.squared, summary(line2)$r.squared, summary(line3)$r.squared,
                    summary(line4)$r.squared, summary(line5)$r.squared))

                  regions[modno, ] <- c(5, Xvar[x], Xvar[y], Xvar[z], Xvar[a], 0, totRSS, rsq)
                  modno <- modno + 1

                }
            }
        }
    }

    if (noregions == 5) {
        return(regions)
    }

    # six regions
    for (x in 2:(noverts - 10)) {
        for (y in (x + 2):(noverts - 8)) {
            for (z in (y + 2):(noverts - 6)) {
                for (a in (z + 2):(noverts - 4)) {
                  for (b in (a + 2):(noverts - 2)) {

                    line1 <- lm(Yvar[1:x] ~ Xvar[1:x])
                    RSSline1 <- sum(line1$residuals^2)
                    line2 <- lm(Yvar[(x + 1):y] ~ Xvar[(x + 1):y])
                    RSSline2 <- sum(line2$residuals^2)
                    line3 <- lm(Yvar[(y + 1):z] ~ Xvar[(y + 1):z])
                    RSSline3 <- sum(line3$residuals^2)
                    line4 <- lm(Yvar[(z + 1):a] ~ Xvar[(z + 1):a])
                    RSSline4 <- sum(line4$residuals^2)
                    line5 <- lm(Yvar[(a + 1):b] ~ Xvar[(a + 1):b])
                    RSSline5 <- sum(line5$residuals^2)
                    line6 <- lm(Yvar[(b + 1):noverts] ~ Xvar[(b + 1):noverts])
                    RSSline6 <- sum(line6$residuals^2)

                    totRSS <- RSSline1 + RSSline2 + RSSline3 + RSSline4 + RSSline5 + RSSline6
                    rsq <- mean(c(summary(line1)$r.squared, summary(line2)$r.squared, summary(line3)$r.squared,
                      summary(line4)$r.squared, summary(line5)$r.squared, summary(line6)$r.squared))

                    regions[modno, ] <- c(6, Xvar[x], Xvar[y], Xvar[z], Xvar[a], Xvar[b], totRSS,
                      rsq)
                    modno <- modno + 1

                  }
                }
            }
        }
    }

    if (noregions == 6) {
        return(regions)
    }

}
