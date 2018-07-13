data("alligator")
Xvar<-alligator[,1]
nvert<-length(Xvar)
data<-alligator[,2:ncol(alligator)]
data<-scale(data)
pco.gower<-svdPCO(data, "gower")#PCO using svd
PCOscores<-pco.gower$scores[,1:ncol(pco.gower$scores)]
noregions<-2 #Set the maximum number of regions which will be calculated
regiondata<-compileregions(Xvar,PCOscores[,1:5],noregions)
models<-modelselect(regiondata,noregions,2)
support<-model_support(models,nvert, 2)

test_that("regions alligator test", {
  expect_equal(support$Model_support[1,1],2)
  expect_equal(support$Model_support[1,2],10)
  expect_equal(support$Model_support[1,7],0.08875614)
})
