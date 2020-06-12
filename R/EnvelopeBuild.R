
EnvelopeBuild<-function(bStar, A, y, x, mu, P, alpha, wt, family = "binomial",
                         link = "logit", Gridtype = 2L, n = 1L, sortgrid = FALSE){
  
  return(.EnvelopeBuild_cpp(bStar, A, y, x, mu, P, alpha, wt, family = family,
                         link = link, Gridtype = Gridtype, n = n, sortgrid))
}
  
