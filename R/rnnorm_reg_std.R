
rnnorm_reg_std<-function(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = 1L){
  
  return(.rnnorm_reg_std_cpp(n, y, x, mu, P, alpha, wt, f2, Envelope, family, link, progbar = progbar))

  }
  
