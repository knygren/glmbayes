
glmb_Standardize_Model<-function(y, x, P, bstar, A1){
  
  return(.glmb_Standardize_Model_cpp(y, x, P, bstar, A1))

  }
  
