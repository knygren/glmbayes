#' Design matrix and model frame for newdata
#'
#' Derives the model frame and design matrix for a complete version of the 
#' newdata data frame passed to the predict.glmb function
#' @param object an object of class \code{glmb}
#' @param olddata A data frame containing the model variables (and values) from the original model
#' @param newdata A data frame containing the model variables wih possible new values
#' @return A design matrix x_new and a model frame mod_frame_new
#' @example inst/examples/Ex_confint.glmb.R


get_x_matrix<-function(object,olddata,newdata){
  
  nrow1=nrow(olddata)
  nrow2=nrow(newdata)
  
  ## Stack olddata and newdata together
  
  combodata=rbind(olddata,newdata)
  combo_mod_frame=model.frame(formula(object$formula),combodata)
  combo_mod_matrix=model.matrix(formula(object$formula),combo_mod_frame)
  
  x_old=combo_mod_matrix[1:(nrow1),]
  x_new=combo_mod_matrix[(nrow1+1):(nrow1+nrow2),]
  
  ## assign attributes of combo matrix to sub-matrices 
  
  attr(x_old,which="assign")<-attr(combo_mod_matrix,which="assign")
  attr(x_old,which="contrasts")<-attr(combo_mod_matrix,which="contrasts")
  
  attr(x_new,which="assign")<-attr(combo_mod_matrix,which="assign")
  attr(x_new,which="contrasts")<-attr(combo_mod_matrix,which="contrasts")
  
  mod_frame_old=combo_mod_frame[1:(nrow1),]
  mod_frame_new=combo_mod_frame[(nrow1+1):(nrow1+nrow2),]
  
  x=object$glm$x
  
  if(!dim(x_old)[1]==dim(x)[1]) stop("Number of rows in final constucted x matrix for olddata does not match original")
  if(!dim(x_old)[2]==dim(x)[2]) stop("Number of columns in final constucted x matrix for olddata does not match original")
  
  for(i in 1:dim(x_old)[2]){
    if(!colnames(x_old)[i]==colnames(x)[i]) stop("Column names in final matrix do not match original")
  }
  
  x_new2=model.matrix(formula(object$glm$formula),mod_frame_new)
  

  return(list(x_new=x_new2,mod_frame_new=mod_frame_new))
}

