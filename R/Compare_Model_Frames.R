#' Compare Model Frames
#'
#' Compares model frames and stops if the dimensions, names, or attributes don't match 
#' @param oldframe The target model frame.
#' @param newframe The model frame being compared to the target
#' @param Check_Rows Logical indicating whether Number of Rows Should be checked (the default). If FALSE,
#' differing number of rows are allowed
#' @return 0 if the two model frames pass the comparison. If not, the function stops.
#' @example inst/examples/Ex_confint.glmb.R
#' @noRd

Compare_Model_Frames<-function(oldframe,newframe,Check_Rows=TRUE){
  
  # Check that dimensions match
  if(Check_Rows==TRUE)  if(!nrow(newframe)==nrow(oldframe)) {
    print("New Model Frame Does Not Have The Same Number of Rows")
    return(FALSE)
  }
  if(!ncol(newframe)==ncol(oldframe))  {
    print("New Model Frame Does Not Have The Same Number of Columns")
    return(FALSE)
  }
  # Check that Column Names Match
  
  for(i in 1:length(names(oldframe))) {
    if(!names(newframe)[i]==names(oldframe)[i])
      {
      print("New Model Frame Has Different Column Names")  
      return(FALSE)
    }
    }
  
  # Check that other attributes match
  
  for(i in 1:length(names(oldframe)))
  {
    if(!is.null(attr.all.equal(newframe[, names(newframe)[i]], oldframe[, names(oldframe)[i]])))
    { 
      not_equal=attr.all.equal(newframe[, names(newframe)[i]], oldframe[, names(oldframe)[i]])
      if(isTRUE(grepl("Attributes", not_equal))) 
      {
        print("Attributes are Different.")
        print(not_equal)
        return(FALSE)
      }
      
    }
  }
  return(TRUE)
}
