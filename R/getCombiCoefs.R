getCombiCoefs <- function(model){
  classes = attr(model$terms,"dataClasses")
  factors = ifelse(classes[2:length(classes)]!="numeric",T,F)
  f = i = var = 1
  result = data.frame(i=1:length(coef(model)),var=NA)
  for(factor in factors){
    if(factor){
      n = length(unlist(model$xlevels[f]))
      for(j in 1:(n-1)){
        result[i,"var"] = var
        i = i + 1
      }
      var = var + 1
      f = f + 1
    }else{
      result[i,"var"] = var
      var = var + 1
      i = i + 1
    }
  }
  return(result)
}