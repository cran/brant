print.testresult <- function(model,X2,df.v,by.var) {
  p.values = pchisq(X2,df.v,lower.tail=FALSE) 
  if(by.var){
    var.names = unlist(strsplit(as.character(formula(model))[3],split=" \\+ "))
  }else{
    var.names = names(coef(model))
  }
  # longest name
  longest.char = max(nchar(var.names))
  n.tabs = ceiling(longest.char/7)
  n.tabs = ifelse(n.tabs<2,2,n.tabs)
  cat(paste0(rep("-",28+8*n.tabs),collapse = ""),"\n")
  cat(paste0("Test for",paste0(rep("\t",n.tabs-1),collapse = ""),"X2\tdf\tprobability"),"\n")
  cat(paste0(rep("-",28+8*n.tabs),collapse = ""),"\n")
  cat(paste0("Omnibus",paste0(rep("\t",n.tabs),collapse = ""),round(X2[1],digits=2),"\t",df.v[1],"\t",round(p.values[1],digits=2)))
  cat("\n")
  for(i in 1:length(var.names)){
    name = var.names[i]
    tabs.sub = ceiling(nchar(name)/7)-1
    cat(paste0(name,paste0(rep("\t",n.tabs-tabs.sub),collapse = ""),round(X2[i+1],digits=2),"\t",df.v[i+1],"\t",round(p.values[i+1],digits=2),"\n"))
  }
  cat(paste0(rep("-",28+8*n.tabs),collapse = ""),"\n\n")
  cat("H0: Parallel Regression Assumption holds")
}