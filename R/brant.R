brant <- function(model,by.var=F){
  m_model <- model$call
  if (is.matrix(eval.parent(m_model$data))) 
    m_model$data <- as.data.frame(data)
  m_model[-which(names(m_model) %in% c("", "formula", "data"))] <- NULL
  m_model[[1L]] <- quote(stats::model.frame)
  m_model <- eval.parent(m_model)
  Terms <- attr(m_model, "terms")
  x <- model.matrix(Terms, m_model)
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)
  x <- x[, -xint, drop = FALSE]
  y <- as.numeric(model.response(m_model))
  x.variables = names(m_model)[-1]
  temp.data = data.frame(m_model, y)
  if(grepl(":",paste0(colnames(x), collapse = "")) & by.var){
    by.var = FALSE
    warning("by.var = TRUE currently not supported for interactions, setting by.var to FALSE")
  }
  
  x.factors = c()
  for (name in x.variables) {
    if (!is.numeric(m_model[, name])) {
      x.factors = c(x.factors, name)
    }
  }
  if (length(x.factors) > 0) {
    tab = table(data.frame(temp.data$y, m_model[,x.factors]))
    count0 = sum(tab == 0)
  }else {
    count0 = 0
  }
  
  
  J = max(y,na.rm=T)
  K = length(coef(model))
  for(m in 1:(J-1)){
    temp.data[[paste0("z",m)]] = ifelse(y>m,1,0)
  }
  binary.models = list()
  beta.hat = matrix(NA,nrow=J-1,ncol=K+1,byrow=T)
  var.hat = list()
  for(m in 1:(J-1)){
    mod = glm(paste0("z",m," ~ ",as.character(formula(model)[3])),data=temp.data, family="binomial")
    binary.models[[paste0("model",m)]] = mod
    beta.hat[m,] = coef(mod)
    var.hat[[m]] = vcov(mod)
  }
  
  X = cbind(1, x)
  tau = matrix(model$zeta,nrow=1,ncol=J-1,byrow=T)
  pi.hat = matrix(NA,nrow=length(model$model[,1]),ncol=J-1,byrow=T)
  for(m in 1:(J-1)){
    pi.hat[,m] = binary.models[[m]]$fitted.values
  }
  
  
  varBeta = matrix(NA,nrow = (J-1)*K, ncol = (J-1)*K)
  for(m in 1:(J-2)){
    for(l in (m+1):(J-1)){
      Wml = Matrix::Diagonal(x=pi.hat[,l] - pi.hat[,m]*pi.hat[,l])
      Wm = Matrix::Diagonal(x=pi.hat[,m] - pi.hat[,m]*pi.hat[,m])
      Wl = Matrix::Diagonal(x=pi.hat[,l] - pi.hat[,l]*pi.hat[,l])
      Xt = t(X)
      varBeta[((m-1)*K+1):(m*K),((l-1)*K+1):(l*K)] = as.matrix((solve(Xt %*% Wm %*% X)%*%(Xt %*% Wml %*% X)%*%solve(Xt %*% Wl %*% X))[-1,-1])
      varBeta[((l-1)*K+1):(l*K),((m-1)*K+1):(m*K)] = varBeta[((m-1)*K+1):(m*K),((l-1)*K+1):(l*K)]
    }
  }
  
  betaStar = c()
  for(m in 1:(J-1)){
    betaStar = c(betaStar,beta.hat[m,-1])
  }
  for(m in 1:(J-1)){
    varBeta[((m-1)*K+1):(m*K),((m-1)*K+1):(m*K)] = var.hat[[m]][-1,-1]
  }
  
  I = diag(1,K)
  E0 = diag(0,K)
  for(i in 1:(J-2)){
    for(j in 1:(J-1)){
      if(j == 1){
        temp = I
      }else if(j == i+1){
        temp = cbind(temp,-I)
      }else{
        temp = cbind(temp,E0)
      }
    }
    if(i==1){
      D = temp
    }else{
      D = rbind(D,temp)
    }
  }
  X2 = t(D%*%betaStar) %*% solve(D %*% varBeta %*% t(D)) %*% (D %*% betaStar)
  df.v = (J-2)*K
  
  if(by.var){
    combinations = getCombiCoefs(model)
    for(v in unique(combinations$var)){
      k = subset(combinations,var==v)$i
      s = c()
      df.v.temp = 0
      for(e in k){
        s = c(s,seq(from=e,to=K*(J-1),by=K))
        df.v.temp = df.v.temp + J-2
      }
      s = sort(s)
      Ds = D[,s]
      if (!is.null(dim(Ds))){
        Ds = Ds[which(!apply(Ds == 0, 1, all)), ]
      }
      if(!is.null(dim(Ds)))
        X2 = c(X2,t(Ds%*%betaStar[s]) %*% solve(Ds %*% varBeta[s,s] %*% t(Ds)) %*% (Ds %*% betaStar[s]))
      else
        X2 = c(X2,t(Ds%*%betaStar[s]) %*% solve(Ds %*% varBeta[s,s] %*% t(t(Ds))) %*% (Ds %*% betaStar[s]))
      df.v = c(df.v,df.v.temp)
    }
  }else{
    for(k in 1:K){
      s = seq(from=k,to=K*(J-1),by=K)
      Ds = D[,s]
      if (!is.null(dim(Ds))){
        Ds = Ds[which(!apply(Ds == 0, 1, all)), ]
      }
      if(!is.null(dim(Ds)))
        X2 = c(X2,t(Ds%*%betaStar[s]) %*% solve(Ds %*% varBeta[s,s] %*% t(Ds)) %*% (Ds %*% betaStar[s]))
      else
        X2 = c(X2,t(Ds%*%betaStar[s]) %*% solve(Ds %*% varBeta[s,s] %*% t(t(Ds))) %*% (Ds %*% betaStar[s]))
      df.v = c(df.v,J-2)
    }
  }
  
  result.matrix = print.testresult(model,X2,df.v,by.var)
  if(count0!=0){
    warning(paste0(count0," combinations in table(dv,ivs) do not occur. Because of that, the test results might be invalid."))
  }
  invisible(result.matrix)
}
