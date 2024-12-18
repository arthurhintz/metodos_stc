source("funcoes.R")

#MPS UMOC ----------------------------------------------------------------------

est_MPS_MW<-function(y,method=c("NR","BFGS","CG","NM","SANN","BFGSR")){
  
  #Start Point an Limits For the Parameters
  p<-1
  inf<-rep(0.1,length.out=p)
  sup<-rep(10,length.out=p)
  k1<-2*p
  nsp<-5
  mint<-100
  
  mod.gosolnp<-try(Rsolnp::gosolnp(fun=MPS_Mw,n.restarts=nsp,LB=inf,UB=sup,
                                   n.sim=mint,y=y,m.optim=-1,control=list(trace=F)),T)
  #-----------------------------------------------------------------------------
  tmp<-test.fun.gosolnp(mod.gosolnp)
  if(tmp==TRUE){
    start.teta<-try(mod.gosolnp$pars,T)
    #------Estimação------------------------------------------------------------
    mod<-try(optim(par=start.teta,fn=MPS_Mw,method=method,hessian=TRUE,
                   control=list(fnscale=-1),y=y,m.optim=1.0),T)
    #------Testes estimativas e variância---------------------------------------
    tmp2<-test.fun(mod)
    if(class(tmp2)=="numeric"){coef<-tmp2}
    else{coef<-rep(NA,k1)}
  }#Fim IF gosolnp
  else{coef<-rep(NA,k1)}
  return(coef)
}

#Valores Positivo-----
is.positive<-function(a){
  #todas as posições do vetor são positivas
  k<-length(a)
  tmp<-sum(a>0)
  return(k==tmp)
}

#Testa optim----
test.fun<-function(object){
  if(class(object)=="list"){
    if(object$convergence == 0){
      parameters<-try(object$par,T)
      hess<-try(object$hessian,T)
      var.coef<-try(diag(solve(-hess)),T)
      if(is.numeric(parameters)==TRUE){
        if(is.numeric(var.coef)==TRUE){
          if(is.positive(var.coef)==TRUE){
            z<-c(parameters,var.coef)
            return(z)
          }else{return(FALSE)}
        }else{return(FALSE)}
      }else{return(FALSE)}
    }else{return(FALSE)}
  }else{return(FALSE)}
}

#Testa Gosolnp----
test.fun.gosolnp<-function(object){
  if(class(object)=="list"){
    parameters<-try(object$pars,T)
    if(is.numeric(parameters)==TRUE){return(TRUE)}
    else{return(FALSE)}
  }
}

#test.gosolnp----
# test.gosolnp<-function(object){
#   if(class(object)=="list"){
#     parameters<-try(object$pars,T)
#     hess<-try(object$hessian,T)
#     #Não tem o menos pois a função gosolnp
#     #encontra o ponto de mínimo e fazemos isso
#     #fornecendo -log-vero
#     var.coef<-try(diag(solve(hess)),T)
#     if(is.numeric(parameters)==TRUE){
#       if(is.numeric(var.coef)==TRUE){
#         if(is.positive(var.coef)==TRUE){
#           z<-c(parameters,var.coef)
#           return(z)
#         }else{return(FALSE)}
#       }else{return(FALSE)}
#     }else{return(FALSE)}
#   }else{return(FALSE)}
# }

#Retira as linhas NA----
which.NA<-function(x){
  x<-data.frame(x)
  lines<-GLDEX::which.na(x[,1])
  tmp<-length(lines)
  if(tmp>0){
    y<-x[-lines,]
  }
  else{y<-x}
  return(y)
}
# Paralelismo----
progresso<-function(iterations){
  iterations <- RS  # used for the foreach loop  
  
  pb <- progress::progress_bar$new(
    format = "[:bar] :elapsed | Faltam: :eta",
    total = iterations,    # 100 
    width = 60)
  
  progress_letter <- rep(LETTERS[1:1], 1)  # token reported in progress bar
  
  # allowing progress bar to be used in foreach 
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  } 
  
  opts <- list(progress = progress)
  return(opts)
}
#-------------------------------------------------------------------------------





