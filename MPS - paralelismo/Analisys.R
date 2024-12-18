#-------------------------------------------------------------------------------
#Compila os resultados, calculando média, viés, eqm, assimetria, curtose, IC e TH
#-------------------------------------------------------------------------------
SIM.COMP<-function(coef,par,names.pars=NULL,conf=c(0.95),nsign=c(0.05),asymptotic=TRUE,arquivo=NULL)
{
  coef=output
  par=theta0
  names.pars=names_pars
  conf=c(0.90,0.95,0.99)
  nsign=c(0.01,0.05,0.1)
  arquivo=arquivo
  library(knitr)
  library(xtable)
  #-----------------------------------------------------------------------------
  #nsign<-nível de significância do Teste de Wald
  theta0<-par
  thetaH0<-par#as.vector(rep(0,length(theta0)))
  re<-length(coef[,1])
  ne<-length(coef[1,])/2
  coeff<-coef[,1:ne]
  var.coef<-coef[,(ne+1):length(coef[1,])]
  contIC<-matrix(0,ncol=length(theta0),nrow=length(conf))
  WTest<-matrix(0,ncol=length(theta0),nrow=length(nsign))
  q.norm<-qnorm((1-conf)/2,lower.tail = FALSE)
  #-------Calculation of accuracy measures--------------------------------------
  vmean<-mean(coeff)
  vbias<-vmean-theta0
  vsd<-sd(coeff)
  vmse<-vbias^2+vsd^2
  vac<-moments::skewness(coeff)
  vk<-moments::kurtosis(coeff)
  #-----------IC e TH-----------------------------------------------------------
  for(i in 1:re){
    zstat<-as.numeric(abs((coeff-thetaH0)/sqrt(var.coef)))
    pvalues<-2*(1-pnorm(zstat))
    for(v in 1:length(conf)){
      LS<-coeff+q.norm[v]*sqrt(var.coef)
      LI<-coeff-q.norm[v]*sqrt(var.coef)
      contIC[v]<-contIC[v]+(LI<=theta0 & theta0<=LS)
    }
    for(v in 1:length(nsign)){
      WTest[v]<-WTest[v]+(pvalues<nsign[v])
    }
  }
  tc<-contIC/re
  Wth<-WTest/re
  results<-matrix(0,nrow=6+length(conf)+length(nsign),ncol=length(theta0))
  results<-round(rbind(vmean,vbias,vsd,vmse,vac,vk,tc,Wth),4)
  if(!is.null(names.pars)){colnames(results)<-names.pars}
  P<-c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11")
  if(is.null(names.pars)){colnames(results)<-P[1:ne]}
  NIC<-paste0("IC",conf*100,"%")
  NTH<-paste0("TH",nsign*100,"%")
  rownames(results)<-c("Mean","Bias","SD","MSE","AC","K",NIC,NTH)
  if(!is.null(arquivo)){
    #salva a Tabela de Resultados no arquivo------------------------------------
    #write(c("#---------------------------"),file=arquivo,append=TRUE,ncolumns=1)
    #print(xtable(results,type="latex",digits = 4),file=arquivo,append=TRUE)
    write(c("#---------------------------"),file=arquivo,append=TRUE,ncolumns=1)
    print(xtable(t(results),type="latex",digits = 4),file=arquivo,append=TRUE)
    #---------------------------------------------------------------------------
  }else{
    return(kable(results))
  }
}
#-------------------------------------------------------------------------------
#Retira as linhas da matrix que possuem NA
which.NA<-function(x)
{
  library(GLDEX)
  x<-data.frame(x)
  lines<-which.na(x[,1])
  tmp<-length(lines)
  if(tmp>0){
    y<-x[-lines,]
  }
  else{y<-x}
  return(y)
}
