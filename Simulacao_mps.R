#------Pacotes Necessários------------------------------------------------------
{suppressPackageStartupMessages(library(foreach))
  
  #------Carrega as Funções Necessárias-------------------------------------------
  source("funcoes.R")
  source("Analisys.R")
  source("MW_MPS.R")
  
  #-------Geração de Cenários-----------------------------------------------------
  mu <- c(0.2, 0.5, 1, 1.5)
  ns <- data.frame(seq(50, 1000, by = 50))
  df <-merge(ns,mu)
  colnames(df)<-c("n","mu")
  rm(ns,mu)
  
  #------------Initial Parameters for Estimation----------------------------------
  RS<-7000 #Réplicas de Monte Carlo Simuladas
  RA<-5200 #Réplicas de Monte Carlo Analisadas
  RR<-5000 #Réplicas de Monte Carlo Requeridas
}


cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores) #not to overload your computer
doSNOW::registerDoSNOW(cl)

set.seed(1234)
for(k in 1:length(df$mu)){
  
  #-------Alocação de Núcleos-----------------------------------------------------

  #Cenário----
  {
    n <- df$n[k]
    mu <- df$mu[k]
    method<-"BFGS"
    arquivo<-"MPS_MK_simu.txt"
  }
  theta0<-c(mu)
  rm(mu)
  
  #------Geração das amostras-----------------------------------------------------
  sample<-bigstatsr::FBM(nrow=RS,ncol=n)
  #for(j in 1:RS){
  lixo<-foreach(j=1:RS, .packages=c("foreach"),.combine=rbind) %dopar%{
    source("funcoes.R")
    sample[j,]<-as.numeric(rmax(n=n, mu = theta0[1]))
  }
  rm(lixo)
  #-----Estimação dos Parêmetros e Cálculo da Variância---------------------------
  opts<-progresso(iterations = RS)
  estimate <- bigstatsr::FBM(nrow=RS,ncol=2*length(theta0))
  time<-system.time( 
    foreach(j=1:RS, .packages=c("foreach"),.combine=rbind,.options.snow = opts) %dopar%{
      source("MW_MPS.R")
      output0<-as.data.frame(estimate[])
      output0[output0==0]<-NA
      output0<-which.NA(output0)
      l_out<-length(output0$V1)
      if(l_out>=RA){
      }
      if(l_out<=RA-1){
        estimate[j,]<-try(suppressWarnings(est_MPS_MW(y=sample[j,],method=method)))
      }
    }
  )
  rm(l_out, output0)
  #Salvar as estimativas de Monte Carlo
  output<-estimate[]
  output<-which.NA(output)
  output<-output[1:RR,]
  if(n<100){
    write(c("#---------------------------"),file=arquivo,append=TRUE,ncolumns=1)
    save(output,file=paste0("MPS_MW_n0",n,"mu=",theta0[1],".RData"))
    write(c("#n=",n,"mu",theta0[1]),file=arquivo,append=TRUE, ncolumns=12)
    names_pars<-c("mu")
    #Cálculo Estatísticas e tabela latex
    SIM.COMP(coef=output,par=theta0,names.pars=names_pars,conf=c(0.90,0.95,0.99),nsign=c(0.01,0.05,0.1),
             arquivo=arquivo)
  }
  if(n>=100){
    write(c("#---------------------------"),file=arquivo,append=TRUE,ncolumns=1)
    save(output,file=paste0("MPS_MW_n",n,"mu=",theta0[1],".RData"))
    write(c("#n=",n,"mu",theta0[1]),file=arquivo,append=TRUE, ncolumns=12)
    names_pars<-c("mu")
    #Cálculo Estatísticas e tabela latex
    SIM.COMP(coef=output,par=theta0,names.pars=names_pars,conf=c(0.90,0.95,0.99),nsign=c(0.01,0.05,0.1),
             arquivo=arquivo)
  }
  rm(theta0, estimate,output)

}

#-------End Parallelism---------------------------------------------------------
foreach::registerDoSEQ()
parallel::stopCluster(cl)
