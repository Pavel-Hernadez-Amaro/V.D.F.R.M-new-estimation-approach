library(refund)
library(fda)
library(mgcv)
library(SOP)
library(expm)
library(writexl)
library(Tplyr)
library(ggpubr)
library(tidyverse)
library(MASS)
library(gridExtra)
options(error=NULL)


for (index in 1:4) {
  
  print(c("index =",index))
  
  for (jj in 1:3) {
    
    
    print(c("jj =",jj))
    
    if (index%%4==1) {
      
      caso="Normal_NegBin"
      if (jj%%3==1) {
        
        size="_N_100.RData"
        
      }
      if (jj%%3==2) {
        
        size="_N_200.RData"
        
      }
      if (jj%%3==0) {
        
        size="_N_500.RData"
        
      }
    }
    
    if (index%%4==2) {
      
      caso="Normal_Uniform"
      if (jj%%3==1) {
        
        size="_N_100.RData"
        
      }
      if (jj%%3==2) {
        
        size="_N_200.RData"
        
      }
      if (jj%%3==0) {
        
        size="_N_500.RData"
        
      }
    }
    
    if (index%%4==3) {
      
      caso="Poisson_NegBin"
      if (jj%%3==1) {
        
        size="_N_100.RData"
        
      }
      if (jj%%3==2) {
        
        size="_N_200.RData"
        
      }
      if (jj%%3==0) {
        
        size="_N_500.RData"
        
      }
    }
    
    if (index%%4==0) {
      
      caso="Poisson_Uniform"
      if (jj%%3==1) {
        
        size="_N_100.RData"
        
      }
      if (jj%%3==2) {
        
        size="_N_200.RData"
        
      }
      if (jj%%3==0) {
        
        size="_N_500.RData"
        
      }
    }
    

    # nam_X_aux <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Last Simulations/Single Beta/MissD",caso,caso ,sep="/")
    
    nam_X_aux <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/ULTIMAS SIMULACIONES/Variable Domain",caso ,sep="/")
    
    nam_X=paste(nam_X_aux,size,sep="")
    
    aux_size=size
    aux_caso=caso
    aux_index=index
    aux_jj=jj
    
    load(nam_X)
    
    size=aux_size
    caso=aux_caso
    index=aux_index
    jj=aux_jj
    
    Results_Y=Y_std=Y_means=matrix(0,nrow = 4,ncol=8)
    
    if ((index==1 || index==2) & (jj==1 || jj==2)) {
      if (jj==1) {
        Y_means=colMeans(Y_ERRORES)
        Y_matrix=matrix(data = Y_means,nrow = 4,byrow = 1)
        
        for (scenario in 1:8) {
          
          VDFR=scenario
          SOF=1*8+scenario
          FF_VDFR=2*8+scenario
          SB_VDFR=3*8+scenario
          
          Y_std[1,scenario]=sd(Y_ERRORES[,VDFR])
          Y_std[2,scenario]=sd(Y_ERRORES[,SOF])
          Y_std[3,scenario]=sd(Y_ERRORES[,FF_VDFR])
          Y_std[4,scenario]=sd(Y_ERRORES[,SB_VDFR])
          
        }
        
        
      }
      
      
      
      if (jj==2 & index==1) {
        
        Y_means=colMeans(Normal_NegBin_N_200_results$Y_ERRORES)
        
        Y_matrix=matrix(data = Y_means,nrow = 4,byrow = 1)
        
        for (scenario in 1:8) {
          
          VDFR=scenario
          SOF=1*8+scenario
          FF_VDFR=2*8+scenario
          SB_VDFR=3*8+scenario
          
          Y_std[1,scenario]=sd(Normal_NegBin_N_200_results$Y_ERRORES[,VDFR])
          Y_std[2,scenario]=sd(Normal_NegBin_N_200_results$Y_ERRORES[,SOF])
          Y_std[3,scenario]=sd(Normal_NegBin_N_200_results$Y_ERRORES[,FF_VDFR])
          Y_std[4,scenario]=sd(Normal_NegBin_N_200_results$Y_ERRORES[,SB_VDFR])
          Y_values=rbind(as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,c(VDFR)]),as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,c(SOF)]),as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,c(FF_VDFR)]),as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,c(SB_VDFR)]))
          Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",R),rep("SOF",R),rep("FF_VDFR",R),rep("SB_VDFR",R))))
          
        }
        
        
        
        
      }
      
      if (jj==2 & index==2) {
        
        Y_means=colMeans(Normal_Uniform_N_200_results$Y_ERRORES)
        
        Y_matrix=matrix(data = Y_means,nrow = 4,byrow = 1)
        
        for (scenario in 1:8) {
          
          VDFR=scenario
          SOF=1*8+scenario
          FF_VDFR=2*8+scenario
          SB_VDFR=3*8+scenario
          
          Y_std[1,scenario]=sd(Normal_Uniform_N_200_results$Y_ERRORES[,VDFR])
          Y_std[2,scenario]=sd(Normal_Uniform_N_200_results$Y_ERRORES[,SOF])
          Y_std[3,scenario]=sd(Normal_Uniform_N_200_results$Y_ERRORES[,FF_VDFR])
          Y_std[4,scenario]=sd(Normal_Uniform_N_200_results$Y_ERRORES[,SB_VDFR])
          
      }
        
        
        
      }
    }else if(index==3 || index==4){
      
      outliers_Y=outliers_Y_length=list()
      
      q_ric_Y=NULL

      for (scenario in 1:8) {
        
        for (iind in 1:4) {
          
          q_Y=quantile(ERROR_Y[,iind,scenario])
          
          q_ric_Y[iind]=q_Y[4]-q_Y[2] #sqrt(diag(var(Y_ERRORES)))
          
          outliers_Y[[iind]]=which(ERROR_Y[,iind,scenario]>q_Y[4]+1.5*q_ric_Y[iind] | ERROR_Y[,iind,scenario]<q_Y[2]-1.5*q_ric_Y[iind])
          
          outliers_Y_length[[iind]] = length(unique(outliers_Y[[iind]]))
          
          if (is_empty(outliers_Y[[iind]])) { 
            outliers_Y[[iind]]=dim(ERROR_Y)[1]+1
          }
          
        Y_means[iind,scenario]=mean(ERROR_Y[-outliers_Y[[iind]],iind,scenario])
        
        Y_std[iind,scenario]=sd(ERROR_Y[-outliers_Y[[iind]],iind,scenario])
        
        if (Y_means[iind,scenario]>6) {
          
          if (iind==2 & scenario==7) {
            outliers_Y[[iind]] = which(ERROR_Y[,iind,scenario]>12)# MAS DE 10
            outliers_Y_length[[iind]] = length(unique(outliers_Y[[iind]]))
            
            Y_means[iind,scenario]=mean(ERROR_Y[-outliers_Y[[iind]],iind,scenario])
            Y_std[iind,scenario]=sd(ERROR_Y[-outliers_Y[[iind]],iind,scenario])
            
            }else{
          outliers_Y[[iind]] = which(ERROR_Y[,iind,scenario]>5)# MAS DE 5
          outliers_Y_length[[iind]] = length(unique(outliers_Y[[iind]]))
          
          Y_means[iind,scenario]=mean(ERROR_Y[-outliers_Y[[iind]],iind,scenario])
          Y_std[iind,scenario]=sd(ERROR_Y[-outliers_Y[[iind]],iind,scenario])
          }}
        
        }
      }

      Y_matrix=Y_means
      
    }else{
      
      print("ENTRO")
      
      Y_means=colMeans(ERROR_Y)
      Y_matrix=Y_means
      
      for (scenario in 1:8) {
        for (iind in 1:4) {
          
          Y_std[iind,scenario]=sd(ERROR_Y[,iind,scenario])
          
          q_Y=quantile(ERROR_Y[,iind,scenario])
          
          q_ric_Y[iind]=q_Y[4]-q_Y[2] #sqrt(diag(var(Y_ERRORES)))
          
          outliers_Y[[iind]]=which(ERROR_Y[,iind,scenario]>q_Y[4]+1.5*q_ric_Y[iind] | ERROR_Y[,iind,scenario]<q_Y[2]-1.5*q_ric_Y[iind])
          
          outliers_Y_length[[iind]] = length(unique(outliers_Y[[iind]]))
          
          if (is_empty(outliers_Y[[iind]])) { 
            outliers_Y[[iind]]=dim(ERROR_Y)[1]+1
          }
          
        }
        

      }
    }
    
    
    for (i in 1:dim(Y_matrix)[1]) {
      for (j in 1:dim(Y_matrix)[2]) {
        
        
        aux=paste(round(Y_matrix[i,j], digits = 3),round(Y_std[i,j], digits = 3),sep=" (")
        Results_Y[i,j]=paste(aux,NULL,sep = ")")
      }}
    
    aux=paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/ULTIMAS SIMULACIONES/Variable Domain/Excel/",caso,substr(size,1,6),"_Y_mean",".xlsx",sep="")
    
    write_xlsx(as.data.frame(Results_Y),path = aux )
    
    
  }
}