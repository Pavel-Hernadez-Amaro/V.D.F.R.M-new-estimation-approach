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

Results_B=Results_Y=Y_mean=Y_std=matrix(0,nrow = 2,ncol=4)
B_mean=B_std=matrix(0,nrow = 2,ncol=4)

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
    
    nam_X_aux <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/ULTIMAS SIMULACIONES/Single Beta/Variable Domain",caso,caso ,sep="/")
    
    nam_X=paste(nam_X_aux,size,sep="")
    
    aux_size=size
    
    aux_caso=caso
    
    load(nam_X)
    
    size=aux_size
    
    caso=aux_caso
    
    Results_Y=Y_std=matrix(0,nrow = 2,ncol=4)
    
    for (iind in 1:(dim(Y_ERRORES)[2])) {
      
      nam_RIC_Y=paste("RIC_case_Y",iind,sep = "_")
      
      q_Y=quantile(Y_ERRORES[,iind])
      
      q_ric_Y[iind]=assign(nam_RIC_Y,q_Y[4]-q_Y[2]) #sqrt(diag(var(Y_ERRORES)))
      
      outliers_Y[[iind]]=which(Y_ERRORES[,iind]>q_Y[4]+1.5*q_ric_Y[iind] | Y_ERRORES[,iind]<q_Y[2]-1.5*q_ric_Y[iind])
      
      
      
      nam_RIC_B=paste("RIC_case_B",iind,sep = "_")
      
      q_B=quantile(B_ERRORES[,iind])
      
      q_ric_B[iind]=assign(nam_RIC_B,q_B[4]-q_B[2]) #sqrt(diag(var(Y_ERRORES)))
      
      outliers_B[[iind]]=which(B_ERRORES[,iind]>q_B[4]+1.5*q_ric_B[iind] | B_ERRORES[,iind]<q_B[2]-1.5*q_ric_B[iind])
      
      
      
    }
    
    for (j_ind in 1:4) {
      
      print(c("j_ind =",j_ind))
      
      SB_VDFR=j_ind
      SOF=1*4+j_ind
      
      
      if (length(as.double(outliers_Y[[SB_VDFR]]))==0 & length(as.double(outliers_Y[[SOF]]))==0) {
        
        Y_mean[1,j_ind]=mean(Y_ERRORES[,SB_VDFR])
        Y_mean[2,j_ind]=mean(Y_ERRORES[,SOF])
        
        Y_std[1,j_ind]=sd(Y_ERRORES[,SB_VDFR])
        Y_std[2,j_ind]=sd(Y_ERRORES[,SOF])
        
        
        
        
      }
      else if(length(as.double(outliers_Y[[SB_VDFR]]))==0 & length(as.double(outliers_Y[[SOF]]))!=0){        
        
        Y_mean[1,j_ind]=mean(Y_ERRORES[,SB_VDFR])
        Y_mean[2,j_ind]=mean(Y_ERRORES[-outliers_Y[[SOF]],SOF])
        
        Y_std[1,j_ind]=sd(Y_ERRORES[,SB_VDFR])
        Y_std[2,j_ind]=sd(Y_ERRORES[-outliers_Y[[SOF]],SOF])
        
      }else if(length(as.double(outliers_Y[[SB_VDFR]]))!=0 & length(as.double(outliers_Y[[SOF]]))==0){
        
        Y_mean[1,j_ind]=mean(Y_ERRORES[-outliers_Y[[SB_VDFR]],SB_VDFR])
        Y_mean[2,j_ind]=mean(Y_ERRORES[,SOF])
        
        Y_std[1,j_ind]=sd(Y_ERRORES[-outliers_Y[[SB_VDFR]],SB_VDFR])
        Y_std[2,j_ind]=sd(Y_ERRORES[,SOF])
        
      }else{
        
        
        Y_mean[1,j_ind]=mean(Y_ERRORES[-outliers_Y[[SB_VDFR]],SB_VDFR])
        Y_mean[2,j_ind]=mean(Y_ERRORES[-outliers_Y[[SOF]],SOF])
        
        Y_std[1,j_ind]=sd(Y_ERRORES[-outliers_Y[[SB_VDFR]],SB_VDFR])
        Y_std[2,j_ind]=sd(Y_ERRORES[-outliers_Y[[SOF]],SOF])
        
      }
      
      
      if (length(as.double(outliers_B[[SB_VDFR]]))==0 & length(as.double(outliers_B[[SOF]]))==0) {
        
        B_mean[1,j_ind]=mean(B_ERRORES[,SB_VDFR])
        B_mean[2,j_ind]=mean(B_ERRORES[,SOF])
        
        B_std[1,j_ind]=sd(B_ERRORES[,SB_VDFR])
        B_std[2,j_ind]=sd(B_ERRORES[,SOF])
        
      }else if(length(as.double(outliers_B[[SB_VDFR]]))==0 & length(as.double(outliers_B[[SOF]]))!=0){        
        
        B_mean[1,j_ind]=mean(B_ERRORES[,SB_VDFR])
        B_mean[2,j_ind]=mean(B_ERRORES[-outliers_B[[SOF]],SOF])
        
        B_std[1,j_ind]=sd(B_ERRORES[,SB_VDFR])
        B_std[2,j_ind]=sd(B_ERRORES[-outliers_B[[SOF]],SOF])
        
      }else if(length(as.double(outliers_B[[SB_VDFR]]))!=0 & length(as.double(outliers_B[[SOF]]))==0){
        
        B_mean[1,j_ind]=mean(B_ERRORES[-outliers_B[[SB_VDFR]],SB_VDFR])
        B_mean[2,j_ind]=mean(B_ERRORES[,SOF])
        
        B_std[1,j_ind]=sd(B_ERRORES[-outliers_B[[SB_VDFR]],SB_VDFR])
        B_std[2,j_ind]=sd(B_ERRORES[,SOF])
        
      }else{
        
        
        B_mean[1,j_ind]=mean(B_ERRORES[-outliers_B[[SB_VDFR]],SB_VDFR])
        B_mean[2,j_ind]=mean(B_ERRORES[-outliers_B[[SOF]],SOF])
        
        B_std[1,j_ind]=sd(B_ERRORES[-outliers_B[[SB_VDFR]],SB_VDFR])
        B_std[2,j_ind]=sd(B_ERRORES[-outliers_B[[SOF]],SOF])
        
      }
      
      
    }
    
    for (i in 1:dim(Y_mean)[1]) {
      for (j in 1:dim(Y_mean)[2]) {
        
        aux=paste(round(Y_mean[i,j], digits = 3),round(Y_std[i,j], digits = 3),sep=" (")
        Results_Y[i,j]=paste(aux,NULL,sep = ")")
      }}
    
    aux_Y=paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/ULTIMAS SIMULACIONES/Single Beta/Variable Domain/Excel/",caso,substr(size,1,6),"_Y_mean",".xlsx",sep="")
    
    write_xlsx(as.data.frame(Results_Y),path = aux_Y )
    
    for (i in 1:dim(B_mean)[1]) {
      for (j in 1:dim(B_mean)[2]) {
        
        
        aux=paste(round(B_mean[i,j], digits = 4),round(B_std[i,j], digits = 3),sep=" (")
        Results_B[i,j]=paste(aux,NULL,sep = ")")
        
        
      }}
    
    aux_B=paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/ULTIMAS SIMULACIONES/Single Beta/Variable Domain/Excel/",caso,substr(size,1,6),"_B_mean",".xlsx",sep="")
    
    write_xlsx(as.data.frame(Results_B),path = aux_B )
    
    
  }}
