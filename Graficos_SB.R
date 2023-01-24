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

outliers_Y=outliers_Y_length=list()

q_ric_Y=NULL

B_plot_Poisson_Uniform_aux=rep(list(list()),3)
Y_plot_Poisson_Uniform_aux=rep(list(list()),3)
B_plot_Poisson_Uniform_arrange_aux=list()
Y_plot_Poisson_Uniform_arrange_aux=list()

B_plot_Poisson_NegBin_aux=rep(list(list()),3)
Y_plot_Poisson_NegBin_aux=rep(list(list()),3)
B_plot_Poisson_NegBin_arrange_aux=list()
Y_plot_Poisson_NegBin_arrange_aux=list()

B_plot_Normal_Uniform_aux=rep(list(list()),3)
Y_plot_Normal_Uniform_aux=rep(list(list()),3)
B_plot_Normal_Uniform_arrange_aux=list()
Y_plot_Normal_Uniform_arrange_aux=list()

B_plot_Normal_NegBin_aux=rep(list(list()),3)
Y_plot_Normal_NegBin_aux=rep(list(list()),3)
B_plot_Normal_NegBin_arrange_aux=list()
Y_plot_Normal_NegBin_arrange_aux=list()


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
    
    temp_B=aux_temp_B=paste("B_plot_",caso,"_aux",sep="")
    temp_Y=aux_temp_Y=paste("Y_plot_",caso,"_aux",sep="")
    
    temp_BB=aux_temp_BB=paste("B_plot_",caso,"_arrange","_aux",sep="")
    temp_YY=aux_temp_YY=paste("Y_plot_",caso,"_arrange","_aux",sep="")
    
    aux_size=size
    aux_jj=jj
    aux_caso=caso
    aux_index=index
    
    nam_X_aux <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/ULTIMAS SIMULACIONES/Single Beta/Variable Domain",caso,caso ,sep="/")
    
    nam_X=paste(nam_X_aux,size,sep="")
    
    load(nam_X)
    
    caso=aux_caso
    size=aux_size
    jj=aux_jj
    index=aux_index
    
    # print(c("ESTE ES EL JJ QUE SE USA",jj))
    
    temp_B=aux_temp_B
    temp_Y=aux_temp_Y
    temp_BB=aux_temp_BB
    temp_YY=aux_temp_YY

    # assign(temp_B,NULL)
    # assign(temp_Y,NULL)
    # assign(temp_BB,NULL)
    # assign(temp_YY,NULL)
    
    # temp_1=temp_2=NULL
    
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
        

        Y_values=rbind(as.matrix(Y_ERRORES[,c(SB_VDFR)]),as.matrix(Y_ERRORES[,c(SOF)]))
        Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("SB_VDFR",R),rep("SOF",R))))
        
        temp_1=get(temp_Y)
        
        temp_1[[jj]][[j_ind]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_Y,temp_1)
        
        
        
      }else if(length(as.double(outliers_Y[[SB_VDFR]]))==0 & length(as.double(outliers_Y[[SOF]]))!=0){        
        

        Y_values=rbind(as.matrix(Y_ERRORES[,c(SB_VDFR)]),as.matrix(Y_ERRORES[-outliers_Y[[SOF]],c(SOF)]))
        Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("SB_VDFR",R),rep("SOF",R-length(unique(outliers_Y[[SOF]]))))))
        
        temp_1=get(temp_Y)
        
        temp_1[[jj]][[j_ind]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_Y,temp_1)
      }else if(length(as.double(outliers_Y[[SB_VDFR]]))!=0 & length(as.double(outliers_Y[[SOF]]))==0){
        

        Y_values=rbind(as.matrix(Y_ERRORES[-outliers_Y[[SB_VDFR]],c(SB_VDFR)]),as.matrix(Y_ERRORES[,c(SOF)]))
        Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("SB_VDFR",R-length(unique(outliers_Y[[SB_VDFR]]))),rep("SOF",R))))
        
        temp_1=get(temp_Y)
        
        temp_1[[jj]][[j_ind]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_Y,temp_1)
        
        
        
      }else{
        
        
        Y_values=rbind(as.matrix(Y_ERRORES[-outliers_Y[[SB_VDFR]],c(SB_VDFR)]),as.matrix(Y_ERRORES[-outliers_Y[[SOF]],c(SOF)]))
        Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("SB_VDFR",R-length(unique(outliers_Y[[SB_VDFR]]))),rep("SOF",R-length(unique(outliers_Y[[SOF]]))))))
        
        temp_1=get(temp_Y)
        
        temp_1[[jj]][[j_ind]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_Y,temp_1)
        
        
      }
      
      
      if (length(as.double(outliers_B[[SB_VDFR]]))==0 & length(as.double(outliers_B[[SOF]]))==0) {
        
        B_values=rbind(as.matrix(B_ERRORES[,c(SB_VDFR)]),as.matrix(B_ERRORES[,c(SOF)]))
        B_ERRORES_DF=data.frame(values=B_values,method=as.factor(c(rep("SB_VDFR",R),rep("SOF",R))))
        
        temp_2=get(temp_B)
        
        temp_2[[jj]][[j_ind]]=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_B,temp_2)
        
        
        
      }else if(length(as.double(outliers_B[[SB_VDFR]]))==0 & length(as.double(outliers_B[[SOF]]))!=0){        
        

        B_values=rbind(as.matrix(B_ERRORES[,c(SB_VDFR)]),as.matrix(B_ERRORES[-outliers_B[[SOF]],c(SOF)]))
        B_ERRORES_DF=data.frame(values=B_values,method=as.factor(c(rep("SB_VDFR",R),rep("SOF",R-length(unique(outliers_B[[SOF]]))))))
        
        temp_2=get(temp_B)
        
        temp_2[[jj]][[j_ind]]=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_B,temp_2)
      }else if(length(as.double(outliers_B[[SB_VDFR]]))!=0 & length(as.double(outliers_B[[SOF]]))==0){
        
    
        B_values=rbind(as.matrix(B_ERRORES[-outliers_B[[SB_VDFR]],c(SB_VDFR)]),as.matrix(B_ERRORES[,c(SOF)]))
        B_ERRORES_DF=data.frame(values=B_values,method=as.factor(c(rep("SB_VDFR",R-length(unique(outliers_B[[SB_VDFR]]))),rep("SOF",R))))
        
        temp_2=get(temp_B)
        
        temp_2[[jj]][[j_ind]]=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_B,temp_2)
        
        
        
      }else{
        
        
        B_values=rbind(as.matrix(B_ERRORES[-outliers_B[[SB_VDFR]],c(SB_VDFR)]),as.matrix(B_ERRORES[-outliers_B[[SOF]],c(SOF)]))
        B_ERRORES_DF=data.frame(values=B_values,method=as.factor(c(rep("SB_VDFR",R-length(unique(outliers_B[[SB_VDFR]]))),rep("SOF",R-length(unique(outliers_B[[SOF]]))))))
        
        temp_2=get(temp_B)
        
        temp_2[[jj]][[j_ind]]=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_B,temp_2)
        
        
      }
      
      
    }
    
    
  }
  
  
  
  
  for (j_aux in 1:2) {
    
    # temp_22=get(temp_BB)
    # 
    # temp_22[[j_aux]] = ggarrange(get(temp_B)[[1]][[j_aux]],get(temp_B)[[1]][[j_aux+2]],get(temp_B)[[2]][[j_aux]],get(temp_B)[[2]][[j_aux+2]],get(temp_B)[[3]][[j_aux]],get(temp_B)[[3]][[j_aux+2]],
    #                              ncol = 2, nrow = 3)
    # assign(temp_BB,temp_22)
    
    
    temp_11=get(temp_YY)
    
    temp_11[[j_aux]] = ggarrange(get(temp_Y)[[1]][[j_aux]],get(temp_Y)[[1]][[j_aux+2]],get(temp_Y)[[2]][[j_aux]],get(temp_Y)[[2]][[j_aux+2]],get(temp_Y)[[3]][[j_aux]],get(temp_Y)[[3]][[j_aux+2]],
                                 ncol = 2, nrow = 3)
    assign(temp_YY,temp_11)
    
    
  }
  
  
}

B_plot_Normal_Uniform_arrange_aux[[2]]
Y_plot_Normal_Uniform_arrange_aux[[2]]


B_plot_Poisson_Uniform_arrange_aux[[2]]
Y_plot_Poisson_Uniform_arrange_aux[[2]] 


B_plot_Poisson_NegBin_arrange_aux[[2]]
Y_plot_Poisson_NegBin_arrange_aux[[2]]

B_plot_Normal_NegBin_arrange_aux[[2]]
Y_plot_Normal_NegBin_arrange_aux[[2]]

