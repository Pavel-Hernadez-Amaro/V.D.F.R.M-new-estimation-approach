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

Results_Y=Y_std=matrix(0,nrow = 4,ncol=8)

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

for (index in 3:4) {
  
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
    
    # nam_X_aux <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Last Simulations/Single Beta/MissD",caso,caso ,sep="/")
    
    nam_X_aux <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/ULTIMAS SIMULACIONES/Variable Domain",caso ,sep="/")
    
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
    
    Results_Y=Y_std=matrix(0,nrow = 4,ncol=8)
    
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
          
          Y_values=rbind(as.matrix(Y_ERRORES[,c(VDFR)]),as.matrix(Y_ERRORES[,c(SOF)]),as.matrix(Y_ERRORES[,c(FF_VDFR)]),as.matrix(Y_ERRORES[,c(SB_VDFR)]))
          Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",R),rep("SOF",R),rep("FF_VDFR",R),rep("SB_VDFR",R))))
          
          temp_1=get(temp_Y)
          
          temp_1[[jj]][[scenario]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
            geom_violin()+
            ylab("Error") +
            scale_x_discrete(labels = NULL, breaks = NULL)+
            xlab("")+
            theme_bw() +
            stat_summary(fun=median, geom="point", size=2, color="black")
          
          assign(temp_Y,temp_1)
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
          
          temp_1=get(temp_Y)
          
          temp_1[[jj]][[scenario]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
            geom_violin()+
            ylab("Error") +
            scale_x_discrete(labels = NULL, breaks = NULL)+
            xlab("")+
            theme_bw() +
            stat_summary(fun=median, geom="point", size=2, color="black")
          
          assign(temp_Y,temp_1)
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
          
          Y_values=rbind(as.matrix(Normal_Uniform_N_200_results$Y_ERRORES[,c(VDFR)]),as.matrix(Normal_Uniform_N_200_results$Y_ERRORES[,c(SOF)]),as.matrix(Normal_Uniform_N_200_results$Y_ERRORES[,c(FF_VDFR)]),as.matrix(Normal_Uniform_N_200_results$Y_ERRORES[,c(SB_VDFR)]))
          Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",R),rep("SOF",R),rep("FF_VDFR",R),rep("SB_VDFR",R))))
          
          temp_1=get(temp_Y)
          
          temp_1[[jj]][[scenario]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
            geom_violin()+
            ylab("Error") +
            scale_x_discrete(labels = NULL, breaks = NULL)+
            xlab("")+
            theme_bw() +
            stat_summary(fun=median, geom="point", size=2, color="black")
          
          assign(temp_Y,temp_1)
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
        }
        
        if (index==3) {
        
        if (scenario==1 || scenario==5) {

          outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>4)# MAS DE 5
          outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
          
          outliers_Y[[4]] = which(ERROR_Y[,4,scenario]>4)# MAS DE 5
          outliers_Y_length[[4]] = length(unique(outliers_Y[[4]]))
          
        }

        if (scenario==3) {
          
          outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>8)# MAS DE 5
          outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
          
          outliers_Y[[4]] = which(ERROR_Y[,4,scenario]>7.5)# MAS DE 5
          outliers_Y_length[[4]] = length(unique(outliers_Y[[4]]))
          
        }

        if (scenario==7) {
          
          if (jj==1 || jj==2) {
          
          outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>7.5)# MAS DE 5
          outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
          
          outliers_Y[[4]] = which(ERROR_Y[,4,scenario]>7.5)# MAS DE 5
          outliers_Y_length[[4]] = length(unique(outliers_Y[[4]]))
          }
          
          if (jj==3) {
            
            outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>15)# MAS DE 5
            outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
            
            outliers_Y[[4]] = which(ERROR_Y[,4,scenario]>7.5)# MAS DE 5
            outliers_Y_length[[4]] = length(unique(outliers_Y[[4]]))
          }
          
        }

        if (scenario==2 || scenario==4 || scenario==6 || scenario==8) {
          
          outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>4)# MAS DE 5
          outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
          

        }
}
 
        if (index==4) {
          
          if (scenario==1 || scenario==5) {
            
            outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>4)# MAS DE 5
            outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
            
            outliers_Y[[4]] = which(ERROR_Y[,4,scenario]>4)# MAS DE 5
            outliers_Y_length[[4]] = length(unique(outliers_Y[[4]]))
            
          }
          
          if (scenario==3 || scenario==7) {
            
            outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>5)# MAS DE 5
            outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
            
            outliers_Y[[4]] = which(ERROR_Y[,4,scenario]>5)# MAS DE 5
            outliers_Y_length[[4]] = length(unique(outliers_Y[[4]]))
            
          }
          
          if (scenario==4 || scenario==8) {
            
            outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>5)# MAS DE 5
            outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
            

          }
          
          if (scenario==2 || scenario==6) {
            
            outliers_Y[[2]] = which(ERROR_Y[,2,scenario]>4)# MAS DE 5
            outliers_Y_length[[2]] = length(unique(outliers_Y[[2]]))
            
            
          }
        }
        
        Y_values=rbind(as.matrix(ERROR_Y[-outliers_Y[[1]],1,scenario]),as.matrix(ERROR_Y[-outliers_Y[[2]],2,scenario]),as.matrix(ERROR_Y[-outliers_Y[[3]],3,scenario]),as.matrix(ERROR_Y[-outliers_Y[[4]],4,scenario]))
        Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",R-outliers_Y_length[[1]]),rep("SOF",R-outliers_Y_length[[2]]),rep("FF-VDFR",R-outliers_Y_length[[3]]),rep("SB-VDFR",R-outliers_Y_length[[4]]))))
        
        temp_1=get(temp_Y)
        
        temp_1[[jj]][[scenario]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_Y,temp_1)
      }
      
    }else{
      
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
        
        Y_values=rbind(as.matrix(ERROR_Y[-outliers_Y[[1]],1,scenario]),as.matrix(ERROR_Y[-outliers_Y[[2]],2,scenario]),as.matrix(ERROR_Y[-outliers_Y[[3]],3,scenario]),as.matrix(ERROR_Y[-outliers_Y[[4]],4,scenario]))
        Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",R-outliers_Y_length[[1]]),rep("SOF",R-outliers_Y_length[[2]]),rep("FF-VDFR",R-outliers_Y_length[[3]]),rep("SB-VDFR",R-outliers_Y_length[[4]]))))
        
        temp_1=get(temp_Y)
        
        temp_1[[jj]][[scenario]]=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
          geom_violin()+
          ylab("Error") +
          scale_x_discrete(labels = NULL, breaks = NULL)+
          xlab("")+
          theme_bw() +
          stat_summary(fun=median, geom="point", size=2, color="black")
        
        assign(temp_Y,temp_1)
        
      }
    }
    
    
    for (i in 1:dim(Y_matrix)[1]) {
      for (j in 1:dim(Y_matrix)[2]) {
        
        
        aux=paste(round(Y_matrix[i,j], digits = 3),round(Y_std[i,j], digits = 3),sep=" (")
        Results_Y[i,j]=paste(aux,NULL,sep = ")")
      }}
    
    aux=paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Last Simulations/Variable Domain/Excel/",caso,substr(size,1,6),"_Y_mean",".xlsx",sep="")
    
    # write_xlsx(as.data.frame(Results_Y),path = aux )
    
    
  }
  
  for (j_aux in 1:4) {
    
    # temp_22=get(temp_BB)
    # 
    # temp_22[[j_aux]] = ggarrange(get(temp_B)[[1]][[j_aux]],get(temp_B)[[1]][[j_aux+4]],get(temp_B)[[2]][[j_aux]],get(temp_B)[[2]][[j_aux+4]],get(temp_B)[[3]][[j_aux]],get(temp_B)[[3]][[j_aux+4]],
    #                              ncol = 2, nrow = 3)
    # assign(temp_BB,temp_22)
    # 
    
    temp_11=get(temp_YY)
    
    temp_11[[j_aux]] = ggarrange(get(temp_Y)[[1]][[j_aux]],get(temp_Y)[[1]][[j_aux+4]],get(temp_Y)[[2]][[j_aux]],get(temp_Y)[[2]][[j_aux+4]],get(temp_Y)[[3]][[j_aux]],get(temp_Y)[[3]][[j_aux+4]],
                                 ncol = 2, nrow = 3)
    assign(temp_YY,temp_11)
    
    
  }
  
  
}

# B_plot_Normal_Uniform_arrange_aux[[2]]
Y_plot_Normal_Uniform_arrange_aux[[4]]


# B_plot_Poisson_Uniform_arrange_aux[[2]]
Y_plot_Poisson_Uniform_arrange_aux[[2]] 


# B_plot_Poisson_NegBin_arrange_aux[[2]]
Y_plot_Poisson_NegBin_arrange_aux[[3]]

# B_plot_Normal_NegBin_arrange_aux[[2]]
Y_plot_Normal_NegBin_arrange_aux[[4]]


