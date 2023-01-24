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
library(abind)
options(error=NULL)

library(foreach)
library(doParallel)

used_cores=detectCores()-4

registerDoParallel(cores = 4)

########### Here we generate the data and set some fixed parameters

N=500 # NUMBER OF SUBJECTS

J=100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

k=10  # SIZE OF THE TEST DATA SET 

fold=N/k

R=100 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY 

case=8 # THIS IS THE TOTAL NUMBER OF SCENARIOS SIMULATED 

registrated_points=function(x) {
  J.i <- sum(!is.na(x))
  y <- approx(x=seq(0,1,length=J.i), y=x[1:J.i],
              xout=seq(0,1,length=J))$y
  y
}

time_no_vc=time_SOP=time_Gellar=time_Goldsmith=array(dim = c(k,R,case))


# # HERE WE GENERATE TO EMPTY ARRAYS THAT WILL BE THE TRUE FUNCTIONAL COEFFICIENTS
Beta_h=Beta_h_no_vc=array(dim =c((k-1)*fold,J,R))
Beta_h_test=array(dim =c(fold,J,R))
# #

# EMPTY LISTS THAT WILL BE THE ESTIMATED FUNCTIONAL COEFFICIENTS USING THE Gellar APPROACH
Beta_G=list() 
Beta_G_te=list()
#

# EMPTY ARRAYS WHERE THE DOMAINS OF EVERY SUBJECTS IN EVERY ITERATION WILL BE STORE
M_it=array(dim=c(R,N,case))
M_it_test=array(dim=c(R,fold,case))
M_it_train=array(dim=c(R,(k-1)*fold,case))
#

# EMPTY ARRAYS WHERE THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
y_h_ellos=array(dim=c(R,fold,case))
y_h_ellos_te=array(dim=c(R,fold,case))
y_h_sop=array(dim=c(R,fold,case))
y_h_no_vc=array(dim=c(R,fold,case))
y_h_adaptive=array(dim=c(R,fold,case))
#


# THIS ARE THE N/K ERRORS IN EACH GROUPS
error_group_FF_VDFR=error_group_SB=error_group_Carmen=error_group_VDFR=array(dim=c(R,k,case))

error_group_FF_VDFR_nu=error_group_SB_nu=error_group_Carmen_nu=error_group_VDFR_nu=array(dim=c(R,k,case))

# THIS ARE THE N/K RESPONSE VARIABLE (y and nu) IN EACH GROUPS

Y_group_FF_VDFR=Y_group_SB=Y_group_Carmen=Y_group_VDFR=array(dim=c(R,N,case))

Y_group_FF_VDFR_nu=Y_group_SB_nu=Y_group_Carmen_nu=Y_group_VDFR_nu=array(dim=c(R,N,case))

ERROR_B=array(dim = c(R,2,case))

Y_HAT=NU_HAT=Y_hat=nu_hat= array(dim = c(R,N,4*case))
  
ERROR_Y=ERROR_Y_nu=Times=array(dim = c(R,4,case))

start=proc.time() 

for (iter_out in 1:case) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE 
  
  print(c("case = ",iter_out))
  
  fold_results=foreach (iter = 1:R, .combine = rbind, .packages = c("refund","MASS","abind")) %dopar% {
    
    print(c("iter = ",iter))
    
    set.seed(1000+iter) 
    
  # M = round(runif(N,10,J),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
    
    M = rnegbin(N,24,1) # Para 1000 poner (240,2)
    
    if (max(M)>J) {
      # print("si")
      M[which(M>J)]=J
    }
    
    
    if (min(M)<=10) {
      # print("si")
      # M[which(M<=10)]=round(runif(length(which(M<=10)),31,max(M)))
      M[which(M<=10)]=10
    }
    
    T=max(M)
    
    t=1:T
    
    M=sort(M) # WE SORT THE DATA WITHOUT LOSS OF GENERALITY
    
    M_it[iter,,iter_out]=M # WE STORE THE DOMAINS FOR EVERY ITERATION AND SCENARIO
    
    
    ############ HERE WE GENERATE THE FUNCTIONAL DATA
    
    X_se=matrix(NA,N,T) # NOISY
    
    X_s=matrix(NA,N,T) # NOT NOISY
    
    for (i in 1:N) {
      
      u=rnorm(1)
      
      temp=matrix(NA,10,T)
      
      for (k in 1:10) {
        
        v_i1=rnorm(1,0,4/k^2)
        v_i2=rnorm(1,0,4/k^2)
        
        temp[k,1:M[i]]=v_i1*sin(2*pi*k*(1:M[i])/100)+v_i2*cos(2*pi*k*(1:M[i])/100)
      }
      
      B=apply(temp,2,sum)
      
      B=B+u
      
      X_s[i,]=B 
      X_se[i,]=(B)+rnorm(T,0,1) # WE ADD NOISE
      
    }
    
    X_reg = t(apply(X_se, 1,registrated_points)) # THESE ARE THE REGISTRATED POINTS
    
    # SOME SAVE CHECKS FOR UNWANTED NAs
    for (i in 1:T) {
      if (length(which(is.na(X_s[,i])))==N) {
        print(c("iteracion",l,"columna",i))
        stop("Hay columnas con NA")  
      }}
    
    
    
    
    ###### HERE WE GENERATE THE TRUE Beta COEFFICIENT AND THE RESPONSE VARIABLE
    
    Beta=array(dim = c(N,T,4))  
    nu=y=rep(0,N)
    
    for (i in 1:N) {
      
      # TRUE FUNCTIONAL COEFFICIENTS
      
      Beta[i,1:(M[i]),1]=((10*t[1:(M[i])]/M[i])-5)/10
      Beta[i,1:(M[i]),2]=((1-(2*M[i]/T))*(5-40*((t[1:(M[i])]/M[i])-0.5)^2))/10
      Beta[i,1:(M[i]),3]=(5-10*((M[i]-t[1:(M[i])])/T))/10
      Beta[i,1:(M[i]),4]=(sin(2*pi*M[i]/T)*(5-10*((M[i]-t[1:(M[i])])/T)))/10
      #
      
      
      # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)  
      
      if (iter_out==8) {
        nu[i]=sum(X_se[i,]*Beta[i,,4],na.rm = 1)/(M[i]) # NOISY
      }
      if (iter_out<=4) {
        nu[i]=sum(X_s[i,]*Beta[i,,iter_out],na.rm = 1)/(M[i]) #NOT NOISY
      }
      if(iter_out>4 & iter_out<8){
        nu[i]=sum(X_se[i,]*Beta[i,,iter_out%%4],na.rm = 1)/M[i] # NOISY
      }
      
    }
    
    
    # y=nu+rnorm(N,sd = 1) # ADDING NOISE TO THE GAUSSIAN MODEL
    
    y=rpois(N,exp(nu)) # POISSON MODEL ##### CAHNGE THE ERROR EQUATIONS IN LINE 477
    
    ############## HERE ENDS THE GENERATION OF THE DATA  
    
    ############ HERE BEGINS THE ESTIMATION OF THE MODEL
    
    # ASSIGNING THE TRAIN AND TEST SETS 
    
    for (group in 1:k) {
      
      print(c("group = ",group))
      
      current_group=(fold*(group-1)+1):(group*fold)
      
      X_test=X_se[current_group,]
      X_train=X_se[-current_group,]
      
      X_reg_train=X_reg[-current_group,]
      X_reg_test=X_reg[current_group,]
      
      M_test=M[current_group]
      M_train=M[-current_group]
      
      M_it_test[iter,,iter_out]=M_test
      M_it_train[iter,,iter_out]=M_train
      
      y_test=y[current_group]
      y_train=y[-current_group]
      
      nu_test=nu[current_group]
      nu_train=nu[-current_group]
      
      # NUMBER OF BASIS FOR EVERY MARGINAL FUNCTION IN OUR APPROACH
      c1=25
      c2=25
      c3=25
      
      # MODEL ESTIMATION BY OUR APPROACH
      
      start_SOP=proc.time()
      
      BB=Data2B_simpson(X_train, M_train, nbasis=c(c1,c2,c3),sub = 25, lim =c(min(M),max(M)),) # HERE WE TRASFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL
      E=B2XZG(BB$B,c=c(c2,c3)) # HERE THE FUNCTION B2XZG() GENERATE THE NECESSARY MATRICES FOR THE MIXED MODEL
      res=XZG2theta(E$X, E$Z, E$G, E$T, y_train,family = poisson()) # HERE THE MODEL COEFICIENT ARE ESTIMATED.
      
      end_SOP=proc.time()
      time_SOP[group,iter,iter_out]=end_SOP[1]-start_SOP[1]
      #
      
      start_no_vc=proc.time() 
      
      BB_EPOC_no_vc=Data2B_simpson_no_vc(X_train, M_train, nbasis=c(c1,c2),sub = 25)
      E_EPOC_no_vc=B2XZG_1d(BB_EPOC_no_vc$B,c=c(c2))
      res_EPOC_no_vc=XZG2theta(X = E_EPOC_no_vc$X, Z = E_EPOC_no_vc$Z, G = E_EPOC_no_vc$G, T = E_EPOC_no_vc$T, y = y_train, family = poisson())
      
      end_no_vc=proc.time()
      time_no_vc[group,iter,iter_out]=end_no_vc[1]-start_no_vc[1]
      #

      # MODEL ESTIMATION USING GELLAR APPROACH
      
      start_Gellar=proc.time()
      
      fit_Gellar <- pfr(y_train ~ lf.vd(X_train,k=89),family = poisson())
      
      end_Gellar=proc.time()
      time_Gellar[group,iter,iter_out]=end_Gellar[1]-start_Gellar[1]
      #
      
      # MODEL ESTIMATION USING GOLDSMITH APPROACH
      
      start_Goldsmith=proc.time()
      
      fit_Gellar_te <- pfr(y_train ~ lf(X_reg_train, bs = "ps",k=25,presmooth = "bspline",presmooth.opts = list(nbasis=25)),family = poisson())#, offset = log((M_train)/365)
      
      end_Goldsmith=proc.time()
      time_Goldsmith[group,iter,iter_out]=end_Goldsmith[1]-start_Goldsmith[1]
      
      #################### HERE ENDS THE ESTIMATION OF THE MODELS
      
      #################### HERE WE CALCULATE THE ESTIMATION ERRORS
      
      error_2_ellos=error_2_te_ellos=matrix(0,nrow = (k-1)*fold,ncol = 1)
      error_2_sop=matrix(0,nrow = (k-1)*fold,ncol = 1)
      
      error_2_ellos_test=error_2_te_ellos_test=matrix(0,nrow = fold,ncol = 1)
      error_2_sop_test=matrix(0,nrow = fold,ncol = 1)
      
      ERROR_2_sop=0
      ERROR_2_ellos=ERROR_2_te_ellos=0
      
      # ESTIMATED COEFFICIENTS FOR THE GELLAR AND GOLDSMITH APPROACHES 
      
      Beta_G[[iter]]=coef(fit_Gellar, n=c(length(t),length(unique(M))))$value
      Beta_G_te[[iter]]=coef(fit_Gellar_te, n=T)$value
      
      for (j in 1:((k-1)*fold)) { # THESE ARE THE ESTIMATED COEFFICIENTS BY MY APPROACH
        
        prod=as.matrix(kronecker(BB$B_Phi[[j]]$B,t(BB$B_T$B[j,])))
        
        Beta_h[j,1:M_it_train[iter,j,iter_out],iter]=as.vector(prod %*% res$theta) # HERE WE MULTIPLY THE BIDIMENSIONAL BASIS FOR THE COEFFICIENT
        
        # HERE WE SELECT THE CORRECT INDEX FOR THE COEFFICIENT OF THE GELLAR AND GOLDSMITH APPROACHES
        
        ind=which(unique(M_it[iter,,iter_out])== M_it_train[iter,j,iter_out])
        ind_t=max(M_it[iter,,iter_out]) # HERE IS CRITICAL THAT OBSERVATIONS ARE CONSECUTIVES IN 1:MAX(M)
        Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
        Beta_refund=Beta_refund[1:M_it_train[iter,j,iter_out]]
        #
        
        
        # HERE WE CALCULATE THE ERROR OF THE ESTIMATED FUNCTIONAL COEFFICIENTS
        
        if (iter_out<=4) {
          
          True_Beta=Beta[-(current_group),1:M_it_train[iter,j,iter_out],iter_out]
          
          error_2_sop[j,1]=sum((True_Beta[j,]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
          
          error_2_ellos[j,1]=sum((True_Beta[j,]-Beta_refund)^2)
          
        }
        
        if (iter_out>4 & iter_out<8) {
          
          True_Beta=Beta[(-current_group),1:M_it_train[iter,j,iter_out],iter_out%%4]
          
          error_2_sop[j,1]=sum((True_Beta[j,]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
          
          error_2_ellos[j,1]=sum((True_Beta[j,]-Beta_refund)^2)
          
        }
        
        if (iter_out==8) {
          
          True_Beta=Beta[(-current_group),1:M_it_train[iter,j,iter_out],4]
          
          error_2_sop[j,1]=sum((True_Beta[j,]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
          
          error_2_ellos[j,1]=sum((True_Beta[j,]-Beta_refund)^2)

        }
        
        
        
      }
      
      ERROR_2_ellos=sum(error_2_ellos[,1])/(J*(J+1))
      
      ERROR_2_sop=sum(error_2_sop[,1])/(J*(J+1))
      
      
      #############
      
      B_ERROR_2_ellos=ERROR_2_ellos
      
      B_ERROR_2_sop=ERROR_2_sop
      
      ########### END OF ESTIMATION ERRORS
      
      ############# HERE WE CALCULATE THE ERROR OF THE ESTIMATED RESPONSE VARIABLE
      
      BB_test=Data2B_simpson(X_test, M_test, nbasis=c(c1,c2,c3),sub = 25,lim = c(min(M),max(M))) # WE GENERATE THE BASIS OF THE TEST DATA SET
      
      y_h_sop[iter,,iter_out] = BB_test$B %*% res$theta # ESTIMATED REPSONSE VARIABLE USING OUR APPROACH
      
      
      BB_test_no_vc=Data2B_simpson_no_vc(X_test, M_test, nbasis=c(c1,c2),sub = 25,lim = c(min(M),max(M))) # WE GENERATE THE BASIS OF THE TEST DATA SET
      
      y_h_no_vc[iter,,iter_out] = BB_test_no_vc$B %*% res_EPOC_no_vc$theta 
      
      
      # ESTIMATED REPSONSE VARIABLE USING HE GELLAR AND GOLDSMITH APPROACHES
      
      for (j in 1:fold) {
        
        prod=as.matrix(kronecker(BB_test$B_Phi[[j]]$B,t(BB_test$B_T$B[j,])))
        
        Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter]=as.vector(prod %*% res$theta) # HERE WE MULTIPLY THE BIDIMENSIONAL BASIS FOR THE COEFFICIENT
        
        ind=which(unique(M_it[iter,,iter_out])== M_it_test[iter,j,iter_out])
        ind_t=max(M_it[iter,,iter_out]) 
        
        Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
        Beta_refund=Beta_refund[1:M_it_test[iter,j,iter_out]]
        
        Beta_refund_te=Beta_G_te[[iter]][1:M_it_test[iter,j,iter_out]]
        
        y_h_ellos[iter,j,iter_out] = sum(X_test[j,1:M_it_test[iter,j,iter_out]]*Beta_refund,na.rm = 1)/M_it_test[iter,j,iter_out]
        
        y_h_ellos_te[iter,j,iter_out] = sum(X_test[j,1:M_it_test[iter,j,iter_out]]*Beta_refund_te,na.rm = 1)/M_it_test[iter,j,iter_out]
        
        
        # THIS SECTION OF CODE IS FOR THE CASE WHEN THE USER WANTS TO CALCULATE THE ERRORS FOR THE ESTIMATION OF THE FUNCTIONAL COEFFICIENT ONLY IN THE TEST SET
        
        #
        # if (iter_out<=4) {
        #
        #   error_2_sop_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   error_2_sop_ad_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_h_ad_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   # error_2[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen[1:M[j]])^2)
        #   #
        #   # error_2_te[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen_te[1:M[j]])^2)
        #
        #   error_2_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_refund)^2)
        #
        #   error_2_te_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_refund_te)^2)
        #
        # }
        #
        # if (iter_out>4 & iter_out<8) {
        #
        #   error_2_sop_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   error_2_sop_ad_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_h_ad_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   # error_2[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen[1:M[j]])^2)
        #   #
        #   # error_2_te[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen_te[1:M[j]])^2)
        #
        #   error_2_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_refund)^2)
        #
        #   error_2_te_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_refund_te)^2)
        #
        # }
        #
        # if (iter_out==8) {
        #
        #   error_2_sop_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   error_2_sop_ad_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_h_ad_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   # error_2[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen[1:M[j]])^2)
        #   #
        #   # error_2_te[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen_te[1:M[j]])^2)
        #
        #   error_2_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_refund)^2)
        #
        #   error_2_te_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_refund_te)^2)
        #
        # }
        #
        # ERROR_2_ellos_test=sum(error_2_ellos_test[,1])/(J*(J+1))
        #
        # ERROR_2_te_ellos_test=sum(error_2_te_ellos_test[,1])/(J*(J+1))
        #
        # ERROR_2_sop_test=sum(error_2_sop_test[,1])/(J*(J+1))
        #
        # # ERROR_2_sop_ad_test=sum(error_2_sop_ad_test[,1])/(J*(J+1))
        #
        # #############
        #
        # B_ERROR_2_ellos_test[iter,iter_out]=ERROR_2_ellos_test
        #
        # B_ERROR_2_te_ellos_test[iter,iter_out]=ERROR_2_te_ellos_test
        #
        # B_ERROR_2_sop_test[iter,iter_out]=ERROR_2_sop_test
        #
        # # B_ERROR_2_sop_ad_test[iter,iter_out]=ERROR_2_sop_ad_test
        
        
      } # HERE ENDS THE INNNER FOR ITERATION (RUNS ON K)
      
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF NORMAL RESPONSE
      
      # error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos[iter,,iter_out])^2)/fold)
      # error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_sop[iter,,iter_out])^2)/fold)
      # error_group_SB[iter,group,iter_out]=sqrt(sum((y_test-y_h_no_vc[iter,,iter_out])^2)/fold)
      # error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos_te[iter,,iter_out])^2)/fold)
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF POISSON RESPONSE
      
      error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos[iter,,iter_out])))^2)/fold)
      error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos_te[iter,,iter_out])))^2)/fold)
      error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_sop[iter,,iter_out])))^2)/fold)
      error_group_SB[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_no_vc[iter,,iter_out])))^2)/fold)
      
      error_group_VDFR_nu[iter,group,iter_out]=sqrt(sum((nu_test-((y_h_ellos[iter,,iter_out])))^2)/fold)
      error_group_Carmen_nu[iter,group,iter_out]=sqrt(sum((nu_test-((y_h_ellos_te[iter,,iter_out])))^2)/fold)
      error_group_FF_VDFR_nu[iter,group,iter_out]=sqrt(sum((nu_test-((y_h_sop[iter,,iter_out])))^2)/fold)
      error_group_SB_nu[iter,group,iter_out]=sqrt(sum((nu_test-((y_h_no_vc[iter,,iter_out])))^2)/fold)
      
      Y_group_VDFR[iter,current_group,iter_out]=exp(y_h_ellos[iter,,iter_out])
      Y_group_Carmen[iter,current_group,iter_out]=exp(y_h_ellos_te[iter,,iter_out])
      Y_group_FF_VDFR[iter,current_group,iter_out]=exp(y_h_sop[iter,,iter_out])
      Y_group_SB[iter,current_group,iter_out]=exp(y_h_no_vc[iter,,iter_out])
      
      Y_group_VDFR_nu[iter,current_group,iter_out]=y_h_ellos[iter,,iter_out]
      Y_group_Carmen_nu[iter,current_group,iter_out]=y_h_ellos_te[iter,,iter_out]
      Y_group_FF_VDFR_nu[iter,current_group,iter_out]=y_h_sop[iter,,iter_out]
      Y_group_SB_nu[iter,current_group,iter_out]=y_h_no_vc[iter,,iter_out]
      
      
      
    } # HERE END THE FOR OF group = 1:fold
    
    Y_ERROR_2_ellos = mean(error_group_VDFR[iter,,iter_out])
    Y_ERROR_2_te_ellos = mean(error_group_Carmen[iter,,iter_out])
    Y_ERROR_2_sop = mean(error_group_FF_VDFR[iter,,iter_out])
    Y_ERROR_2_no_vc = mean(error_group_SB[iter,,iter_out])
    
    Y_ERROR_2_ellos_nu = mean(error_group_VDFR_nu[iter,,iter_out])
    Y_ERROR_2_te_ellos_nu = mean(error_group_Carmen_nu[iter,,iter_out])
    Y_ERROR_2_sop_nu = mean(error_group_FF_VDFR_nu[iter,,iter_out])
    Y_ERROR_2_no_vc_nu = mean(error_group_SB_nu[iter,,iter_out])
    
    ERROR_Y= cbind(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop,Y_ERROR_2_no_vc)
    
    ERROR_Y_nu= cbind(Y_ERROR_2_ellos_nu,Y_ERROR_2_te_ellos_nu,Y_ERROR_2_sop_nu,Y_ERROR_2_no_vc_nu)

    aux=rbind(Y_group_VDFR[iter,,iter_out],Y_group_Carmen[iter,,iter_out],Y_group_FF_VDFR[iter,,iter_out],Y_group_SB[iter,,iter_out])
    
    aux_2=rbind(Y_group_VDFR_nu[iter,,iter_out],Y_group_Carmen_nu[iter,,iter_out],Y_group_FF_VDFR_nu[iter,,iter_out],Y_group_SB_nu[iter,,iter_out])
    
    
    case=iter_out
    Gellar=case
    Gellar_te=1*8+case
    FF_VDFR=2*8+case
    SB_VDFR=3*8+case

    
    Y_hat[iter,,c(Gellar,Gellar_te,FF_VDFR,SB_VDFR)]=aux
    
    nu_hat[iter,,c(Gellar,Gellar_te,FF_VDFR,SB_VDFR)]=aux_2
    
    ERROR_B=cbind(B_ERROR_2_ellos, B_ERROR_2_sop)
    
    time_Gellar_group= mean(time_Gellar[,iter,iter_out])
    
    time_Goldsmith_group= mean(time_Goldsmith[,iter,iter_out])
    
    time_FF_VDFR_group= mean(time_SOP[,iter,iter_out])
    
    time_SB_VDFR_group= mean(time_no_vc[,iter,iter_out])
    
    TIMES = cbind(time_Gellar_group,time_Goldsmith_group,time_FF_VDFR_group,time_SB_VDFR_group)
    
    return(list(ERROR_B=ERROR_B, ERROR_Y=ERROR_Y, ERROR_Y_nu=ERROR_Y_nu, TIMES=TIMES, Y_hat=Y_hat, nu_hat=nu_hat))
    
  } # HERE ENDS THE MIDDLE FOR ITERATION (RUNS ON R) 
  
  case=iter_out
  Gellar=case
  Gellar_te=1*8+case
  FF_VDFR=2*8+case
  SB_VDFR=3*8+case
  
  for (indx in 1:R) {
    
    
    ERROR_B[indx,,iter_out]=fold_results[,1][[indx]]
    ERROR_Y[indx,,iter_out]=fold_results[,2][[indx]]
    ERROR_Y_nu[indx,,iter_out]=fold_results[,3][[indx]]
    Times[indx,,iter_out]=fold_results[,4][[indx]]
    Y_HAT[indx,,c(Gellar,Gellar_te,FF_VDFR,SB_VDFR)]=fold_results[,5][[indx]][indx,,c(Gellar,Gellar_te,FF_VDFR,SB_VDFR)]
    NU_HAT[indx,,c(Gellar,Gellar_te,FF_VDFR,SB_VDFR)]=fold_results[,6][[indx]][indx,,c(Gellar,Gellar_te,FF_VDFR,SB_VDFR)]
    
  }
  
  
} # HERE ENDS THE OUTTER FOR ITERATION (RUNS ON iter_out)

end=proc.time()
time=end[3]-start[3]

time/60/60


################# BOXPLOTS

B_ERRORES_test=data.frame(B_ERROR_2_ellos_test,B_ERROR_2_te_ellos_test,B_ERROR_2_sop_test,B_ERROR_2_sop_ad_test) #B_ERROR_2,B_ERROR_2_te,

B_ERRORES=data.frame(B_ERROR_2_ellos,B_ERROR_2_te_ellos,B_ERROR_2_sop,B_ERROR_2_sop_ad) #B_ERROR_2,B_ERROR_2_te,

Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop,Y_ERROR_2_no_vc) #,Y_ERROR_2_sop_manual) #Y_ERROR_2,Y_ERROR_2_te,

