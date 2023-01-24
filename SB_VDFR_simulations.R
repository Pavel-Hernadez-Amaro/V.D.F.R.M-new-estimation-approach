library(refund)
library(fda)
library(mgcv)
library(SOP)
library(expm)
library(writexl)
library(Tplyr)
library(ggpubr)
library(tidyverse)
library(gridExtra)
options(error=NULL)

# library(foreach)
# library(doParallel)
########### Here we generate the data and set some fixed parameters

N=200 # NUMBER OF SUBJECTS

J=100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

k=10  # SIZE OF THE TEST DATA SET 

fold=N/k

R=100 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY 

case=4 # THIS IS THE TOTAL NUMBER OF SCENARIOS SIMULATED (THE NUMBER OF BETAS *2)

registrated_points=function(x) {
  J.i <- sum(!is.na(x))
  y <- approx(x=seq(0,1,length=J.i), y=x[1:J.i],
              xout=seq(0,1,length=J))$y
  y
}

time_SB_VDFR=time_Goldsmith=array(dim = c(R,fold,case))
Times=array(dim=c(R,2,case))

# # HERE WE GENERATE TO EMPTY ARRAYS THAT WILL BE THE ESTIMATED FUNCTIONAL COEFFICIENTS
Beta_SB_VDFR=array(dim =c(k,J,R,case)) # THE DIMENSIONS ARE DIFFERENTS TO THE FF-VDFR SIMULATIONS BECAUSE 
                                       # WE DO NOT NEED ONE BETA PER SUBJECT BUT ONLY ONE BETA


Beta_SB_VDFR_2=Beta_Goldsmith_2=array(dim =c(J,R,case))

Beta_SB_VDFR_test=array(dim =c(1,J,R,case))             
Beta_Goldsmith=array(dim =c(k,J,R,case))

##

# EMPTY ARRAYS WHERE THE DOMAINS OF EVERY SUBJECTS IN EVERY ITERATION WILL BE STORE
M_it=array(dim=c(R,N,case))
M_it_test=array(dim=c(R,fold,case))
M_it_train=array(dim=c(R,(k-1)*fold,case))
#

# EMPTY ARRAYS WHERE THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
y_Goldsmith=array(dim=c(R,fold,case))
y_SB_VDFR=array(dim=c(R,fold,case))

Y_group_SB_nu=Y_group_Goldsmith_nu=array(dim=c(R,N,case))
Y_group_SB=Y_group_Goldsmith=array(dim=c(R,N,case))
#



# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED COEFFICIENT FUNCTION WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO

error_SB_VDFR=error_Goldsmith=matrix(0,nrow =k,ncol = 1)


B_ERROR_SB_VDFR=B_ERROR_Goldsmith=array(dim = c(R,case))

B_ERROR_SB_VDFR_2=B_ERROR_Goldsmith_2=array(dim = c(R,case))

B_ERROR_Goldsmith_test=array(dim = c(R,case))
B_ERROR_SB_VDFR_test=array(dim = c(R,case))
#


# THIS ARE THE N/K ERRORS IN EACH GROUPS

error_group_Goldsmith_nu=error_group_SB_VDFR_nu=error_group_SB_VDFR=error_group_Goldsmith=array(dim=c(R,k,case))


# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO

Y_ERROR_Goldsmith_nu=array(dim = c(R,case))
Y_ERROR_SB_VDFR_nu=array(dim = c(R,case))

Y_ERROR_Goldsmith=array(dim = c(R,case))
Y_ERROR_SB_VDFR=array(dim = c(R,case))
#

start=proc.time() 

for (iter_out in 1:case) {
  
  print(c("case = ",iter_out))
  
  for (iter in 1:R) {
    
print(c("iter = ",iter))
  
set.seed(1000+iter) 

# M = round(runif(N,10,J),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS

M = rnegbin(N,24,1)

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

M=sort(M)

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

# HERE WE REGISTRATE THE OBSERVED POINTS OF THE CURVES

X_reg = t(apply(X_se, 1,registrated_points)) # THESE ARE THE REGISTRATED POINTS

# SOME SAVE CHECKS FOR UNWANTED NAs
for (i in 1:T) {
  if (length(which(is.na(X_s[,i])))==N) {
    print(c("iteracion",l,"columna",i))
    stop("Hay columnas con NA")  
  }}




###### HERE WE GENERATE THE TRUE Beta COEFFICIENT AND THE RESPONSE VARIABLE

Beta=array(dim = c(2,T))  
nu=y=rep(0,N)

# TRUE FUNCTIONAL COEFFICIENTS

Beta[1,1:T]=((10*t/T)-5)/50
Beta[2,1:T]=-1*(5-40*((t/T)-0.5)^2)/50
#

dim_beta=dim(Beta)[1]

for (i in 1:N) {
  
  # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)  
  
  if (iter_out==dim_beta*2) {
    nu[i]=sum(X_se[i,]*Beta[2,],na.rm = 1)/M[i] # NOISY
  }
  if (iter_out<=dim_beta) {
    nu[i]=sum(X_s[i,]*Beta[iter_out,],na.rm = 1)/M[i] #NOT NOISY
  }
  if(iter_out>dim_beta & iter_out<dim_beta*2){
    nu[i]=sum(X_se[i,]*Beta[iter_out%%dim_beta,],na.rm = 1)/M[i] # NOISY
  }
  
}


# y=nu+rnorm(N,sd = 1) # ADDING NOISE TO THE GAUSSIAN MODEL

y=rpois(N,exp(nu)) # POISSON MODEL

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
  
##########
  
c1=25
c2=25
c3=25

# MODEL ESTIMATION BY OUR APPROACH

start_SB_VDFR=proc.time() 

BB_SB_VDFR=Data2B_simpson_no_vc(X_train, M_train, nbasis=c(c1,c2),sub = 25, lim =c(min(M),max(M)))
E_SB_VDFR=B2XZG_1d(BB_SB_VDFR$B,c=c(c2))
res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y_train, family = poisson())

end_SB_VDFR=proc.time()
time_SB_VDFR[iter,group,iter_out]=end_SB_VDFR[1]-start_SB_VDFR[1]


start_Goldsmith=proc.time()

fit_Goldsmith <- pfr(y_train ~ lf(X_reg_train, bs = "ps",k=25,presmooth = "bspline",presmooth.opts = list(nbasis=25)), family = poisson())#, offset = log((M_train)/365)      

end_Goldsmith=proc.time()
time_Goldsmith[iter,group,iter_out]=end_Goldsmith[1]-start_Goldsmith[1]

#################### HERE ENDS THE ESTIMATION OF THE MODELS

#################### HERE WE CALCULATE THE ESTIMATION ERRORS

  # HERE WE CALCULATE THE ERROR OF THE ESTIMATED FUNCTIONAL COEFFICIENTS
  
  if (iter_out==dim_beta*2) {
    True_Beta=Beta[dim_beta,]
  }
  if (iter_out<=dim_beta) {
    True_Beta=Beta[iter_out,]
  }
  if(iter_out>dim_beta & iter_out<dim_beta*2){
    True_Beta=Beta[iter_out%%dim_beta,]
  }
  

prod_SB=BB_SB_VDFR$Phi$B

  if (group==10) {
    
    Beta_Goldsmith[group,1:M_it_train[iter,(k-1)*fold,iter_out],iter,iter_out]=coef(fit_Goldsmith, n=M_it_train[iter,(k-1)*fold,iter_out])$value
    Beta_SB_VDFR[group,1:M_it_train[iter,(k-1)*fold,iter_out],iter,iter_out]=as.vector(prod_SB %*% res_SB_VDFR$theta) 
    #
    
  }else{
    
    Beta_Goldsmith[group,1:T,iter,iter_out]=coef(fit_Goldsmith, n=T)$value
    Beta_SB_VDFR[group,1:T,iter,iter_out]=as.vector(prod_SB %*% res_SB_VDFR$theta) 
    }
#
    

error_SB_VDFR[group,1]=sum((True_Beta-Beta_SB_VDFR[group,1:T,iter,iter_out])^2,na.rm = 1)/J
error_Goldsmith[group,1]=sum((True_Beta-Beta_Goldsmith[group,1:T,iter,iter_out])^2,na.rm = 1)/J

############# HERE WE CALCULATE THE ERROR OF THE ESTIMATED RESPONSE VARIABLE


BB_test_SB_VDFR=Data2B_simpson_no_vc(X_test, M_test, nbasis=c(c1,c2),sub = 25,lim = c(min(M),max(M))) # WE GENERATE THE BASIS OF THE TEST DATA SET

y_SB_VDFR[iter,,iter_out] = BB_test_SB_VDFR$B %*% res_SB_VDFR$theta # ESTIMATED REPSONSE VARIABLE USING SB_VDFR

#

for (j in 1:fold) {

y_Goldsmith[iter,j,iter_out] = sum(X_test[j,1:M_it_test[iter,j,iter_out]]*Beta_Goldsmith[group,1:M_it_test[iter,j,iter_out],iter,iter_out],na.rm = 1)/M_it_test[iter,j,iter_out]
  
}

#

# ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF NORMAL RESPONSE

# error_group_SB_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_SB_VDFR[iter,,iter_out])^2)/fold)
# error_group_Goldsmith[iter,group,iter_out]=sqrt(sum((y_test-y_Goldsmith[iter,,iter_out])^2)/fold)

# ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF POISSON RESPONSE

error_group_Goldsmith[iter,group,iter_out]=sqrt(sum((y_test-exp(y_Goldsmith[iter,,iter_out])))^2)/fold
error_group_SB_VDFR[iter,group,iter_out]=sqrt(sum((y_test-exp(y_SB_VDFR[iter,,iter_out])))^2)/fold

error_group_Goldsmith_nu[iter,group,iter_out]=sqrt(sum((nu_test-((y_Goldsmith[iter,,iter_out])))^2)/fold)
error_group_SB_VDFR_nu[iter,group,iter_out]=sqrt(sum((nu_test-((y_SB_VDFR[iter,,iter_out])))^2)/fold)

Y_group_Goldsmith[iter,current_group,iter_out]=exp(y_Goldsmith[iter,,iter_out])
Y_group_SB[iter,current_group,iter_out]=exp(y_SB_VDFR[iter,,iter_out])

Y_group_Goldsmith_nu[iter,current_group,iter_out]=y_Goldsmith[iter,,iter_out]
Y_group_SB_nu[iter,current_group,iter_out]=y_SB_VDFR[iter,,iter_out]


  } # HERE END THE FOR OF group = 1:fold
  

Times[iter,,iter_out]=cbind(mean(time_SB_VDFR[iter,,iter_out]),mean(time_Goldsmith[iter,,iter_out]))

#### THIS SECTION OF CODE IS FOR CALCULATE THE AMSE ERROR NOT USING K-fold: SUFIX _2 TO ALL ERRORS AND COEFFICIENTS THIS IS THE ONE INCLUDED IN THE PAPER

BB_SB_VDFR=Data2B_simpson_no_vc(X_se, M, nbasis=c(c1,c2),sub = 25, lim =c(min(M),max(M)))
E_SB_VDFR=B2XZG_1d(BB_SB_VDFR$B,c=c(c2)) 
res_SB_VDFR=XZG2theta(X = E_SB_VDFR$X, Z = E_SB_VDFR$Z, G = E_SB_VDFR$G, T = E_SB_VDFR$T, y = y, family = poisson())

fit_Goldsmith = pfr(y ~ lf(X_reg, bs = "ps",k=25,presmooth = "bspline",presmooth.opts = list(nbasis=25)), family = poisson())

# HERE WE CALCULATE THE ERROR OF THE ESTIMATED FUNCTIONAL COEFFICIENTS

if (iter_out==dim_beta*2) {
  True_Beta=Beta[dim_beta,]
}
if (iter_out<=dim_beta) {
  True_Beta=Beta[iter_out,]
}
if(iter_out>dim_beta & iter_out<dim_beta*2){
  True_Beta=Beta[iter_out%%dim_beta,]
}

prod_SB=BB_SB_VDFR$Phi$B

Beta_Goldsmith_2[1:T,iter,iter_out]=coef(fit_Goldsmith, n=T)$value
Beta_SB_VDFR_2[1:T,iter,iter_out]=as.vector(prod_SB %*% res_SB_VDFR$theta) 


ERROR_Goldsmith_2=sum((True_Beta-Beta_Goldsmith_2[1:T,iter,iter_out])^2,na.rm = 1)/T
ERROR_SB_VDFR_2=sum((True_Beta-Beta_SB_VDFR_2[1:T,iter,iter_out])^2,na.rm = 1)/T

#############

B_ERROR_Goldsmith_2[iter,iter_out]=ERROR_Goldsmith_2
B_ERROR_SB_VDFR_2[iter,iter_out]=ERROR_SB_VDFR_2

#### END OF THE AMSE ERROR USING K-foldç

ERROR_Goldsmith=mean(error_Goldsmith[,1])
ERROR_SB_VDFR=mean(error_SB_VDFR[,1])

#############

  B_ERROR_Goldsmith[iter,iter_out]=ERROR_Goldsmith
  B_ERROR_SB_VDFR[iter,iter_out]=ERROR_SB_VDFR

  Y_ERROR_Goldsmith[iter,iter_out]= mean(error_group_Goldsmith[iter,,iter_out])
  Y_ERROR_SB_VDFR[iter,iter_out]= mean(error_group_SB_VDFR[iter,,iter_out])

  Y_ERROR_Goldsmith_nu[iter,iter_out]= mean(error_group_Goldsmith_nu[iter,,iter_out])
  Y_ERROR_SB_VDFR_nu[iter,iter_out]= mean(error_group_SB_VDFR_nu[iter,,iter_out])
  
}
}

end=proc.time()
time=end[1]-start[1]

time/60/60


B_ERRORES=data.frame(B_ERROR_SB_VDFR,B_ERROR_Goldsmith)

B_ERRORES_2=data.frame(B_ERROR_SB_VDFR_2,B_ERROR_Goldsmith_2)

Y_ERRORES=data.frame(Y_ERROR_SB_VDFR,Y_ERROR_Goldsmith)
NU_ERRORES=data.frame(Y_ERROR_SB_VDFR_nu,Y_ERROR_Goldsmith_nu)
