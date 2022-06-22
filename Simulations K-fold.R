library(refund)
library(fda)
library(mgcv)
library(SOP)
library(expm)

# ALL SIMULATED DATA AND RESULTS ARE FOUND IN https://drive.google.com/file/d/1P_VlMo4agsdR7MhE9p1THNcSKEqiqvpR/view?usp=sharing 
# AND ARE READY PLOT

########### Here we generate the data and set some fixed parameters

N=100 # NUMBER OF SUBJECTS

J=100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

k=10  # SIZE OF THE TEST DATA SET 

R=100 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY 

case=8 # THIS IS THE TOTAL NUMBER OF SCENARIOS SIMULATED 

#ndb=15 # THIS IS THE NUMBER OF LAMBDAS FOR THE ADPATIVE CASE

test=sort(sample(1:N,k)) # HERE WE SET THE TEST SET INDEX

train=1:N 
train=train[-test] # HERE WE SET THE TRAIN SET INDEX


# # HERE WE GENERATE TO EMPTY ARRAYS THAT WILL BE THE TRUE FUNCTIONAL COEFFICIENTS
Beta_h=Beta_h_ad=array(dim =c(N-k,J,R))
# Beta_h_test=Beta_h_ad_test=array(dim =c(k,J,R))
# #

# EMPTY LISTS THAT WILL BE THE ESTIMATED FUNCTIONAL COEFFICIENTS USING THE Gellar APPROACH
Beta_G=list() 
Beta_G_te=list()
#

# EMPTY ARRAYS WHERE THE DOMAINS OF EVERY SUBJECTS IN EVERY ITERATION WILL BE STORE
M_it=array(dim=c(R,N,case))
M_it_test=array(dim=c(R,k,case))
M_it_train=array(dim=c(R,N-k,case))
#

# EMPTY ARRAYS WHERE THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
y_h_ellos=array(dim=c(R,k,case))
y_h_ellos_te=array(dim=c(R,k,case))
y_h_sop=array(dim=c(R,k,case))
y_h_adaptive=array(dim=c(R,k,case))
#

# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED COEFFICIENT FUNCTION WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
B_ERROR_2_ellos=B_ERROR_2_te_ellos=array(dim = c(R,case))
B_ERROR_2_sop=B_ERROR_2_sop_ad=array(dim = c(R,case))
B_ERROR_2_ellos_test=B_ERROR_2_te_ellos_test=array(dim = c(R,case))
B_ERROR_2_sop_test=B_ERROR_2_sop_ad_test=array(dim = c(R,case))
#

# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
Y_ERROR_2_ellos=Y_ERROR_2_te_ellos=array(dim = c(R,case))
Y_ERROR_2_sop=Y_ERROR_2_sop_ad=Y_ERROR_2_sop_manual=array(dim = c(R,case))
#

start=proc.time() 

for (iter_out in 1:case) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE 
  
  print(c("case =",iter_out))
  
  for (iter in 1:R) {
    
    
    print(iter)
    
    set.seed(1000+iter) 
    
    M = round(runif(N,10,J),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
    
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
      X_se[i,]=B+rnorm(T,0,1) # WE ADD NOISE
      
    }
    
    # SOME SAVE CHECKS FOR UNWANTED NAs
    for (i in 1:T) {
      if (length(which(is.na(X_s[,i])))==N) {
        print(c("iteracion",l,"columna",i))
        stop("Hay columnas con NA")  
      }}
    
    
    # ASSIGNING THE TRAIN AND TEST SETS 
    X_test=X_se[test,]
    X_train=X_se[-test,]
    
    M_test=M[test]
    M_train=M[-test]
    
    M_it_test[iter,,iter_out]=M_test
    M_it_train[iter,,iter_out]=M_train
    
    ###### HERE WE GENERATE THE TRUE Beta COEFFICIENT AND THE RESPONSE VARIABLE

    Beta=array(dim = c(N,T,4))  
    nu=y=rep(0,N)
    
    for (i in 1:N) {
      
      # TRUE FUNCTIONAL COEFFICIENTS

      Beta[i,1:(M[i]),1]=(10*t[1:(M[i])]/M[i])-5
      Beta[i,1:(M[i]),2]=(1-(2*M[i]/T))*(5-40*((t[1:(M[i])]/M[i])-0.5)^2)
      Beta[i,1:(M[i]),3]=5-10*((M[i]-t[1:(M[i])])/T)
      Beta[i,1:(M[i]),4]=sin(2*pi*M[i]/T)*(5-10*((M[i]-t[1:(M[i])])/T))
      #
      
      # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)  
      
      if (iter_out==8) {
        nu[i]=sum(X_se[i,]*Beta[i,,4],na.rm = 1)/M[i] # NOISY
      }
      if (iter_out<=4) {
        nu[i]=sum(X_s[i,]*Beta[i,,iter_out],na.rm = 1)/M[i] #NOT NOISY
      }
      if(iter_out>4 & iter_out<8){
        nu[i]=sum(X_se[i,]*Beta[i,,iter_out%%4],na.rm = 1)/M[i] # NOISY
      }
      
    }
    
    
    y=nu+rnorm(N,sd = 1) # ADDING NOISE TO THE MODEL
    
    y_test=y[test]
    y_train=y[-test]
    
 
    ############## HERE ENDS THE GENERATION OF THE DATA 
    
    ############ HERE BEGINS THE ESTIMATION OF THE MODEL
    
    # NUMBER OF BASIS FOR EVERY MARGINAL FUNCTION IN OUR APPROACH
    c1=25
    c2=25
    c3=25
    
    # MODEL ESTIMATION BY OUR APPROACH
    BB=Data2B_simpson(X_train, M_train, nbasis=c(c1,c2,c3),sub = 25) # HERE WE TRASFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL
    E=B2XZG(BB$B,c=c(c2,c3)) # HERE THE FUNCTION B2XZG() GENERATE THE NECESSARY MATRICES FOR THE MIXED MODEL 
    res=XZG2theta(E$X, E$Z, E$G, E$T, y_train) # HERE THE MODEL COEFICIENT ARE ESTIMATED.
    
    # ADAPTIVE
    # BB_ad=Data2B_simpson_ad(X_train, M_train, nbasis=c(c1,c2,c3),sub = 25, ndb = ndb)
    # E_ad=B2XZG(BB_ad$B,c=c(c2,c3)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
    # res_ad=XZG2theta(E_ad$X, E_ad$Z, E_ad$G, E_ad$T, y_train) # AQUÍ AJUSTO EL MODELO MIXTO Y RECUPERO LOS COEFICIENTES ORIGINALES
    
    # MODEL ESTIMATION USING GELLAR APPROACH
    fit_Gellar <- pfr(y_train ~ lf.vd(X_train,k=89))
    
    # MODEL ESTIMATION USING GOLDSMITH APPROACH
    fit_Gellar_te <- pfr(y_train ~ lf.vd(X_train,transform = "standardized",basistype = "te",k=c(11,8)))
    
    #################### HERE ENDS THE ESTIMATION OF THE MODELS
    
    #################### HERE WE CALCULATE THE ESTIMATION ERRORS
    
    error_2_ellos=error_2_te_ellos=matrix(0,nrow = N-k,ncol = 1)
    error_2_sop=error_2_sop_ad=matrix(0,nrow = N-k,ncol = 1)
    
    error_2_ellos_test=error_2_te_ellos_test=matrix(0,nrow = k,ncol = 1)
    error_2_sop_test=error_2_sop_ad_test=matrix(0,nrow = k,ncol = 1)
    
    ERROR_2_sop=ERROR_2_sop_ad=0
    ERROR_2_ellos=ERROR_2_te_ellos=0
    
    
    # ESTIMATED COEFFICIENTS FOR THE GELLAR AND GOLDSMITH APPROACHES 
    
    Beta_G[[iter]]=coef(fit_Gellar, n=c(length(t),length(unique(M))))$value
    Beta_G_te[[iter]]=coef(fit_Gellar_te, n=c(length(t),length(unique(M))))$value
    
    
    for (j in 1:(N-k)) { # THESE ARE THE ESTIMATED COEFFICIENTS BY MY APPROACH 

      prod=as.matrix(kronecker(BB$B_Phi[[j]]$B,t(BB$B_T$B[j,])))
      
      Beta_h[j,1:M_it_train[iter,j,iter_out],iter]=as.vector(prod %*% res$theta) # HERE WE MULTIPLY THE BIDIMENSIONAL BASIS FOR THE COEFFICIENT
      
      # prod=as.matrix(kronecker(BB_ad$B_Phi[[j]]$B,t(BB_ad$B_T$B[j,])))
      # Beta_h_ad[j,1:M_it_train[iter,j,iter_out],iter]=as.vector(prod %*% res_ad$theta) # Aquí estamos multiplicando la base bidimesional por los theta
      
      
      # HERE WE SELECT THE CORRECT INDEX FOR THE COEFFICIENT OF THE GELLAR AND GOLDSMITH APPROACHES 
      
      ind=which(unique(M_it[iter,,iter_out])== M_it_train[iter,j,iter_out])
      ind_t=max(M_it[iter,,iter_out]) # HERE IS CRITICAL THAT OBSERVATIONS ARE CONSECUTIVES IN 1:MAX(M)
      Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
      Beta_refund=Beta_refund[1:M_it_train[iter,j,iter_out]]
      Beta_refund_te=Beta_G_te[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
      Beta_refund_te=Beta_refund_te[1:M_it_train[iter,j,iter_out]]
      #
      
      # HERE WE CALCULATE THE ERROR OF THE ESTIMATED FUNCTIONAL COEFFICIENTS
      
      if (iter_out<=4) {
        
        error_2_sop[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
        
        # error_2_sop_ad[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out]-Beta_h_ad[j,1:M_it_train[iter,j,iter_out],iter])^2)
        
        error_2_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out]-Beta_refund)^2)
        
        error_2_te_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out]-Beta_refund_te)^2)
        
      }
      
      if (iter_out>4 & iter_out<8) {
        
        error_2_sop[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out%%4]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
        
        # error_2_sop_ad[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out%%4]-Beta_h_ad[j,1:M_it_train[iter,j,iter_out],iter])^2)
        
        error_2_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out%%4]-Beta_refund)^2)
        
        error_2_te_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out%%4]-Beta_refund_te)^2)
        
      }
      
      if (iter_out==8) {
        
        error_2_sop[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],4]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
        
        # error_2_sop_ad[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],4]-Beta_h_ad[j,1:M_it_train[iter,j,iter_out],iter])^2)
        
        error_2_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],4]-Beta_refund)^2)
        
        error_2_te_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],4]-Beta_refund_te)^2)
        
      }
      
      
     
    }
    
    ERROR_2_ellos=sum(error_2_ellos[,1])/(J*(J+1))
    
    ERROR_2_te_ellos=sum(error_2_te_ellos[,1])/(J*(J+1))
    
    ERROR_2_sop=sum(error_2_sop[,1])/(J*(J+1))
    
    # ERROR_2_sop_ad=sum(error_2_sop_ad[,1])/(J*(J+1))
    
    #############
    
    B_ERROR_2_ellos[iter,iter_out]=ERROR_2_ellos
    
    B_ERROR_2_te_ellos[iter,iter_out]=ERROR_2_te_ellos
    
    B_ERROR_2_sop[iter,iter_out]=ERROR_2_sop
    
    # B_ERROR_2_sop_ad[iter,iter_out]=ERROR_2_sop_ad
    
    ############# HERE WE CALCULATE THE ERROR OF THE ESTIMATED RESPONSE VARIABLE
    
    BB_test=Data2B_simpson(X_test, M_test, nbasis=c(c1,c2,c3),sub = 25) # WE GENERATE THE BASIS OF THE TEST DATA SET
    
    # BB_ad_test=Data2B_simpson_ad(X_test, M_test, nbasis=c(c1,c2,c3),sub = 25,ndb = ndb)
    
    y_h_sop[iter,,iter_out] = BB_test$B %*% res$theta # ESTIMATED REPSONSE VARIABLE USING OUR APPROACH
    
    # y_h_adaptive[iter,,iter_out] = BB_ad_test$B %*% res_ad$theta
    
    
    # ESTIMATED REPSONSE VARIABLE USING HE GELLAR AND GOLDSMITH APPROACHES 
    
     for (j in 1:k) {

      ind=which(unique(M_it[iter,,iter_out])== M_it_test[iter,j,iter_out])
      ind_t=max(M_it[iter,,iter_out]) # OJO AQUÍ ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
    
      Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
      Beta_refund=Beta_refund[1:M_it_test[iter,j,iter_out]]
    
      Beta_refund_te=Beta_G_te[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
      Beta_refund_te=Beta_refund_te[1:M_it_test[iter,j,iter_out]]
  
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
    
    # ESTIMATION ERRORS FOR THE TEST DATA SET
    Y_ERROR_2_ellos[iter,iter_out]=sqrt(sum((y_test-y_h_ellos[iter,,iter_out])^2)/k)
    Y_ERROR_2_te_ellos[iter,iter_out]=sqrt(sum((y_test-y_h_ellos_te[iter,,iter_out])^2)/k)
    Y_ERROR_2_sop[iter,iter_out]=sqrt(sum((y_test-y_h_sop[iter,,iter_out])^2)/k)
    # Y_ERROR_2_sop_ad[iter,iter_out]=sqrt(sum((y_test-y_h_adaptive[iter,,iter_out])^2)/k)
    
    # GENERAL ESTIMATION ERRORS (NOT USED IN THE k-FOLD APPROACH)
    # Y_ERROR_2_ellos[iter,iter_out]=sqrt(sum((y-fit_Gellar$fitted.values)^2)/N)
    # Y_ERROR_2_te_ellos[iter,iter_out]=sqrt(sum((y-fit_Gellar_te$fitted.values)^2)/N)
    # Y_ERROR_2_sop[iter,iter_out]=sqrt(sum((y-res$fit$fitted.values)^2)/N)
    # Y_ERROR_2_sop_ad[iter,iter_out]=sqrt(sum((y-res_ad$fit$fitted.values)^2)/N)

  } # HERE ENDS THE MIDDLE FOR ITERATION (RUNS ON R) 
  
} # HERE ENDS THE OUTTER FOR ITERATION (RUNS ON iter_out)

end=proc.time()
time=end[1]-start[1]

time/60/60 # SIMULATION TIME IN HOURS


################# BOXPLOTS 

B_ERRORES_test=data.frame(B_ERROR_2_ellos_test,B_ERROR_2_te_ellos_test,B_ERROR_2_sop_test,B_ERROR_2_sop_ad_test) #B_ERROR_2,B_ERROR_2_te,

B_ERRORES=data.frame(B_ERROR_2_ellos,B_ERROR_2_te_ellos,B_ERROR_2_sop,B_ERROR_2_sop_ad) #B_ERROR_2,B_ERROR_2_te,

Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop,Y_ERROR_2_sop_ad) #,Y_ERROR_2_sop_manual) #Y_ERROR_2,Y_ERROR_2_te,

# SELECT THE SCENARIO TO PLOT 
case=5

Gellar=case
Gellar_te=1*8+case
New_approach=2*8+case
Adaptive=3*8+case

boxplot(B_ERRORES[,c(Gellar,New_approach)],pars  =  list(xaxt = "n"),xlab="")
axis(1, at=c(1,2),gap.axis = 0.75, labels = c("Gellar","New Approach"))
#

boxplot(Y_ERRORES[,c(Gellar,Gellar_te,New_approach)],pars  =  list(xaxt = "n"),xlab="")
axis(1, at=c(1,2,3),gap.axis = 0.75, labels = c("Gellar","Goldsmith","New Approach"))

# boxplot(B_ERRORES_test[,c(Gellar,New_approach)],pars  =  list(xaxt = "n"),xlab="")
# axis(1, at=c(1,2),gap.axis = 0.75, labels = c("Gellar","SOP"))


############# IN CASE WE NEED SOME METRICS OF THE ERRORS

case=8

Gellar=case
Gellar_te=1*8+case
New_approach=2*8+case
Adaptive=3*8+case

# mean(Y_ERRORES[,c(Gellar)])
# mean(Y_ERRORES[,c(Gellar_te)])
# mean(Y_ERRORES[,c(New_approach)])
# mean(Y_ERRORES[,c(Adaptive)])

# median(Y_ERRORES[,c(Gellar)])
# median(Y_ERRORES[,c(Gellar_te)])
# median(Y_ERRORES[,c(New_approach)])
# median(Y_ERRORES[,c(Adaptive)])

# sd(Y_ERRORES[,c(Gellar)])
# sd(Y_ERRORES[,c(Gellar_te)])
# sd(Y_ERRORES[,c(New_approach)])
# sd(Y_ERRORES[,c(Adaptive)])

# mean(B_ERRORES[,c(Gellar)])
# mean(B_ERRORES[,c(New_approach)])
# mean(B_ERRORES[,c(Adaptive)])

# median(B_ERRORES[,c(Gellar)])
# median(B_ERRORES[,c(New_approach)])
# median(B_ERRORES[,c(Adaptive)])

# sd(B_ERRORES[,c(Gellar)])
# sd(B_ERRORES[,c(New_approach)])
# sd(B_ERRORES[,c(Adaptive)])
#
