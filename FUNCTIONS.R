# NECESSARY FUNCTIONS IN OUR NEW METHODOLOGY

############## FUNCTION TO CREATE THE CORRECT B-SPLINE BASIS

bspline <-function(X., XL., XR., NDX., BDEG.){
  dx <- (XR. - XL.)/NDX.
  knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
  B <- spline.des(knots, X., BDEG.+1, 0*X.)$design
  res <- list(B = B, knots = knots)
  res
}

############## THE FUNCTION OF THE PACKAGE SOP AND SOPExamples WITH SOME NECESSARY ADJUSTMENS

sop.fit=function (y, X, Z, weights = NULL, G = NULL, vcstart = NULL, 
                  etastart = NULL, mustart = NULL, offset = NULL, family = gaussian(), 
                  control = sop.control()) 
{
  deviance <- function(C, G, w, sigma2, ssr, edf) {
    log_det_C <- determinant(C)$modulus
    log_det_G <- determinant(G)$modulus
    deviance <- log_det_C + log_det_G + sum(log(sigma2 * 
                                                  1/w)) + ssr/sigma2 + edf
    deviance
  }
  control <- do.call("sop.control", control)
  if (missing(X)) 
    stop("Missing argument: 'X' must be provided")
  if (missing(y)) 
    stop("Missing argument: 'y' must be provided")
  if (missing(Z)) 
    stop("Missing argument: 'Z' must be provided")
  if (!is.null(vcstart)) {
    if (length(vcstart) != (length(G) + 1)) {
      stop("The length of 'vcstart' should be equal to the length of 'G' + 1)")
    }
  }
  trace <- control$trace
  ynames <- if (is.matrix(y)) {
    rownames(y)
  }
  else {
    names(y)
  }
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- NCOL(X)
  EMPTY <- nvars == 0
  nc <- ncol(X)
  ncomp <- length(G)
  nq <- ncol(Z)
  if (!is.null(vcstart)) {
    la <- vcstart
  }
  else {
    la <- rep(1, len = ncomp + 1)
  }
  devold <- 1e+10
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
  }
  prior.weights <- weights
  if (is.null(offset)) {
    offset <- rep.int(0, nobs)
  }
  if (is.null(G)) {
    stop("Missing argument: 'G' must be provided")
  }
  Xnames <- colnames(X)
  na.ind <- is.na(y)
  y.tmp <- y
  y[na.ind] <- 1
  weights <- weights * (!na.ind)
  X <- as.matrix(X)
  start <- NULL
  variance <- family$variance
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) {
    stop("illegal `family' argument")
  }
  valideta <- family$valideta
  if (is.null(valideta)) {
    valideta <- function(eta) TRUE
  }
  validmu <- family$validmu
  if (is.null(validmu)) {
    validmu <- function(mu) TRUE
  }
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (NCOL(y) > 1) {
    stop("y must be univariate")
  }
  eta <- if (!is.null(etastart)) {
    etastart
  }
  else if (!is.null(start)) {
    if (length(start) != nvars) {
      stop(gettextf("Length of start should equal %d and correspond to initial coefs.", 
                    nvars))
    }
    else {
      coefold <- start
      offset + as.vector(if (nvars == 1) 
        c(X, Z) * start
        else c(X, Z) %*% start)
    }
  }
  else {
    family$linkfun(mustart)
  }
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta))) {
    stop("Can't find valid starting values: please specify some")
  }
  for (i in 1:control$maxit) {
    deriv <- family$mu.eta(eta)
    z <- (eta - offset) + (y - mu)/deriv
    w <- as.vector(deriv^2/family$variance(mu))
    w <- w * weights
    z[!weights] <- 0
    mat <- construct.matrices(X, Z, z, w,GLAM = FALSE)
    if (trace) 
      start1 <- proc.time()[3]
    for (it in 1:control$maxit) {
      Ginv <- 0
      for (ii in 1:ncomp) {
        Ginv <- Ginv + 1/la[ii + 1] * G[[ii]]
      }
      GG <- 1/Ginv
      V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., 
                           mat$ZtZ.)
      D <- diag(c(rep(0, nc), Ginv))
      H <- (1/la[1]) * V + D
      Hinv <- try(solve(H))
      if (inherits(Hinv, "try-error")) {
        Hinv <- ginv(H)
      }
      b <- (1/la[1]) * Hinv %*% mat$u
      b.fixed <- b[1:nc]
      b.random <- b[-(1:nc)]
      aux <- GG - diag(Hinv[-(1:nc), -(1:nc)])
      ed.sum <- 0
      ied <- taus <- NULL
      for (i.fun in 1:ncomp) {
        d1.d <- (1/la[i.fun + 1]) * G[[i.fun]]
        ed1 <- sum(aux * d1.d)
        ed1 <- ifelse(ed1 <= 1e-10, 1e-10, ed1)
        tau1 <- sum(b.random^2 * G[[i.fun]])/ed1
        tau1 <- ifelse(tau1 <= 1e-10, 1e-10, tau1)
        taus <- c(taus, tau1)
        ied <- c(ied, ed1)
      }
      ssr <- mat$yty. - t(c(b.fixed, b.random)) %*% (2 * 
                                                       mat$u - V %*% b)
      dev <- deviance(H, diag(GG), w[w != 0], la[1], ssr, 
                      sum(b.random^2 * Ginv))[1]
      if (family$family == "gaussian" | family$family == 
          "Gamma" | family$family == "quasipoisson") {
        sig2 <- as.numeric((ssr/(length(y[weights != 
                                            0]) - sum(ied) - nc)))
      }
      else {
        sig2 <- 1
      }
      lanew <- c(sig2, taus)
      dla <- abs(devold - dev)
      if (trace) {
        cat(sprintf("%1$3d %2$10.6f", it, dev))
        cat(sprintf("%8.3f", ied), "\n")
      }
      if (dla < control$epsilon) 
        break
      la <- lanew
      if (la[1]<1e-9) {
        la[1]=1e-9
      }
      devold <- dev
    }
    if (trace) {
      end1 <- proc.time()[3]
      cat("Timings:\nSOP", (end1 - start1), "seconds\n")
    }
    eta.old <- eta
    eta <- X %*% b.fixed + Z %*% b.random + offset
    mu <- linkinv(eta)
    tol <- sum((eta - eta.old)^2)/sum(eta^2)
    if (tol < control$epsilon | (family$family == "gaussian" & 
                                 family$link == "identity")) 
      break
  }
  end <- proc.time()[3]
  mu.eta <- family$mu.eta
  mu.eta.val <- mu.eta(eta)
  linear.predictor <- eta
  mu <- linkinv(eta)
  names(mu) <- ynames
  names(linear.predictor) <- ynames
  residuals <- family$dev.resids(y, mu, weights)
  s <- attr(residuals, "sign")
  if (is.null(s)) {
    s <- sign(y - mu)
  }
  residuals <- sqrt(pmax(residuals, 0)) * s
  names(residuals) <- ynames
  residuals[na.ind] <- NA
  names(b.fixed) <- Xnames
  names(b.random) <- colnames(Z)
  names(la) <- c("ssr", names(G))
  names(ied) <- c(names(G))
  dev.residuals <- family$dev.resids(y, mu, weights)
  dev.residuals[na.ind] <- NA
  deviance <- sum(dev.residuals, na.rm = TRUE)
  null.deviance <- glm(y ~ offset(offset), family = family, 
                       weights = prior.weights)$deviance
  out <- list(tol.ol = tol, it.ol = i, tol.il = dla, it.in = it, 
              vc = la, edf = ied)
  fit <- list()
  fit$b.fixed <- b.fixed
  fit$b.random <- b.random
  fit$fitted.values <- mu
  fit$linear.predictor <- linear.predictor
  fit$residuals <- residuals
  fit$X <- X
  fit$Z <- Z
  fit$G <- G
  fit$y <- y.tmp
  fit$weights <- weights
  fit$family <- family
  fit$out <- out
  fit$deviance <- deviance
  fit$null.deviance <- null.deviance
  fit$Vp <- Hinv
  class(fit) <- "sop"
  invisible(fit)
}


#################

fit.SOP = function (y, X, Z, Lambda, diagonal = FALSE, family = gaussian(), 
                    offset = NULL, weights = NULL, maxit = 200, thr = 0.001, 
                    trace = TRUE) 
{
  if (trace) 
    start <- proc.time()[3]
  if (is.null(offset)) 
    offset <- rep(0, length(y))
  if (is.null(weights)) 
    weights <- rep(1, length(y))
  if (is.null(X)) {
    stop("The design matrix associated with the fixed part of the model can not be NULL")
  }
  if (diagonal) {
    g <- construct.capital.lambda(Lambda)
  }
  else {
    g <- construct.capital.lambda.matrices(Lambda)
  }
  la = c(1, rep(1, length = length(g)))
  devold = 1e+10
  np <- c(ncol(X), ncol(Z))
  mustart <- etastart <- NULL
  nobs <- length(y)
  eval(family$initialize)
  mu <- mustart
  eta <- family$linkfun(mustart)
  if (trace) {
    cat("Effective dimensions\n")
    cat("-------------------------\n")
    cat(sprintf("%1$3s %2$12s", "It.", "Deviance"), sep = "")
    cat(sprintf("%12s", names(g)), sep = "")
    cat("\n")
  }
  for (it in 1:maxit) {
    deriv <- family$mu.eta(eta)
    z <- (eta - offset) + (y - mu)/deriv
    w <- as.vector(deriv^2/family$variance(mu))
    w <- w * weights
    mat <- construct.matrices(X, Z, z, w, FALSE)
    for (it in 1:maxit) {
      if (diagonal) {
        Ginv <- vector(length = length(g[[1]]))
        for (i in 1:length(g)) {
          Ginv <- Ginv + (1/la[i + 1]) * g[[i]]
        }
        G <- 1/Ginv
        V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., 
                             mat$ZtZ.)
        D <- Matrix::bdiag(diag(rep(0, np[1]), ncol = np[1]), 
                           diag(Ginv))
      }
      else {
        Ginv = 0
        for (i in 1:length(g)) {
          Ginv <- Ginv + (1/la[i + 1]) * g[[i]]
        }
        Ginv <- Matrix(Ginv)
        G <- try(solve(Ginv))
        if (class(G) == "try-error") {
          G <- ginv(Ginv)
        }
        V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., 
                             mat$ZtZ.)
        D <- Matrix::bdiag(diag(rep(0, np[1]), nrow = np[1], 
                                ncol = np[1]), Ginv)
      }
      H <- (1/la[1]) * V + D
      Hinv <- try(solve(H))
      if (class(Hinv) == "try-error") {
        Hinv <- ginv(as.matrix(H))
      }
      b <- (1/la[1]) * Hinv %*% mat$u
      b.fixed <- b[1:np[1]]
      b.random <- b[-(1:np[1])]
      if (diagonal) {
        aux <- G - diag(Hinv[-(1:np[1]), -(1:np[1])])
        ssv <- ed <- tau <- vector(mode = "list", length = length(g))
        for (i in 1:length(g)) {
          g.inv.d <- (1/la[i + 1]) * g[[i]]
          ed[[i]] <- sum(aux * g.inv.d)
          ed[[i]] <- ifelse(ed[[i]] <= 1e-10, 1e-10, 
                            ed[[i]])
          ssv[[i]] <- sum(b.random^2 * g[[i]])
          tau[[i]] <- ssv[[i]]/ed[[i]]
          tau[[i]] <- ifelse(tau[[i]] <= 1e-10, 1e-10, 
                             tau[[i]])
        }
        ssr = mat$yty. - t(c(b.fixed, b.random)) %*% 
          (2 * mat$u - V %*% b)
        dev <- deviance(H, diag(G), w[w != 0], la[1], 
                        ssr, sum(b.random^2 * Ginv))[1]
      }
      else {
        aux <- G - Hinv[-(1:np[1]), -(1:np[1])]
        updates <- lapply(1:length(g), function(i, g, 
                                                la, b.random, aux) {
          g.inv.d <- (1/la[i + 1]) * g[[i]]
          ed <- sum(colSums(t(aux) * g.inv.d))
          tau <- as.vector((t(b.random) %*% g[[i]] %*% 
                              b.random)/ed)
          tau <- ifelse(tau <= 1e-10, 1e-10, tau)
          res <- list(ed = ed, tau = tau)
          res
        }, g = g, la = la, b.random = b.random, aux = aux)
        ed <- as.list(unlist(updates)[seq(1, 2 * length(g), 
                                          by = 2)])
        tau <- as.list(unlist(updates)[seq(2, 2 * length(g), 
                                           by = 2)])
        ssr = mat$yty. - t(c(b.fixed, b.random)) %*% 
          (2 * mat$u - V %*% b)
        dev <- deviance(H, G, w[w != 0], la[1], ssr, 
                        (b.random) %*% Ginv %*% b.random)[1]
      }
      if (family$family == "gaussian" | family$family == 
          "quasipoisson") {
        phi <- as.numeric((ssr/(length(y[w != 0]) - sum(unlist(ed)) - 
                                  np[1])))
      }
      else {
        phi <- 1
      }
      lanew = c(phi, unlist(tau))
      dla = abs(devold - dev)
      if (trace) {
        cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
        cat(sprintf("%12.3f", unlist(ed)), sep = "")
        cat("\n")
      }
      if (dla < thr) 
        break
      la <- lanew
      if (la[1]<1e-9) {
        la[1]=1e-9
      }
      devold <- dev
    }
    eta.old <- eta
    eta <- X %*% b.fixed + Z %*% b.random + offset
    mu <- family$linkinv(eta)
    tol <- sum((eta - eta.old)^2)/sum(eta^2)
    if (tol < 1e-06 | family$family == "gaussian") 
      break
  }
  if (trace) {
    end <- proc.time()[3]
    cat("All process", end - start, "seconds\n")
  }
  res <- list()
  res$dat <- list(y = y, X = X, Z = Z)
  res$family = family
  res$ed <- unlist(ed)
  res$tot_ed <- sum(unlist(ed))
  res$vc <- unlist(tau)
  res$phi <- phi
  res$coeff <- c(b.fixed, b.random)
  res$eta <- eta
  res$mu <- mu
  res$deviance <- dev
  res$convergence <- ifelse(it < maxit, TRUE, FALSE)
  res$Vp <- Hinv
  res
}

################################## INNER_PRODUCT NECESSARY FUNCTIONS:

fdchk <- function(fdobj) {
  
  #  check the class of FDOBJ and extract coefficient matrix
  
  if (inherits(fdobj, "fd")) {
    coef  <- fdobj$coefs
  } else {
    if (inherits(fdobj, "basisfd")) {
      coef  <- diag(rep(1,fdobj$nbasis - length(fdobj$dropind)))
      fdobj <- fd(coef, fdobj)
    } else { 
      stop("FDOBJ is not an FD object.")
    }
  }
  
  #  extract the number of replications and basis object
  
  coefd <- dim(as.matrix(coef))
  if (length(coefd) > 2) stop("Functional data object must be univariate")
  nrep     <- coefd[2]
  basisobj <- fdobj$basis
  
  return(list(nrep, fdobj))
  
}

knotmultchk <- function(basisobj, knotmult) {
  type <- basisobj$type
  if (type == "bspline") {
    # Look for knot multiplicities in first basis
    params  <- basisobj$params
    nparams <- length(params)
    norder  <- basisobj$nbasis - nparams
    if (norder == 1) {
      knotmult <- c(knotmult, params)
    } else {
      if (nparams > 1) {
        for (i in 2:nparams) 
          if (params[i] == params[i-1]) knotmult <- c(knotmult, params[i])
      }
    }
  }
  return(knotmult)
}

# SIMPSON NUMERIC INTEGRATION VERSION 1.0 SEE BURDEN AND FAIRES (2016)

Simpson <- function(fdobj1, fdobj2=NULL, fdobj3=NULL, Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),rng = range1, sub=25, wtfd = 0){
  
  #  Check FDOBJ1 and get no. replications and basis object
  
  result1   <- fdchk(fdobj1)
  nrep1     <- result1[[1]]
  fdobj1    <- result1[[2]]
  coef1     <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1     <- basisobj1$type
  range1    <- basisobj1$rangeval
  
  #  Default FDOBJ2 to a constant function, using a basis that matches
  #  that of FDOBJ1 if possible.
  
  if (is.null(fdobj2)) {
    tempfd    <- fdobj1
    tempbasis <- tempfd$basis
    temptype  <- tempbasis$type
    temprng   <- tempbasis$rangeval
    if (temptype == "bspline") {
      basis2 <- create.bspline.basis(temprng, 1, 1)
    } else {
      if (temptype == "fourier") basis2 <- create.fourier.basis(temprng, 1)
      else                       basis2 <- create.constant.basis(temprng)
    }
    fdobj2 <- fd(1,basis2)
  }
  
  #  Check FDOBJ2 and get no. replications and basis object
  
  result2   <- fdchk(fdobj2)
  nrep2     <- result2[[1]]
  fdobj2    <- result2[[2]]
  coef2     <- fdobj2$coefs
  basisobj2 <- fdobj2$basis
  type2     <- basisobj2$type
  range2    <- basisobj2$rangeval
  
  
  if (is.null(fdobj3)) {
    return(inprod(fdobj1, fdobj2, Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                  rng = range1, wtfd = 0))
  }
  
  # check ranges
  
  if (rng[1] < range1[1] || rng[2] > range1[2]) stop(
    "Limits of integration are inadmissible.")

  #  check LFDOBJ1 and LFDOBJ2
  
  Lfdobj1 <- int2Lfd(Lfdobj1)
  Lfdobj2 <- int2Lfd(Lfdobj2)
  # Lfdobj3 <- int2Lfd(Lfdobj3)
  
  #  set iter
  
  iter <- 0
  
  # The default case, no multiplicities.
  
  rngvec <- rng
  
  #  check for any knot multiplicities in either argument
  
  knotmult <- numeric(0)

  if (type1 == "bspline") knotmult <- knotmultchk(basisobj1, knotmult)
  if (type2 == "bspline") knotmult <- knotmultchk(basisobj2, knotmult)

  #  Modify RNGVEC defining subinvervals if there are any
  #  knot multiplicities.
  
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < rng[2]]
    rngvec   <- c(rng[1], knotmult, rng[2])
  }
  
  #  check for either coefficient array being zero
  
  if ((all(c(coef1) == 0) || all(c(coef2) == 0)))
    return(matrix(0,nrep1,nrep2*length(fdobj3)))
  
   n=2*sub # THIS IS THE INTERVALS FOR THE SIMPSON METHOD

  nrng <- length(rngvec)
  for (irng  in  2:nrng) {  
    rngi <- c(rngvec[irng-1],rngvec[irng])
    
    #  change range so as to avoid being exactly on multiple knot values
    
    if (irng > 2   ) rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) rngi[2] <- rngi[2] - 1e-10

    width = (rngi[2] - rngi[1])/n
    
    # the first iteration uses just the endpoints
    fx1 <- eval.fd(rngi, fdobj1, Lfdobj1)
    fx2 <- eval.fd(rngi, fdobj2, Lfdobj2)
    fx3 <- matrix(fdobj3,nrow = dim(fx2)[1],ncol = length(fdobj3),byrow = TRUE)
    fx2 <- Rten2(fx2,fx3)
    
    XI0 = matrix(crossprod(fx1,fx2),nrep1,nrep2*length(fdobj3))
    
    XI1=0
    XI2=0
    
    for (i in 1:(n-1)) {

### THESE LINES OF CODE ARE THE CORE OF THE INNER PRODUCT. NOTICE HOW WE PERFORM THE INTEGRATION ONLY IN THE t VARIABLE BUT THEN WE CREATE THE BIDEMENSIONAL BASIS IN EVERY ITERATION
      
      x = rngi[1] + i*width
      fx1 <- eval.fd(x, fdobj1, Lfdobj1)
      fx2 <- eval.fd(x, fdobj2, Lfdobj2)
      fx3 <- matrix(fdobj3,nrow = dim(fx2)[1],ncol = length(fdobj3),byrow = TRUE)
      fx2 <- Rten2(fx2,fx3)
      
      Fx=matrix(crossprod(fx1,fx2),nrep1,nrep2*length(fdobj3))
      
      if (i%%2==0) {
        
        XI2= XI2 + Fx
      }else{
        XI1= XI1 + Fx
      }
      
    }
    
    XI= width*(XI0 + 2*XI2 + 4*XI1)/3
  }
  
  if(length(dim(XI) == 2)) {
    #  coerce inprodmat to be nonsparse
    return(as.matrix(XI))
  } else {
    #  allow inprodmat to be sparse if it already is
    return(XI)
  }
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL:

Data2B_simpson=function(X, M, nbasis=c(30,30,30), bdeg=c(3,3,3),sub=25){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  c3=nbasis[3]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M)!=N) {
    stop("length of 'M' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M[i]) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }

    ############### HERE WE CREATE THE BASIS FOR THE DATA
    
    rng[i,]= c(1, M[i])
    
    XL=rng[i,1]-0.001
    XR=rng[i,2]+0.001
    
    c=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(1:M[i], XL, XR, c, bdeg[1]))
    
    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT 
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d(aux,2,c1)
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,1:M[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,1:M[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)
    
    c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_Phi <- paste("Phi", i, sep = "_")
    L_Phi[[i]]=assign(nam_Phi, bspline(1:M[i], XL, XR, c_t, bdeg[2]))
    
  }
  
  ####### HERE WE CREATE THE MARGINAL BASIS FOR THE T VARIABLE In B(t,T)
  
  xlim_T <- c(min(M), max(M))
  XL_T=xlim_T[1]-0.001
  XR_T=xlim_T[2]+0.001
  
  c_T=c3-bdeg[3] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  B_T=bspline(M, XL_T, XR_T, c_T, bdeg[3])
  

  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT  
                  # IS NOT MECESSARY TO DO THIS FOR THE T MARGINAL BASIS
  for (i in 1:N) {
    
    # DATA BASIS
    
    breaks=L_X[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_X_aux <- paste("B_aux", i, sep = "_")
    L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))
    
    # t MARGINAL BASIS
    
    breaks=L_Phi[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_Phi_aux <- paste("Phi_aux", i, sep = "_")
    L_Phi_aux[[i]]=assign(nam_Phi_aux, create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2))))
    
  }
  
  # PERFORMING THE INNER PRODUCT
  
  L_aux=list()
  
  for (i in 1:N) {
    
    PROD=Simpson(L_X_aux[[i]],L_Phi_aux[[i]],B_T$B[i,],rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
    }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, B_T=B_T, B_Phi=L_Phi)
  
}

################### SAME AS THE PREVIOUS FUNCTION BUT USING THE ADAPTIVE APPORACH, SEE RODRÍGUEZ ET AL (2019)

Data2B_simpson_ad=function(X, M, nbasis=c(30,30,30), bdeg=c(3,3,3), sub=25, ndb=25){
  
  require(fda)
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  c3=nbasis[3]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M)!=N) {
    stop("length of 'M' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M[i]) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }

    
    ###############
    
    rng[i,]= c(1, M[i])
    
    XL=rng[i,1]-0.001
    XR=rng[i,2]+0.001
    
    c=c1-bdeg[1]
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(1:M[i], XL, XR, c, bdeg[1]))
    
    ######### 
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d_ad(aux,2,c1,ndb = ndb)
    
    aux_3=XZG2theta_1d_ad(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,1:M[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,1:M[i]])-L_y[[i]]))
    
    ###############
    
    c_t=c2-bdeg[2]
    
    nam_Phi <- paste("Phi", i, sep = "_")
    L_Phi[[i]]=assign(nam_Phi, bspline(1:M[i], XL, XR, c_t, bdeg[2]))
    
  }
  
  #######
  
  xlim_T <- c(min(M), max(M))
  XL_T=xlim_T[1]-0.001
  XR_T=xlim_T[2]+0.001
  
  c_T=c3-bdeg[3]
  
  B_T=bspline(M, XL_T, XR_T, c_T, bdeg[3])
  
  #################
  
  
  for (i in 1:N) {
    
    #
    
    breaks=L_X[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_X_aux <- paste("B_aux", i, sep = "_")
    L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))
    
    #
    
    breaks=L_Phi[[i]]$knots
    
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    
    nam_Phi_aux <- paste("Phi_aux", i, sep = "_")
    L_Phi_aux[[i]]=assign(nam_Phi_aux, create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2))))
    
  }
    L_aux=list()
  
  for (i in 1:N) {
    
    PROD=Simpson(L_X_aux[[i]],L_Phi_aux[[i]],B_T$B[i,],rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, B_T=B_T, B_Phi=L_Phi)
  
}


################### THE NEXT FUNCTIONS GENERATE THE NECESSARY MATRICES FOR THE MIXED MODEL:

# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 2D CASE, SEE LEE AND DURBÁN REGUERA (2010)

B2XZG <-function (B, pord=c(2,2),c=c(10,10)) {
    
    c1=c[1]
    c2=c[2]
    
    c1c2 = ncol(B)
    if (c1c2!=c1*c2) {
      stop("c1 * c2 must me equal to the number of colums of B")
    }
    
    D_1 = diff(diag(c1), differences=pord[1])
    D_2 = diff(diag(c2), differences=pord[2])
    
    P1.svd = svd(crossprod(D_1))
    P2.svd = svd(crossprod(D_2))
    
    U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
    U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
    d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues
    
    U_2s = (P2.svd$u)[,1:(c2-pord[2])] # eigenvectors
    U_2n=((P2.svd$u)[,-(1:(c2-pord[2]))])
    d2 = (P2.svd$d)[1:(c2-pord[2])]  # eigenvalues
    
    
    T_n=kronecker(U_1n,U_2n)
    
    AUX_1=kronecker(U_1n,U_2s)
    AUX_2=kronecker(U_1s,U_2n)
    AUX_3=kronecker(U_1s,U_2s)
    
    T_s=cbind(AUX_1,AUX_2,AUX_3)
    
    
    Z = B%*%T_s
    X = B%*%T_n
    
    ####
    
    d_1s=diag(P1.svd$d)[1:(c1-pord[1]),1:(c1-pord[1])]
    d_2s=diag(P2.svd$d)[1:(c2-pord[2]),1:(c2-pord[2])]
    
    T_1=kronecker(diag(pord[1]),d_2s)
    T_2=matrix(0,nrow = pord[2]*(c1-pord[1]),ncol = pord[2]*(c1-pord[1]))
    T_3=kronecker(diag(c1-pord[1]),d_2s)
    
    T_21=cbind(T_1,matrix(0,nrow=dim(T_1)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(T_1)[2]))
    T_22=cbind(matrix(0,nrow=dim(T_2)[1],ncol=dim(T_1)[2]),T_2,matrix(0,nrow=dim(T_2)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(T_1)[2]-dim(T_2)[2]))
    T_23=cbind(matrix(0,nrow=((c2-pord[2])*(c1-pord[1])),ncol=(c1*c2-pord[1]*pord[2])-dim(T_3)[2]),T_3)
    
    H_1=matrix(0,nrow = pord[1]*(c2-pord[2]),ncol = pord[1]*(c2-pord[2]))
    H_2=kronecker(d_1s,diag(pord[2]))
    H_3=kronecker(d_1s,diag(c2-pord[2]))
    
    H_11=cbind(H_1,matrix(0,nrow=dim(H_1)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(H_1)[2]))
    H_12=cbind(matrix(0,nrow=dim(H_2)[1],ncol=dim(H_1)[2]),H_2,matrix(0,nrow=dim(H_2)[1],ncol=(c1*c2-pord[1]*pord[2])-dim(H_1)[2]-dim(H_2)[2]))
    H_13=cbind(matrix(0,nrow=((c2-pord[2])*(c1-pord[1])),ncol=(c1*c2-pord[1]*pord[2])-dim(H_3)[2]),H_3)
    
    L_2=rbind(T_21,T_22,T_23)
    L_1=rbind(H_11,H_12,H_13)
    
    t_2=diag(L_1)
    t_1=diag(L_2)
    
    G=list(t_1,t_2)
    names(G)=c("t_1","t_2")
    
    T=cbind(T_n,T_s)
    
    ####
    
    list(X = X, Z = Z, G=G, T=T, d1 = d1, d2=d2, D_1 = D_1, D_2=D_2, U_1n=U_1n, U_1s=U_1s, U_2n=U_2n, U_2s=U_2s, T_n=T_n, T_s=T_s, t_1=t_1, t_2=t_2)
  }

# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 1D CASE, SEE LEE AND DURBÁN REGUERA (2010)
B2XZG_1d <-function (B, pord=c(2),c=c(10)) {
    
    c1=c[1]
    
    
    
    
    D_1 = diff(diag(c1), differences=pord[1])
    
    P1.svd = svd(crossprod(D_1))
    
    U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
    U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
    d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues
    
    T_n=(U_1n)
    
    T_s=U_1s
    
    
    Z = B%*%T_s
    X = B%*%T_n
    
    ####
    G=list(d1)
    # print(G)
    
    # G=list(t_1,t_2)
    # names(G)=c("t_1","t_2")
    # 
    T=cbind(T_n,T_s)
    
    ####
    
    list(X = X, Z = Z, G=G, T=T, d1 = d1, D_1 = D_1, U_1n=U_1n, U_1s=U_1s, T_n=T_n, T_s=T_s)
  }

# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 1D CASE WITH THE ADAPTIVE APPROACH, SEE LEE AND DURBÁN REGUERA (2010)
B2XZG_1d_ad <-function (B, pord=c(2),c=c(10), ndb=25) {
  
  c1=c[1]
  
  D_1 = diff(diag(c1), differences=pord[1])
  
  P1.svd = svd(crossprod(D_1))
  
  U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
  U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
  d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues
  
  T_n=(U_1n)
  
  T_s=U_1s
  
  Z = B%*%T_s
  X = B%*%T_n
  
  bdeg=3
  
  v.comp_1 <- 1:ncol(U_1s)/ncol(U_1s)
  C_1 <- t(SOPExamples:::bbase(v.comp_1, min(v.comp_1), max(v.comp_1), ndx = ndb-bdeg, bdeg = bdeg)$B)
  # One random (smooth) component
  Lambda_1 <- vector("list", length = 1) 
  # Precision matrices associated with that component
  Lambda_1[[1]] <- vector("list", length = nrow(C_1))
  for (i in 1:nrow(C_1)) {
    Lambda_1[[1]][[i]] <-  t(P1.svd$u[,1:(c1-pord[1])]) %*% t(D_1) %*% diag(C_1[i,]) %*% D_1 %*% P1.svd$u[,1:(c1-pord[1])]
    }
  
  
  ####
  G=Lambda_1
  T=cbind(T_n,T_s)
  ####
  
  list(X = X, Z = Z, G=G, T=T, d1 = d1, D_1 = D_1, U_1n=U_1n, U_1s=U_1s, T_n=T_n, T_s=T_s)
}


##################### FUNCTIONS FOR ESTIMATING THE COEFFICIENTS OF THE MODEL

# THIS FUNCTION FIT THE MODEL USING OUR APPROACH 
# AND RECOVER THE ORIGINAL THETA COEFFICIENTS THAT ARE NECESSARY FOR CALCULATE THE ESTIMATED FUNCTIONAL COEFFICIENT BETA FOR THE 2D CASE.
XZG2theta=function(X, Z, G, T, y, family=gaussian(), offset = NULL){
  require(SOP)
  
  if (dim(X)[1]!=dim(Z)[1]) {
    stop("'X' and 'Z'must have same numbers of rows")
  }
  
  if ( dim(Z)[2]!=length(G[[1]]) || dim(Z)[2]!=length(G[[2]]) ) {
    stop("The number of columns of 'Z' must be equal to the length of 'G'")
  }
  
  w = as.vector(rep(1,dim(X)[1]))
  
  fit=sop.fit(X = X, Z = Z, G = G, 
              y = y, family = family,
              control = list(trace = FALSE),offset = offset)
  
  if (dim(fit$Vp)[1]-dim(T)[1]==0) {
  
    theta_aux=c(fit$b.fixed,fit$b.random) 
    
  }else{
  
    aux=dim(fit$Vp)[1]-dim(T)[1]
    
    theta_aux=c(fit$b.fixed[-(1:aux)],fit$b.random) # ESTE [1:4] FUE AGREGADO PARA QUE AL AGREGAR OTRAS VARIABLES LINEALES SOLO COGIERA LA PARTE FUNCIONAL 
    
  }
  
  theta=T %*% theta_aux
  
  
  
  if (dim(fit$Vp)[1]==dim(T)[1]) {
    
    covar_theta=T %*% fit$Vp %*% t(T)
    std_error_theta=sqrt(diag(T %*% fit$Vp %*% t(T)))
    std_error_non_functional=NULL
    
  }else{
    
    covar_theta=T %*% fit$Vp[-(1:aux),-(1:aux)] %*% t(T)
    std_error_theta=sqrt(diag(T %*% fit$Vp[-(1:aux),-(1:aux)] %*% t(T)))
    std_error_non_functional=sqrt(diag(fit$Vp[1:aux,1:aux]))
    WALD= (fit$b.fixed[(1:aux)]/std_error_non_functional)  
    p_values=2*pnorm(abs(WALD),lower.tail = FALSE)
    }
  
  
  
  list (fit=fit, theta=theta, std_error_theta=std_error_theta, std_error_non_functional=std_error_non_functional, covar_theta=covar_theta, p_values=p_values)
}

# THIS FUNCTION FIT THE MODEL USING OUR APPROACH 
# AND RECOVER THE ORIGINAL THETA COEFFICIENTS THAT ARE NECESSARY FOR CALCULATE THE ESTIMATED FUNCTIONAL COEFFICIENT BETA FOR THE 1D CASE.
XZG2theta_1d=function(X, Z, G, T, y, family=gaussian()){
  require(SOP)
  
  fit=sop.fit(X = X, Z = Z, G = G, 
                y = y, family = family, weights = rep(1,dim(X)[1]), 
                control = list(trace = FALSE))
  
  theta_aux=c(fit$b.fixed,fit$b.random)
  
  theta=T %*% theta_aux
  
  list (fit=fit,theta=theta)
}

# THIS FUNCTION FIT THE MODEL USING OUR APPROACH 
# AND RECOVER THE ORIGINAL THETA COEFFICIENTS THAT ARE NECESSARY FOR CALCULATE THE ESTIMATED FUNCTIONAL COEFFICIENT BETA FOR THE 1D CASE 
# WITH THE ADAPTIVE APPROACH.
XZG2theta_1d_ad=function(X, Z, G, T, y, family=gaussian()){
  require(SOP)
  require(SOPExamples)
  
  fit=fit.SOP(X = X, Z = Z, Lambda = G, 
              y = y, family = family ,trace = 0, diagonal = 0)
  
  theta_aux=fit$coeff
  
  theta=T %*% theta_aux
  
  list (fit=fit,theta=theta)
}
