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

inprod=function (fdobj1, fdobj2 = NULL, Lfdobj1 = int2Lfd(0), Lfdobj2 = int2Lfd(0), rng = range1, wtfd = 0) {
  
  result1 <- fdchk(fdobj1)
  nrep1 <- result1[[1]]
  fdobj1 <- result1[[2]]
  coef1 <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1 <- basisobj1$type
  range1 <- basisobj1$rangeval
  if (is.null(fdobj2)) {
    tempfd <- fdobj1
    tempbasis <- tempfd$basis
    temptype <- tempbasis$type
    temprng <- tempbasis$rangeval
    if (temptype == "bspline") {
      basis2 <- create.bspline.basis(temprng, 1, 1)
    }
    else {
      if (temptype == "fourier") 
        basis2 <- create.fourier.basis(temprng, 1)
      else basis2 <- create.constant.basis(temprng)
    }
    fdobj2 <- fd(1, basis2)
  }
  result2 <- fdchk(fdobj2)
  nrep2 <- result2[[1]]
  fdobj2 <- result2[[2]]
  coef2 <- fdobj2$coefs
  basisobj2 <- fdobj2$basis
  type2 <- basisobj2$type
  range2 <- basisobj2$rangeval
  if (rng[1] < range1[1] || rng[2] > range1[2]) 
    stop("Limits of integration are inadmissible.")
  # if (is.fd(fdobj1) && is.fd(fdobj2) && type1 == "bspline" && 
  #     type2 == "bspline" && is.eqbasis(basisobj1, basisobj2) && 
  #     is.integer(Lfdobj1) && is.integer(Lfdobj2) && length(basisobj1$dropind) == 
  #     0 && length(basisobj1$dropind) == 0 && wtfd == 0 && all(rng == 
  #                                                             range1)) {
  #   inprodmat <- inprod.bspline(fdobj1, fdobj2, Lfdobj1$nderiv, 
  #                               Lfdobj2$nderiv)
  #   return(inprodmat)
  # }
  Lfdobj1 <- int2Lfd(Lfdobj1)
  Lfdobj2 <- int2Lfd(Lfdobj2)
  iter <- 0
  rngvec <- rng
  knotmult <- numeric(0)
  if (type1 == "bspline") 
    knotmult <- knotmultchk(basisobj1, knotmult)
  if (type2 == "bspline") 
    knotmult <- knotmultchk(basisobj2, knotmult)
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < 
                           rng[2]]
    rngvec <- c(rng[1], knotmult, rng[2])
  }
  if ((all(c(coef1) == 0) || all(c(coef2) == 0))) 
    return(matrix(0, nrep1, nrep2))
  JMAX <- 15
  JMIN <- 5
  EPS <- 1e-04
  inprodmat <- matrix(0, nrep1, nrep2)
  nrng <- length(rngvec)
  for (irng in 2:nrng) {
    rngi <- c(rngvec[irng - 1], rngvec[irng])
    if (irng > 2) 
      rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) 
      rngi[2] <- rngi[2] - 1e-10
    iter <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1, JMAXP)
    h[2] <- 0.25
    s <- array(0, c(JMAXP, nrep1, nrep2))
    sdim <- length(dim(s))
    fx1 <- eval.fd(rngi, fdobj1, Lfdobj1)
    fx2 <- eval.fd(rngi, fdobj2, Lfdobj2)
    if (!is.numeric(wtfd)) {
      wtd <- eval.fd(rngi, wtfd, 0)
      fx2 <- matrix(wtd, dim(wtd)[1], dim(fx2)[2]) * fx2
    }
    s[1, , ] <- width * matrix(crossprod(fx1, fx2), nrep1, 
                               nrep2)/2
    tnm <- 0.5
    for (iter in 2:JMAX) {
      tnm <- tnm * 2
      if (iter == 2) {
        x <- mean(rngi)
      }
      else {
        del <- width/tnm
        x <- seq(rngi[1] + del/2, rngi[2] - del/2, del)
      }
      fx1 <- eval.fd(x, fdobj1, Lfdobj1)
      fx2 <- eval.fd(x, fdobj2, Lfdobj2)
      if (!is.numeric(wtfd)) {
        wtd <- eval.fd(wtfd, x, 0)
        fx2 <- matrix(wtd, dim(wtd)[1], dim(fx2)[2]) * 
          fx2
      }
      chs <- width * matrix(crossprod(fx1, fx2), nrep1, 
                            nrep2)/tnm
      s[iter, , ] <- (s[iter - 1, , ] + chs)/2
      if (iter >= 5) {
        ind <- (iter - 4):iter
        ya <- s[ind, , ]
        ya <- array(ya, c(5, nrep1, nrep2))
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y <- ya[ns, , ]
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5 - m)) {
            ho <- xa[i]
            hp <- xa[i + m]
            w <- (cs[i + 1, , ] - ds[i, , ])/(ho - hp)
            ds[i, , ] <- hp * w
            cs[i, , ] <- ho * w
          }
          if (2 * ns < 5 - m) {
            dy <- cs[ns + 1, , ]
          }
          else {
            dy <- ds[ns, , ]
            ns <- ns - 1
          }
          y <- y + dy
        }
        ss <- y
        errval <- max(abs(dy))
        ssqval <- max(abs(ss))
        if (all(ssqval > 0)) {
          crit <- errval/ssqval
        }
        else {
          crit <- errval
        }
        if (crit < EPS && iter >= JMIN) 
          break
      }
      s[iter + 1, , ] <- s[iter, , ]
      h[iter + 1] <- 0.25 * h[iter]
      if (iter == JMAX) 
        warning("Failure to converge.")
    }
    inprodmat <- inprodmat + ss
  }
  if (length(dim(inprodmat) == 2)) {
    return(as.matrix(inprodmat))
  }
  else {
    return(inprodmat)
  }
}

# SIMPSON NUMERIC INTEGRATION VERSION 1.0 SEE BURDEN AND FAIRES (2016)

Simpson <- function(fdobj1, fdobj2=NULL, fdobj3=NULL, Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),rng = rng, sub=25, wtfd = 0){
  
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
                  rng = rng, wtfd = 0))
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

Data2B_simpson=function(X, M, nbasis=c(30,30,30), bdeg=c(3,3,3),sub=25, lim=NULL){
  
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
  
  if (is.null(lim)) {
  
    B_T=bspline(M, XL_T, XR_T, c_T, bdeg[3])
      
  }else{
    
    B_T=bspline(M, lim[1]-0.001, lim[2]+0.001, c_T, bdeg[3])
  }
  
  
  # matplot(B_T$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
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
  
  
  for (i in 1:N) {
    
    PROD=Simpson(L_X_aux[[i]],L_Phi_aux[[i]],B_T$B[i,],rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
    }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, B_T=B_T, B_Phi=L_Phi)
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT FOR ALL SUBJECTS:

Data2B_simpson_no_vc_old=function(X, M, nbasis=c(30,31), bdeg=c(3,3),sub=25, lim=NULL){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
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
    
  }
  
  rng_t= c(1, max(M))
  
  XL_t=rng_t[1]-0.001
  XR_t=rng_t[2]+0.001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(1:max(M), XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
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
    
  }
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n=length(breaks)
  
  
  Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
  
  for (i in 1:N) {
    

    PROD=Simpson(L_X_aux[[i]],Phi_aux,rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT FOR ALL SUBJECTS:

Data2B_simpson_no_vc=function(X, M, nbasis=c(30,31), bdeg=c(3,3),sub=25, lim=NULL){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
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
    
  }
  
  rng_t= c(1, max(M))
  
  XL_t=rng_t[1]-0.001
  XR_t=rng_t[2]+0.001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(1:max(M), XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT

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
  
  }
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n_knots=length(breaks)
  
  
  # Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n_knots-2):(n_knots+2)))
  
  for (i in 1:N) {
  

    # if (i!=N) {

      ind_knots=max(which(Phi$knots<=(M[i]+0.001)))

    # }
    #else{
    #   ind_knots=length(Phi$knots)
    # }

   # IDEA 1: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE)

    # breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_1=c(breaks[1:5],rep(0,n_knots-n),breaks[6:(n-5)],breaks[(n-4):n])
    # 
    # n_breaks_no_vc_1=length(breaks_no_vc_1)
    # 
    # Phi_aux_1=create.bspline.basis(breaks=breaks_no_vc_1,norder=bdeg[2]+1,dropin=c(1:5,(n_breaks_no_vc_1-2):(n_breaks_no_vc_1+2)))
      
   # IDEA 2: TOMAR LOS NODOS DEL INTERVALO

    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_2=c(breaks[1:2],rep(0,n_knots-n),breaks[3:(n-2)],breaks[(n-1):n])
    # 
    # n_breaks_no_vc_2=length(breaks_no_vc_2)
    # 
    # Phi_aux_2=create.bspline.basis(breaks=breaks_no_vc_2,norder=bdeg[2]+1,dropin=c(1:7,(n_breaks_no_vc_2):(n_breaks_no_vc_2+2)))
      
    # IDEA 3: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE) PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    tag=3
    Phi_aux_3=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
      
      
    # IDEA 4: TOMAR LOS NODOS DEL INTERVALO PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
      
      # breaks=Phi$knots[4:(ind_knots)]
      # dife=diff(breaks)[1]
      # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
      # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
      # n=length(breaks)
      # 
      # tag=4
      # 
      # Phi_aux_4=create.bspline.basis(breaks=breaks, norder=bdeg[2]+1,dropin=c(1:2,(n):(n+2)))
      
    
    PROD=Simpson(L_X_aux[[i]],Phi_aux_3,rng = c(1,M[i]),sub = sub)/M[i]
    
    if (tag==3 || tag==4) {
      PROD=cbind(PROD,matrix(0, nrow = dim(PROD)[1], ncol=nbasis[2]-dim(PROD)[2]))
    }
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
}

# matplot(eval.basis(1:max(M),Phi_aux_2),type="l")

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
    
  for (i in 1:N) {
    
    PROD=Simpson(L_X_aux[[i]],L_Phi_aux[[i]],B_T$B[i,],rng = c(1,M[i]),sub = sub)/M[i]
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, B_T=B_T, B_Phi=L_Phi)
  
}

####################### TRANSFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL BUT WITH THE SAME FUNCTIONAL COEFFICIENT AND PARTIALLY OBSERVED CURVES:

Data2B_simpson_no_vc_MissD=function(X, M, nbasis=c(30,31), bdeg=c(3,3),sub=25, lim=NULL){
  
  require(fda)
  
  ## SETTING SOME MATRICES AND PARAMETERS
  
  # X is the matrix of Data # dim(X) == N x max(M_b)
  # M is the matrix with columns M_a and M_b of numbers of observations dim(M) == N x 2 
  
  N=dim(X)[1]
  
  c1=nbasis[1]
  c2=nbasis[2]
  
  M_a=M[,1]
  M_b=M[,2]
  
  error=K=NULL
  
  rng=matrix(0,ncol = 2,nrow = N)
  
  L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  
  A=matrix(0,nrow = N, ncol = N*c1)
  
  if (length(M_a)!=N | length(M_b)!=N) {
    stop("length of 'M_a' and 'M_b' must be equal to 'N'")
  }
  
  for (i in 1:N) {
    
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M_b[i]-M_a[i]+1) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }
    
    ############### HERE WE CREATE THE BASIS FOR THE DATA
    
    rng[i,]= c(M_a[i], M_b[i])
    
    XL=rng[i,1]-0.001
    XR=rng[i,2]+0.001
    
    c=c1-bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
    
    nam_X <- paste("B", i, sep = "_")
    L_X[[i]]=assign(nam_X, bspline(M_a[i]:M_b[i], XL, XR, c, bdeg[1]))
    
    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT 
    
    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)
    
    aux=L_X[[i]]$B
    
    aux_2=B2XZG_1d(aux,2,c1)
    
    aux_3=XZG2theta_1d(X = aux_2$X,Z = aux_2$Z,G = aux_2$G,T = aux_2$T,y=X[i,M_a[i]:M_b[i]] )
    
    A[i,((c1*(i-1))+1):(i*c1)]=aux_3$theta
    
    L_theta[[i]]=aux_3$theta
    
    nam_y <- paste("y_h", i, sep = "_")
    
    L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])
    
    error[i]=mean(abs((X[i,M_a[i]:M_b[i]])-L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA
    
    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)
    
  }
  
  rng_t= c(1, max(M_b))
  
  XL_t=rng_t[1]-0.001
  XR_t=rng_t[2]+0.001
  
  c_t=c2-bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  
  Phi=bspline(1:max(M_b), XL_t, XR_t, c_t, bdeg[2])
  
  # matplot(Phi$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT
  
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
    
  }
  
  # PERFORMING THE INNER PRODUCT
  
  # t MARGINAL BASIS
  
  breaks=Phi$knots #[1:ind_knots]
  dife=diff(breaks)[1]
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
  n_knots=length(breaks)
  
  
  # Phi_aux=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n_knots-2):(n_knots+2)))
  
  for (i in 1:N) {
    
    
    # if (i!=N) {
    
    ind_knots=max(which(Phi$knots<=(M_b[i]+0.001)))
    
    # }
    #else{
    #   ind_knots=length(Phi$knots)
    # }
    
    # IDEA 1: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE)
    
    # breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_1=c(breaks[1:5],rep(0,n_knots-n),breaks[6:(n-5)],breaks[(n-4):n])
    # 
    # n_breaks_no_vc_1=length(breaks_no_vc_1)
    # 
    # Phi_aux_1=create.bspline.basis(breaks=breaks_no_vc_1,norder=bdeg[2]+1,dropin=c(1:5,(n_breaks_no_vc_1-2):(n_breaks_no_vc_1+2)))
    
    # IDEA 2: TOMAR LOS NODOS DEL INTERVALO
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # breaks_no_vc_2=c(breaks[1:2],rep(0,n_knots-n),breaks[3:(n-2)],breaks[(n-1):n])
    # 
    # n_breaks_no_vc_2=length(breaks_no_vc_2)
    # 
    # Phi_aux_2=create.bspline.basis(breaks=breaks_no_vc_2,norder=bdeg[2]+1,dropin=c(1:7,(n_breaks_no_vc_2):(n_breaks_no_vc_2+2)))
    
    # IDEA 3: TOMAR LOS NODOS DEL INTERVALO +-3(GRADO DEL BSPLINE) PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    breaks=Phi$knots[1:(ind_knots+bdeg[2])]
    dife=diff(breaks)[1]
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    n=length(breaks)
    tag=3
    Phi_aux_3=create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2)))
    
    
    # IDEA 4: TOMAR LOS NODOS DEL INTERVALO PERO RELLENAR CON CEROS LA BASE Y NO LOS NODOS
    
    # breaks=Phi$knots[4:(ind_knots)]
    # dife=diff(breaks)[1]
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # breaks=c(breaks[1]-dife,breaks, breaks[length(breaks)]+dife)
    # n=length(breaks)
    # 
    # tag=4
    # 
    # Phi_aux_4=create.bspline.basis(breaks=breaks, norder=bdeg[2]+1,dropin=c(1:2,(n):(n+2)))
    
    
    PROD=Simpson(L_X_aux[[i]],Phi_aux_3,rng = c(M_a[i],M_b[i]),sub = sub)/(M_b[i]-M_a[i])
    
    if (tag==3 || tag==4) {
      PROD=cbind(PROD,matrix(0, nrow = dim(PROD)[1], ncol=nbasis[2]-dim(PROD)[2]))
    }
    
    K=rbind(K,PROD)
  }
  
  res=A%*%K
  
  list(B=res, A=A, K=K, x_h=L_y, error=error, B_X=L_X, theta_x=L_theta, Phi=Phi)
  
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
  
  # if ( dim(Z)[2]!=length(G[[1]]) || dim(Z)[2]!=length(G[[2]]) ) {
  #   stop("The number of columns of 'Z' must be equal to the length of 'G'")
  # }
  
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
    p_values=NULL
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


# AUXILIARY FUNCTIONS FROM SOP PACKAGE MIGHT BE NEEDED

deviance2 <- function(C, w, sigma2, ssr, edf) {
  log_det_C <- determinant(C)$modulus
  deviance <- log_det_C + sum(log(sigma2*1/w)) + ssr/sigma2 + edf
  deviance  
}
deviance <- function(C, G, w, sigma2, ssr, edf) {
  log_det_C <- determinant(C)$modulus
  log_det_G <- determinant(G)$modulus
  deviance <- log_det_C + log_det_G + sum(log(sigma2*1/w)) + ssr/sigma2 + edf
  deviance  
}
###################################################### 
insert.na <- function(vec, na.ind) {
  aux <- rep(NA, l = length(na.ind))
  aux[!na.ind] <- vec
  aux
}
insert.matrix.na <- function(mat, na.ind) {
  aux <- matrix(NA, ncol = ncol(mat), nrow = length(na.ind))
  aux[!na.ind,] <- mat
  aux
}
add_coeff_to_terms <- function(object, nterms, b.fixed, b.random) {
  nlin <- nterms[1]
  nre <- nterms[2]
  nfun <- nterms[3]
  
  dim.fixed <- c(if(nlin > 0) object$lin$dim, if(nfun > 0) unlist(lapply(object$f, function(x) sum(x$Xmat$dim$fixed))))
  dim.random <- c(if(nre > 0) object$random$dim, if(nfun > 0) unlist(lapply(object$f, function(x) x$Xmat$dim$random)))
  
  e.f <- cumsum(dim.fixed)
  s.f <- e.f - dim.fixed + 1
  
  e.r <- cumsum(dim.random)
  s.r <- e.r - dim.random + 1
  
  if (nlin > 0) {
    object$lin$coeff <- b.fixed[2:(sum(object$lin$dim)+1)]
  }
  
  if (nre > 0) {
    object$random$coeff <- b.random[1:sum(object$random$dim)]
  }
  
  if (nfun > 0) {
    for (i.fun in 1:nfun) {
      if(ncol(object$f[[i.fun]]$Xmat$X) == 0) {
        vcum.f <- NULL
      } else {                      
        vcum.f <- s.f[i.fun + nlin]:e.f[i.fun + nlin]
      }
      vcum.r <- s.r[i.fun + nre]:e.r[i.fun + nre]
    }
    object$f[[i.fun]]$Xmat$coeff.fixed <- b.fixed[vcum.f + 1]
    object$f[[i.fun]]$Xmat$coeff.random <- b.random[vcum.r]
  }
  object
}
######################################################
interpret.f.formula <- function(formula) {
  env <- environment(formula) 
  if(inherits(formula, "character"))  {          
    formula <- as.formula(formula)
  }
  
  tf <- terms.formula(formula, specials = c("f"))
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  res <- eval(parse(text = terms[1]), envir = env)
  res
}
######################################################
construct.1D.pspline <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character")) {          
    formula <- as.formula(formula)
  }	
  
  res <- interpret.f.formula(formula)
  x1 <- data[ ,res$vars]
  type <- res$type
  
  if(!is.null(res$type) && res$type == "adaptive") {
    MM1 <- MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 5)
    smooth.comp <- paste("ad(", res$vars,")", sep = "")
  } else {
    MM1 <- MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 4)
    smooth.comp <- paste("f(", res$vars,")", sep = "")
  }  
  X <- MM1$X[,-1,drop = FALSE]
  Z <- MM1$Z
  d <- MM1$d
  B <- MM1$B
  c1 <- ncol(B)    
  
  x.fixed <-  ""
  for(i in 0:(res$pord[1]-1)){
    if(i == 1) 
      x.fixed <- c(x.fixed, res$vars)
    else if( i > 1)
      x.fixed <- c(x.fixed, paste(res$vars ,"^", i, sep = ""))
  }
  names.fixed <- x.fixed	
  names.random <- paste(smooth.comp, c(res$vars), sep = "|")    
  
  dim.random <- (c1 -res$pord[1])
  dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
  
  names(dim$fixed) <- names.fixed[-1]
  names(dim$random) <- paste(smooth.comp, "Global")
  
  
  v.comp <- 1:ncol(Z)/ncol(Z)
  if (!is.null(res$type) && res$type == "adaptive") {
    C <-  bbase(v.comp, min(v.comp), max(v.comp), res$nseg.sp, res$degree.sp)$B
  } else {
    C <- matrix(MM1$d, ncol = 1)
  }
  G.comp <- t(C)
  g <- vector(mode = "list", length = nrow(G.comp)) 
  for (i in 1:nrow(G.comp)) {
    g[[i]] <- G.comp[i,]
  }
  
  names(g) <- rep(names.random, length(g))
  colnames(X) <- names.fixed[-1]
  colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
  
  attr(dim$fixed, "random") <- rep(FALSE, length(dim$fixed))
  attr(dim$fixed, "smooth") <- rep(TRUE, length(dim$fixed))
  
  attr(dim$random, "random") <- rep(TRUE, length(dim$random)) 
  attr(dim$random, "smooth") <- rep(TRUE, length(dim$random))
  
  terms <- list()
  terms$MM <- MM1
  
  attr(terms, "term") <- smooth.comp    
  init.var <- rep(1, length(g))
  
  # Center the matrices
  cmX <- colMeans(X)
  cmZ <- colMeans(Z)
  
  X <- sweep(X, 2, cmX)
  Z <- sweep(Z, 2, cmZ)
  
  res <- list(X = X, Z = Z, g = g, init.var = init.var, dim = dim, terms = terms, edflabel = names.random, cm = list(X = cmX, Z = cmZ))
  res
}
######################################################
######################################################
construct.2D.pspline <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character")) {
    formula <- as.formula(formula)
  }
  
  res <- interpret.f.formula(formula)    
  x1 <- data[ ,res$vars[1]]
  x2 <- data[ ,res$vars[2]]    
  
  type <- res$type
  
  MM1 <- MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 4)
  MM2 <- MM.basis(x2, min(x2), max(x2), res$nseg[2], res$degree[2], res$pord[2], 4)    
  
  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B    
  c1 <- ncol(B1); c2 <- ncol(B2)
  
  # Very dirty: names of the fixed part (null) part of the 2D P-spline:
  # Depend on the penatly order assume for each covariate.
  # For intance, if we assume pord = 2 for x1 (i.e. X1 = [1|x1] and pord = 3 for x2 (X2 = [1|x2|x2^2]), then we have that
  # X = [1 | x1 | x2 | x1:x2| x2^2 | x1:x2^2]  
  
  x.fixed <- y.fixed <- ""
  for(i in 0:(res$pord[1]-1)){
    if(i == 1) 
      x.fixed <- c(x.fixed, res$vars[1])
    else if( i > 1)
      x.fixed <- c(x.fixed, paste(res$vars[1], "^", i, sep = ""))
  }
  for(i in 0:(res$pord[2]-1)){
    if(i == 1) 
      y.fixed <- c(y.fixed, res$vars[2])
    else if( i > 1)
      y.fixed <- c(y.fixed, paste(res$vars[2], "^", i, sep = ""))
  }
  xy.fixed <- NULL
  for(i in 1:length(y.fixed)) {
    xy.fixed <- c(xy.fixed, paste(y.fixed[i], x.fixed, sep= ""))
  }
  xy.fixed <- xy.fixed[xy.fixed != ""]
  names.fixed <- xy.fixed
  
  smooth.comp <- paste("f(", res$vars[1],",", res$vars[2],")", sep = "")
  names.random <- paste(smooth.comp, c(res$vars[1], res$vars[2]), sep = "|")
  
  X <- Rten2(X2, X1)
  
  # Delete the intercept
  X <- X[,-1,drop = FALSE]
  Z <- cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2, Z1))
  
  dim.random <- c((c1 -res$pord[1])*res$pord[2] , (c2 - res$pord[2])*res$pord[1] , (c1 - res$pord[1])*(c2 - res$pord[2]))  	
  dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
  names(dim$fixed) <- names.fixed
  names(dim$random) <- paste(smooth.comp, "Global")
  
  # Variance/Covariance components: two variances, one for each covariate
  g1u <- rep(1, res$pord[2])%x%d1
  g2u <- d2%x%rep(1, res$pord[1])
  g1b <- rep(1,c2 - res$pord[2])%x%d1
  g2b <- d2%x%rep(1,c1 - res$pord[1])  
  g <- list()	
  g[[1]] <- c(g1u, rep(0, dim.random[2]), g1b)
  g[[2]] <- c(rep(0, dim.random[1]), g2u, g2b)    
  names(g) <- names.random
  
  colnames(X) <- names.fixed
  colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
  
  attr(dim$fixed, "random") <- rep(FALSE, length(dim$fixed))
  attr(dim$fixed, "smooth") <- rep(TRUE, length(dim$fixed))
  
  attr(dim$random, "random") <- rep(TRUE, length(dim$random)) 
  attr(dim$random, "smooth") <- rep(TRUE, length(dim$random))
  
  terms <- list()
  terms$MM <- list(MM1 = MM1, MM2 = MM2)    
  attr(terms, "term") <- smooth.comp
  
  # Initialize variance components
  init.var <- rep(1, length(g))
  
  # Center the matrices
  cmX <- colMeans(X)
  cmZ <- colMeans(Z)
  
  X <- sweep(X, 2, cmX)
  Z <- sweep(Z, 2, cmZ)
  
  res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var, terms = terms, cm = list(X = cmX, Z = cmZ))	
  res
}
######################################################
######################################################
construct.3D.pspline <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character"))  {          
    formula <- as.formula(formula)
  }
  
  res <- interpret.f.formula(formula)
  x1 <- data[ ,res$vars[1]]
  x2 <- data[ ,res$vars[2]]    
  x3 <- data[ ,res$vars[3]]    
  
  type <- res$type
  intercept <- TRUE
  
  MM1 = MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 4 , intercept)
  MM2 = MM.basis(x2, min(x2), max(x2), res$nseg[2], res$degree[2], res$pord[2], 4 , intercept)    
  MM3 = MM.basis(x3, min(x3), max(x3), res$nseg[3], res$degree[3], res$pord[3], 4 , intercept)
  
  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B    
  X3 <- MM3$X; Z3 <- MM3$Z; d3 <- MM3$d; B3 <- MM3$B
  
  c1 <- ncol(B1)
  c2 <- ncol(B2)
  c3 = ncol(B3)      
  
  x.fixed <- y.fixed <- z.fixed <- ""
  
  for(i in 0:(res$pord[1]-1)){
    if(i == 1) 
      x.fixed <- c(x.fixed, res$vars[1])
    else if( i > 1)
      x.fixed <- c(x.fixed, paste(res$vars[1], "^", i, sep = ""))
  }
  for(i in 0:(res$pord[2]-1)){
    if(i == 1) 
      y.fixed <- c(y.fixed, res$vars[2])
    else if( i > 1)
      y.fixed <- c(y.fixed, paste(res$vars[2], "^", i, sep = ""))
  }
  for(i in 0:(res$pord[3]-1)){
    if(i == 1) 
      z.fixed <- c(z.fixed, res$vars[3])
    else if( i > 1)
      z.fixed <- c(z.fixed, paste(res$vars[3], "^", i, sep = ""))
  }	
  xy.fixed <- NULL
  for(i in 1:length(x.fixed)) {
    for(j in 1:length(y.fixed)) {
      xy.fixed <- c(xy.fixed, paste(z.fixed, y.fixed[j], x.fixed[i], sep= ""))
    }
  }
  xy.fixed <- xy.fixed[xy.fixed != ""]	
  names.fixed <- xy.fixed
  
  smooth.comp <- paste("f(", res$vars[1],",", res$vars[2],",", res$vars[3],")", sep = "")
  names.random <- paste(smooth.comp, c(res$vars[1], res$vars[2], res$vars[3]), sep = "|")
  
  rx12 <- Rten2(X1,X2)			
  X <- Rten2(rx12, X3)
  # Delete the intercept
  X <- X[,-1,drop = FALSE]
  colnames(X) <- names.fixed
  
  Z <- cbind(Rten2(Rten2(Z1,X2), X3), 
             Rten2(Rten2(X1,Z2), X3),
             Rten2(rx12, Z3),
             Rten2(Rten2(Z1,Z2), X3),
             Rten2(Rten2(Z1,X2), Z3),
             Rten2(Rten2(X1,Z2), Z3),
             Rten2(Rten2(Z1,Z2), Z3))
  
  dim.random <- c((c1 -res$pord[1])*res$pord[2] ,
                  (c2 - res$pord[2])*res$pord[1] ,
                  (c1 - res$pord[1])*(c2 - res$pord[2]))		
  dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
  names(dim$fixed)  <- names.fixed
  names(dim$random) <- paste(smooth.comp, "Global")
  np <-  c(prod(res$pord),
           (c1-res$pord[1])*res$pord[2]*res$pord[3],
           (c2-res$pord[2])*res$pord[1]*res$pord[3],
           (c3-res$pord[3])*res$pord[1]*res$pord[2],
           (c1-res$pord[1])*(c2-res$pord[2])*res$pord[3],
           (c1-res$pord[1])*(c3-res$pord[3])*res$pord[2],
           (c2-res$pord[2])*(c3-res$pord[3])*res$pord[1],
           (c1-res$pord[1])*(c2-res$pord[2])*(c3-res$pord[3])) 
  
  # Variance/Covariance components: two variances, one for each covariate
  r1<-rep(1,res$pord[1])
  r2<-rep(1,res$pord[2])
  r3<-rep(1,res$pord[3])
  d1u <- d1%x%r1%x%r2
  d2u <- r1%x%d2%x%r3
  d3u <- r1%x%r2%x%d3
  
  r01<-rep(1,c1-res$pord[1])
  r02<-rep(1,c2-res$pord[2])
  r03<-rep(1,c3-res$pord[3])
  d11b <- d1%x%r02%x%r3
  d12b <- d1%x%r2%x%r03
  
  d21b <- r01%x%d2%x%r3
  d22b <- r1%x%d2%x%r03
  
  d31b <- r01%x%r2%x%d3
  d32b <- r1%x%r02%x%d3
  
  d1t <- d1%x%r02%x%r03
  d2t <- r01%x%d2%x%r03
  d3t <- r01%x%r02%x%d3
  
  # Number of parameters in each part: 8, 44, 44, 44, 242, 242, 242, 1331
  
  D <- diag(c(rep(0,np[1]), rep(1,sum(np[-1]))))
  
  G1inv.n <- c(d1u, rep(0, sum(np[3:4])), d11b, d12b, rep(0, np[7]), d1t)
  G2inv.n <- c(rep(0, np[2]), d2u, rep(0, np[4]), d21b, rep(0, np[6]), d22b, d2t)
  G3inv.n <- c(rep(0, sum(np[2:3])), d3u, rep(0, np[5]), d31b, d32b, d3t)
  
  g <- list(G1inv.n, G2inv.n, G3inv.n)
  names(g) <- names.random
  colnames(X) <- names.fixed
  colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
  
  attr(dim$fixed, "random") <- rep(FALSE, length(dim$fixed))
  attr(dim$fixed, "smooth") <- rep(TRUE, length(dim$fixed))
  
  attr(dim$random, "random") <- rep(TRUE, length(dim$random)) 
  attr(dim$random, "smooth") <- rep(TRUE, length(dim$random))
  
  terms <- list()
  terms$MM <- list(MM1 = MM1, MM2 = MM2, MM3 = MM3)
  attr(terms, "term") <- smooth.comp
  
  # Initialize variance components
  init.var <- rep(1, length(g))
  
  # Center the matrices
  cmX <- colMeans(X)
  cmZ <- colMeans(Z)
  
  X <- sweep(X, 2, cmX)
  Z <- sweep(Z, 2, cmZ)
  
  res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var, terms = terms, cm = list(X = cmX, Z = cmZ))	
  res
}
######################################################
######################################################
construct.random.part <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character"))          
    formula <- as.formula(formula)
  
  mf <- model.frame(formula, data=data, drop.unused.levels = TRUE, na.action = NULL)
  mt <- terms(mf)    
  f.terms <- attr(mt, "term.labels")[attr(mt,"dataClasses") == "factor"]
  Z <- model.matrix(mt, data = mf, contrasts.arg = lapply(mf[,f.terms, drop = FALSE], contrasts, contrasts=FALSE))
  Z[is.na(Z)] <- 0
  
  attr(mt, "contrast") <- attr(Z,"contrast")
  attr(mt, "xlev") <- .getXlevels(mt, mf)
  
  
  dim <- table(attr(Z,"assign"))[-1]
  
  e <- cumsum(dim)
  s <- e - dim + 1
  
  g <- list()
  for(i in 1:length(dim)) {
    g[[i]] <- rep(0,sum(dim))
    g[[i]][s[i]:e[i]] <- 1
  }
  names(g) <- names(dim) <- attr(mt,"term.labels")
  attr(dim, "random") <- rep(TRUE, length(dim)) 
  attr(dim, "smooth") <- rep(FALSE, length(dim))
  
  # Initialize variance components
  init.var <- rep(0.01, length(g))
  
  res <- list(Z = Z[,-1, drop = FALSE], dim = dim, g = g, init.var = init.var, terms = mt)
  res
}
########################################################
########################################################	
construct.fixed.part <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character"))          
    formula <- as.formula(formula)
  
  mf <- model.frame(formula, data, drop.unused.levels = TRUE)
  mt <- terms(mf)   
  X <- model.matrix(mt, mf)
  
  dim <- table(attr(X,"assign"))[-1]
  names(dim) <- attr(mt, "term.labels")
  
  attr(mt, "contrast") <- attr(X,"contrast")
  attr(mt, "xlev") <- .getXlevels(mt, mf)
  
  var.aux <- attr(mt, "term.labels")[attr(mt, "order") == 1]
  dataClasses <- attr(mt, "dataClasses")
  int <- attr(mt, "term.labels")[attr(mt, "order") != 1]
  for (i in int){
    var.aux <- c(var.aux, i)
    if(any(dataClasses[strsplit(i, ":")[[1]]] == "factor")) {
      dataClasses <- c(dataClasses, "factor")     
    } else {
      dataClasses <- c(dataClasses, "numeric")
      names(dataClasses) <- var.aux
    }
  }
  attr(mt, "dataClasses") <- dataClasses
  
  attr(dim, "random") <-  attr(dim, "smooth") <- rep(FALSE, length(dim)) 	
  res <- list(X = X[,-1, drop = FALSE], dim = dim, terms = mt)
  res	
}
########################################################
########################################################
########################################################	
rae <- function(x) {
  args <- match.call()
  res <- args$x
  res    	
}
######################################################## 
########################################################  
f <- function (..., nseg = 10, pord = 2, degree = 3) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  
  if (length(nseg)<d) nseg=rep(nseg,d)
  if (length(pord)<d) pord=rep(pord,d)
  if (length(degree)<d) degree=rep(degree,d)
  
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (term[1] == ".") {
    stop("f(.) not yet supported.")
  }	
  if (d > 1) { 
    for (i in 2:d) {
      term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      if (term[i] == ".") { 
        stop("f(.) not yet supported.")
      }
    }
  }
  for (i in 1:d){
    term[i] <- attr(terms(reformulate(term[i])), "term.labels")
  } 
  nseg.new <- round(nseg)
  if (all.equal(nseg.new,nseg) != TRUE) {
    warning("argument nseg of f() should be integer and has been rounded")
  }
  nseg <- nseg.new
  pord.new <- round(pord)
  if (all.equal(pord.new,pord) != TRUE) {
    warning("argument pord of f() should be integer and has been rounded")
  }
  pord <- pord.new
  
  if (length(unique(term)) != d) { 
    stop("Repeated variables as arguments of a smooth are not permitted")
  }	
  full.call <- paste("f(", term[1], sep = "")
  if (d > 1) { 
    for (i in 2:d) {
      full.call <- paste(full.call, ",", term[i], sep = "")
    }
  }	
  label <- paste(full.call, ")", sep = "")
  ret <- list(vars = term, nseg = nseg, pord = pord, degree = degree, dim = d, label = label)
  ret
}
########################################################
########################################################  
ad <- function (..., nseg = 10, pord = 2 , degree = 3, nseg.sp = 5, degree.sp = 3) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  if (d > 1) {
    stop("Adaptive B-splines for more than one dimension are not yet supported.")
  }  
  if (length(nseg) < d) nseg <- rep(nseg, d)
  if (length(pord) < d) pord <- rep(pord, d)
  if (length(degree) < d) degree <- rep(degree, d)
  
  if (length(nseg.sp) < d) nseg.sp <- rep(nseg.sp, d)
  if (length(degree.sp) < d) degree.sp <- rep(degree.sp, d)  
  
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (term[1] == ".") {
    stop("f(.) not yet supported.")
  } 
  
  for (i in 1:d){
    term[i] <- attr(terms(reformulate(term[i])), "term.labels")
  }
  
  nseg.new <- round(nseg)
  if (all.equal(nseg.new, nseg) != TRUE) {
    warning("argument nseg of ad() should be integer and has been rounded")
  }
  nseg <- nseg.new
  
  nseg.sp.new <- round(nseg.sp)
  if (all.equal(nseg.sp.new, nseg.sp) != TRUE) {
    warning("argument nseg.sp of ad() should be integer and has been rounded")
  }
  nseg.sp <- nseg.sp.new
  
  pord.new <- round(pord)
  if (all.equal(pord.new, pord) != TRUE) {
    warning("argument pord of ad() should be integer and has been rounded")
  }
  pord <- pord.new
  
  label <- paste0("ad(", term[1], ")")
  ret <- list(vars = term, nseg = nseg, pord = pord, degree = degree, nseg.sp = nseg.sp, degree.sp = degree.sp, dim = d, label = label, type = "adaptive")
  ret
}
#######################################################################
#######################################################################
tpower <- function(x, t, p) {
  # Function for truncated p-th power function
  return((x - t) ^ p * (x > t))
}

#######################################################################
#######################################################################
spline.bbase<-function (knots, X., bdeg, eps = 1e-05) {
  if(is.null(attributes(knots)$dx))
    dx <- mean(diff(knots))
  else dx <- attributes(knots)$dx
  P <- outer(X., knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1) * 
                                         dx^bdeg)
  B <- (-1)^(bdeg + 1) * P %*% t(D)
  B[B < eps] = 0
  B
}

##################################################################
##################################################################
bbase <- function(x, xl = min(x), xr = max(x), ndx = 10,   bdeg = 3, eps = 1e-5) {
  # Function for B-spline basis
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)  
  attr(knots,"dx")<-dx
  B <- spline.bbase(knots, x, bdeg, eps)  
  res <- list(B = B, knots = knots)
  res 
}
#######################################################################
#######################################################################
MM.basis <- function (x, xl, xr, ndx, bdeg, pord, decom = 1, intercept = TRUE) {
  #print("MM.basis")  
  Bb = bbase(x,xl,xr,ndx,bdeg)
  knots <- Bb$knots
  B = Bb$B
  m = ncol(B)
  n = nrow(B)
  D = diff(diag(m), differences=pord)
  P.svd = svd(crossprod(D))
  U.Z = (P.svd$u)[,1:(m-pord)]  # eigenvectors
  d = (P.svd$d)[1:(m-pord)]     # eigenvalues
  Z = B%*%U.Z
  U.X = NULL
  if(decom == 1) {
    U.X = ((P.svd$u)[,-(1:(m-pord))])
    X = B%*%U.X
  } else if (decom == 2){
    X = NULL
    for(i in 0:(pord-1)){
      X = cbind(X,x^i)
    }
  } else if(decom == 3) {
    U.X = NULL
    for(i in 0:(pord-1)){
      U.X = cbind(U.X, knots[-c((1:(bdeg - 1)),(length(knots)- (bdeg - 1) + 1):length(knots))]^i)
    }
    X = B%*%U.X
  } else if(decom == 4) { # Wood's 2013
    X = B%*%((P.svd$u)[,-(1:(m-pord))])
    #id.v <- rep(1, nrow(X))
    #D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
  } else if(decom == 5) { # Paul's parameterization
    U.Z <- (t(D)%*%solve(D%*%t(D)))
    Z <- B%*%U.Z
    X = B%*%((P.svd$u)[,-(1:(m-pord))])
    #id.v <- rep(1, nrow(X))
    #D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf					
  } else if(decom == 6) { # martin's approach
    U.Z <- t(D)
    Z <- B%*%U.Z  
    X = B%*%((P.svd$u)[,-(1:(m-pord))])
    #id.v <- rep(1, nrow(X))
    #D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf			
  }
  if (!intercept) X <- X[,-1,drop = FALSE]
  list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
}
#######################################################################
#######################################################################
construct.block <- function(A1,A2,A3,A4) {
  block <- rbind(cbind(A1,A2), cbind(A3,A4))
  return(block)
}
construct.capital.lambda <- function(g) {
  length.eq <- all(sapply(g, function(x) {
    diff(range(unlist(lapply(x, length)))) < .Machine$double.eps ^ 0.5
  }))
  if(length.eq) {
    l <- length(g)
    if(l == 1) {
      if(length(g[[1]]) == 1) {
        res <- g
      } else {
        res <- do.call("c", lapply(g, function(x) x))
      }
    } else {
      dim <- sapply(g, function(x) {
        if(is.list(x))
          unlist(lapply(x, length))[1]
        else
          length(x)
      })		
      end <- cumsum(dim)
      init <- end - dim + 1
      
      res <- do.call("c", lapply(1:length(g), function(x, g, init, end, dim) {
        temp <- g[[x]]
        if(is.list(temp)) {
          lapply(temp, function(y, x, dim) {
            aux <- rep(0, l = sum(dim))
            aux[init[x]:end[x]] <- y
            aux
          }, x = x, dim = dim)
        } else {
          aux <- rep(0, l = sum(dim))
          aux[init[x]:end[x]] <- temp
          list(aux)
        }
      }, g = g, init = init, end = end, dim = dim))	
    }
  } else {
    stop("Error in construct.capital.lambda")
  }	
  res
}
###################################################################
###################################################################
# Model matrices and GLAM
###################################################################
###################################################################
Rten <- function(X) {
  one <- matrix(1, 1, ncol(X))
  kronecker(X,one)*kronecker(one,X)
}
###################################################################
###################################################################
Rten2 <- function(X1,X2) {
  one.1 <- matrix(1,1,ncol(X1))
  one.2 <- matrix(1,1,ncol(X2))
  kronecker(X1,one.2)*kronecker(one.1,X2)
}
###################################################################
###################################################################
H <- function(X,A) {
  d <- dim(A)
  M <- matrix(A, nrow = d[1])
  XM <- X%*%M
  array(XM, c(nrow(XM),d[-1]))
}
###################################################################
###################################################################
Rotate <- function(A) {
  d <- 1:length(dim(A))
  d1 <- c(d[-1],d[1])
  aperm(A, d1)
}
###################################################################
###################################################################
RH <- function(X,A) {
  Rotate(H(X,A))
}
###################################################################
###################################################################
A1.form <- function(l, w = NULL){
  d <- length(l)
  n <- rev(sapply(l,nrow))
  c <- rev(sapply(l,ncol))
  if (is.null(w)) {
    W <- array(1, n)
  } else {
    W <- array(w, n)
  }
  tmp <- RH(t(Rten(l[[d]])), W)
  for (i in (d-1):1) {
    tmp <- RH(t(Rten(l[[i]])),tmp)
  }
  dim(tmp)<- rep(c, rep(2,d))
  Fast1 <- aperm(tmp, as.vector(matrix(1:(d*2), byrow = TRUE, ncol = 2)))
  Fast <- if(prod(c)) matrix(Fast1, nrow = prod(c)) else Fast1
  return(Fast)
}
###################################################################
###################################################################
A2.form <- function(l1, l2, w = NULL) {
  d1 <- length(l1)
  d2 <- length(l2)
  if(!(d1 == d2)) {
    stop("l1 and l2 should have the same dimension")
  }
  n <- rev(sapply(l1, nrow))
  d <- rev(sapply(l1, ncol))
  c <- rev(sapply(l2, ncol))
  
  if (is.null(w)) {
    W <- array(1, n)
  } else {
    W <- array(w, n)
  }
  tmp <- RH(t(Rten2(l2[[d1]], l1[[d1]])), W)
  for (i in (d1-1):1) {
    tmp <- RH(t(Rten2(l2[[i]], l1[[i]])),tmp)
  }
  dim(tmp)<- as.vector(rbind(d,c))
  Fast1 <- aperm(tmp, as.vector(matrix(1:(d1*2), byrow = TRUE, ncol = 2)))
  Fast <- if(prod(d)) matrix(Fast1, nrow = prod(d)) else Fast1
  return(Fast)
}
###################################################################
###################################################################
XtX <- function(X, w = NULL) {
  A1.form(X, w)
}
###################################################################
###################################################################
XtZ <- function(X, Z, w = NULL) {
  d <- length(Z)
  res <- NULL
  for (i in 1:d) {
    res <- cbind(res, A2.form(X, Z[[i]], w))
  }
  res
}
###################################################################
###################################################################
ZtZ <- function(Z, w = NULL) {
  d <- length(Z)
  upper <- list()
  for(i in 1:d) {
    upper[[i]] <- list()
    upper[[i]][[i]] <- A1.form(Z[[i]], w)
  }
  # Obtain the elements of the matrix
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      upper[[i]][[j]] <- A2.form(Z[[i]], Z[[j]], w)
    }
  }
  # Create the matrix
  res <- NULL
  for (i in 1:d) {
    if( i == 1) {
      res <- do.call("cbind", upper[[1]])
    } else {
      tmp <- do.call("cbind", upper[[i]])
      for(j in (i-1):1) {
        if(length(upper[[j]][[i]]))
          tmp <- cbind(t(upper[[j]][[i]]), tmp)
      }
      if(nrow(tmp))
        res <- rbind(res, tmp)
    }
  }
  res
}
###################################################################
###################################################################
Xty <- function(X, y, w = NULL) {
  d <- length(X)
  n <- rev(sapply(X, nrow))
  if(is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w*y, n)
  }
  tmp <- RH(t(X[[d]]),Y)
  for(i in (d-1):1) {
    tmp <- RH(t(X[[i]]), tmp)
  }
  as.vector(tmp)
}
###################################################################
###################################################################
Zty <- function(Z, y, w = NULL) {
  d <- length(Z)
  n <- rev(sapply(Z[[1]], nrow))
  if(is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w*y, n)
  }
  res <- NULL
  for(i in 1:d) {
    k <- length(Z[[i]])
    tmp <- RH(t(Z[[i]][[k]]),Y)
    for(j in (k-1):1) {
      tmp <- RH(t(Z[[i]][[j]]), tmp)
    }
    res <- c(res, as.vector(tmp))
  }
  res
}
###################################################################
###################################################################
Xtheta <- function(X, theta) {
  d <- length(X)
  n <- rev(sapply(X, ncol))
  Theta <- array(theta, n)
  tmp <- RH(X[[d]], Theta)
  for(i in (d-1):1) {
    tmp <- RH(X[[i]], tmp)
  }
  as.vector(tmp)
}
###################################################################
###################################################################
Ztheta <- function(Z, theta, np) {
  d <- length(Z)
  for(i in 1:d) {
    if (i == 1) {
      res <- Xtheta(Z[[i]], theta[1:(np[1])])
    } else {
      init <- sum(np[1:(i-1)])
      fin  <- np[i]
      if(fin) res <- res + Xtheta(Z[[i]], theta[(init+1):(init+fin)])
    }
  }
  res
}
###################################################################
###################################################################
construct.matrices <- function(X, Z, z, w, GLAM) {
  
  # print(dim(X))
  # print(dim(w))
  # 
  
  if(GLAM) {
    XtX. <- XtX(X,w) 
    XtZ. <- XtZ(X,Z,w)
    ZtX. <- t(XtZ.)
    ZtZ. <- ZtZ(Z,w)
    Zty. = Zty(Z,z,w)
    Xty. = Xty(X,z,w)
    yty. <- sum((z^2)*w)
    ZtXtZ = rbind(XtZ., ZtZ.)
    u <- c(Xty.,Zty.)
  } else {
    XtW. = t(X*w)
    XtX. = XtW.%*%X
    XtZ. = XtW.%*%Z
    ZtX. = t(XtZ.)
    ZtW. =  t(Z*w)
    ZtZ. = ZtW.%*%Z
    Xty. = XtW.%*%z
    Zty. = ZtW.%*%z
    yty. <- sum((z^2)*w)
    ZtXtZ = rbind(XtZ., ZtZ.)
    u <- c(Xty.,Zty.)
  }
  res <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
}

