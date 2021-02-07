#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Christophe Dutang
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the
#   Free Software Foundation, Inc.,
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#
#############################################################################
### maximum likelihood estimation for censored or non-censored data
###
###         R functions
###
### many ideas are taken from the fitdistr function of the MASS package and
### the mle function of the stat package.





#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Regis Pouillot, Jean-Baptiste Denis, Christophe Dutang
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the
#   Free Software Foundation, Inc.,
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#
#############################################################################
### fit parametric distributions for non-censored data
###
###         R functions
###

fitdist <- function (data, distr, method = c("mle", "mme", "qme", "mge", "mse"), start=NULL,
                     fix.arg=NULL, discrete, keepdata = TRUE, keepdata.nb=100, ...)
{
  ##HelprFuncSTART
  ##HelprFuncSTART
  ##HelprFuncSTART
  mledist <- function (data, distr, start=NULL, fix.arg=NULL, optim.method="default",
      lower=-Inf, upper=Inf, custom.optim=NULL, weights=NULL, silent=TRUE, gradient=NULL,
      checkstartfix=FALSE, ...)
      # data may correspond to a vector for non censored data or to
      # a dataframe of two columns named left and right for censored data
  {
      if (!is.character(distr))
          stop("distr must be a character string naming a distribution")
      else
          distname <- distr
      ddistname <- paste("d", distname, sep="")
      argddistname <- names(formals(ddistname))

      if (!exists(ddistname, mode="function"))
          stop(paste("The ", ddistname, " function must be defined"))
      if(is.null(custom.optim))
        optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))

      start.arg <- start #to avoid confusion with the start() function of stats pkg (check is done lines 87-100)
      if(is.vector(start.arg)) #backward compatibility
        start.arg <- as.list(start.arg)

      txt1 <- "data must be a numeric vector of length greater than 1 for non censored data"
      txt2 <- "or a dataframe with two columns named left and right and more than one line for censored data"
      if(!is.null(weights))
      {
        if(any(weights < 0))
          stop("weights should be a vector of integers greater than 0")
        if(!is.allint.w(weights))
          stop("weights should be a vector of (strictly) positive integers")
        if(length(weights) != NROW(data))
          stop("weights should be a vector with a length equal to the observation number")
        warning("weights are not taken into account in the default initial values")
      }

      if (is.vector(data)) {
          cens <- FALSE
          if (!(is.numeric(data) & length(data)>1))
              stop(paste(txt1, txt2))
      }
      else {
          cens <- TRUE
          censdata <- data
          if (!(is.vector(censdata$left) & is.vector(censdata$right) & length(censdata[, 1])>1))
              stop(paste(txt1, txt2))
          pdistname<-paste("p", distname, sep="")
          if (!exists(pdistname, mode="function"))
              stop(paste("The ", pdistname, " function must be defined to apply maximum likelihood to censored data"))

      }

      if (cens) {
          #format data for calculation of starting values and fitting process
          dataformat <- cens2pseudo(censdata)
          data <- dataformat$pseudo
          rcens <- dataformat$rcens; lcens <- dataformat$lcens
          icens <- dataformat$icens; ncens <- dataformat$ncens

          irow <- cens2idxrow(censdata)
          irow.rcens <- irow$rcens; irow.lcens <- irow$lcens
          irow.icens <- irow$icens; irow.ncens <- irow$ncens
      }

      if(!checkstartfix) #pre-check has not been done by fitdist() or bootdist()
      {
        # manage starting/fixed values: may raise errors or return two named list
        arg_startfix <- manageparam(start.arg=start, fix.arg=fix.arg, obs=data,
                                    distname=distname)

        #check inconsistent parameters
        hasnodefaultval <- sapply(formals(ddistname), is.name)
        arg_startfix <- checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg,
                                       argddistname, hasnodefaultval)
        #arg_startfix contains two names list (no longer NULL nor function)

        #set fix.arg.fun
        if(is.function(fix.arg))
          fix.arg.fun <- fix.arg
        else
          fix.arg.fun <- NULL
      }else #pre-check has been done by fitdist<cens>() or bootdist<cens>()
      {
        arg_startfix <- list(start.arg=start, fix.arg=fix.arg)
        fix.arg.fun <- NULL
      }

      #unlist starting values as needed in optim()
      vstart <- unlist(arg_startfix$start.arg)
      #sanity check
      if(is.null(vstart))
        stop("Starting values could not be NULL with checkstartfix=TRUE")

      #erase user value
      #(cannot coerce to vector as there might be different modes: numeric, character...)
      fix.arg <- arg_startfix$fix.arg


      ############# closed-form formula for uniform distribution ##########
      if(distname == "unif")
      {
        if(length(fix.arg) >= 2)
        {
          stop("'fix.arg' sets all distribution parameters without any parameter to estimate.")
        }else if(length(fix.arg) == 1)
        {
          if(names(fix.arg) == "min")
            par <- c(max=max(data))
          else if(names(fix.arg) == "max")
            par <- c(min=min(data))
          else
            stop("'fix.arg' must specify names which are arguments to 'distr'.")
        }else
          par <- c(min=min(data), max=max(data))

          myarg <- c(list(data), as.list(par), as.list(fix.arg))
          loglikval <- sum(log(do.call(dunif, myarg)))
          res <- list(estimate = par[!names(par) %in% names(fix.arg)], convergence = 0,
                      loglik = loglikval, meth = "closed formula",
                      hessian = NA, optim.function= NA, fix.arg = fix.arg)

          return(res)
      }


      ############# MLE fit using optim or custom.optim ##########

      # definition of the function to minimize : - log likelihood
      # for non censored data
      if (!cens && is.null(weights)) {
          # the argument names are:
          # - par for parameters (like in optim function)
          # - fix.arg for optional fixed parameters
          # - obs for observations (previously dat but conflicts with genoud data.type.int argument)
          # - ddistnam for distribution name
          if ("log" %in% argddistname){
              fnobj <- function(par, fix.arg, obs, ddistnam){
                  -sum(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg), log=TRUE) ) )
              }
          }
          else{
          fnobj <- function(par, fix.arg, obs, ddistnam) {
              -sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
              }
          }
      }
      else if(cens && is.null(weights)) #censored data
      {
        argpdistname<-names(formals(pdistname))
          if (("log" %in% argddistname) & ("log.p" %in% argpdistname))
              fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam)
                  -sum(do.call(ddistnam, c(list(ncens), as.list(par), as.list(fix.arg), list(log=TRUE)))) -
                  sum(do.call(pdistnam, c(list(lcens), as.list(par), as.list(fix.arg), list(log=TRUE)))) -
                  sum(do.call(pdistnam, c(list(rcens), as.list(par), as.list(fix.arg), list(lower.tail=FALSE), list(log=TRUE)))) -
                  sum(log(do.call(pdistnam, c(list(icens$right), as.list(par), as.list(fix.arg))) - # without log=TRUE here
                  do.call(pdistnam, c(list(icens$left), as.list(par), as.list(fix.arg))) )) # without log=TRUE here
          else
              fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam)
                  -sum(log(do.call(ddistnam, c(list(ncens), as.list(par), as.list(fix.arg))))) -
                  sum(log(do.call(pdistnam, c(list(lcens), as.list(par), as.list(fix.arg))))) -
                  sum(log(1-do.call(pdistnam, c(list(rcens), as.list(par), as.list(fix.arg))))) -
                  sum(log(do.call(pdistnam, c(list(icens$right), as.list(par), as.list(fix.arg))) -
                  do.call(pdistnam, c(list(icens$left), as.list(par), as.list(fix.arg))) ))
      }else if(!cens && !is.null(weights))
      {
          fnobj <- function(par, fix.arg, obs, ddistnam) {
            -sum(weights * log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
          }
      }else if(cens && !is.null(weights))
      {
        fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam)
        {
          p1 <- log(do.call(ddistnam, c(list(ncens), as.list(par), as.list(fix.arg))))
          p2 <- log(do.call(pdistnam, c(list(lcens), as.list(par), as.list(fix.arg))))
          p3 <- log(1-do.call(pdistnam, c(list(rcens), as.list(par), as.list(fix.arg))))
          p4 <- log(do.call(pdistnam, c(list(icens$right), as.list(par), as.list(fix.arg))) -
                      do.call(pdistnam, c(list(icens$left), as.list(par), as.list(fix.arg))) )
          -sum(weights[irow.ncens] * p1) -
            sum(weights[irow.lcens] * p2) -
            sum(weights[irow.rcens] * p3) -
            sum(weights[irow.icens] * p4)
        }
      }

      #get warning value
      owarn <- getOption("warn")

      # Try to minimize the minus (log-)likelihood using the base R optim function
      if(is.null(custom.optim))
      {
          hasbound <- any(is.finite(lower) | is.finite(upper))

          # Choice of the optimization method
          if (optim.method == "default")
          {
            meth <- ifelse(length(vstart) > 1, "Nelder-Mead", "BFGS")
          }else
            meth <- optim.method

          if(meth == "BFGS" && hasbound && is.null(gradient))
          {
            meth <- "L-BFGS-B"
            txt1 <- "The BFGS method cannot be used with bounds without provided the gradient."
            txt2 <- "The method is changed to L-BFGS-B."
            warning(paste(txt1, txt2))
          }

          options(warn=ifelse(silent, -1, 0))
          #select optim or constrOptim
          if(hasbound) #finite bounds are provided
          {
            if(!is.null(gradient))
            {
              opt.fun <- "constrOptim"
            }else #gradient == NULL
            {
              if(meth == "Nelder-Mead")
                opt.fun <- "constrOptim"
              else if(meth %in% c("L-BFGS-B", "Brent"))
                opt.fun <- "optim"
              else
              {
                txt1 <- paste("The method", meth, "cannot be used by constrOptim() nor optim() without gradient and bounds.")
                txt2 <- "Only optimization methods L-BFGS-B, Brent and Nelder-Mead can be used in such case."
                stop(paste(txt1, txt2))
              }
            }
            if(opt.fun == "constrOptim")
            {
              #recycle parameters
              npar <- length(vstart) #as in optim() line 34
              lower <- as.double(rep_len(lower, npar)) #as in optim() line 64
              upper <- as.double(rep_len(upper, npar))

              # constraints are : Mat %*% theta >= Bnd, i.e.
              # +1 * theta[i] >= lower[i];
              # -1 * theta[i] >= -upper[i]

              #select rows from the identity matrix
              haslow <- is.finite(lower)
              Mat <- diag(npar)[haslow, ]
              #select rows from the opposite of the identity matrix
              hasupp <- is.finite(upper)
              Mat <- rbind(Mat, -diag(npar)[hasupp, ])
              colnames(Mat) <- names(vstart)
              rownames(Mat) <- paste0("constr", 1:NROW(Mat))

              #select the bounds
              Bnd <- c(lower[is.finite(lower)], -upper[is.finite(upper)])
              names(Bnd) <- paste0("constr", 1:length(Bnd))

              initconstr <- Mat %*% vstart - Bnd
              if(any(initconstr < 0))
                stop("Starting values must be in the feasible region.")

              if(!cens)
              {
                opttryerror <- try(opt <- constrOptim(theta=vstart, f=fnobj, ui=Mat, ci=Bnd, grad=gradient,
                      fix.arg=fix.arg, obs=data, ddistnam=ddistname, hessian=!is.null(gradient), method=meth,
                      ...), silent=TRUE)
              }
              else #cens == TRUE
                opttryerror <- try(opt <- constrOptim(theta=vstart, f=fnobjcens, ui=Mat, ci=Bnd, grad=gradient,
                      ddistnam=ddistname, rcens=rcens, lcens=lcens, icens=icens, ncens=ncens, pdistnam=pdistname,
                      fix.arg=fix.arg, hessian=!is.null(gradient), method=meth, ...), silent=TRUE)
              if(!inherits(opttryerror, "try-error"))
                if(length(opt$counts) == 1) #appears when the initial point is a solution
                  opt$counts <- c(opt$counts, NA)


            }else #opt.fun == "optim"
            {
              if(!cens)
                opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                                ddistnam=ddistname, hessian=TRUE, method=meth, lower=lower, upper=upper,
                                                ...), silent=TRUE)
              else #cens == TRUE
                opttryerror <- try(opt <- optim(par=vstart, fn=fnobjcens, fix.arg=fix.arg, gr=gradient,
                                                rcens=rcens, lcens=lcens, icens=icens, ncens=ncens, ddistnam=ddistname,
                                                pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, upper=upper,
                                                ...), silent=TRUE)
            }

          }else #hasbound == FALSE
          {
            opt.fun <- "optim"
            if(!cens)
              opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                              ddistnam=ddistname, hessian=TRUE, method=meth, lower=lower, upper=upper,
                                              ...), silent=TRUE)
            else #cens == TRUE
              opttryerror <- try(opt <- optim(par=vstart, fn=fnobjcens, fix.arg=fix.arg, gr=gradient,
                                              rcens=rcens, lcens=lcens, icens=icens, ncens=ncens, ddistnam=ddistname,
                                              pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, upper=upper,
                                              ...), silent=TRUE)
          }
          options(warn=owarn)

          if (inherits(opttryerror, "try-error"))
          {
              warnings("The function optim encountered an error and stopped.")
              if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))
              return(list(estimate = rep(NA, length(vstart)), convergence = 100, loglik = NA,
                          hessian = NA, optim.function=opt.fun, fix.arg = fix.arg,
                          optim.method=meth, fix.arg.fun = fix.arg.fun, counts=c(NA, NA)))
          }

          if (opt$convergence>0) {
              warnings("The function optim failed to converge, with the error code ",
                       opt$convergence)
          }
          if(is.null(names(opt$par)))
            names(opt$par) <- names(vstart)
          res <- list(estimate = opt$par, convergence = opt$convergence, value=opt$value,
                      hessian = opt$hessian, optim.function=opt.fun, optim.method=meth,
                      fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights = weights,
                      counts=opt$counts, optim.message=opt$message, loglik = -opt$value)
      }
      else # Try to minimize the minus (log-)likelihood using a user-supplied optim function
      {
          options(warn=ifelse(silent, -1, 0))
          if (!cens)
              opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data,
                  ddistnam=ddistname, par=vstart, ...), silent=TRUE)
          else
              opttryerror <-try(opt<-custom.optim(fn=fnobjcens, fix.arg=fix.arg, rcens=rcens,
                  lcens=lcens, icens=icens, ncens=ncens, ddistnam=ddistname, pdistnam=pdistname,
                  par=vstart, ...), silent=TRUE)
          options(warn=owarn)

          if (inherits(opttryerror, "try-error"))
          {
              warnings("The customized optimization function encountered an error and stopped.")
              if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))
              return(list(estimate = rep(NA, length(vstart)), convergence = 100, loglik = NA,
                          hessian = NA, optim.function=custom.optim, fix.arg = fix.arg,
                          fix.arg.fun = fix.arg.fun, counts=c(NA, NA)))
          }

          if (opt$convergence>0) {
              warnings("The customized optimization function failed to converge, with the error code ",
                       opt$convergence)
          }
          if(is.null(names(opt$par)))
            names(opt$par) <- names(vstart)
          argdot <- list(...)
          method.cust <- argdot$method
          res <- list(estimate = opt$par, convergence = opt$convergence, value=opt$value,
                      hessian = opt$hessian, optim.function = custom.optim, optim.method = method.cust,
                      fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights = weights,
                      counts=opt$counts, optim.message=opt$message, loglik = -opt$value)
      }

      return(res)
  }

  ## old function with previous name for censored data
  mledistcens <- function(censdata, distr, start=NULL, optim.method="default", lower=-Inf, upper=Inf)
  {
      stop("The function \"mledistcens\" is no more used. Now the same function \"mledist\" must be used both for censored and non censored data.")
  }


  # testdpqfun function returns a vector of TRUE when the d,p,q functions exist
  # and behave like in R (e.g. d,p,qexp()), otherwise a list of messages.

  # INPUTS
  #distr : the distribution name
  #fun: a character vector with letters among d, p, q
  #start.arg: an initial value list
  #fix.arg: a fixed value list
  #discrete: a logical whether the distribution is discrete

  # OUTPUTS
  # a vector of logical TRUE or a vector of text messages
  testdpqfun <- function(distr, fun=c("d","p","q"), start.arg,
                         fix.arg=NULL, discrete=FALSE)
  {
    stopifnot(all(is.character(fun)))
    fun <- fun[fun %in% c("d","p","q")]
    stopifnot(length(fun) > 0)

    if(is.vector(start.arg))
      start.arg <- as.list(start.arg)
    if(is.function(fix.arg))
      stop("fix.arg should be either a named list or NULL but not a function")

    op <- options() #get current options
    #print(getOption("warn"))
    options(warn=-1)

    res <- NULL
    if("d" %in% fun)
      res <- rbind(res, test1fun(paste0("d", distr), start.arg, fix.arg))
    if("p" %in% fun)
      res <- rbind(res, test1fun(paste0("p", distr), start.arg, fix.arg))
    if("q" %in% fun)
      res <- rbind(res, test1fun(paste0("q", distr), start.arg, fix.arg))

    options(op)     # reset (all) initial options
    res
  }

  test1fun <- function(fn, start.arg, fix.arg, dpqr)
  {
    res <- data.frame(ok=FALSE, txt="")
    stopifnot(is.list(start.arg))
    if(!is.null(fix.arg))
      stopifnot(is.list(fix.arg))

    #does the function exist?
    if(!exists(fn, mode="function"))
    {
      res$txt <- paste("The", fn, "function must be defined")
      return(res)
    }

    #naming convention
    if(missing(dpqr))
      dpqr <- substr(fn, 1, 1)
    firstarg_theo <- switch(dpqr, "d"="x", "p"="q", "q"="p", "r"="n")
    firstarg_found <- names(formals(fn))[1]
    if(firstarg_found != firstarg_theo)
    {
      t0 <- paste("The", fn, "function should have its first argument named:", firstarg_theo)
      res$txt <- paste(t0, "as in base R")
      return(res)
    }

    #zero-component vector
    res0 <- try(do.call(fn, c(list(numeric(0)), start.arg, fix.arg)), silent=TRUE)
    t0 <- paste("The", fn, "function should return a zero-length vector when input has length zero and not raise an error")
    t1 <- paste("The", fn, "function should return a zero-length vector when input has length zero")
    if(class(res0) == "try-error")
    {
      res$txt <- t0
      return(res)
    }
    if(length(res0) != 0)
    {
      res$txt <- t1
      return(res)
    }

    #inconsistent value
    x <- c(0, 1, Inf, NaN, -1)
    res1 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
    t2 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent values and not raise an error")
    if(class(res1) == "try-error")
    {
      res$txt <- t2
      return(res)
    }

    #missing value
    x <- c(0, 1, NA)
    res2 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
    t4 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not raise an error")
    t5 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not remove missing values")
    if(class(res2) == "try-error")
    {
      res$txt <- t4
      return(res)
    }
    if(length(res2) != length(x))
    {
      res$txt <- t5
      return(res)
    }

    #inconsistent parameter
    x <- 0:1
    start.arg <- lapply(start.arg, function(x) -x)
    res3 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
    t6 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent parameters and not raise an error")
    if(class(res3) == "try-error")
    {
      res$txt <- t6
      return(res)
    }

    #wrong parameter name
    x <- 0:1
    names(start.arg) <- paste0(names(start.arg), "_")
    res4 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
    t8 <- paste("The", fn, "function should raise an error when names are incorrectly named")
    if(class(res4) != "try-error")
    {
      res$txt <- t8
      return(res)
    }
    return(data.frame(ok=TRUE, txt=""))
  }


  # INPUTS
  # argdistname : argument names of the distribution from names(formals())

  # OUTPUTS
  # parameter names (as a vector) of the distribution (excluding non parameter argument)

  computegetparam <- function(argdistname)
  {
    #remove first argument, that should be "x", "p", "q", or "n", see ?dgamma, pgamma, qgamma
    argdistname <- argdistname[-1]
    nonparaminR <- c("x", "p", "q", "n") #defensive programming
    #remove other arguments, see ?dgamma, pgamma, qgamma, dbeta
    nonparaminR <- c(nonparaminR, "log", "log.p", "lower.tail", "ncp")
    nonparaminActuar <- c("limit", "order", "t")
    nonparaminGamlssdist <- "fast"
    nonparamspecial <- c("...", "..1", "..2")
    #see ?dnig, dhyperb, dskewlap, dgig,...
    nonparaminGenHyperbolic <- c("param", "KOmega", "ibfTol", "nmax", "method", "intTol",
                                 "valueOnly", "nInterpol", "uniTol", "subdivisions", "logPars")
    #see ?dsn
    nonparamsn <- "dp"

    plist <- setdiff(argdistname, nonparaminR)
    plist <- setdiff(plist, nonparaminActuar)
    plist <- setdiff(plist, nonparaminGamlssdist)
    plist <- setdiff(plist, nonparamspecial)
    plist <- setdiff(plist, nonparaminGenHyperbolic)
    plist <- setdiff(plist, nonparamsn)

    plist
  }


  # checkparam function checks start.arg and fix.arg that parameters are named correctly

  # INPUTS
  # start.arg : a named list
  # fix.arg : NULL or a named list
  # argdistname : argument names of the distribution
  # hasnodefaultval : vector of logical indicating no default value of argument

  # OUTPUTS
  # a named list with components: ok (TRUE or FALSE), txt (NULL or the error message),
  # start.arg : a named list of starting values for optimization
  # or a function to compute them from data
  checkparamlist <- function(start.arg, fix.arg, argdistname, hasnodefaultval)
  {
    errtxt <- list(t3="'start' must specify names which are arguments to 'distr'.",
            t4="'fix.arg' must specify names which are arguments to 'distr'.",
            t5="A distribution parameter cannot be specified both in 'start' and 'fix.arg'.",
            t6="'start' should not have NA or NaN values.",
            t7="'fix.arg' should not have NA or NaN values.",
            t8="Some parameter names have no starting/fixed value: ",
            t9="Some parameter names have no starting/fixed value but have a default value: ")

    vstart <- unlist(start.arg)
    #check unexpected names
    m <- match(names(vstart), argdistname)
    if (any(is.na(m)))
      stop(errtxt$t3)
    #check NA/NaN values
    if(any(is.na(vstart)) || any(is.nan(vstart)))
      stop(errtxt$t6)
    if(!is.null(fix.arg))
    {
      vfix <- unlist(fix.arg)
      #check unexpected names
      mfix <- match(names(vfix), argdistname)
      if (any(is.na(mfix)))
        stop(errtxt$t4)

      # check that some parameters are not both in fix.arg and start
      minter <- match(names(vstart), names(vfix))
      if (any(!is.na(minter)))
        stop(errtxt$t5)

      #check NA/NaN values
      if(any(is.na(vfix)) || any(is.nan(vfix)))
        stop(errtxt$t7)
      allparname <- names(c(vstart, vfix))
    }else
      allparname <- names(vstart)

    theoparam <- computegetparam(argdistname)
    #special case where both scale and rate are allowed, see ?dgamma
    if("scale" %in% theoparam && "rate" %in% theoparam)
    {
      errt8 <- any(!allparname %in% theoparam) || length(allparname) != length(theoparam)-1
    #special case where both prob and mu are allowed, see ?dnbinom
    }else if(length(theoparam) == 3 && all(c("size", "prob", "mu") %in% theoparam))
    {
      errt8 <- any(!allparname %in% theoparam) || length(allparname) != length(theoparam)-1
    }else
      errt8 <- any(!theoparam %in% allparname)
    #only make a warning if unset arguments have a default value
    if(errt8)
    {
      unsetarg <- theoparam[!theoparam %in% allparname]
      if(any(hasnodefaultval[unsetarg]))
        stop(paste0(errtxt$t8, unsetarg, "."))
      else
        warning(paste0(errtxt$t9, unsetarg, "."))
    }

    list("start.arg"=start.arg, "fix.arg"=fix.arg)
  }

  # start.arg.default function returns initial values of parameters generally using moments or quantiles

  # INPUTS
  #x : data vector or matrix
  #distr : the distribution name

  # OUTPUTS
  # a named list or raises an error
  start.arg.default <- function(x, distr)
  {
    if (distr == "norm") {
      n <- length(x)
      sd0 <- sqrt((n - 1)/n) * sd(x)
      mx <- mean(x)
      start <- list(mean=mx, sd=sd0)
    }else if (distr == "lnorm") {
      if (any(x <= 0))
        stop("values must be positive to fit a lognormal distribution")
      n <- length(x)
      lx <- log(x)
      sd0 <- sqrt((n - 1)/n) * sd(lx)
      ml <- mean(lx)
      start <- list(meanlog=ml, sdlog=sd0)
    }else if (distr == "pois") {
      start <- list(lambda=mean(x))
    }else if (distr == "exp") {
      if (any(x < 0))
        stop("values must be positive to fit an exponential  distribution")
      start <- list(rate=1/mean(x))
    }else if (distr == "gamma") {
      if (any(x < 0))
        stop("values must be positive to fit an gamma  distribution")
      n <- length(x)
      m <- mean(x)
      v <- (n - 1)/n*var(x)
      start <- list(shape=m^2/v, rate=m/v)
    }else if (distr == "nbinom") {
      n <- length(x)
      m <- mean(x)
      v <- (n - 1)/n*var(x)
      size <- ifelse(v > m, m^2/(v - m), 100)
      start <- list(size = size, mu = m)
    }else if (distr == "geom" ) {
      m <- mean(x)
      prob <- ifelse(m>0, 1/(1+m), 1)
      start <- list(prob=prob)
    }else if (distr == "beta") {
      if (any(x < 0) | any(x > 1))
        stop("values must be in [0-1] to fit a beta distribution")
      n <- length(x)
      m <- mean(x)
      v <- (n - 1)/n*var(x)
      aux <- m*(1-m)/v - 1
      start <- list(shape1=m*aux, shape2=(1-m)*aux)
    }else if (distr == "weibull") {
      if (any(x < 0))
        stop("values must be positive to fit an Weibull  distribution")
      m <- mean(log(x))
      v <- var(log(x))
      shape <- 1.2/sqrt(v)
      scale <- exp(m + 0.572/shape)
      start <- list(shape = shape, scale = scale)
    }else if (distr == "logis") {
      n <- length(x)
      m <- mean(x)
      v <- (n - 1)/n*var(x)
      start <- list(location=m, scale=sqrt(3*v)/pi)
    }else if (distr == "cauchy") {
      start <- list(location=median(x), scale=IQR(x)/2)
    }else if (distr == "unif"){
      start <- list(min=0, max=1)
    }else if (distr == "invgamma")
    {
      if (any(x < 0))
        stop("values must be positive to fit an inverse gamma  distribution")
      #http://en.wikipedia.org/wiki/Inverse-gamma_distribution
      m1 <- mean(x)
      m2 <- mean(x^2)
      shape <- (2*m2-m1^2)/(m2-m1^2)
      scale <- m1*m2/(m2-m1^2)
      start <- list(shape=shape, scale=scale)
    }else if (distr == "llogis")
    {
      if (any(x < 0))
        stop("values must be positive to fit a log-logistic  distribution")
      p25 <- as.numeric(quantile(x, 0.25))
      p75 <- as.numeric(quantile(x, 0.75))
      shape <- 2*log(3)/(log(p75)-log(p25))
      scale <- exp(log(p75)+log(p25))/2
      start <- list(shape=shape, scale=scale)
    }else if (distr == "invweibull")
    {
      if (any(x < 0))
        stop("values must be positive to fit an inverse Weibull distribution")
      g <- log(log(4))/(log(log(4/3)))
      p25 <- as.numeric(quantile(x, 0.25))
      p75 <- as.numeric(quantile(x, 0.75))
      shape <- exp((g*log(p75)-log(p25))/(g-1))
      scale <-log(log(4))/(log(shape)-log(p25))
      start <- list(shape=shape, scale=max(scale, 1e-9))
    }else if (distr == "pareto1")
    {
      if (any(x < 0))
        stop("values must be positive to fit a Pareto distribution")
      #http://www.math.umt.edu/gideon/pareto.pdf
      x1 <- min(x)
      m1 <- mean(x)
      n <- length(x)
      shape <- (n*m1-x1)/(n*(m1-x1))
      min <- x1*(n*shape - 1)/(n*shape)
      start <- list(shape=shape, min=min)
    }else if (distr == "pareto")
    {
      if (any(x < 0))
        stop("values must be positive to fit a Pareto distribution")
      m1 <- mean(x)
      m2 <- mean(x^2)
      scale <- (m1*m2)/(m2-2*m1^2)
      shape <- 2*(m2-m1^2)/(m2-2*m1^2)
      start <- list(shape=shape, scale=scale)
    }else if (distr == "lgamma")
    {
      if (any(x < 0))
        stop("values must be positive to fit a log-gamma distribution")
      #p228 of Klugmann and Hogg (1984)
      m1 <- mean(log(x))
      m2 <- mean(log(x)^2)
      alpha <- m1^2/(m2-m1^2)
      lambda <- m1/(m2-m1^2)
      start <- list(shapelog=alpha, ratelog=lambda)
    }else if (distr == "trgamma") {
      if (any(x < 0))
        stop("values must be positive to fit an trans-gamma  distribution")
      #same as gamma with shape2=tau=1
      n <- length(x)
      m <- mean(x)
      v <- (n - 1)/n*var(x)
      start <- list(shape1=m^2/v, shape2=1, rate=m/v)
    }else if (distr == "invtrgamma") {
      if (any(x < 0))
        stop("values must be positive to fit an inverse trans-gamma  distribution")
      #same as gamma with shape2=tau=1
      n <- length(1/x)
      m <- mean(1/x)
      v <- (n - 1)/n*var(1/x)
      start <- list(shape1=m^2/v, shape2=1, rate=m/v)
    }else
      stop(paste0("Unknown starting values for distribution ", distr, "."))

    return(start)
  }


  # checkparam function checks start.arg and fix.arg that parameters are named correctly

  # INPUTS
  # start.arg : starting values for optimization or the function to compute them from data or NULL
  # fix.arg : fixed values of paramaters or the function to compute them from data or NULL
  # obs : the full dataset
  # distname : name of the distribution


  # OUTPUTS
  # two named list with untested components
  manageparam <- function(start.arg, fix.arg, obs, distname)
  {
    #if clause with 3 different cases:
    #start.arg : NULL | named list | a function

    if(is.null(start.arg))
    {
      trystart <- try(start.arg.default(obs, distname), silent = TRUE)
      if(class(trystart) == "try-error")
      {
        cat("Error in computing default starting values.\n")
        stop(trystart)
      }
      lstart <- trystart
      #lstart should be a named list but check it
      hasnoname <- is.null(names(lstart)) || !is.list(lstart)
      if(hasnoname)
        stop("Starting values must be a named list, error in default starting value.")

    }else if(is.list(start.arg))
    {
      hasnoname <- is.null(names(start.arg))
      if(hasnoname)
        stop("Starting values must be a named list (or a function returning a named list).")
      lstart <- start.arg
    }else if(is.function(start.arg))
    {
      trystart <- try(start.arg(obs), silent = TRUE)
      if(class(trystart) == "try-error")
      {
        cat("Error in computing starting values with your function.\n")
        stop(trystart)
      }
      lstart <- trystart
      hasnoname <- is.null(names(lstart)) || !is.list(lstart)
      if(hasnoname)
        stop("Starting values must be a named list, your function does not return that.")
    }else
      stop("Wrong type of argument for start")

    #if clause with 3 different cases:
    #fix.arg : NULL | named list | a function
    if(is.null(fix.arg))
    {
      lfix <- NULL
    }else if(is.list(fix.arg))
    {
      hasnoname <- is.null(names(fix.arg))
      if(hasnoname)
        stop("Fixed parameter values must be a named list (or a function returning a named list).")
      lfix <- fix.arg
    }else if(is.function(fix.arg))
    {
      tryfix <- try(fix.arg(obs), silent = TRUE)
      if(class(tryfix) == "try-error")
      {
        cat("Error in computing fixed parameter values with your function.\n")
        stop(tryfix)
      }
      lfix <- tryfix
      hasnoname <- is.null(names(lfix)) || !is.list(lfix)
      if(hasnoname)
        stop("Fixed parameter values must be a named list, your function does not return that.")
    }else
      stop("Wrong type of argument for fix.arg")

    #eliminate arguments both in lstart and lfix (when start.arg was NULL)
    if(is.null(start.arg) && !is.null(lfix))
    {
      lstart <- lstart[!names(lstart) %in% names(lfix)]
      if(length(lstart) == 0)
        stop("Don't need to use fitdist() if all parameters have fixed values")
    }

    list("start.arg"=lstart, "fix.arg"=lfix)
  }
##HelprFuncEND
##HelprFuncEND
##HelprFuncEND

    #check argument distr
    if (!is.character(distr))
        distname <- substring(as.character(match.call()$distr), 2)
    else
        distname <- distr
    ddistname <- paste("d", distname, sep="")
    if (!exists(ddistname, mode="function"))
        stop(paste("The ", ddistname, " function must be defined"))

    #pdistname <- paste("p", distname, sep="")
    #if (!exists(pdistname, mode="function"))
    #    stop(paste("The ", pdistname, " function must be defined"))
    #check argument discrete
    if(missing(discrete))
    {
      if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois")))
        discrete <- TRUE
      else
        discrete <- FALSE
    }
    if(!is.logical(discrete))
      stop("wrong argument 'discrete'.")
    if(!is.logical(keepdata) || !is.numeric(keepdata.nb) || keepdata.nb < 2)
      stop("wrong arguments 'keepdata' and 'keepdata.nb'")
    #check argument method
    if(any(method == "mom"))
        warning("the name \"mom\" for matching moments is NO MORE used and is replaced by \"mme\"")

    method <- match.arg(method, c("mle", "mme", "qme", "mge", "mse"))
    if(method %in% c("mle", "mme", "mge", "mse"))
      dpq2test <- c("d", "p")
    else
      dpq2test <- c("d", "p", "q")
    #check argument data
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
        stop("data must be a numeric vector of length greater than 1")

    #encapsulate three dots arguments
    my3dots <- list(...)
    if (length(my3dots) == 0)
      my3dots <- NULL
    n <- length(data)

    # manage starting/fixed values: may raise errors or return two named list
    arg_startfix <- manageparam(start.arg=start, fix.arg=fix.arg, obs=data,
                                 distname=distname)

    #check inconsistent parameters
    argddistname <- names(formals(ddistname))
    hasnodefaultval <- sapply(formals(ddistname), is.name)
    arg_startfix <- checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg,
                                   argddistname, hasnodefaultval)
    #arg_startfix contains two names list (no longer NULL nor function)
    #store fix.arg.fun if supplied by the user
    if(is.function(fix.arg))
      fix.arg.fun <- fix.arg
    else
      fix.arg.fun <- NULL

    # check d, p, q, functions of distname
    resdpq <- testdpqfun(distname, dpq2test, start.arg=arg_startfix$start.arg,
               fix.arg=arg_startfix$fix.arg, discrete=discrete)
    if(any(!resdpq$ok))
    {
      for(x in resdpq[!resdpq$ok, "txt"])
        warning(x)
    }


    # Fit with mledist, qmedist, mgedist or mmedist
    if (method == "mme")
    {
        if (!is.element(distname, c("norm", "lnorm", "pois", "exp", "gamma",
                "nbinom", "geom", "beta", "unif", "logis")))
            if (!"order" %in% names(my3dots))
                stop("moment matching estimation needs an 'order' argument")

        mme <- mmedist(data, distname, start=arg_startfix$start.arg,
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)

        sd <- NULL
        correl <- varcovar <- NULL

        estimate <- mme$estimate
        loglik <- mme$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        convergence <- mme$convergence
        fix.arg <- mme$fix.arg
        weights <- mme$weights
    }else if (method == "mle")
    {
        mle <- mledist(data, distname, start=arg_startfix$start.arg,
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)
        if (mle$convergence>0)
           stop("the function mle failed to estimate the parameters,
                with the error code ", mle$convergence, "\n")
        estimate <- mle$estimate
        if(!is.null(mle$hessian)){
            #check for NA values and invertible Hessian
            if(all(!is.na(mle$hessian)) && qr(mle$hessian)$rank == NCOL(mle$hessian)){
                varcovar <- solve(mle$hessian)
                sd <- sqrt(diag(varcovar))
                correl <- cov2cor(varcovar)
            }else{
                varcovar <- NA
                sd <- NA
                correl <- NA
            }
        }else{
            varcovar <- NA
            sd <- NA
            correl <- NA
        }
        loglik <- mle$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        convergence <- mle$convergence
        fix.arg <- mle$fix.arg
        weights <- mle$weights
    }else if (method == "qme")
    {
        if (!"probs" %in% names(my3dots))
            stop("quantile matching estimation needs an 'probs' argument")

        qme <- qmedist(data, distname, start=arg_startfix$start.arg,
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)

        estimate <- qme$estimate
        sd <- NULL
        loglik <- qme$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        correl <- varcovar <- NULL

        convergence <- qme$convergence
        fix.arg <- qme$fix.arg
        weights <- qme$weights
    }else if (method == "mge")
    {
        if (!"gof" %in% names(my3dots))
            warning("maximum GOF estimation has a default 'gof' argument set to 'CvM'")

        mge <- mgedist(data, distname, start=arg_startfix$start.arg,
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)

        estimate <- mge$estimate
        sd <- NULL
        loglik <- mge$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        correl <- varcovar <- NULL

        convergence <- mge$convergence
        fix.arg <- mge$fix.arg
        weights <- NULL
    }else if (method == "mse")
    {
      mse <- msedist(data, distname, start=arg_startfix$start.arg,
                     fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)

      estimate <- mse$estimate
      sd <- NULL
      loglik <- mse$loglik
      npar <- length(estimate)
      aic <- -2*loglik+2*npar
      bic <- -2*loglik+log(n)*npar
      correl <- varcovar <- NULL

      convergence <- mse$convergence
      fix.arg <- mse$fix.arg
      weights <- mse$weights
    }else
    {
        stop("match.arg() did not work correctly")
    }

    #needed for bootstrap
    if (!is.null(fix.arg))
      fix.arg <- as.list(fix.arg)

    if(keepdata)
    {
      reslist <- list(estimate = estimate, method = method, sd = sd, cor = correl,
                  vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=data,
                  distname = distname, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun,
                  dots = my3dots, convergence = convergence, discrete = discrete,
                  weights = weights)
    }else #just keep a sample set of all observations
    {
      n2keep <- min(keepdata.nb, n)-2
      imin <- which.min(data)
      imax <- which.max(data)
      subdata <- data[sample((1:n)[-c(imin, imax)], size=n2keep, replace=FALSE)]
      subdata <- c(subdata, data[c(imin, imax)])

      reslist <- list(estimate = estimate, method = method, sd = sd, cor = correl,
                  vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=subdata,
                  distname = distname, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun,
                  dots = my3dots, convergence = convergence, discrete = discrete,
                  weights = weights)
    }


    return(structure(reslist, class = "fitdist"))

}
