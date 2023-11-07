## =====================
   version <- function()
## =====================
{
   print('This is Version 1.0 of estimateN.r (1 Aug 2007)')
}

## ============================
   ranstart <- function(nparam)
## ============================
{
        ##-----------------------------------------------------
        ## Set random starting values. nparam is the number of 
        ## parameters. Random starting values are logit of 
        ## uniform(0,1) variables.
        ##-----------------------------------------------------
    u = runif(nparam)
    log(u/(1-u))
}


## ======================================================================
   binbbinprob <- function( kmax, pbin, pbbin, theta, alphabbin)
## ======================================================================
{
        ##----------------------------------------------------------
        ## Calculates binomial + beta-binomial mixture probabilities
        ##  pbin = binomial p
        ## pbbin = beta-binomial p
        ## theta = beta-binomial theta
        ## alphabbin = mixture proportion for beta-bin
        ##----------------------------------------------------------
    prob = numeric(kmax)

        ## ----------------------------------
        ## Begin by finding the logarithm of
        ## the probability of getting a zero
        ## ----------------------------------
    r <- c(0:(kmax-1))
    logpz <- sum(log(1-pbbin+r*theta) - log(1+r*theta))
    prob[1] <- alphabbin * exp(logpz) + (1-alphabbin)*(1-pbin)^kmax
    cumprb <- prob[1]

        ## ----------------------------
        ## Find the logarithm of other 
        ## probabilities from 1 to kmax
        ## ----------------------------
    eps <- 1.0E-20
    
    r1 = c(1:kmax)
    r2 = r1 - 1
    term1 = cumsum(log(pbbin+r2*theta))
    term2 = cumsum(log(1-pbbin+r2*theta))
    term3 = sum(log(1+r2*theta))
    logcomb <- lgamma(kmax+1) - lgamma(r1+1) - lgamma(kmax-r1+1)
    logprb = term1 + c(rev(term2[1:(kmax-1)]), 0) - term3
    logbin <- r1*log(pbin+eps) + (kmax-r1)*log(1-pbin+eps)
    prob[2:(kmax+1)] <- alphabbin * exp(logcomb+logprb) + (1-alphabbin)*exp(logcomb+logbin)
    prob
}


## ======================================
   unpack <- function(param, modelnumber)
## ======================================
{
    small <- 1.0E-20

        ##--------------------------
        ## Binomial model 
        ##    param[2] = logit(pbin)
        ##--------------------------
    if (modelnumber == 1) {
       pbin = 1 / (1 + exp(-param[2]))
       pbbin = small 
       theta = small
       alphabbin = small
    }

        ##--------------------------
        ## Beta-binomial model 
        ##    param[2] = logit(pbin)
        ##    param[3] = log(theta)
        ##--------------------------
    if (modelnumber == 2) {
       pbbin = 1 / (1 + exp(-param[2]))
       theta = exp(param[3])
       pbin = small 
       alphabbin = 1 - small
    }

        ##-----------------------------------------------------
        ## Two-binomial model 
        ##    param[2] = logit(p1)
        ##    param[3] = logit(p2)
        ##    param[4] = logit(alpha1) = proportion component 1
        ## NB Force p1 <= p2, to avoid redundancy
        ##-----------------------------------------------------
    if (modelnumber == 3) {
#       if (param[2] > param[3]) {
#          tmp = param[3]
#          param[3] = param[2]
#          param[2] = tmp
#       }
       pbbin = 1 / (1 + exp(-param[2]))
       pbin = 1 / (1 + exp(-param[3]))
       alphabbin = 1 / (1 + exp(-param[4]))
       theta = small
    }

        ##-----------------------------------------------------
        ## Binomial beta-binomial model 
        ##    param[2] = logit(pbin)
        ##    param[3] = logit(pbbin)
        ##    param[4] = log(theta)
        ##    param[5] = logit(alphabbin)
        ##-----------------------------------------------------
    if (modelnumber == 4) {
       pbin = 1 / (1 + exp(-param[2]))
       pbbin = 1 / (1 + exp(-param[3]))
       theta = exp(param[4])
       alphabbin = 1 / (1 + exp(-param[5]))
    }

    c(pbin, pbbin, theta, alphabbin)
}

## ===========================================================
   loglik <- function(param, modelnumber, freq, kobsmax, kmax)
## ===========================================================
{
        ##--------------------------------------------------------
        ## Calculates the log-likelihood for the appropriate model 
        ##--------------------------------------------------------
    n0 = exp(param[1])    # Size of zero class
    sumf = sum(freq)      # Total observed frequency
    ntot = sumf + n0      # Current estimate of popn size

        ##-----------------------------
        ## Set parameters appropriately
        ##-----------------------------
    parvals = unpack(param, modelnumber)
    pbin = parvals[1]
    pbbin = parvals[2]
    theta = parvals[3]
    alphabbin = parvals[4]

        ##-------------------------------------------
        ## Calculate probabilities and log-likelihood
        ##-------------------------------------------
    prob = binbbinprob(kmax, pbin, pbbin, theta, alphabbin)
    logl = -lgamma(ntot+1)+lgamma(n0+1)-n0*log(prob[1]) +
           sum( lgamma(freq[1:kobsmax]+1) - freq[1:kobsmax] * log(prob[2:(kobsmax+1)]) )
    logl=Re(logl)
}


## ===================================================================
   profloglik <- function(param, modelnumber, freq, kobsmax, kmax, n0)
## ===================================================================
{
        ##----------------------------------------------------------------
        ## Calculates the profile log-likelihood for the appropriate model 
        ##----------------------------------------------------------------
    sumf = sum(freq)      # Total observed frequency
    ntot = sumf + n0      # Current estimate of popn size

        ##-----------------------------
        ## Set parameters appropriately
        ##-----------------------------
    parvals = unpack(c(0,param), modelnumber)
    pbin = parvals[1]
    pbbin = parvals[2]
    theta = parvals[3]
    alphabbin = parvals[4]

        ##-------------------------------------------
        ## Calculate probabilities and log-likelihood
        ##-------------------------------------------
    prob = binbbinprob(kmax, pbin, pbbin, theta, alphabbin)
    logl = -lgamma(ntot+1)+lgamma(n0+1)-n0*log(prob[1]) +
           sum( lgamma(freq[1:kobsmax]+1) - freq[1:kobsmax] * log(prob[2:(kobsmax+1)]) )
    logl=Re(logl)
}


## ==================================================================
   ploglik <- function(logn0val, loglmin, starter, modelnumber, freq, 
                       kobsmax, kmax, confpercent)
## ==================================================================
{
        ## Function used for calculating profile log-likelihood
        ## intervals. Calculates the difference between the 
        ## profile-log-likelihood at logn0val and at the overall
        ## mle, and subtracts the appropriate point of the 
        ## chi-squared distribution. So this function takes the
        ## value zero at the confidence limits.

    nmodparam = c(2,3,4,5)     # number of parameters for each model
    nstartprof = 10
    n0 = exp(logn0val)
    minlogl = 9.99E29 
    for (j in 1: nstartprof) {
        if (modelnumber==1) options(warn=-1)
        rndmod <- optim(starter, profloglik,, modelnumber, freq, kobsmax, kmax, n0)
        options(warn=0)
        if (rndmod$value < minlogl) {
           minlogl <- rndmod$value
           profmodel <- rndmod
        }
#              cat(j, rndmod$value, rndmod$convergence, "\n")
        starter = ranstart(nmodparam[modelnumber]-1)
    }
    profmodel$val - loglmin - qchisq(confpercent/100, 1) / 2
}



## ============================================================================
   printest <- function(loglmin, param, modelnumber, freq, kobsmax, kmax, vcov, 
                        noses, lower, upper, confpercent)
## ============================================================================
{
        ##--------------------------
        ## Produces formatted output
        ##--------------------------
    n0 = exp(param[1])    # Size of zero class
    sumf = sum(freq)      # Total observed frequency
    ntot = sumf + n0      # Current estimate of popn size

        ##-----------------------------
        ## Set parameters appropriately
        ##-----------------------------
    parvals = unpack(param, modelnumber)
    pbin = parvals[1]
    pbbin = parvals[2]
    theta = parvals[3]
    alphabbin = parvals[4]

    logitpbin = log(pbin/(1-pbin))
    logitpbbin = log(pbbin/(1-pbbin))
    logitalphabbin = log(alphabbin/(1-alphabbin))
    logn0 = param[1]
    logtheta = log(theta) 

    pr.pbin = format(round(pbin, digits=4))
    pr.pbbin = format(round(pbbin, digits=4))
    pr.alphabbin = format(round(alphabbin, digits=4))
    pr.n0 = format(round(n0, digits=1))
    pr.ntot = format(round(ntot, digits=1))
    pr.theta = format(round(theta, digits=3))

    pr.logitpbin = format(round(logitpbin, digits=3))
    pr.logitpbbin = format(round(logitpbbin, digits=3))
    pr.logitalphabbin = format(round(logitalphabbin, digits=3))
    pr.logn0 = format(round(logn0, digits=2))
    pr.logtheta = format(round(logtheta, digits=3))

  if (noses) {
    if (modelnumber == 1) {
       cat("\n*** BINOMIAL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
          cat("\nWARNING Hessian matrix is ill-conditioned or not +ve definite\n")
          cat("No SEs available\n")
       cat("\nParameter estimates\n")
       cat("Binomial p =",pr.pbin,"(logit =",pr.logitpbin,")\n")
    }

    if (modelnumber == 2) {
       cat("\n*** BETA-BINOMIAL MODEL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
          cat("\nWARNING Hessian matrix is ill-conditioned or not +ve definite\n")
          cat("No SEs available\n")
       cat("\nParameter estimates\n")
       cat("Beta-binomial p =",pr.pbbin,"(logit =",pr.logitpbbin,")\n")
       cat("Beta-binomial theta=",pr.theta,"(log =",pr.logtheta,")\n")
    }

    if (modelnumber == 3) {
       cat("\n*** TWO-BINOMIAL MODEL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
          cat("\nWARNING Hessian matrix is ill-conditioned or not +ve definite\n")
          cat("No SEs available\n")
       cat("\nParameter estimates\n")
       cat("Binomial p1 =",pr.pbin,"(logit =",pr.logitpbin,")\n")
       cat("Binomial p2 =",pr.pbbin,"(logit =",pr.logitpbbin,")\n")
       cat("Proportion p1 =",pr.alphabbin,"(logit =",pr.logitalphabbin,")\n")
    }

    if (modelnumber == 4) {
       cat("\n*** BINOMIAL/BETA-BINOMIAL MODEL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
          cat("\nWARNING Hessian matrix is ill-conditioned or not +ve definite\n")
          cat("No SEs available\n")
       cat("\nParameter estimates\n")
       cat("Binomial p =",pr.pbin,"(logit =",pr.logitpbin,")\n")
       cat("Beta-binomial p =",pr.pbbin,"(logit =",pr.logitpbbin,")\n")
       cat("Beta-binomial theta=",pr.theta,"(log =",pr.logtheta,")\n")
       cat("Proportion beta-bin =",pr.alphabbin,"(logit =",pr.logitalphabbin,")\n")
    }
    cat("Number of unseen individuals =", pr.n0,"(log =",pr.logn0,")\n")
} else {
    if (modelnumber == 1) {
       cat("\n*** BINOMIAL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
       cat("\nParameter estimates\n")
       cat("Binomial p =",pr.pbin,"(logit =",pr.logitpbin,", SE = ",
            format(round(sqrt(vcov[2,2]),digits=4)),")\n")
    }

    if (modelnumber == 2) {
       cat("\n*** BETA-BINOMIAL MODEL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
       cat("\nParameter estimates\n")
       cat("Beta-binomial p =",pr.pbbin,"(logit =",pr.logitpbbin,", SE = ",
            format(round(sqrt(vcov[2,2]),digits=4)),")\n")
       cat("Beta-binomial theta=",pr.theta,"(log =",pr.logtheta,", SE = ",
            format(round(sqrt(vcov[3,3]),digits=4)),")\n")
    }

    if (modelnumber == 3) {
       cat("\n*** TWO-BINOMIAL MODEL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
       cat("\nParameter estimates\n")
       cat("Binomial p1 =",pr.pbin,"(logit =",pr.logitpbin,", SE = ",
            format(round(sqrt(vcov[3,3]),digits=4)),")\n")
       cat("Binomial p2 =",pr.pbbin,"(logit =",pr.logitpbbin,", SE = ",
            format(round(sqrt(vcov[2,2]),digits=4)),")\n")
       cat("Proportion p1 =",pr.alphabbin,"(logit =",pr.logitalphabbin,", SE = ",
            format(round(sqrt(vcov[4,4]),digits=4)),")\n")
    }

    if (modelnumber == 4) {
       cat("\n*** BINOMIAL/BETA-BINOMIAL MODEL ***\n")
       cat("Minimised negative log-likelihood value =", format(round(loglmin, digits=4)),"\n")
       cat("\nParameter estimates\n")
       cat("Binomial p =",pr.pbin,"(logit =",pr.logitpbin,", SE = ",
            format(round(sqrt(vcov[2,2]),digits=4)),")\n")
       cat("Beta-binomial p =",pr.pbbin,"(logit =",pr.logitpbbin,", SE = ",
            format(round(sqrt(vcov[3,3]),digits=4)),")\n")
       cat("Beta-binomial theta=",pr.theta,"(log =",pr.logtheta,", SE = ",
            format(round(sqrt(vcov[4,4]),digits=4)),")\n")
       cat("Proportion beta-bin =",pr.alphabbin,"(logit =",pr.logitalphabbin,", SE = ",
            format(round(sqrt(vcov[5,5]),digits=4)),")\n")
    }
    cat("Number of unseen individuals =", pr.n0,"(log =",pr.logn0,", SE = ",
            format(round(sqrt(vcov[1,1]),digits=3)),")\n")
}

    cat("Number of individuals seen = ",format(round(sum(freq),digits=0)),"\n")
    cat("\nEstimated total population (N)=", pr.ntot,"\n")
    cat(confpercent,"% profile limits for N\n")
    if (lower < 0 & upper < 0) {
       cat("lower =",format(round(abs(lower),digits=1)), "(number seen)   upper = ", 
           "Infinite (probably)\n\n")
    }
    if (lower < 0 & upper > 0) {
       cat("lower =",format(round(abs(lower),digits=1)), "(number seen)   upper = ", 
           format(round(abs(upper),digits=1)),"\n\n")
    }
    if (lower > 0 & upper < 0) {
       cat("lower =",format(round(abs(lower),digits=1)), "upper = ", 
           "Infinite (probably)\n\n")
    }
    if (lower > 0 & upper > 0) {
       cat("lower =",format(round(abs(lower),digits=1)), "upper = ", 
           format(round(abs(upper),digits=1)),"\n\n")
    }


    if (!noses) {
       cat("Correlation matrix of TRANSFORMED parameters\n")
       cat("(First parameter is log(n0). Other parameters\n")
       cat(" are in the order given above)\n")
       correlations <- cov2cor(vcov)
       print(correlations)
       cat("\n")
    }

      ### Fitted values
       fitprob <- binbbinprob( kmax, pbin, pbbin, theta, alphabbin)
       tmp = fitprob[2:(length(freq)+1)]
       Fitted = sum(freq) * tmp / (1-fitprob[1])
       print(cbind(freq, Fitted), digits=3)
       if ( (length(freq)+1) < length(fitprob)) {
          tailfreq = sum(fitprob[(length(freq)+1):length(fitprob)]) / (1-fitprob[1])
          tailfreq = sum(freq) * tailfreq
          cat("\nTail frequency = ",format(round(tailfreq, digits=2)),"\n\n")
       }

}


## ============================================================================
   fitonemodel <- function(freq, kmax, model="bin", initparam=NULL,
                      nrandstart=50, plotprofile=FALSE, logn0values=NULL,
                      confpercent=95, nrandprof=20, print=1) {
## ============================================================================
   ## Arguments are as follows:
   ## 
   ##   freq   an array of frequencies. freq[k] is the number of animals
   ##          captured on k occasions
   ##
   ##   kmax   the maximum number of capture occasions. Note that the 
   ##          dimension of array freq may be less than kmax, as terminal
   ##          zero frequencies may be omitted. For example, the taxicab
   ##          frequencies may be given as  
   ##              freq = c(142, 81, 49, 7, 3, 1, 0, 0, 0, 0)
   ##          or just as 
   ##              freq = c(142, 81, 49, 7, 3, 1)
   ##
   ##   model  one of "bin", "betabin", "twobin", "binbbin"
   ##
   ##   initparam   array of initial parameter values. If unset, the program
   ##               chooses default values if nrandstart=0. Otherwise uses
   ##               random starting values.
   ## 
   ##   nrandstart  number of random starting values. If zero, program uses
   ##               fixed starting values, either those supplied as argument
   ##               initparam, or default values chosen by the program. 
   ##               Otherwise, best results from random starting points are
   ##               presented by the program.
   ##   plotprofile whether to plot a profile log-likelihood wrt n0 (default
   ##               FALSE)
   ##   logn0values values at which profile is to be calculated for plot
   ##   confpercent level for confidence intervals
   ##   nrandprof   number of random start points to use for profile
   ##               likelihood calculations (reset to 1 if <= 1)

        ## --------------------------------
        ## Set appropriate value of kobsmax
        ## --------------------------------
    kobsmax = max( c(1:length(freq)) * (freq > 0))
    if (kobsmax > kmax) {
       stop('Error: kobsmax > kmax')
    }

        ## -------------------------------------------------------------
        ## Determine model number, 1=bin, 2=betabin, 3=twobin, 4=binbbin
        ## -------------------------------------------------------------
    modelnumber = sum( c(1:4) * (model == c("bin", "betabin", "twobin", "binbbin")))
    if (modelnumber == 0) {
       stop('Error: unknown model')
    }

        ## --------------------------------------------
        ## Decide whether random starts and set initial 
        ## parameter values if appropriate
        ## --------------------------------------------
    nmodparam = c(2,3,4,5)     # number of parameters for each model
    if ( (nmodparam[modelnumber]) > kmax) {
       stop ('Error: not enough data points')
    }

        ## ----------------------------
        ## Reset nrandprof if necessary 
        ## ----------------------------
    if (nrandprof < 1) {
       nrandprof = 1
    }

        ## -------------------------
        ## Set/check starting values 
        ## -------------------------
    randstart = FALSE
    if (length(initparam) == 0) {
       if (nrandstart == 0) {
          initparam = rep(0, times=nmodparam[modelnumber])
       }  else {
          randstart = TRUE
       }
    } else if (length(initparam) != nmodparam[modelnumber]) {
             stop('Error: initparam has wrong number of elements')
    }


        ## -------------------------
        ## Fit the appropriate model
        ## -------------------------

    minlogl = 9.99E29


    if (randstart) {

           ##=========================================================
           ## Calculate a preliminary coverage-based estimator for N0
           ## (number unseen) and use this to set a preliminary value
           ## for logN0 called logN0guess. This is used to centre
           ## the random starting points for logN0. Two forms of 
           ## coverage estimator are used; the Good-Turing estimator
           ## for model="bin" and an estimator recommended by 
           ## Ashbridge and Goudie (2000, Communications in Statistics
           ## - Simulation and Computation, 29, 1215-1237) for all
           ## other models.
           ##=========================================================
#ä       if (model=="bin") {
              ## Good-Turing estimator
          ntotal = sum( c(1:kobsmax) * freq[1:kobsmax])
          goodunseen = sum(freq)*freq[1] / (ntotal - freq[1])
          logN0guess = log(goodunseen)
       # } else {
       #        ## Ashbridge-Goudie estimator
       #    KA = sum(freq)
       #    CA = 1 - freq[1]/KA + 2/(kmax-1)*freq[2]/KA - 
       #         6/(kmax-1)/(kmax-2) * freq[3]/KA
       #    CA = min(CA,1)
       #    AGestimator = sum(freq[1:kobsmax] / (1 - (1 - c(1:kobsmax) * CA/kmax)^kmax))
       #    AGunseen = AGestimator - sum(freq)
       #    logN0guess = log(AGunseen)
       # }

#llvals = numeric(randstart)
       for (j in 1: nrandstart) {
           initparam = ranstart(nmodparam[modelnumber])

              ## Centre the estimator of logN0 around logN0guess
              ## Also, prevent very large values of logN0
           initparam[1] = initparam[1] + logN0guess
           if (initparam[1] > 15) initparam[1] = 15

              ## Do the optimisation from this starting point

           rndmod <- optim(initparam, loglik,, modelnumber, freq, kobsmax, kmax)
#      llvals[j] = rndmod$val


              ## Update if this random start is an improvement on those 
              ## found previously. Ignore occasional negative values of
              ## the log-likelihood which can occur for very extreme parameter
              ## values, as a result of numerical problems.
           if (rndmod$value < minlogl & rndmod$value > 0) {
              minlogl <- rndmod$value
              fitmodel <- rndmod
#      tmpmod <- optim(fitmodel$par, loglik,, modelnumber, freq, kobsmax, kmax)
#      llvals[j] = -tmpmod$val
#     cat(j, rndmod$val, tmpmod$val, "\n")
           }
#     cat(j, rndmod$val, rndmod$par, "\n")
       }

             ## Refit the model from the best current estimate
             ## This typically improves the fit slightly
       fitmodel <- optim(fitmodel$par, loglik,, modelnumber, freq, kobsmax, kmax)
#       cat("Final", fitmodel$val, fitmodel$par, "\n")
    } else {

             ## Specific user-supplied starting values, not random
           fitmodel <- optim(initparam, loglik,, modelnumber, freq, kobsmax, kmax)
    }

        ## -----------------------------------------------
        ## Find the approximate variance-covariance matrix
        ## -----------------------------------------------
    fithess <- optim(fitmodel$par, loglik,, modelnumber, freq, kobsmax, kmax, hessian=TRUE)
    conditionnumber = kappa(fithess$hessian)
    determinant = det(fithess$hessian)

    if (conditionnumber > 1.0E7 | determinant <= 0) {
       vcov = diag(length(fitmodel$par))
       noses = TRUE
    } else {
       vcov = solve(fithess$hessian)
       noses = FALSE
    }

        ## ------------------------------------------------------
        ## Now draw the log-likelihood profile for n0 if required
        ## ------------------------------------------------------
    if (plotprofile | (length(logn0values) != 0)) {
       n0hat = exp(fitmodel$par[1])
       if (length(logn0values) == 0) {
          logn0values = seq(log(n0hat/100), log(n0hat*100), length.out=25)
       }
       profn0 = logn0values
       for (k in 1:length(logn0values)) {
           n0 = exp(logn0values[k])
           minlogl = 9.99E29 
               ## First try for starting values is the mle
           starter = fitmodel$par[2:length(fitmodel$par)]
#print(starter)
           for (j in 1: nrandprof) {
               if (modelnumber==1) options(warn=-1)
               rndmod <- optim(starter, profloglik,, modelnumber, freq, kobsmax, kmax, n0)
               options(warn=0)
               if (rndmod$value < minlogl) {
                  minlogl <- rndmod$value
                  profmodel <- rndmod
               }
#              cat(j, rndmod$value, rndmod$convergence, "\n")
                   ## Then switch to random starting values
               starter = ranstart(nmodparam[modelnumber]-1)
           }
               ## Refit the model from the best current estimate
               ## This typically improves the fit slightly
           profmodel <- optim(profmodel$par, profloglik,, modelnumber, freq, kobsmax, kmax, n0)
           profn0[k] = profmodel$val
       }
       plot(logn0values, profn0)
       lines(logn0values, profn0)
    }

## TEMPORARY PAUSE
## return(llvals)

        ## ==================================
        ## Now calculate the profile interval
        ## ==================================
     logn0hat = fitmodel$par[1]
     lowlim = exp(logn0hat - log(100))
     upplim = exp(logn0hat + log(100))

     loglmin = fitmodel$val
     starter = fitmodel$par[2:length(fitmodel$par)]
     uppbound = ploglik(logn0hat+log(100), loglmin, starter, modelnumber, freq, 
                        kobsmax, kmax, confpercent)
#     cat("uppbound = ",uppbound,"\n")
     starter = fitmodel$par[2:length(fitmodel$par)]
     lowbound = ploglik(logn0hat-log(100), loglmin, starter, modelnumber, freq, 
                        kobsmax, kmax, confpercent)
#     cat("lowbound = ",lowbound,"\n")

       ## Set the lower bound
     if (lowbound > 0) {
        starter = fitmodel$par[2:length(fitmodel$par)]
        findlow = uniroot(ploglik,,loglmin, starter, modelnumber, freq, 
                          kobsmax, kmax, confpercent, lower=log(lowlim), upper=logn0hat)
        lower = findlow$root
        lower = exp(lower) + sum(freq)
     } else {
        lower = -sum(freq)
     }

       ## Set the upper bound
     if (uppbound > 0) {
        starter = fitmodel$par[2:length(fitmodel$par)]
        findupp = uniroot(ploglik,, loglmin, starter, modelnumber, freq, 
                          kobsmax, kmax, confpercent, lower=logn0hat, upper=log(upplim))
        upper = findupp$root
        upper = exp(upper) + sum(freq)
     } else {
        upper = -1
     }
    if (print==1) {
       printest(fitmodel$val, fitmodel$par, modelnumber, freq, kobsmax, kmax, 
                vcov, noses, lower, upper, confpercent)
    }
    if (print==2) {
       printsimres(fitmodel$val, fitmodel$par, modelnumber, freq, kobsmax, kmax, 
                vcov, noses, lower, upper, confpercent)
    }
}

## ====================================
   fitallmodels <- function(freq, kmax)
## ====================================
{
        ## Quick program for fitting all models with default settings 
        ## of most parameters
   if (kmax < 5) {
      stop('Error: kmax must be >= 5')
   }
   fitonemodel(freq, kmax, model="bin", nrandstart=5, plotprofile=FALSE, nrandprof=10) 
   fitonemodel(freq, kmax, model="betabin", nrandstart=5, plotprofile=FALSE, nrandprof=10) 
   fitonemodel(freq, kmax, model="twobin", nrandstart=100, plotprofile=FALSE, nrandprof=10) 
   fitonemodel(freq, kmax, model="binbbin", nrandstart=200, plotprofile=FALSE, nrandprof=10) 
}





## ============================================================================
   printsimres <- function(loglmin, param, modelnumber, freq, kobsmax, kmax, vcov, 
                        noses, lower, upper, confpercent)
## ============================================================================
{
        ##--------------------------
        ## Produces formatted output
        ##--------------------------
    n0 = exp(param[1])    # Size of zero class
    sumf = sum(freq)      # Total observed frequency
    ntot = sumf + n0      # Current estimate of popn size

        ##-----------------------------
        ## Set parameters appropriately
        ##-----------------------------
    parvals = unpack(param, modelnumber)
    pbin = parvals[1]
    pbbin = parvals[2]
    theta = parvals[3]
    alphabbin = parvals[4]

    logitpbin = log(pbin/(1-pbin))
    logitpbbin = log(pbbin/(1-pbbin))
    logitalphabbin = log(alphabbin/(1-alphabbin))
    logn0 = param[1]
    logtheta = log(theta) 

    pr.pbin = format(round(pbin, digits=4))
    pr.pbbin = format(round(pbbin, digits=4))
    pr.alphabbin = format(round(alphabbin, digits=4))
    pr.n0 = format(round(n0, digits=1))
    pr.ntot = format(round(ntot, digits=1))
    pr.theta = format(round(theta, digits=3))

    pr.logitpbin = format(round(logitpbin, digits=3))
    pr.logitpbbin = format(round(logitpbbin, digits=3))
    pr.logitalphabbin = format(round(logitalphabbin, digits=3))
    pr.logn0 = format(round(logn0, digits=2))
    pr.logtheta = format(round(logtheta, digits=3))

    cat(format(round(loglmin, digits=4))," ")
    cat(pr.ntot," ")
    if (lower < 0 & upper < 0) {
       cat(format(round(-abs(lower),digits=3)), "-1 ")
    }
    if (lower < 0 & upper > 0) {
       cat(format(round(-abs(lower),digits=3)), 
           format(round(abs(upper),digits=3))," ")
    }
    if (lower > 0 & upper < 0) {
       cat(format(round(abs(lower),digits=3)),
           "-1 ")
    }
    if (lower > 0 & upper > 0) {
       cat(format(round(abs(lower),digits=3)), 
           format(round(abs(upper),digits=3))," ")
    }
    cat(pr.pbin, pr.pbbin, pr.theta, pr.alphabbin,"\n")

}



## ============================================================
   simdata <- function( kmax, N, pbin, pbbin, theta, alphabbin)
## ============================================================
{
        ## Simulate a set of frequencies. First set the probabilities
   prob = binbbinprob( kmax, pbin, pbbin, theta, alphabbin)
   as.vector(rmultinom(1, N, prob))
}
