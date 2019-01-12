####################################################################################################################################################
    # Phase I Xbar Chart Based on Champ and Jones(2004); see also Yao et al. (2017), and Yao and Chakraborti (2018)
####################################################################################################################################################

#use too many functions from package mvtnorm
require(mvtnorm)
#require(adehabitatLT)
#dchi <- adehabitatLT::dchi

#focus on these 3 functions from package mvtnorm


#pmvt <- mvtnorm::pmvt
#qmvt <- mvtnorm::qmvt
#pmvnorm <- mvtnorm::pmvnorm


####################################################################################################################################################

#source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/MKLswitch.R')

####################################################################################################################################################
    #parts of getting the charting constant l
####################################################################################################################################################

#########################################################################################################
    # Calculating Owen's table [92]
#########################################################################################################

# pdf of W

W.f <- function(w, n) { # n is sample size
  
  integrand <- function(x, w) {
    
    (pnorm(w + x) - pnorm(x)) ^ (n - 2) * dnorm(w + x) * dnorm(x)
    
  }
  
  n * (n - 1) * integrate(integrand, -Inf, Inf, w = w)$value
  
}

# moment of W

W.moment <- function(r, n) { # r is order of moment and n is sample size
  
  integrand <- function(w, r, n) {
    
    w ^ r * W.f(w, n)
    
  }
  
  integrand <- Vectorize(integrand, 'w')
  
  integrate(integrand, 0, Inf, r = r, n = n)$value
  
}

# the result is the 1st moment
W.moment(r = 1, n = 5)
# result
# E(W) = 2.325926


#########################################################################################################
    # calculating the approximate degrees of freedom and the scalar in H. A. David pp 243
#########################################################################################################

# pars finding function for the approximate method
# k is number of batches and n is sample size

pars.root.finding <- function(k, n, lower = 1e-6) {
  
  root.finding <- function(nu, k, n, ratio) {
    
    ratio / k - nu / 2 / pi * beta(nu / 2, 1 / 2) ^ 2 + 1
    
  }
  
  M <- W.moment(1, n)
  M2 <- W.moment(2, n)
  
  V <- M2 - M ^ 2
  
  ratio <- V / M ^ 2
  
  cat('dn2/Vn = ', 1 / ratio, '\n')
  
  nu <- uniroot(root.finding, c(lower, k * n), k = k, n = n, ratio = ratio)$root
  
  C <- sqrt(V / k + M ^ 2)
  
  c(nu, C)
  
}

#########################################################################################################

d2.f <- function(n) {

    integrand <- function(x, n) {
    
        1 - (1 - pnorm(x)) ^ n - pnorm(x) ^ n
    
    }
    
    integrate(integrand, -Inf, Inf, n = n)$value

}

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function

ub.f <- function(n, nu, var.option = 'c4') {

        if (var.option == 'c4') {
            corr.par <- c4.f(nu)
        } else if (var.option == 'd2') {
            corr.par <- d2.f(n)
        }
        else {
            corr.par <- 1
        }
        
        corr.par
        
}

PH1.corr.f <- function(k, off.diag = - 1 / (k - 1)){
                                                                                   #correlation matrix
    crr <- diag(k)
    crr[which(crr == 0)] <- off.diag

    crr

}

####################################################################################################################################################
    #get l using the multivariate t cdf
####################################################################################################################################################

PH1.get.cc.mvt <- function(
                    n.dim
                    ,n.batch
                    ,n.samplesize
                    ,nu
                    ,sigma
                    ,corr.par
                    #,alternative = '2-sided'
                    ,var.option = 'c4'
                    ,maxiter = 10000
) {

}{

    alternative = '2-sided'                                                   #turn off the alternative

    pu <- 1 - FAP

    #corr.par <- ub.f(n.samplesize, nu, var.option)

    L <- ifelse(
            alternative == '2-sided',
            qmvt(pu, df = nu, sigma = sigma, tail = 'both.tails', maxiter = maxiter)$quantile,
            qmvt(pu, df = nu, sigma = sigma, maxiter = maxiter)$quantile
    )

    c.i <- L * corr.par * sqrt((n.batch - 1) / n.batch)             #get c

    list(l = L, c.i = c.i)

}

####################################################################################################################################################
    #get L by multivariate Normal
####################################################################################################################################################

PH1.joint.pdf.mvn.chisq <- function(
                                Y
                                ,c.i
                                ,n.dim
                                ,n.batch
                                ,n.samplesize
                                ,nu
                                ,C
                                ,sigma
                                ,corr.par
                                #,alternative = '2-sided'
                                ,var.option = 'c4')
{

    alternative = '2-sided'                                                   #turn off the alternative

    s <- length(Y)

    #corr.par <- ub.f(n.samplesize, nu, var.option)

    L <- c.i / sqrt((n.batch - 1) / n.batch * nu) * sqrt(Y) / corr.par * C

    dpp <- lapply(
            1:s,
            function(i){

                LL <- rep(L[i], n.dim)

                ifelse(
                    alternative == '2-sided',
                    pmvnorm(lower = -LL, upper = LL, sigma = sigma),
                    pmvnorm(lower = rep(-Inf, n.dim), upper = LL, sigma = sigma)
                )

            }

    )

    dpp <- unlist(dpp)

    dpp * dchisq(Y, df = nu)


}

PH1.root.mvn.F <- function(
                    c.i
                    , n.dim
                    , n.batch
                    , n.samplesize
                    , nu
                    , C
                    , sigma
                    , corr.par
                    , pu
                    #, alternative = '2-sided'
                    , var.option = 'c4'
                    , subdivisions = 2000
                    , rel.tol = 1e-2)
{
    alternative = '2-sided'                                                   #turn off the alternative
                                                          
    pp <- integrate(
            PH1.joint.pdf.mvn.chisq
            ,lower = 0
            ,upper = Inf
            ,c.i = c.i
            ,n.dim = n.dim
            ,n.batch = n.batch
            ,n.samplesize = n.samplesize
            ,nu = nu
            ,C = C
            ,sigma = sigma
            ,corr.par = corr.par
            #,alternative = '2-sided'
            ,var.option = 'c4'
            ,subdivisions = subdivisions
            ,rel.tol = rel.tol
        )$value

    pu - pp


}


PH1.get.cc.mvn <- function(
                 n.dim
                 ,n.batch
                 ,n.samplesize
                 ,nu
                 ,C = 1
                 ,FAP = 0.1
                 ,sigma = PH1.corr.f(k = n.dim, off.diag = off.diag)
                 ,corr.par = ub.f(n.samplesize, nu, var.option)
                 #,alternative = '2-sided'
                 ,var.option = 'c4'
                 ,interval = c(1, 7)
                 ,maxiter = 10000
                 ,subdivisions = 2000
                 ,tol = 1e-2

){
    alternative = '2-sided'                                                   #turn off the alternative

    #corr.par <- ub.f(n.samplesize, nu, var.option)

    pu <- 1 - FAP

    c.i <- uniroot(
            PH1.root.mvn.F
            ,n.dim = n.dim
            ,n.batch = n.batch
            ,n.samplesize = n.samplesize
            ,nu = nu
            ,C = C
            ,sigma = sigma
            ,corr.par = corr.par
            ,pu = pu
            #, alternative = '2-sided'
            ,var.option = var.option
            ,subdivisions = subdivisions
            ,tol = tol
            ,rel.tol = tol
            ,maxiter = maxiter
    )$root
    
    L <- c.i / corr.par * sqrt(n.batch / (n.batch - 1)) * C

    list(l = L, c.i = c.i)


}

#get.cc.mvn(20, 80)

####################################################################################################################################################

PH1.get.cc <- function(
            n.dim
            ,n.batch
            ,n.samplesize
            ,nu
            ,FAP = 0.1
            ,corr.P = PH1.corr.f(k = n.dim, off.diag = off.diag)
            #,alternative = '2-sided'
            ,var.option = 'c4'
            ,method = 'direct'
            ,var.est = 'MS'
            ,lower.MR = 1e-6
            ,maxiter = 10000            
            ,indirect.interval = c(1, 7)
            ,indirect.subdivisions = 100L
            ,indirect.tol = .Machine$double.eps^0.25


){
    alternative = '2-sided'                                                   #turn off the alternative

    is.int <- ifelse(nu == round(nu), 1, 0)

    
    corr.par <- ub.f(n.samplesize, nu, var.option)
    
    C <- 1
    
    if (var.est == 'MS') {
    
        pars.MR <- pars.root.finding(n.batch, n.samplesize, lower = lower.MR)
        
        nu <- pars.MR[1]
        C <- pars.MR[2]
        
    }
    
    if (method == 'direct') {
    
        if (is.int == 1) {
            PH1.get.cc.mvt(
            
                 n.dim = n.dim
                ,n.batch = n.batch
                ,n.samplesize = n.samplesize
                ,nu = nu
                ,sigma = sigma
                ,corr.par = corr.P
                #,alternative = '2-sided'
                ,var.option = var.option
                ,maxiter = maxiter
            
            )
        } else if (is.int == 0) {
        
             stop('Nu is not an integer. Please use the indirect method instead.')
        
        }
    
    } else if (method == 'indirect') {
    
        if (is.int == 1) {
        
            cat('Nu is an integer. Using the indirect method may slow the computation process down.', '\n')
        
        }
        
        PH1.get.cc.mvn(
        
             n.dim = n.dim
            ,n.batch = n.batch
            ,n.samplesize = n.samplesize
            ,nu = nu
            ,C = C
            ,FAP = FAP
            ,sigma = sigma
            ,corr.par = corr.par
            #,alternative = '2-sided'
            ,var.option = 'c4'
            ,interval = indirect.interval
            ,maxiter = maxiter
            ,subdivisions = indirect.subdivisions
            ,tol = indirect.tol

        )
    
    } else {
    
        stop('Unknown method. Please select the direct method or the indirect method.')
    
    }
    

}



####################################################################################################################################################
    #above need to test, below need to revise
####################################################################################################################################################



















####################################################################################################################################################
    #reverse the process
####################################################################################################################################################
#only support the direct method
#PH1.get.FAP0 <- function(
#            c.i
#            ,k
#            ,nu
#            ,off.diag = -1/(k - 1)
#            ,alternative = '2-sided'
#            ,var.option = TRUE
#            ,maxiter = 10000
#            ,indirect.subdivisions = 100L
#            ,indirect.tol = .Machine$double.eps^0.25
#
#){
#
#    #method <- 'direct'
#
#    #Phase1 <- TRUE
#
#    is.int <- ifelse(nu == round(nu), 1, 0)
#
#    #if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
#
#    corr.P <- PH1.corr.f(k = k, off.diag = off.diag)
#
#    corr.par <- corr.par.f(nu, var.option)
#
#    L <- c.i / corr.par * sqrt(k / (k - 1))
#
#
#    if (method == 'direct' & is.int == 1) {                       #using multivariate T to obtain L and K
#
#        ll <- rep(L, k)
#
#        ifelse(
#            alternative == '2-sided',
#            1 - pmvt(lower = -ll, upper = ll, df = nu, corr = corr.P, algorithm = TVPACK, abseps = 1e-12),
#            1 - pmvt(lower = -Inf, upper = ll, df = nu, corr = corr.P, algorithm = TVPACK, abseps = 1e-12)
#        )
#
#    }
#
#
#
#
#}

####################################################################################################################################################
    #Build Control Chart
####################################################################################################################################################

PH1XBAR <- function(
			X
            ,c.i = NULL
			,FAP = 0.1
			,off.diag = -1/(k - 1)
			#,alternative = '2-sided'
            ,model = 'ANOVA-based'
            ,var.option = TRUE
			,plot.option = TRUE
			,maxiter = 10000
			,method = 'direct'
			,indirect.interval = c(1, 7)
			,indirect.subdivisions = 100L
			,indirect.tol = .Machine$double.eps^0.25
) {

    alternative = '2-sided'                                                   #turn off the alternative


	k <- dim(X)[1]
	n <- dim(X)[2]

	X.bar <- rowMeans(X)

	X.bar.bar <- mean(X)

    if (model == 'ANOVA-based') {

        nu <- k - 1

        corr.par <- ub.f(n, nu, var.option)

        sigma.v <- sqrt(var(X.bar)) / corr.par

    } else if (model == 'basic') {

        nu <- k * (n - 1)

        corr.par <- ub.f(n, nu, var.option)

        sigma.v <- sqrt(sum(apply(X, 1, var)) / k) / corr.par / sqrt(n)

    } else {

        stop("need to specify whether it is based on the ANOVA-based model or others")

    }

    if (is.null(c.i)) {
        c.i <- PH1.get.cc(
                k = k
                ,nu = nu
                ,FAP = FAP
                ,off.diag = off.diag
                #,alternative = alternative
                ,var.option = var.option
                ,maxiter = maxiter
                ,method = method
                ,indirect.interval = indirect.interval
                ,indirect.subdivisions = indirect.subdivisions
                ,indirect.tol = indirect.tol
         )$c.i
    } else {

            c.i <- c.i

    }

	LCL <- X.bar.bar - c.i * sigma.v
	UCL <- X.bar.bar + c.i * sigma.v

	if (plot.option == TRUE){

		plot(c(1, k), c(LCL, UCL), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Mean', type = 'n')

		axis(side = 1, at = 1:k)

		points(1:k, X.bar, type = 'o')
		points(c(-1, k + 2), c(LCL, LCL), type = 'l', col = 'red')
		points(c(-1, k + 2), c(UCL, UCL), type = 'l', col = 'red')
		points(c(-1, k + 2), c(X.bar.bar, X.bar.bar), type = 'l', col = 'blue')
		text(round(k * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
		text(round(k * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)


	}

	list(CL = X.bar.bar, sigma = sigma.v, PH1.cc = c.i, k = k, nu = nu, LCL = LCL, UCL = UCL, CS = X.bar)

}


####################################################################################################################################################
    #Example about table 6.3 (Montgomery, 2013)
####################################################################################################################################################
PH1XBAR.data <- function(){
	matrix(c(
		74.03	,	74.002	,	74.019	,	73.992	,	74.008	,
		73.995	,	73.992	,	74.001	,	74.011	,	74.004	,
		73.988	,	74.024	,	74.021	,	74.005	,	74.002	,
		74.002	,	73.996	,	73.993	,	74.015	,	74.009	,
		73.992	,	74.007	,	74.015	,	73.989	,	74.014	,
		74.009	,	73.994	,	73.997	,	73.985	,	73.993	,
		73.995	,	74.006	,	73.994	,	74	,	74.005	,
		73.985	,	74.003	,	73.993	,	74.015	,	73.988	,
		74.008	,	73.995	,	74.009	,	74.005	,	74.004	,
		73.998	,	74	,	73.99	,	74.007	,	73.995	,
		73.994	,	73.998	,	73.994	,	73.995	,	73.99	,
		74.004	,	74	,	74.007	,	74	,	73.996	,
		73.983	,	74.002	,	73.998	,	73.997	,	74.012	,
		74.006	,	73.967	,	73.994	,	74	,	73.984	,
		74.012	,	74.014	,	73.998	,	73.999	,	74.007	,
		74	,	73.984	,	74.005	,	73.998	,	73.996	,
		73.994	,	74.012	,	73.986	,	74.005	,	74.007	,
		74.006	,	74.01	,	74.018	,	74.003	,	74	,
		73.984	,	74.002	,	74.003	,	74.005	,	73.997	,
		74	,	74.01	,	74.013	,	74.02	,	74.003	,
		73.982	,	74.001	,	74.015	,	74.005	,	73.996	,
		74.004	,	73.999	,	73.99	,	74.006	,	74.009	,
		74.01	,	73.989	,	73.99	,	74.009	,	74.014	,
		74.015	,	74.008	,	73.993	,	74	,	74.01	,
		73.982	,	73.984	,	73.995	,	74.017	,	74.013

	), ncol = 5, byrow = T)
}





####################################################################################################################################################
    # Phase II Xbar chart. Please see Yao and Chakraborti (2018)
####################################################################################################################################################

PH2.inner.normal <- function(Z, Y, c.ii, k, n, nu = k - 1, var.option = 'c4'){

    corr.par <- ub.f(n, nu, var.option)

    qn <- Z / sqrt(k) + c.ii / corr.par / sqrt(nu) * sqrt(Y)

    pnorm(qn)

}

PH2.CFAR.intgrand <- function(u, v, c.ii, k, nu = k - 1, var.option = TRUE){

    Z <- qnorm(u)
    Y <- qchisq(v, nu)

    p1 <- PH2.inner.normal(Z, Y, c.ii, k, nu, var.option)
    p2 <- PH2.inner.normal(Z, Y, -c.ii, k, nu, var.option)

    1 - p1 + p2
}

PH2.CARL.intgrand <- function(u, v, c.ii, k, nu = k - 1, var.option = TRUE) 1 / PH2.CFAR.intgrand(u, v, c.ii, k, nu, var.option)


####################################################################################################################################################


PH2.get.cc.uc <- function(ARL0, k, nu = k - 1, var.option = TRUE, interval = c(1, 10), u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

        PH2.root.finding.uc <- function(c.ii, k, nu = k - 1, var.option = TRUE, ARL0, u, v) {

            ARL0 - mean(PH2.CARL.intgrand(u, v, c.ii, k, nu, var.option))

    }

    k <- k

    rt <- uniroot(PH2.root.finding.uc, interval = interval, k = k, nu = nu, var.option = var.option, ARL0 = ARL0, u = u, v = v, tol = tol)$root

    rt


}

####################################################################################################################################################


PH2.root.finding.EPC <- function(p, k, nu = k - 1, c.ii, ARLb, var.option = TRUE, u = runif(100000), v = runif(100000)){

    1 - p - mean(PH2.CARL.intgrand(u, v, c.ii, k, nu, var.option) >= ARLb)

}

PH2.get.cc.EPC <- function(p, k, nu = k - 1, eps = 0.1, ARL0 = 370, var.option = TRUE, interval = c(1, 10), u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

    ARLb <- (1 - eps) * ARL0

    rt <- uniroot(
            PH2.root.finding.EPC,
            interval = interval,
            p = p,
            k = k,
            nu = nu,
            ARLb = ARLb,
            var.option = var.option,
            u = u,
            v = v,
            tol = tol
        )$root

    rt

}

#rnum <- 100000
#u <- runif(rnum)
#v <- runif(rnum)

#debug(PH2.get.c.ii.EPC)
#PH2.get.c.ii.EPC(p = 0.05, k = 100, eps = 0.1, ARL0 = 500, var.option = TRUE, interval = c(1, 5), u = u, v = v)


PH2.get.ARLb.EPC <- function(p, k, nu = k - 1, c.ii, ARL0 = NULL, var.option = TRUE, interval = c(1, 1000), u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

    rt <- uniroot(
            PH2.root.finding.EPC,
            interval = interval,
            p = p,
            k = k,
            nu = nu,
            c.ii = c.ii,
            var.option = var.option,
            u = u,
            v = v,
            tol = tol
        )$root

    if (is.null(ARL0)) {
        list(ARLb = rt, eps = NULL, ARL0 = NULL)
    } else {
        list(ARLb = rt, eps = 1 - rt / ARL0, ARL0 = ARL0)
    }


}
#PH2.get.ARLb.EPC(p = 0.1, k = 100, c.ii = 3, var.option = TRUE, u = u, v = v)
#PH2.get.ARLb.EPC(p = 0.1, k = 100, c.ii = 3, ARL0 = 370, var.option = TRUE, u = u, v = v)

PH2.get.k.EPC <- function(p, c.ii, eps = 0.1, ARL0 = 370, var.option = TRUE, interval = c(1000, 3000), model = 'ANOVA-based', n = 10, u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

    if (model == 'ANOVA-based') {

        PH2.root.finding.k.EPC <- function(p, k, c.ii, ARLb, var.option = TRUE, n = 10, u, v){

            nu <- k - 1

            PH2.root.finding.EPC(p, k, nu, c.ii, ARLb, var.option, u, v)

        }


    } else if (model == 'basic') {

        PH2.root.finding.k.EPC <- function(p, k, c.ii, ARLb, var.option = TRUE, n = 10, u, v){

            nu <- (n - 1) * k

            PH2.root.finding.EPC(p, k, nu, c.ii, ARLb, var.option, u, v)

        }

    } else {

        stop("need to specify whether it is based on the ANOVA-based model or others")
    }


    ARLb <- (1 - eps) * ARL0

    rt <- uniroot(
            PH2.root.finding.k.EPC,
            interval = interval,
            p = p,
            c.ii = c.ii,
            ARLb = ARLb,
            var.option = var.option,
            n = n,
            u = u,
            v = v,
            tol = tol
        )$root

    list(k = rt, PH1.sample.size = rt * n)

}

####################################################################################################################################################


PH2XBAR <- function(
            X
            ,PH1.info = list(X = NULL, mu = NULL, sigma = NULL, k = NULL, n = NULL, model = "ANOVA-based")
            ,c.ii.info = list(c.ii = NULL, method = 'UC', ARL0 = 370, p = 0.05, eps = 0.1, interval.c.ii.UC = c(1, 3.2), interval.c.ii.EPC = c(1, 10), UC.tol = .Machine$double.eps^0.25, EPC.tol = .Machine$double.eps^0.25)
            #,alternative = '2-sided'
            ,var.option = TRUE
            ,plot.option = TRUE
            ,maxsim = 100000
) {

    alternative = '2-sided'                                                   #turn off the alternative

    k.ii <- dim(X)[1]

    X.bar <- rowMeans(X)


    if (is.null(PH1.info$X)){

        k <- PH1.info$k

        if (PH1.info$model == 'ANOVA-based') {

            nu <- k - 1

        } else if (PH1.info$model == 'basic') {

            if (is.null(PH1.info$n)) {

                stop('Subgroup size needs to be defined')

            } else {

                n <- PH1.info$n
                nu <- k * (n - 1)

            }

        } else {

            stop('The model in Phase 1 needs to be defined')

        }

        corr.par <- ub.f(n, nu, var.option)

        X.bar.bar <- PH1.info$mu
        sigma.v <- PH1.info$sigma / corr.par

    } else {

        PH1.chart <- PH1XBAR(X = PH1.info$X, c.i = 1, model = PH1.info$model, plot.option = FALSE)

        k <- PH1.chart$k
        nu <- PH1.chart$nu

        X.bar.bar <- PH1.chart$CL
        sigma.v <- PH1.chart$sigma

    }

    u <- runif(maxsim)
    v <- runif(maxsim)

    UC <- NULL
    EPC <- NULL

    n.c.ii <- 0

    UC.flag <- 0
    EPC.flag <- 0

    if (is.null(c.ii.info$c.ii)) {

        if (c.ii.info$method == 'both') {

            UC.flag <- 1
            EPC.flag <- 1

        } else if (c.ii.info$method == 'UC') {

            UC.flag <- 1

        } else if (c.ii.info$method == 'EPC'){

            EPC.flag <- 1

        } else {

            stop("Please check the method of obtaining charting constants in Phase 2")
        }

        if (UC.flag == 1) {

            n.c.ii <- n.c.ii + 1

            UC <- PH2.get.cc.uc(
                    ARL0 = c.ii.info$ARL0
                    , k = k
                    , nu = nu
                    , var.option = var.option
                    , interval = c.ii.info$interval.c.ii.UC
                    , u = u
                    , v = v
                    , tol = c.ii.info$UC.tol
                )

        }

        if (EPC.flag == 1){

            n.c.ii <- n.c.ii + 1

            EPC <- PH2.get.cc.EPC(
                    p = c.ii.info$p
                    , k = k
                    , nu = nu
                    , eps = c.ii.info$eps
                    , ARL0 = c.ii.info$ARL0
                    , var.option = var.option
                    , interval = c.ii.info$interval.c.ii.EPC
                    , u = u
                    , v = v
                    , tol = c.ii.info$EPC.tol
                )

        }

        c.ii <- list(UC = UC, EPC = EPC)

    } else {

            n.c.ii <- length(c.ii.info$c.ii)

            c.ii <- c.ii.info$c.ii

    }

    LCL <- X.bar.bar - unlist(c.ii) * sigma.v
    UCL <- X.bar.bar + unlist(c.ii) * sigma.v

    if (plot.option == TRUE){

        plot(c(1, k.ii), c(min(c(LCL, X.bar)), max(c(UCL, X.bar))), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Mean', type = 'n')

        points(1:k.ii, X.bar, type = 'o')

        axis(side = 1, at = 1:k.ii)
        points(c(-1, k.ii + 2), c(X.bar.bar, X.bar.bar), type = 'l', col = 'blue')

        for (ii in 1:n.c.ii){

            points(c(-1, k.ii + 2), c(LCL[ii], LCL[ii]), type = 'l', col = 'red')
            points(c(-1, k.ii + 2), c(UCL[ii], UCL[ii]), type = 'l', col = 'red')

            if (is.null(c.ii.info$c.ii)) {

                lab.in <- paste('(', names(c.ii)[ii], ')', sep = '')

            } else {

                lab.in <- NULL

            }

            ulab <- paste('UCL', lab.in, ' = ', round(UCL[ii], 4), sep = '')
            llab <- paste('LCL', lab.in, ' = ', round(LCL[ii], 4), sep = '')

            text(round(k.ii * 0.8), UCL[ii], ulab, pos = 1)
            text(round(k.ii * 0.8), LCL[ii], llab, pos = 3)

        }




    }

    list(CL = X.bar.bar, sigma = sigma.v, k = k, nu = nu, PH2.cc = unlist(c.ii), LCL = LCL, UCL = UCL, CS = X.bar)

}

####################################################################################################################################################



PH2XBAR.data <- function(){

    matrix(
        c(
            240, 243, 250, 253, 248,
            238, 242, 245, 251, 247,
            239, 242, 246, 250, 248,
            235, 237, 246, 249, 246,
            240, 241, 246, 247, 249,
            240, 243, 244, 248, 245,
            240, 243, 244, 249, 246,
            245, 250, 250, 247, 248,
            238, 240, 245, 248, 246,
            240, 242, 246, 249, 248,
            240, 243, 246, 250, 248,
            241, 245, 243, 247, 245,
            247, 245, 255, 250, 249,
            237, 239, 243, 247, 246,
            242, 244, 245, 248, 245,
            237, 239, 242, 247, 245,
            242, 244, 246, 251, 248,
            243, 245, 247, 252, 249,
            243, 245, 248, 251, 250,
            244, 246, 246, 250, 246,
            241, 239, 244, 250, 246,
            242, 245, 248, 251, 249,
            242, 245, 248, 243, 246,
            241, 244, 245, 249, 247,
            236, 239, 241, 246, 242,
            243, 246, 247, 252, 247,
            241, 243, 245, 248, 246,
            239, 240, 242, 243, 244,
            239, 240, 250, 252, 250,
            241, 243, 249, 255, 253
        )
        ,nrow = 30
        ,ncol = 5
        ,byrow = TRUE
    )

}

#PH2.get.k.EPC(p = 0.05, c.ii = 3, eps = 0.2, ARL0 = 370, var.option = TRUE, interval = c(100, 400), type = 'ANOVA-based', subgroup.size = 10, u = u, v = v)


#
#m <- dim(x)[1]
#n <- dim(x)[2]
#
#x.bar.bar <- mean(x)
#
#sigma.v <- sqrt(sum(apply(x, 1, var)) / m) / c4.f(m * (n - 1))
#
#k <- get.cc(m, m * (n - 1), 0.05)$K
#
#LCL <- x.bar.bar - k * sigma.v / sqrt(n)
#UCL <- x.bar.bar + k * sigma.v / sqrt(n)
#
#
#x.bar <- rowMeans(x)
#
#plot(c(1, m), c(LCL * 49997/50000, UCL * 50003/50000), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Means', type = 'n')
#
#axis(side = 1, at = 1:25)
#
#points(1:m, x.bar, type = 'o')
#points(c(-1, m + 2), c(LCL, LCL), type = 'l', col = 'red')
#points(c(-1, m + 2), c(UCL, UCL), type = 'l', col = 'red')
#points(c(-1, m + 2), c(x.bar.bar, x.bar.bar), type = 'l', col = 'blue')
#text(20, UCL * 50001/50000, paste('UCL = ', round(UCL, 4)))
#text(20, LCL * 49999/50000, paste('LCL = ', round(LCL, 4)))
#
#pos <- c(rep(c(3, 1), 5), 1, 3, 3, 1, 3, 1, 1, 3, 1, 3, 1, 3, 1, 3, 1) #rep(c(1, 3), 8))[-26]
#
#text(1:m, x.bar, round(x.bar, 4), cex = 0.7, pos = pos)
#
#
####################################################################################################################################################
    #Example about how to get K and L by FAP0
####################################################################################################################################################
#get.cc(
#            20
#            ,80
#            ,FAP = 0.1
#            ,off.diag = NULL
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'direct'
#)
#get.cc(
#            20
#            ,80
#            ,FAP = 0.1
#            ,off.diag = NULL
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'indirect'
#)

####################################################################################################################################################
    #Example about how to get FAP by K
####################################################################################################################################################
#
#k <- c(
#22.84794
#,10.52163
#,7.386245
#,6.057107
#,5.340985
#,4.898076
#,4.598652
#,4.383316
#,3.842428
#,3.619871
#,3.498807
#,3.422752
#,3.280788
#,3.182391
#,3.151007
#,3.135567
#,3.126384
#,3.120294
#
#
#
#
#
#)
#
#m.seq <- c(
#    3
#    ,4
#    ,5
#    ,6
#    ,7
#    ,8
#    ,9
#    ,10
#    ,15
#    ,20
#    ,25
#    ,30
#    ,50
#    ,100
#    ,150
#    ,200
#    ,250
#    ,300
#
#)
#
#
#i <- 0
#
#record <- rep(NA, 18)
#
#for (m in m.seq){
#
#    i <- i + 1
#
#    nu <- m - 1
#
#    record[i] <- get.FAP0(k[i], m, nu)
#
#}
#
#as.matrix(record)




