#' Run the BPLS regression model
#'
#' Posterior inference of the Bayesian partial least squares regression model using a Gibbs sampler. There are three types of models available depending on the assumed prior structure on the model parameters (see details).
#' 
#' @param X Matrix of predictor variables.
#' @param Y Vector or matrix of responses.
#' @param Xtest Matrix of predictor variables to predict for. 
#' @param Prior 	List of hyperparameters specifying the parameter prior distributions. If left \code{NULL}, a generic set of priors will be generated.
#' @param Qs Upper limit on the number of latent components. If \code{NULL} it is chosen automatically.
#' @param N_MCMC 	Number of iterations to run the Markov chain Monte Carlo algorithm.
#' @param BURN 	Number of iteration to be discarded as the burn-in.
#' @param Thin Thinning procedure for the MArkov chain. \code{Thin = 1} results in no thinning. Only use for long chains to reduce memory.
#' @param model.type Type of BPLS model to use; one of \code{standard}, \code{ss} (spike-and-slab), or \code{LASSO} (see details).
#' @param scale. Logical; if \code{TRUE} then the data variables will be scale to have unit variance.
#' @param center. Logical; if \code{TRUE} then the data variables will be zero-centred.
#' @param PredInterval   Coverage of prediction intervals if \code{Xtest} is provided; 0.95 by default. 
#' @details The number of latent variables is inferred using the multiplicative gamma process prior (Bhattacharya and Dunson, 2011).
#' Posterior samples from the fitted model are stored as a list.
#' There are three types of parameter prior structures resulting in three different model types:
#' * **BPLS**: No additional structure assumed; set \code{model.type=standard}. This model mimics the standard partial least squares regression (PLS; Wold, 1973).
#' * **ss-BPLS**: A spike-and-slab variant introducing additonal column-wise sparsity to the loading matrix relating to the response variables \code{Y}; set \code{model.type=ss}. This approach mimics the Two-way Orthogonal PLS regression (O2PLS; Trygg and Wold, 2003).
#' * **L-BPLS**: A LASSO variant introducing additonal element-wise sparsity to the loading matrix relating to the response variables \code{Y}; set \code{model.type=LASSO}. This approach mimics the sparse PLS regression (sPLS; Chun and Keles, 2010).
#'
#' Empirical comparisons in Urbas et al. (2024) suggest that the LASSO variant is the best at point predictions and prediction interval coverage when applied to spectral data.
#' @references Bhattacharya, A. and Dunson, D. B. (2011) Sparse Bayesian infinite factor models, \emph{Biometrika}, 98(2): 291–306
#'
#' Chun, H. and Keles, S. (2010). Sparse partial least squares regression for simultaneous dimension reduction and variable selection. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 72(1):3–25.
#'
#' Trygg, J. and Wold, S. (2003). O2-PLS, a two-block (X–Y) latent variable regression (LVR) method with an integral OSC filter. \emph{Journal of Chemometrics}, 17(1):53–64.
#'
#' Urbas, S., Lovera, P., Daly, R., O'Riordan, A., Berry, D., and Gormley, I. C. (2024). "Predicting milk traits from spectral data using Bayesian probabilistic partial least squares regression." \emph{The Annals of Applied Statistics}, 18(4): 3486-3506. <\doi{10.1214/24-AOAS1947}>
#'
#' Wold, H. (1973). Nonlinear iterative partial least squares (NIPALS) modelling: some current developments. In \emph{Multivariate analysis–III}, pages 383–407. Elsevier.
#' @return A list of:
#' \item{\code{chain}}{A Markov chain of samples from the parameter posterior.}
#' \item{\code{X}}{Original set of predictor variables.}
#' \item{\code{Y}}{Original set of response variables.}
#' \item{\code{Xtest}}{Original set of predictor variables to predict from; if \code{Xtest} is provided.}
#' \item{\code{Ytest}}{Point predictions for new responses; if \code{Xtest} is provided.}
#' \item{\code{Ytest_PI}}{Prediction intervals for new responses (by default 0.95 coverage); if \code{Xtest} is provided.}
#' \item{\code{Ytest_dist}}{Posterior predictive distributions for new responses; if \code{Xtest} is provided.}
#' \item{\code{diag}}{Additional diagnostics for assessing chain convergence.}
#' @export
#' @examples
#' \donttest{
#' # data(milk_MIR)
#' X = milk_MIR$xMIR
#' Y = milk_MIR$yTraits[, c('Casein_content','Fat_content')]
#' 
#' set.seed(1)
#' # fit model to 25% of data and predict on remaining 75%
#' idx = sample(seq(nrow(X)),floor(nrow(X)*0.25),replace = FALSE)
#'
#' Xtrain = X[idx,];Ytrain = Y[idx,]
#' Xtest = X[-idx,];Ytest = Y[-idx,]
#'
#' # fit the model (default MCMC settings can take longer)
#' bplsr_Fit = bplsr(Xtrain,Ytrain)
#'
#' # generate predictions
#' bplsr_pred = bplsr.predict(model = bplsr_Fit, newdata = Xtest)
#'
#' # point predictions
#' head(bplsr_pred$Ytest)
#'
#' # lower and upper limits of prediction interval
#' head(bplsr_pred$Ytest_PI)
#'
#' # plot of predictive posterior distribution for single test sample
#' hist(bplsr_pred$Ytest_dist[1,'Casein_content',], freq = FALSE,
#'      main = 'Posterior predictive density', xlab = 'Casein_content')}
bplsr = function(X,Y, Xtest = NULL, Prior = NULL, Qs = NULL, N_MCMC = 2e4,
						 BURN = ceiling(0.3*N_MCMC), Thin = 1, model.type = 'standard',
						 scale. = TRUE, center. = TRUE, PredInterval = 0.95){

	if(any(is.na(X))){
		stop('The X matrix contains missing values.')
	}
	if(any(is.na(Y))){
		stop('The Y matrix/vector contains missing values.')
	}
	if(is.null(Prior)){
		Prior = list(Asig= 2.5, Bsig=0.1, Apsi = 2.5, Bpsi = 1.5, Alam = 1.5, Blam = 1.5,
				nuW = c(2,3), nuC = c(2,3), alpha = c(2.2,2.2), beta = c(1,1))
	}
	X = as.matrix(X)
	Y = as.matrix(Y)

	muX = colMeans(X)
	muY = colMeans(Y)
	sdX = apply(X,2,sd)
	sdY = apply(Y,2,sd)
	namesX = colnames(X)
	namesY = colnames(Y)
	standards = list(muX = muX, sdX = sdX, muY = muY, sdY = sdY,
						namesX = namesX,namesY = namesY)
	X. = scale(X, center = muX, scale = sdX)
	Y. = scale(Y, center = muY, scale = sdY)

	
	if(model.type == 'ss'){
		Out = bplsrMCMC(X., Y. , Xtest, Prior, N_MCMC, BURN, Thin,
			model.type,mcmc.kernel = ssBPLS_mcmc_kernel,standards,PredInterval)
	} else if(model.type == 'LASSO'){
		Out = bplsrMCMC(X., Y. , Xtest, Prior, N_MCMC, BURN, Thin,
			model.type,mcmc.kernel = LBPLS_mcmc_kernel,standards,PredInterval)
	} else if(model.type == 'standard'){
		Out = bplsrMCMC(X., Y. , Xtest, Prior, N_MCMC, BURN, Thin,
			model.type,mcmc.kernel = BPLS_mcmc_kernel,standards,PredInterval)
	}
	Out$Xtrain = X
	Out$Ytrain = Y
	Out$standards = standards
	return(Out)
	

}

#' Predict from a fitted BPLS regression model
#' 
#' Generates predictions from the fitted BPLS regression model using Monte Carlo simulation.
#'
#' @param model Output of \code{bplsr}.
#' @param newdata Matrix of predictor variables to predict for. 
#' @param PredInterval Intended coverage of prediction intervals (between 0 and 1). Setting the value to 0 only produces point predictions without prediction intervals.
#' @details Predictions of the responses are generated from the posterior predictive distribution, marginalising out the model parameters; see Section 3.5 of Urbas et al. (2024).
#' @references Urbas, S., Lovera, P., Daly, R., O'Riordan, A., Berry, D., and Gormley, I. C. (2024). "Predicting milk traits from spectral data using Bayesian probabilistic partial least squares regression." \emph{The Annals of Applied Statistics}, 18(4): 3486-3506. <\doi{10.1214/24-AOAS1947}>
#' @return A list of:
#' \item{\code{Ytest}}{Point predictions for new responses; if \code{Xtest} is provided.}
#' \item{\code{Ytest_PI}}{Prediction intervals for new responses (by default 0.95 coverage); if \code{Xtest} is provided.}
#' \item{\code{Ytest_dist}}{Posterior predictive distributions for new responses; if \code{Xtest} is provided.}
#' @export
#' @examples
#' \donttest{
#' # data(milk_MIR)
#' X = milk_MIR$xMIR
#' Y = milk_MIR$yTraits[, c('Casein_content','Fat_content')]
#' 
#' set.seed(1)
#' # fit model to 25% of data and predict on remaining 75%
#' idx = sample(seq(nrow(X)),floor(nrow(X)*0.25),replace = FALSE)
#'
#' Xtrain = X[idx,];Ytrain = Y[idx,]
#' Xtest = X[-idx,];Ytest = Y[-idx,]
#'
#' # fit the model (for default MCMC settings leave Qs and N_MCMC blank; can take longer)
#' bplsr_Fit = bplsr(Xtrain,Ytrain, Qs = 10, N_MCMC = 5000)
#'
#' # generate predictions
#' bplsr_pred = bplsr.predict(model = bplsr_Fit, newdata = Xtest)
#'
#' # point predictions
#' head(bplsr_pred$Ytest)
#'
#' # lower and upper limits of prediction interval
#' head(bplsr_pred$Ytest_PI)
#'
#' # plot of predictive posterior distribution for single test sample
#' hist(bplsr_pred$Ytest_dist[1,'Casein_content',], freq = FALSE,
#'      main = 'Posterior predictive density', xlab = 'Casein_content')}
bplsr.predict = function(model, newdata, PredInterval = 0.95){
	Xtest = as.matrix(newdata)
	Xtest. = scale(Xtest, center=model$standards$muX,scale = model$standards$sdX)
	R = nrow(model$chain[[1]]$C)
	EYtmp = matrix(0, nrow = nrow(Xtest.), ncol = R)
	storeY = array(NA, dim = c(nrow(Xtest.),R,length(model$chain)))


	# storeY = rep(list(EYtmp), length(StoreIdx))
	pb <- progress::progress_bar$new(
		format = "Prediction:  [:bar] :percent .. :eta",total = length(model$chain),
		clear = FALSE, width= 40)
	for(iter in seq_along(model$chain)){
		tmp = SampleYpred(Xtest.,model$chain[[iter]])
		EYtmp = EYtmp + tmp$EYpred
		storeY[,,iter] = InvScale(tmp$Ypred,
				 center. = model$standards$muY, scale. = model$standards$sdY)
		pb$tick()
	}
	dimnames(storeY)[[2]] = model$standards$namesY
	EYtest = InvScale(EYtmp/length(model$chain),
				 center. = model$standards$muY, scale. = model$standards$sdY)

	if(PredInterval>0){
		half_alpha = 0.5*(1-PredInterval)
		Ytest_PI = array(NA, dim = c(nrow(Xtest.),R,2))
		for(n in seq(nrow(Xtest.))){
			Ytest_PI[n,,] = t(apply(storeY[n,,],1,quantile, probs = c(half_alpha,1-half_alpha)))
		}

		dimnames(Ytest_PI)[[2]] = model$standards$namesY
		dimnames(Ytest_PI)[[3]] = c(paste0(half_alpha*100,"%"),paste0((1-half_alpha)*100,"%"))
	} else{
		Ytest_PI = NULL
	}
	return(list(Ytest = EYtest[,1:ncol(EYtest)], Ytest_PI = Ytest_PI,Ytest_dist = storeY))
}

bplsr.getESS = function(EYchain){

	if(is.list(EYchain)){
		EYchain = lapply(EYchain, as.vector)
		EY = Reduce(rbind,EYchain)
		ess = coda::effectiveSize(EY)
		return(ess)
	} else {
		ess = coda::effectiveSize(EYchain)
		return(ess)
	}
	
}


bplsrMCMC = function(X,Y, Xtest = NULL, Prior = NULL, N_MCMC = 1e3, BURN = 0.3*N_MCMC,
						Thin = 1, model.type,
						mcmc.kernel = BPLS_mcmc_kernel, standards, PredInterval = 0.95) {	
	N = nrow(X)
	P = ncol(X)
	R = ncol(Y)

	# Initialise storage

	pars = INITbplsr(X,Y, Prior, model.type = model.type)

	pars$model = model.type

	StoreIdx = seq(BURN+1, N_MCMC, Thin)
	Counter = 1

	StorePars = rep(list(ReducePars(pars)),length(StoreIdx))
	pb <- progress::progress_bar$new(format = "MCMC:  [:bar] :percent .. :eta",total = N_MCMC,
										 clear = FALSE, width= 40) 
	for(iter in seq_len(N_MCMC)){
		pars = mcmc.kernel(X, Y, pars, Prior)
		if(iter == StoreIdx[Counter]){
			StorePars[[Counter]] = ReducePars(pars)
			Counter = Counter + 1
		}
		pb$tick()
	}
	if(!is.null(Xtest)){
		Xtest. = scale(Xtest,center = standards$muX,scale = standards$sdX)
		EYtmp = matrix(0, nrow = nrow(Xtest.), ncol = ncol(Y))
		# storeY  = rep(list(EYtmp), length(StoreIdx))
		storeY = array(NA, dim = c(nrow(Xtest.),ncol(Y),length(StoreIdx)))
		storeEY = matrix(NA,nrow = length(StoreIdx),ncol = nrow(Xtest.)*ncol(Y))#rep(list(EYtmp), length(StoreIdx))
		pb <- progress::progress_bar$new(
			format = "Diagnostics:  [:bar] :percent .. :eta",total = length(StoreIdx),
			clear = FALSE, width= 40)
		for(iter in seq_along(StorePars)){
			tmp = SampleYpred(Xpred = Xtest.,small_pars = StorePars[[iter]],SampleY = FALSE)
			storeEY[iter,] = standards$muY + standards$sdY*as.vector(tmp$EYpred)
			pb$tick()
		}

		ess = bplsr.getESS(storeEY)

		tmp_list = bplsr.predict(list(chain = StorePars, standards = standards),
						Xtest,PredInterval = PredInterval)

		Ytest = tmp_list$Ytest
		Ytest_PI = tmp_list$Ytest_PI
		storeY = tmp_list$Ytest_dist
	} else{
		Ytest = NULL
		Ytest_PI = NULL
		storeY = NULL
		storeEY = NULL
		ess = NULL
	}
	return(list(chain = StorePars, Ytest = Ytest, Ytest_PI = Ytest_PI, Ytest_dist = storeY,
	 Xtest = Xtest, diag = list(chain = storeEY, ess = ess)))
}

BPLS_mcmc_kernel=function(X, Y, pars, Prior){
	N = nrow(X)
	P = ncol(X)
	R = ncol(Y)

	S = solve(pars$I_Qs + crossprod(pars$W, pars$W /pars$Sig2)+crossprod(pars$C, pars$C /pars$Psi2))

	pars$Z = (X%*%( pars$W /pars$Sig2) + Y%*%(pars$C /pars$Psi2))%*%S + matrix(rnorm(N*pars$Qs),nrow = N)%*%chol(S)

	# W and C updates

	ZtZ = crossprod(pars$Z)#t(pars$Z) %*%Z
	for(j in 1:P){
		choltDj = chol(diag(1/pars$DELTA[j,]) + ZtZ / pars$Sig2[j])
		pars$W[j,] = backsolve(choltDj,
			backsolve(choltDj,t(pars$Z)%*%X[,j]/ pars$Sig2[j],k=pars$Qs,transpose=1)+rnorm(pars$Qs),
			k=pars$Qs)

	}
	for(j in 1:R){
		choltOj = chol(diag(1/pars$OMEGA[j,]) +ZtZ / pars$Psi2[j])
		pars$C[j,] = backsolve(choltOj,
			backsolve(choltOj,t(pars$Z)%*%Y[,j]/ pars$Psi2[j],transpose=1,k=pars$Qs)+rnorm(pars$Qs),
			k=pars$Qs)
	}

	# Uniqueness updates
	if(pars$isoX){
	    pars$Sig2=rep(1/rgamma(1,Prior$Asig+N*P/2,Prior$Bsig + 0.5 *sum(( X - pars$Z %*% t(pars$W))^2)),P)
	} else {
		pars$Sig2=1/rgamma(P,Prior$Asig+N/2,Prior$Bsig + 0.5 *colSums(( X -pars$Z %*% t(pars$W))^2))
	}
	if(pars$isoY){
		pars$Psi2=rep(1/rgamma(1,Prior$Apsi+N*R/2,Prior$Bpsi + 0.5 *sum(( Y - pars$Z%*% t(pars$C))^2)),R)
	} else {
		pars$Psi2=1/rgamma(R,Prior$Apsi+N/2,Prior$Bpsi + 0.5 *colSums(( Y - pars$Z %*%t(pars$C))^2))
	}

	# Shrinkage

	pars$phiW = rgamma(P*pars$Qs,Prior$nuW[1]+0.5,1)/(Prior$nuW[2]+0.5*pars$W^2 %*%diag(pars$tauWC))
    pars$phiC = rgamma(R*pars$Qs,Prior$nuC[1]+0.5,1)/(Prior$nuC[2]+0.5*pars$C^2 %*%diag(pars$tauWC))

	sum_phiW_W2 = colSums(pars$phiW*pars$W^2)
	sum_phiC_C2 = colSums(pars$phiC*pars$C^2 )
	SumBoth = sum_phiW_W2 + sum_phiC_C2

	# deltaWC[1] = rgamma(1,alphaW[1] + 0.5*P*Qs,betaW[1]+0.5*sum(tauWC/deltaWC[1]*SumBoth))
	for(j in 2:pars$Qs){
		pars$deltaWC[j] = rgamma(1,Prior$alpha[2]+0.5*(P+R)*(pars$Qs-j+1),
		                      Prior$beta[2]+ 0.5*sum( pars$tauWC[j:pars$Qs]/pars$deltaWC[j] * SumBoth[j:pars$Qs]))  
		pars$tauWC = cumprod(pars$deltaWC)
	}

    pars$DELTA = 1/ (pars$phiW %*%diag(pars$tauWC))
    pars$OMEGA = 1/ (pars$phiC %*%diag(pars$tauWC))

    return(pars)
}
ssBPLS_mcmc_kernel=function(X, Y, pars, Prior){
	N = nrow(X)
	P = ncol(X)
	R = ncol(Y)

	S = solve(pars$I_Qs + crossprod(pars$W, pars$W /pars$Sig2)+crossprod(pars$CB, pars$CB /pars$Psi2))

	pars$Z = (X%*%( pars$W /pars$Sig2) + Y%*%(pars$CB /pars$Psi2))%*%S + matrix(rnorm(N*pars$Qs),nrow = N)%*%chol(S)

	# W and C updates

	ZtZ = crossprod(pars$Z)#t(pars$Z) %*%Z
	for(j in 1:P){
		choltDj = chol(diag(1/pars$DELTA[j,]) + ZtZ / pars$Sig2[j])
		pars$W[j,] = backsolve(choltDj,
			backsolve(choltDj,t(pars$Z)%*%X[,j]/ pars$Sig2[j],k=pars$Qs,transpose=1)+rnorm(pars$Qs),
			k=pars$Qs)

	}

	BZt = t(pars$Z) * pars$Bvec#Bcurr%*%t(Zcurr)
	for(j in 1:R){
		choltOj = chol(diag(1/pars$OMEGA[j,]) +ZtZ / pars$Psi2[j])
		pars$C[j,] = backsolve(choltOj,
			backsolve(choltOj,BZt%*%Y[,j]/ pars$Psi2[j],transpose=1,k=pars$Qs)+rnorm(pars$Qs),
			k=pars$Qs)
	}
	pars$CB = makeCB(pars)

	# MH for SS
	LogInvOdds = log(1/pars$probSS - 1)
    for(j in 1:pars$Qs){
		gam = 1-2*pars$Bvec[j]

		tSUM_RESID = t(Y - tcrossprod(pars$Z,pars$CB))# Z%*%t(pars$CB))

		# delta = t(C[,j])%*%PSIinv%*% (C[,j]*ZtZ[j,j] - 2*gam*tSUM_RESID%*%Zcurr[,j])
		delta = crossprod(pars$C[,j] /pars$Psi2, pars$C[,j]*ZtZ[j,j] - 2*gam*tSUM_RESID%*%pars$Z[,j])


		# 1/(1+exp(0.5*delta))
		# exp(-0.5*delta)

		if(runif(1)<as.numeric(exp(-0.5*delta - gam*LogInvOdds))){
			pars$Bvec[j] = pars$Bvec[j]+gam
			if(pars$Bvec[j]==1){
				pars$CB[,j] = pars$C[,j]
			} else {
				pars$CB[,j]
			}
		}
    }
    # update LogInvOdds
    SSsucc = sum(pars$Bvec)
    pars$probSS = rbeta(1, 1+SSsucc, 1+pars$Qs-SSsucc)

    

	# Uniqueness updates
	if(pars$isoX){
	    pars$Sig2=rep(1/rgamma(1,Prior$Asig+N*P/2,Prior$Bsig + 0.5 *sum(( X - pars$Z %*% t(pars$W))^2)),P)
	} else {
		pars$Sig2=1/rgamma(P,Prior$Asig+N/2,Prior$Bsig + 0.5 *colSums(( X -pars$Z %*% t(pars$W))^2))
	}
	if(pars$isoY){
		pars$Psi2=rep(1/rgamma(1,Prior$Apsi+N*R/2,Prior$Bpsi + 0.5 *sum(( Y - pars$Z%*% t(pars$CB))^2)),R)
	} else {
		pars$Psi2=1/rgamma(R,Prior$Apsi+N/2,Prior$Bpsi + 0.5 *colSums(( Y - pars$Z %*%t(pars$CB))^2))
	}

	# Shrinkage

	pars$phiW = rgamma(P*pars$Qs,Prior$nuW[1]+0.5,1)/(Prior$nuW[2]+0.5*pars$W^2 %*%diag(pars$tauWC))
    pars$phiC = rgamma(R*pars$Qs,Prior$nuC[1]+0.5,1)/(Prior$nuC[2]+0.5*pars$C^2 %*%diag(pars$tauWC))

	sum_phiW_W2 = colSums(pars$phiW*pars$W^2)
	sum_phiC_C2 = colSums(pars$phiC*pars$C^2 )
	SumBoth = sum_phiW_W2 + sum_phiC_C2

	# deltaWC[1] = rgamma(1,alphaW[1] + 0.5*P*Qs,betaW[1]+0.5*sum(tauWC/deltaWC[1]*SumBoth))
	for(j in 2:pars$Qs){
		pars$deltaWC[j] = rgamma(1,Prior$alpha[2]+0.5*(P+R)*(pars$Qs-j+1),
		                      Prior$beta[2]+ 0.5*sum( pars$tauWC[j:pars$Qs]/pars$deltaWC[j] * SumBoth[j:pars$Qs]))  
		pars$tauWC = cumprod(pars$deltaWC)
	}

    pars$DELTA = 1/ (pars$phiW %*%diag(pars$tauWC))
    pars$OMEGA = 1/ (pars$phiC %*%diag(pars$tauWC))

    return(pars)
}

LBPLS_mcmc_kernel=function(X, Y, pars, Prior){
	N = nrow(X)
	P = ncol(X)
	R = ncol(Y)

	S = solve(pars$I_Qs + crossprod(pars$W, pars$W /pars$Sig2)+crossprod(pars$C, pars$C /pars$Psi2))

	pars$Z = (X%*%( pars$W /pars$Sig2) + Y%*%(pars$C /pars$Psi2))%*%S + matrix(rnorm(N*pars$Qs),nrow = N)%*%chol(S)

	# W and C updates

	ZtZ = crossprod(pars$Z)#t(pars$Z) %*%Z
	for(j in 1:P){
		choltDj = chol(diag(1/pars$DELTA[j,]) + ZtZ / pars$Sig2[j])
		pars$W[j,] = backsolve(choltDj,
			backsolve(choltDj,t(pars$Z)%*%X[,j]/ pars$Sig2[j],k=pars$Qs,transpose=1)+rnorm(pars$Qs),
			k=pars$Qs)

	}
	for(j in 1:R){
		choltOj = chol(diag(1/pars$OMEGA[j,]) +ZtZ / pars$Psi2[j])
		pars$C[j,] = backsolve(choltOj,
			backsolve(choltOj,t(pars$Z)%*%Y[,j]/ pars$Psi2[j],transpose=1,k=pars$Qs)+rnorm(pars$Qs),
			k=pars$Qs)
	}

	# Uniqueness updates
	if(pars$isoX){
	    pars$Sig2=rep(1/rgamma(1,Prior$Asig+N*P/2,Prior$Bsig + 0.5 *sum(( X - pars$Z %*% t(pars$W))^2)),P)
	} else {
		pars$Sig2=1/rgamma(P,Prior$Asig+N/2,Prior$Bsig + 0.5 *colSums(( X -pars$Z %*% t(pars$W))^2))
	}
	if(pars$isoY){
		pars$Psi2=rep(1/rgamma(1,Prior$Apsi+N*R/2,Prior$Bpsi + 0.5 *sum(( Y - pars$Z%*% t(pars$C))^2)),R)
	} else {
		pars$Psi2=1/rgamma(R,Prior$Apsi+N/2,Prior$Bpsi + 0.5 *colSums(( Y - pars$Z %*%t(pars$C))^2))
	}

	# Shrinkage

	pars$phiW = rgamma(P*pars$Qs,Prior$nuW[1]+0.5,1)/(Prior$nuW[2]+0.5*pars$W^2 %*%diag(pars$tauWC))

	meanPhiC = sqrt(pars$L2_lasso) / ( abs(pars$C)%*%diag(sqrt(pars$tauWC)) )
	pars$phiC = matrix(statmod::rinvgauss(R*pars$Qs,mean = as.vector(meanPhiC), shape  = pars$L2_lasso),
		nrow = R,ncol = pars$Qs)
	L2_lasso = rgamma(1,Prior$Alam + pars$Qs*R,Prior$Blam + 0.5*sum(1/pars$phiC))


	sum_phiW_W2 = colSums(pars$phiW*pars$W^2)
	sum_phiC_C2 = colSums(pars$phiC*pars$C^2 )
	SumBoth = sum_phiW_W2 + sum_phiC_C2

	# deltaWC[1] = rgamma(1,alphaW[1] + 0.5*P*Qs,betaW[1]+0.5*sum(tauWC/deltaWC[1]*SumBoth))
	for(j in 2:pars$Qs){
		pars$deltaWC[j] = rgamma(1,Prior$alpha[2]+0.5*(P+R)*(pars$Qs-j+1),
		                      Prior$beta[2]+ 0.5*sum( pars$tauWC[j:pars$Qs]/pars$deltaWC[j] * SumBoth[j:pars$Qs]))  
		pars$tauWC = cumprod(pars$deltaWC)
	}

    pars$DELTA = 1/ (pars$phiW %*%diag(pars$tauWC))
    pars$OMEGA = 1/ (pars$phiC %*%diag(pars$tauWC))

    return(pars)
}

ReducePars = function(pars){
	return(list(W = pars$W, C = pars$C, Sig2 = pars$Sig2, Psi2 = pars$Psi2))
}

getCoeffMatrix = function(small_pars){
	S = solve(diag(ncol(small_pars$W)) + crossprod(small_pars$W, small_pars$W /small_pars$Sig2))
	a = S%*% t(small_pars$W/small_pars$Sig2)

	return(list(Ca = small_pars$C %*% a, S = S))
}

SampleYpred = function(Xpred,small_pars, EYpred = NULL, S = NULL,
 SampleY = TRUE){
	if(is.null(EYpred)){
		tmp = getCoeffMatrix(small_pars)
		CoeffMat = tmp$Ca
		EYpred =  tcrossprod(Xpred,CoeffMat)
		S = tmp$S
	}
	
	if(SampleY){
		R = ncol(EYpred); N = nrow(EYpred)
		if(R == 1){
			cholVarY = sqrt(tcrossprod(small_pars$C %*% S, small_pars$C) + small_pars$Psi2)
			return(list(EYpred = EYpred,
					Ypred = EYpred +rnorm(N) * as.numeric(cholVarY)))
		} else{
			cholVarY = chol(tcrossprod(small_pars$C %*% S, small_pars$C) + diag(small_pars$Psi2))
		}
		return(list(EYpred = EYpred,
					Ypred = EYpred +matrix(rnorm(N*R),nrow = N) %*% cholVarY))
		
	} else {
		return(list(EYpred = EYpred))
	}
	
}

makeCB = function(pars){
	CB = pars$C
	CB[,pars$Bvec == 0] = 0
	return(CB)
}

bplsrSEM = function(x) {

}

INITbplsr = function(X, Y, Prior,Qs = 0, isoX = FALSE, isoY = FALSE,
		 model.type){
	P = ncol(X)
	R = ncol(Y)
	pcaX = prcomp(X)

	cx <- cumsum(pcaX$sdev)
	cx = cx /tail(cx,1)

	if(Qs==0){
      Qs = min(which(cx>=0.95))
    }
    W = pcaX$rotation[,1:Qs]
    Sig2 = rep(mean(pcaX$sdev[-(1:Qs)])^2,P)
    # SIGMAinv = diag(1/Sig2)
    # Z1 =  X %*%t(solve(diag(Qs) + t(W)%*%SIGMAinv %*%W)%*% t(W)%*%SIGMAinv) 
    Z =  X %*%t(solve(diag(Qs) + crossprod(W, W /Sig2))%*% t(W/Sig2))
    C = t(solve(crossprod(Z))%*%t(Z)%*%Y)
    Psi2 = apply(Y - Z %*% t(C),2,var)

    if(model.type == 'LASSO'){

    	
	    phiW = rgamma(P*Qs,Prior$nuW[1],Prior$nuW[2]);
	    # phiC = rgamma(R*Qs,Prior$nuC[1],Prior$n uC[2])
	    L2_lasso = Prior$Alam/Prior$Blam#rgamma(1,Prior$Alam,Prior$Blam)
	    phiC = 1/rexp(P*Qs,L2_lasso/2)

	    deltaWC = c(1.0, rep(Prior$alpha[2]/Prior$beta[2],Qs-1))
	    tauWC = cumprod(deltaWC)

	    DELTA = matrix(1/(phiW*tauWC), ncol = Qs, byrow = T)
	    OMEGA = matrix(1/(phiC*tauWC), ncol = Qs, byrow = T)

	    return(list(W=W, C=C, Z=Z, Sig2 = Sig2, Psi2 = Psi2,
    			phiW = phiW,phiC = phiC, L2_lasso = L2_lasso,
    			deltaWC = deltaWC, tauWC = tauWC,
    			DELTA = DELTA, OMEGA = OMEGA,
    			isoX = isoX, isoY = isoY, Qs = Qs, I_Qs = diag(Qs)))

    }

    if(model.type == 'ss'){

    	Bvec = rep(1, Qs)
    	CB = C
    	# CB[,Bvec==0] = 0

    	phiW= rgamma(P*Qs,Prior$nuW[1],Prior$nuW[2]); 
    	phiC = rgamma(R*Qs,Prior$nuC[1],Prior$nuC[2])


	    deltaWC = c(1.0, rep(Prior$alpha[2]/Prior$beta[2],Qs-1))
	    tauWC = cumprod(deltaWC)

	    DELTA = matrix(1/(phiW*tauWC), ncol = Qs, byrow = T)
	    OMEGA = matrix(1/(phiC*tauWC), ncol = Qs, byrow = T)


	    return(list(W=W, C=C, Bvec = Bvec, CB = CB, Z=Z, Sig2 = Sig2, Psi2 = Psi2, probSS = 0.5,
	    			phiW = phiW,phiC = phiC,deltaWC = deltaWC, tauWC = tauWC,
	    			DELTA = DELTA, OMEGA = OMEGA,
	    			isoX = isoX, isoY = isoY, Qs = Qs, I_Qs = diag(Qs)))

    } 
    if(model.type == 'standard') {
	    phiW = rgamma(P*Qs,Prior$nuW[1],Prior$nuW[2]) 
	    phiC = rgamma(R*Qs,Prior$nuC[1],Prior$nuC[2])


	    deltaWC = c(1.0, rep(Prior$alpha[2]/Prior$beta[2],Qs-1))
	    tauWC = cumprod(deltaWC)

	    DELTA = matrix(1/(phiW*tauWC), ncol = Qs, byrow = T)
	    OMEGA = matrix(1/(phiC*tauWC), ncol = Qs, byrow = T)


	    return(list(W=W, C=C, Z=Z, Sig2 = Sig2, Psi2 = Psi2,
	    			phiW = phiW,phiC = phiC,deltaWC = deltaWC, tauWC = tauWC,
	    			DELTA = DELTA, OMEGA = OMEGA,
	    			isoX = isoX, isoY = isoY, Qs = Qs, I_Qs = diag(Qs)))

    }




}

InvScale = function(X, center. , scale.){
	Xnew = t(apply(X,1, function(x){x*scale. + center.}))
	colnames(Xnew) = colnames(X)
	return(Xnew)
}

