# just to remember what are the different options
# constraint=TRUE :
# option: 1 sigma (shared sigma)
# option: constrain some beta to be the same of constrain one of the beta to be zero
# error: we can estimate unknown variation no explained by the model in top of known error. Just have to change some parts of thefunction

# Fit a model for which rates depends on a time-serie curve with regime specific parameters estimates.

fit_t_general <- function(tree, data, fun, class.df, input.times, error=NULL, beta=NULL, sigma=NULL, model=c("exponential","linear"), method=c("L-BFGS-B","BB"), upper=Inf, lower=-Inf, control=list(maxit=20000), diagnostic=TRUE, echo=TRUE, constraint=NULL, return.tree=FALSE) {
  
  require(mvMORPH)
  if(!inherits(tree,"simmap")==TRUE) stop("For now only simmap-like mapped trees are allowed.","\n")
  old.tree<-tree
  
  # Parameters
  if(!is.null(names(data))) data <- data[tree$tip.label]
  data<-as.matrix(data)
  method=method[1]
  rownames(data)<-tree$tip.label
  model=model[1]
  # Compute node time from the root to the tips
  times<-max(nodeHeights(tree))-nodeHeights(tree)[match(1:tree$Nnode+length(tree$tip),tree$edge[,1]),1]
  
  names(times)<-1:tree$Nnode+length(tree$tip)
  # Set the root to zero
  times<-max(times)-times
  # Max time
  mtot=max(nodeHeights(tree))
  onestate<-ifelse(dim(tree$mapped.edge)[2]==1,TRUE,FALSE) 

  # Number of species
  n=length(tree$tip.label)
  # Number of traits (for future versions)
  k=1
  tree <- reorderSimmap(tree, order="postorder")

    # Number of maps (selective regimes)
    nstates <- dim(old.tree$mapped.edge)[2]
    if(is.null(constraint)){
           number_maps_beta <- number_maps_sigma <- nstates 
           index.user.sigma <- index.user.beta <- 1:number_maps_sigma
        
    }else{ 
        if(is.null(constraint[["sigma"]])) sigConst <- FALSE else sigConst <- TRUE
        if(is.null(constraint[["beta"]])) betConst <- FALSE else betConst <- TRUE
        
        if(sigConst & betConst){ # i) both beta and sigma constrained
            number_maps_beta <- length(unique(constraint$beta[!is.na(constraint$beta)]))
            number_maps_sigma <- length(unique(constraint$sigma))  
            # set the indices
            index.user.sigma <- constraint$sigma
            index.user.beta <- constraint$beta
            
        }else if(!sigConst & betConst){ # ii) only beta constrained
            number_maps_beta <- length(unique(constraint$beta[!is.na(constraint$beta)]))
            number_maps_sigma <- nstates 
            
            index.user.sigma <- 1:number_maps_sigma
            index.user.beta <- constraint$beta
            
        }else if(sigConst & !betConst){ # iii) only sigma constrained
            number_maps_beta <- nstates 
            number_maps_sigma <- length(unique(constraint$sigma))
            
            index.user.sigma <- constraint$sigma
            index.user.beta <- 1:number_maps_beta
        }

    }
    
  
  # Param likelihood contrasts (we are organizing the tree to be postorder for an efficient traversal)
  #ind=reorder(tree,"postorder",index.only=TRUE)
  phy=tree
  #phy$edge.length<-phy$edge.length[ind]
  #phy$edge<-phy$edge[ind,]
  
  # check for simmap like format
  #if(inherits(tree,"simmap")){
  #  phy$mapped.edge<-phy$mapped.edge[ind,]
  #  phy$maps<-phy$maps[ind]
  #}
  #phy <- reorderSimmap(tree, order="pruningwise")
  
  # Random starting value if not provided
  if(is.null(beta)){
    beta=rep(0, number_maps_beta)
  }
  if(is.null(sigma)){
    sigma=rep(sum(pic(data,old.tree)^2)/n, number_maps_sigma)
  }
  
  if(model=="linear"){
    startval=c(beta,log(sigma))
    nbeta=length(beta)
    nsigma=length(sigma)
  }else if(model=="exponential"){
    startval=c(beta,log(sigma))
    nbeta=length(beta)
    nsigma=length(sigma)
  }
  
  # Error estimation?
  if(!is.null(error)){
    ## Index error
    index_error<-sapply(1:n, function(x){ which(phy$edge[,2]==x)})
    startval=c(startval,0.001)
    # to construct a mixed model (refer to l.172 of the code below)
    if(is.numeric(error)){
      error_meas = error^2
    }else{
      error_meas = numeric(n)       
    }
  }
  
  ##--------------Fonction-generale-DD-Env-------------------------------------------##
  
  BranchtransformMAPS<-function(phy,beta,mtot,times,fun,sigma=NULL,model=NULL,errorValue=NULL){
    #Transformations
    tips <- length(phy$tip.label)
    res <- phy
    
    if(model=="exponential"){  
      # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
      f<-function(x, sigma, beta, funInd){sigma*exp(beta*fun[[funInd]](x,df=class.df,times=input.times))}
      
    }else if(model=="linear"){
      # Clim-lin function
      f<-function(x, sigma, beta, funInd){sigma+beta*fun[[funInd]](x,df=class.df,times=input.times)}
     
    }
    
    # Loops over the edges
    for (i in 1:length(phy$edge.length)) {
      
      age <- times[phy$edge[i, 1] - tips] # retrieve the age at the node
      currentmap<-phy$maps[[i]]           # retrieve the corresponding maps on the edge "i"
      indlength<-length(currentmap)       # How many mapping there are?
      tempedge<-numeric(indlength)        # temporary vector for the mappings
      
      # loop pour traverser les "maps"
      for(betaval in 1:indlength){
          
        if(onestate){
            regimenumber=1
        }else{
            regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[betaval])  # retrieve the regimes within maps
        }

        bet<-beta[regimenumber]           # select the corresponding parameter for beta
        sig<-sigma[regimenumber]          # select the corresponding parameter for sigmz
        bl<-currentmap[[betaval]]         # branch length under the current map
        
        int <- try(integrate(f, lower=age, upper=(age + bl), subdivisions=500, rel.tol = .Machine$double.eps^0.05, sigma=sig, beta=bet, funInd=regimenumber), silent=TRUE)
          
           if(inherits(int ,'try-error')){
             warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
             tempbeta <- NA_real_
           } else {
             tempbeta <- int$value
           }
          
        tempedge[betaval] <- tempbeta 
        # on met à jour age parcequ'on va passer au maps suivant pour la lignée i
        # update "age" because we're moving to the next map for lineage i.
        age<-age+bl
      }
      # update branch length
      res$edge.length[i]<-sum(tempedge)
    }
    
    phy<-res
    
    if(!is.null(errorValue)){
      phy$edge.length[index_error]<-phy$edge.length[index_error] + error_meas + errorValue^2
    }
    
    return(phy)
  }
  
  
  ##---------------------------------------------------------------------------------##
  clikCLIM <- function( param, dat, phylo, mtot, times, fun=fun, model, results=FALSE, return.tree=FALSE) {
    
    if(model=="exponential"){
       
        # create vector of parameters
        beta <- numeric(nstates)
        sigma <- numeric(nstates)
        
        # assign values
        beta[] <- c(param[seq_len(nbeta)])[index.user.beta]
        sigma[] <- c(exp(param[nbeta+seq_len(nsigma)]))[index.user.sigma]
        
        # constrain some values to zero (should be only for beta as sigma=0 is undefined)
        beta[is.na(beta)] <- 0
        
      if(!is.null(error)) errorValue <- param[nbeta+nsigma+1] else errorValue <- NULL
      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      if(return.tree){ return(phylo)}
      if(any(is.na(phylo$edge.length)))  return(1000000)
	  if(any(phylo$edge.length<0)) return(1000000)
      LL<-mvLL(phylo,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))
      
    }else{

        # create vector of parameters
        beta <- numeric(nstates)
        sigma <- numeric(nstates)
        
        # assign values
        beta[] <- c(param[seq_len(nbeta)])[index.user.beta]
        sigma[] <- c(exp(param[nbeta+seq_len(nsigma)]))[index.user.sigma]
        
        # constrain some values to zero (should be only for beta as sigma=0 is undefined)
        beta[is.na(beta)] <- 0
        
      if(!is.null(error)) errorValue <- log(param[nbeta+nsigma+1]) else errorValue <- NULL
      
      # test=sigma+(beta*maxN)
      # if(any(test<=0)){
          #	LL<-list()
        #	LL$logl<-Inf
        #	}else{

      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      if(return.tree){ return(phylo)}
      if(any(is.na(phylo$edge.length))) return(1000000) # instead of checking the parameter values as done previously, I return a high-likelihood value when there are NAs in the branch lengths. Note also that returning Inf value doesn't work with L-BFGS-B algorithm
	  if(any(phylo$edge.length<0)) return(1000000)
      LL<-mvLL(phylo,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))
      
    }
    if(is.na(LL$logl) | is.infinite(LL$logl)){return(1000000)}
    if(results==FALSE){
      return(-LL$logl)
    }else{
      return(list(LL=-LL$logl, mu=LL$theta, s2=sigma))
    }
  }
  
  ##------------------------------------Optimization-------------------------------##
  phyloTrans=NULL
  if(method=="BB"){
    require(BB)
    estim<-spg(par=startval,fn=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)},control=control ,method=3, lower=lower, upper=upper)
  }else if(method=="L-BFGS-B" | method=="Nelder-Mead"){
    estim<-optim(par=startval,fn=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)},control=control, hessian=TRUE, method=method, lower=lower, upper=upper)
    
  }else if(method=="fixed"){
    estim <- list()
    estim$par <- param <- c(beta,log(sigma))
    #estim$par <- c(beta, sigma)
    estim$value <- clikCLIM(param=estim$par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)
    estim$convergence <- 0
    
    ## Return the tree --- just some modifications to the previous code to allow retrieving the tree
    # create vector of parameters
    beta <- numeric(nstates)
    sigma <- numeric(nstates)
    # assign values
    beta[] <- c(param[seq_len(nbeta)])[index.user.beta]
    sigma[] <- c(exp(param[nbeta+seq_len(nsigma)]))[index.user.sigma]
    # constrain some values to zero (should be only for beta as sigma=0 is undefined)
    beta[is.na(beta)] <- 0
    phyloTrans <- BranchtransformMAPS(phy, beta, mtot, times, fun, sigma, model, errorValue=NULL)
  }
  
  # Results
  # Prepar the tables for the results
  beta <- numeric(nstates)
  sigma <- numeric(nstates)
        
  if(model=="exponential"){
      
     # assign values
     beta[] <- c(estim$par[seq_len(nbeta)])[index.user.beta]
     sigma[] <- c(exp(estim$par[nbeta+seq_len(nsigma)]))[index.user.sigma]
     # constrain some values to zero (should be only for beta as sigma=0 is undefined)
     beta[is.na(beta)] <- 0
      
    resultList<-matrix(c(beta, sigma), ncol=nstates, byrow=T)
    colnames(resultList)<-c(colnames(tree$mapped.edge))
    rownames(resultList)<-c("beta","sigma")
      
  }else{
      
     # assign values
     beta[] <- c(estim$par[seq_len(nbeta)])[index.user.beta]
     sigma[] <- c(estim$par[nbeta+seq_len(nsigma)])[index.user.sigma]
     
     # constrain some values to zero (should be only for beta as sigma=0 is undefined)
     beta[is.na(beta)] <- 0
      
    resultList<-matrix(c(beta, exp(sigma)), ncol=nstates, byrow=T)
    colnames(resultList)<-c(colnames(tree$mapped.edge))
    rownames(resultList)<-c("beta","sigma")
  }
  
  if(!is.null(error)) errorValue <- estim$par[nsigma+nbeta+1]^2 else errorValue <- NULL
  
  # LogLikelihood
  LL<--estim$value
  # parameter (anc + sigma)
  if(model=="exponential"){
    nparam=1+length(estim$par)
  }else{
    nparam=1+length(estim$par) #sigma estimated in optimization
  }
  
  # AIC
  AIC<--2*LL+2*nparam
  # AIC corrected
  AICc<-AIC+((2*nparam*(nparam+1))/(Ntip(phy)-nparam-1)) #Hurvich et Tsai, 1989
  
  #ancestral states estimates
  anc<-clikCLIM(param=estim$par, dat=data, phy, mtot, times, fun=fun, model, results=TRUE)$mu
  
  if(return.tree){
  transformed.phylo.ML<-clikCLIM(param=estim$par, dat=data, phy, mtot, times, fun=fun, model, results=TRUE,return.tree=TRUE)
  }
  
  ##---------------------Diagnostics--------------------------------------------##
  
  if(estim$convergence==0 & diagnostic==TRUE){
    cat("\n","successful convergence of the optimizer","\n")
  }else if(estim$convergence==1 & diagnostic==TRUE){
    cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
  }else if(diagnostic==TRUE){
    cat("\n","convergence of the optimizer has not been reached, try simpler model","\n")
  }
  
  # Hessian eigen decomposition to check the derivatives
  if(method=="BB"){
    require(numDeriv)
    #hmat<-hessian(x=estim$par, func=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)$LL})
    #hess<-eigen(hmat)$value
    hess <- 0
  }else if(method=="L-BFGS-B" | method=="Nelder-Mead"){
    hess<-eigen(estim$hessian)$values
  }else{
    hess<-0
  }
  if(any(hess<0)){
    hess.value<-1
    if(diagnostic==TRUE){
      cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
  }else{
    hess.value<-0
    if(diagnostic==TRUE){
      cat("a reliable solution has been reached","\n")}
  }
  
  ##-------------------Print results--------------------------------------------##
  if(echo==TRUE){
    cat("\n")
    cat("Summary results for the",model," model","\n")
    cat("LogLikelihood:","\t",LL,"\n")
    cat("AIC:","\t",AIC,"\n")
    cat("AICc:","\t",AICc,"\n")
    cat(nparam,"parameters")
    cat("\n")
    cat("Estimated rates matrix","\n")
    print(resultList)
    cat("\n")
    cat("Estimated ancestral state","\n")
    cat(anc)
    cat("\n")
    if(!is.null(error)){
      cat("\n")
      cat("Estimated error","\n")
      cat(errorValue)
      cat("\n") 
    }
  }
  
  if(return.tree){
  results<-list(LogLik=LL, AIC=AIC, AICc=AICc, rates=resultList, anc=anc, convergence=estim$convergence, hess.values=hess.value, error=errorValue, param=estim$par, phyloTrans=transformed.phylo.ML)
  } else{
  results<-list(LogLik=LL, AIC=AIC, AICc=AICc, rates=resultList, anc=anc, convergence=estim$convergence, hess.values=hess.value, error=errorValue, param=estim$par, phyloTrans=phyloTrans)
	}
}


