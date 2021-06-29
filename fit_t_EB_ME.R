#this is a wrapper script designed to fit various EB models
# 'including two-rate models and single-rate models (specify 'regime.map' for two-rate models, see details below)
#'
#' @param phylo: a phylogenetic tree or simmap containing all of the lineages in 'data' vector
#' @param data: a named vector of continuous trait values with names corresponding to phylo$tip.labels. 
#' @param regime.map: a simmap object (constructed on the exact same phylogeny passed to the 'phylo' object) which contains a reconstruction of the different rate classes (e.g., tropical and temperate lineages). Currently only tested for two-regime cases.
#' @param error: a named vector (in the same order as 'data') which specifies the standard error of trait measurements, if measurement error is to be accounted for in fit (specify as NULL if not)
#' @param beta: a vector of starting values for the slope parameter estimation
#' @param sigma: a vector of starting values for the rate parameter estimation
#' @param method: optimisation algorithm (see optim())
#' @param upper: upper bound on optimisation algorithm, if "L-BFGS-B" is chosen as 'method'
#' @param lower: lower bound on optimisation algorithm, if "L-BFGS-B" is chosen as 'method'
#' @param control: further commands passed to optim()
#' @param diagnostic: logical specifying whether diagnostic information should printed to the console
#' @param echo: logical specifying whether model fits should printed in the console

library(phytools)

source('/other_scripts/fit_t_general.R')
source('/other_scripts/generalized_functions.R')

fit_t_EB<-function(phylo,data,regime.map=NULL,error=NULL, beta=NULL, sigma=NULL, method=c("L-BFGS-B","BB","Nelder-Mead"), upper=list(beta=0,sigma=Inf), lower=-Inf, control=list(maxit=20000), diagnostic=FALSE, echo=FALSE){
	
	
if(is.null(regime.map)){ 	# single slope version 
			
		#convert phylo to simmap object
		hold<-rep("A",length(phylo$tip.label))
		hold[1]<-"B"
		names(hold)<-phylo$tip.label
		smap<-make.simmap(phylo,hold,message=F)
		new.maps<-list()
		for(i in 1:length(phylo$edge.length)){
			new.maps[[i]]<-phylo$edge.length[i]
			names(new.maps[[i]])<-"A"
			}
		new.mapped.edge<- as.matrix(rowSums(smap$mapped.edge))
		colnames(new.mapped.edge)<-"A"	
		smap$maps<-new.maps
		smap$mapped.edge<-new.mapped.edge
		
		#create function
		new_list_function<-create.function.list.EB1(smap)

		#fit model
		sigma.constraint<-rep(1, dim(smap$mapped.edge)[2])
		beta.constraint<-rep(1, dim(smap$mapped.edge)[2])
		
		up=c(rep(upper$beta,length(beta.constraint)),Inf,0)
				
		out<-fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=NULL,input.times=NULL,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=up, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	

		
}  else if (!is.null(regime.map)) { # multi-slope version where rates are reset to root state at beginning of new regime

		#first check that regime.map and phylo and data are concordant
		if(!all(as.phylo(phylo)$tip.label == as.phylo(regime.map)$tip.label)) { stop("regime map doesn't match phylogeny")}
		if(length(data) != length(as.phylo(regime.map)$tip.label)) { stop("number of lineages in data and regime map don't match")}
		if(! all (names(data) %in% as.phylo(regime.map)$tip.label)) { stop("names of lineages in data and regime map don't match")}
		if(! all (as.phylo(regime.map)$tip.label %in% names(data)) ) { stop("names of lineages in data and regime map don't match")}
		
		
		new_list_function<-create.function.list.EBmulti(regime.map)
				
		#fit model
		sigma.constraint<-rep(1, dim(regime.map$mapped.edge)[2])
		beta.constraint<-seq(1,by=1,length.out=dim(regime.map$mapped.edge)[2])
		
		up=c(rep(upper$beta,length(beta.constraint)),Inf,0)
		
		out<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=NULL,class.df=NULL,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=up, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		
}    
		return(out)
}