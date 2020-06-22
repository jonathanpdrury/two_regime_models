#this is a wrapper script designed to fit various competition models
# 'including two-rate models and single-rate models (specify 'regime.map' for two-rate models, see details below)
#' models with biogeography vs. models without biogeography (specify 'geo.map' for models incorporating biogeography, see details below)
#'
#' @param phylo: a phylogenetic containing all of the lineages in 'data' vector
#' @param data: a named vector of continuous trait values with names corresponding to phylo$tip.labels. 
#' @param error: a named vector (in the same order as 'data') which specifies the standard error of trait measurements
#' @param model: "MC" returns matching competition model fit, "DDexp" returns exponential diversity dependent model fit, "DDlin" returns linear diversity dependent model fit
#' @param pars: starting values for maximum likelihood search (optional); specify as c(log(sigma),rate,log(error)) where 'rate' corresponds to the S parameter in the MC model or the slope parameter in the DD model
#' @param geography.object: an ancestral biogeography reconstruction, (matching 'phylo' object) as a geography object (see documentation for CreateGeoObject in the R-package RPANDA)
#' @param regime.map: a simmap object (constructed on the exact same phylogeny passed to the 'phylo' object) which contains a reconstruction of the different rate classes (e.g., tropical and temperate lineages). Currently only tested for two-regime cases.


source('/other_scripts/PhenotypicModel_PLUSME.R')
source('/other_scripts/PhenotypicADiag.R')
source('/other_scripts/DDlinMulti_geo_ADiag_ME.R')
source('/other_scripts/DDlinMulti_nogeo_ADiag_ME.R')
source('/other_scripts/DDexpMulti_nogeo_ADiag_ME.R')
source('/other_scripts/DDexpMulti_geo_ADiag_ME.R')
source('/other_scripts/MC_twoS_PM_geo_ME.R')
source('/other_scripts/MC_twoS_PM_ME.R')
source('/other_scripts/DDexp_geo_ADiag_ME.R')
source('/other_scripts/DDlin_geo_ADiag_ME.R')
source('/other_scripts/MC_geo_PM_ME.R')
source('/other_scripts/DDexp_nogeo_ADiag_ME.R')
source('/other_scripts/MC_nogeo_ADiag_ME.R')
source('/other_scripts/DDlin_nogeo_ADiag_ME.R')
source('/other_scripts/resortGeoObject.R')
source('/other_scripts/resortSMatrix.R')
source('/other_scripts/CreateSMatrix.R')
source('/other_scripts/CreateClassObject.R')
source('/other_scripts/ReconcileGeoObjectSMatrix.R')

fit_t_comp_ME<-function(phylo,data,error, model=c("MC","DDexp","DDlin"),pars=NULL,geography.object=NULL, regime.map=NULL){

#check to make sure data are univariate, with names matching phylo object
if(length(data)!=length(phylo$tip.label)){stop("length of data does not match length of tree")}
if(is.null(names(data))){stop("data missing taxa names")}
if(!is.null(dim(data))){stop("data needs to be a single trait")}
is_tip <- phylo$edge[,2] <= length(phylo$tip.label)
if(sum(diff(phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp_ME cannot be used with ladderized phylogenies')}


if(is.null(geography.object) & is.null(regime.map)){ #single-slope version for sympatric clades


	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_ME(phylo)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*4 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), S = as.numeric(S), nuisance=as.numeric(nuisance),z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_ME(phylo)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), r = as.numeric(r), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_ME(phylo)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), b = as.numeric(b), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object) & is.null(regime.map)){ #single-slope version with biogeography


	if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
	sgeo<-resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in code
	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_geo_ME(phylo,geo.object=sgeo)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc =(2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), S = as.numeric(S), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_geo_ME(phylo,geo.object=sgeo)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), r = as.numeric(r), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_geo_ME(phylo,geo.object=sgeo)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), b = as.numeric(b), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}

}

if(is.null(geography.object) & !is.null(regime.map)){ #multi-slope version for sympatric clades (i.e., no biogeography)
	
	class.object<-try(CreateClassObject(regime.map))
	if(class(dist.class.object)=="try-error"){
		dist.class.object<-CreateClassObject(regime.simmap.trop.trimmed,rnd=6)
		}
				
	SMatrix<-CreateSMatrix(class.object)
	smat<-resortSMatrix(phylo, SMatrix)
	
	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_twoS_ME(phylo,S.object=smat)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S1<-opt$inferredParams[3]
		S2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), S1 = as.numeric(S1), S2 = as.numeric(S2), nuisance=as.numeric(nuisance),z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_multi_ME(phylo,r.object=smat)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r1<-opt$inferredParams[3]
		r2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), r1 = as.numeric(r1), r2 = as.numeric(r2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_multi_ME(phylo,r.object=smat)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b1<-opt$inferredParams[3]
		b2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), b1 = as.numeric(b1),b2 = as.numeric(b2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object) & !is.null(regime.map)){ #multi-slope version with biogeography

	if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
	sgeo0<-resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in code
	
	class.object<-try(CreateClassObject(regime.map))
	if(class(dist.class.object)=="try-error"){
		dist.class.object<-CreateClassObject(regime.simmap.trop.trimmed,rnd=6)
		}
				
	SMatrix<-CreateSMatrix(class.object)
	smat0<-resortSMatrix(phylo, SMatrix)
	
	int<-try(ReconcileGeoObjectSMatrix(geo.object=sgeo0,S.matrix=smat0))	
	
	#some catches in case there are small rounding issues (happens when events anagenetic in biogeography or regimes happen at a very similar time)				
	if(class(int)=="try-error"){
		int<-try(ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=6))
		}	
	if(class(int)=="try-error"){
		int<-try(ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=7))
		}	
	if(class(int)=="try-error"){
		int<-try(ReconcileGeoObjectSMatrix(geo.object=sgeo,S.matrix=smat,rnd=4))
		}	
	
	sgeo<-int$geo.object
	smat<-int$S.matrix
	
	
	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_twoS_geo_ME(phylo,geo.object=sgeo,S.object=smat)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S1<-opt$inferredParams[3]
		S2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), S1 = as.numeric(S1), S2 = as.numeric(S2), nuisance=as.numeric(nuisance),z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_multi_geo_ME(phylo,geo.object=sgeo,r.object=smat)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r1<-opt$inferredParams[3]
		r2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), r1 = as.numeric(r1), r2 = as.numeric(r2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_multi_geo_ME(phylo,geo.object=sgeo,r.object=smat)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b1<-opt$inferredParams[3]
		b2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), b1 = as.numeric(b1),b2 = as.numeric(b2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}

}

}