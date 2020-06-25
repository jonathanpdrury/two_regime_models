#' Simulates data under a model of matching competition with two competitive regimes, specified in S.matrix
#'
#' @param phylo A phylogenetic tree
#' @param pars A matrix, columns containing values for sig2, S1 (first row of S.matrix), S2 (second row of S.matrix), root.value
#' @param min.Nsegments minimum number of time steps to divide the phylogeny into for the simulation. Default is 2500
#' @param plot should the last run be plotted?
#' @param S.matrix, created using "CreateClassObject" then passed to "CreateSMatrix" (need to wrap this within function in next version)

sim_multiS<-function(phylo,pars, verbose=TRUE, min.Nsegments=2500, plot=FALSE,S.matrix,rnd=6){
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
	#something to check that it isn't sorted
  nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
  root <- length(phylo$tip.label) + 1
  heights<-phytools::nodeHeights(phylo)
  totlen<-max(heights)
  len<-root-1
  nodeDist<-c(as.numeric(sort(max(ape::branching.times(phylo))-ape::branching.times(phylo))),totlen)
  nodeDiff<-diff(nodeDist)
  if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
  old.labels<-as.numeric(names(sort(ape::branching.times(phylo),decreasing=TRUE)))
  old.edge<-phylo$edge
  if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
  	checkmat<-cbind(old.labels,seq(root,len+phylo$Nnode))
  	old.edge<-phylo$edge
  	for(j in 1:phylo$Nnode){phylo$edge[which(old.edge==checkmat[j,1])]<-checkmat[j,2]}
  	}
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:phylo$Nnode){
    other<-phylo$edge[phylo$edge[,1]==i+len, 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+len
      if(b>len){
        counter_three_letters <- counter_three_letters + 1
        int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
        int[3]<-b
      } else {
        int[2]<-phylo$tip.label[[b]]
        int[3]<-0 ##NOTE :I am considering tips to be "0" here and use this below
      }
      mat<-rbind(mat,int)
    }
  }
  nat<-list()
  branchesPresent <- rep(NA, length(nodeDiff))
  for(i in 1:length(nodeDiff)){
    if(i==1){
        nat[[i]]<-mat[mat[,1]==(len+i),2]
      } else {
        IN<-vector()
        P<-mat[as.numeric(mat[,1])<=(len+i),c(2,3)]
        IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(len+i),1])
        nat[[i]]<-IN
      }
    branchesPresent[i] = length(nat[[i]])
  }
	
  masterbranch<-list()
  segmentSize <- rep(NA, phylo$Nnode)
  mappings <- list()
  segsize = sum(nodeDiff)/2500

  for(i in 1:phylo$Nnode){ ##for each node interval

    nati<-nat[[i]]
    mappings[[i]] = rep(NA, branchesPresent[i])

    if(i==1){
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
    }else{
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
      for (j in 1:branchesPresent[i]){
        if(length(which(nati[j] == pastNati)) > 0){
          mappings[[i]][j] = which(nati[j] == pastNati)
        }else{
          mappings[[i]][j] = which(mat[mat[,3]==mat[mat[,2]==nati[j],1],2] == pastNati)
        }
      }
    }
    pastNati <- nati
  }

  ## Looping over the parameters

  out <- matrix( nrow = nrow(pars), ncol = len)

  for (p in 1:nrow(pars)){

    sig2 = pars[p,1]
    S1 = pars[p,2]
    S2 = pars[p,3]
    root.value = pars[p,4]

	smat = S.matrix$S.matrix
	if(!all(colnames(smat[[length(smat)]])==nat[[length(nat)]])){stop("order of S.matrix species is incorrect; perhaps this is a sorted matrix?")}	
    newDist<-S.matrix$times
    newDiff<-S.matrix$spans


    timecount=1
    for(i in 1:phylo$Nnode){
        traitMat <- matrix(nrow = branchesPresent[i], ncol = segmentSize[i]+1)

        if (i == 1){
          traitMat[,1] = root.value
        }else{
          traitMat[,1] = masterbranch[[i-1]][mappings[[i]],(segmentSize[i-1]+1)] #added +1 here since segmentSize[i-1] is penultimate column, not the last one
        }

        tempInd<- 1:branchesPresent[i] # hack to have fast selection of not k, seemed to be faster than a call to which()

        for(k in 1:segmentSize[i]){

	          	for(j in 1:branchesPresent[i]){
	            
	            elmsS1<- if(smat[[timecount]][1,j]==1){which(smat[[timecount]][1,]==1)}else{logical(0)}#these elements are sympatric and in S1
	            elmsS2<- if(smat[[timecount]][2,j]==1){which(smat[[timecount]][2,]==1)}else{logical(0)}#these elements are sympatric and in S2
       		 	tempS1 =  ifelse(length(elmsS1)>0,S1*(mean(traitMat[elmsS1,k])-traitMat[j,k])*segsize,0)
       		 	tempS2 =  ifelse(length(elmsS2)>0,S2*(mean(traitMat[elmsS2,k])-traitMat[j,k])*segsize,0)
      		    temp2 = sqrt(sig2*segsize)
            	traitMat[j,k+1]<-traitMat[j,k] +(tempS1+tempS2)/sum(smat[[timecount]][,j]) +rnorm(1,0,temp2)

				if((k!=segmentSize[i])&&(j==branchesPresent[i])){timecount= ifelse(round(((k*segsize)+nodeDist[i]),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
				###loop for last segment size (to preserve exact branch lengths)
				if(k==segmentSize[i] && nodeDiff[i]%%segsize!=0){
					segsizeT= nodeDiff[i]%%segsize
       		 		tempS1B =  ifelse(length(elmsS1)>0,S1*(mean(traitMat[elmsS1,k])-traitMat[j,k])*segsizeT,0)
       		 		tempS2B =  ifelse(length(elmsS2)>0,S2*(mean(traitMat[elmsS2,k])-traitMat[j,k])*segsizeT,0)
	        		temp2B = sqrt(sig2*segsizeT)
	            	traitMat[j,k+1]<-traitMat[j,k] +(tempS1B+tempS2B)/sum(smat[[timecount]][,j]) +rnorm(1,0,temp2B)
				
				 	if(j==branchesPresent[i]){timecount= ifelse(round((((k-1)*segsize)+nodeDist[i]+segsizeT),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
					}
	        	}	        	
			}
        masterbranch[[i]] = traitMat
      }
    out[p,] = as.vector(masterbranch[[phylo$Nnode]][, tail(segmentSize,n=1)+1])
    if(p==1){colnames(out)<-unlist(nat[[i]])}
  }

  if(plot==TRUE){
    print("plotting last simulated dataset")
    M=seq(0,sum(nodeDiff),length=sum(segmentSize))
    O=list()
    for(i in 1:length(nodeDiff)){
    O[[i]]<-data.frame(seq(nodeDist[[i]],nodeDist[[i+1]],length=(segmentSize[[i]]+1)),as.data.frame(t(masterbranch[[i]])))
    }
    t.plot<-plot(M,1:length(M),col="white", ylim=c(range(sapply(masterbranch,range))), xlab="time", ylab="Value")
    for(i in 1:length(nodeDiff)){
    	for(j in 1:(ncol(O[[i]])-1)){
    		lines(O[[i]][,1],O[[i]][,j+1])
    	}
    }
  }

  return(out)
}


