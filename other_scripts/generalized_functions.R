#out<-make.simmap.BGB(anc.phylo=motmot.tree,subclade.phylo=motmot.tree,ana.events=smap$ana.int,clado.events=smap$clado.int)
#geo.simmap<-out$geo.simmap
#class.object<-out$class.object

##24 July 2018: fixing this to

##this is a generalizable function to flexibly return a class.df function for any geo.simmap, that is, any case with biogeography where ranges are indicated by single capital letters
return.class.df<-function(simmap,class.object){
	states<-colnames(simmap$mapped.edge)
	for(i in 1:length(states)){
		st.id=paste("c(",paste(which(grepl(paste(strsplit(states[i],split="")[[1]],collapse="|"),states)),collapse=","),")",sep="") #this gives the columns to extract from class.df
		#eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]==states[',i,']))',sep="")))
		eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]%in%states[',st.id,']))',sep="")))
	}
	eval(parse(text=paste('return(data.frame(interval=1:length(d1),',paste('d',1:length(states),sep="",collapse=','),'))',sep="")))
}

##this '*_sympatric' version builds a class.df for the case where a simmap only has one state (i.e., the clade is in sympatry)
return.class.df_sympatric<-function(simmap){
	states<-colnames(simmap$mapped.edge)
	d1<-2:length(simmap$tip.label)
	eval(parse(text=paste('return(data.frame(interval=1:length(d1),',paste('d',1:length(states),sep="",collapse=','),'))',sep="")))
}


##this '*_subgroup' version is the version that should be used in all instances where there isn't biogeography (return.class.df should be used for biogeo cases only)
return.class.df_subgroup<-function(simmap,class.object){
	states<-colnames(simmap$mapped.edge)
	for(i in 1:length(states)){
		st.id=i #this gives the columns to extract from class.df
		#eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]==states[',i,']))',sep="")))
		eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]%in%states[',st.id,']))',sep="")))
	}
	eval(parse(text=paste('return(data.frame(interval=1:length(d1),',paste('d',1:length(states),sep="",collapse=','),'))',sep="")))
}


##this is a generalizable function to flexibly return a list of functions (to pass to fit_t_general*) for any geo.simmap

create.function.list_old<-function(geo.simmap, geo.class.object,geo.class.df){
	if(is.null(geo.simmap)){stop('provide geo.simmap')}
	if(is.null(geo.class.object)){stop('provide geo.class.object')}
	if(is.null(geo.class.df)){stop('provide class.df')}
	if(is.null(geo.class.object$times)){stop('$times missing from geo.class.object')}
	states<-colnames(geo.simmap$mapped.edge)
	funlist<-list()
	if(dim(geo.class.df)[2]!=(length(states)+1)){stop('geo.class.df is of incorrect dimensions')}
	
	for(i in 1:length(states)){
		st<-states[i]
		#cols=paste("c(",paste(which(grepl(paste(strsplit(st,split="")[[1]],collapse="|"),states))+1,collapse=","),")",sep="") #this gives the columns to extract from class.df
		#eval(parse(text=paste("funlist[[",i,"]]<-function(x){;values <- geo.class.object$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(as.matrix(geo.class.df[res,",cols,"]));}",sep="")))
		#eval(parse(text=paste("funlist[[",i,"]]<-function(x){;values <- ",deparse(substitute(geo.class.object)),"$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(as.matrix(",deparse(substitute(geo.class.df)),"[res,",cols,"]));}",sep="")))
		eval(parse(text=paste("funlist[[",i,"]]<-function(x){;values <- ",deparse(substitute(geo.class.object)),"$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(",deparse(substitute(geo.class.df)),"[res,",i+1,"]);}",sep="")))
		}
			
	return(funlist)

}

create.function.list<-function(geo.simmap,df, times){
	if(is.null(geo.simmap)){stop('provide geo.simmap')}
	if(is.null(times)){stop('provide vector of times from class.object')}
	if(is.null(df)){stop('provide class.df')}
	#if(is.null(geo.class.object$times)){stop('$times missing from geo.class.object')}
	states<-colnames(geo.simmap$mapped.edge)
	funlist<-list()
	if(dim(df)[2]!=(length(states)+1)){stop('class.df is of incorrect dimensions')}
	
	for(i in 1:length(states)){
		st<-states[i]
		#cols=paste("c(",paste(which(grepl(paste(strsplit(st,split="")[[1]],collapse="|"),states))+1,collapse=","),")",sep="") #this gives the columns to extract from class.df
		#eval(parse(text=paste("funlist[[",i,"]]<-function(x){;values <- geo.class.object$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(as.matrix(geo.class.df[res,",cols,"]));}",sep="")))
		#eval(parse(text=paste("funlist[[",i,"]]<-function(x){;values <- ",deparse(substitute(geo.class.object)),"$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(as.matrix(",deparse(substitute(geo.class.df)),"[res,",cols,"]));}",sep="")))
		#eval(parse(text=paste("funlist[[",i,"]]<-function(x,df,times){;values <- ",deparse(substitute(times)),";res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(",deparse(substitute(df)),"[res,",i+1,"]);}",sep="")))
		eval(parse(text=paste("funlist[[",i,"]]<-function(x,df,times){;values <- times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(df[res,",i+1,"]);}",sep="")))
		}
			
	return(funlist)

}

create.function.list.EB1<-function(geo.simmap){
	if(is.null(geo.simmap)){stop('provide geo.simmap')}
	states<-colnames(geo.simmap$mapped.edge)
	funlist<-list()	
	for(i in 1:length(states)){
		eval(parse(text=paste0("funlist[[",i,"]]<-function(x,df,times){;return(x);}")))
		}
			
	return(funlist)

}

create.function.list.EBmulti<-function(regime.simmap){
	if(is.null(regime.simmap)){stop('provide regime.simmap')}
	states<-colnames(regime.simmap$mapped.edge)
	class.object<-try(CreateClassObject(regime.simmap))
	if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(regime.simmap,rnd=6))}
	if(class(class.object)=="try-error"){class.object<-CreateClassObject(regime.simmap,rnd=7)}
	
	funlist<-list()	
	for(i in 1:length(states)){
		first.time.bin=min(which(lapply(class.object$class.object,function(x) states[i]%in%x[,2])==TRUE))
		if(first.time.bin==1){ #if state is present at the root
		eval(parse(text=paste0("funlist[[",i,"]]<-function(x,df,times){;return(x);}")))
		} else {
		time.since.root = class.object$times[first.time.bin]
		eval(parse(text=paste0("funlist[[",i,"]]<-function(x,df,times){;return(x-",time.since.root,");}")))
		}
		}
			
	return(funlist)

}

create.function.list_sympatric<-function(simmap,class.df){
	if(is.null(class.df)){stop('provide class.df')}
	states<-colnames(simmap$mapped.edge)
	funlist<-list()
	if(dim(class.df)[2]!=(length(states)+1)){stop('class.df is of incorrect dimensions')}
	
	for(i in 1:length(states)){
		st<-states[i]
		#cols=paste("c(",paste(which(grepl(paste(strsplit(st,split="")[[1]],collapse="|"),states))+1,collapse=","),")",sep="") #this gives the columns to extract from class.df
		#eval(parse(text=paste("funlist[[",i,"]]<-function(x){;values <- simmap$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(as.matrix(geo.class.df[res,",cols,"]));}",sep="")))
		#eval(parse(text=paste("funlist[[",i,"]]<-function(x){;values <- ",deparse(substitute(geo.class.object)),"$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(as.matrix(",deparse(substitute(geo.class.df)),"[res,",cols,"]));}",sep="")))
		eval(parse(text=paste("funlist[[",i,"]]<-function(x)		 {;values <- ",deparse(substitute(simmap)),"$times;res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(",deparse(substitute(class.df)),"[res,",i+1,"]);}",sep="")))
		eval(parse(text=paste("funlist[[",i,"]]<-function(x,df,times){;values <- times;								   res <- findInterval(x, values);index <- res==0;res[index==TRUE] <- 1;return(df[res,",i+1,"]);}",sep="")))
		}
			
	return(funlist)

}


###below are other functions from BioGeoBEARS

#code from Nick Matzke's code on BioGeoBEARS wiki
events_txt_list_into_events_table<-function(events_txt_list, trtable=NULL, recalc_abs_ages=TRUE)
	{
	
	if (is.null(events_txt_list))
		{
		errortxt = paste("\nWARNING in events_txt_list_into_events_table(): your events_txt_list has NO events!\n\nThis means your tree has NO d/e/a events across the whole tree.\nThis is *expected* e.g. if you inferred d=e=0 under DEC+J. Input a list of '' or NA to avoid this error.\n\n", sep="")
		cat(errortxt)
		errortxt2 = paste("events_txt_list_into_events_table() is returning NULL which will might cause issues later.\n\n", sep="")
		cat(errortxt2)
		return(NULL)
		}
	
	# Convert NAs to "none"
	events_txt_list[is.na(events_txt_list)] = "none"
	
	# Remove lines with no events or NA:
	noneTF = events_txt_list == "none"
	keepTF = (noneTF == FALSE)
	events_txt_list = events_txt_list[keepTF]
	
	
	# If no anagenetic events, return NULL
	if (length(events_txt_list) == 0)
		{
		events_table = NULL
		return(events_table)
		}


	# Include the trtable, if that is input
	if (length(trtable) > 0)
		{
		trtable_subset = NULL
		}

	
	# Convert the events text back into a table:
	tmptable = NULL
	for (i in 1:length(events_txt_list))
		{
		#print(events_txt_list)
		tmptable_rows = events_txt_into_events_table(events_txt_list[i])
		rownums_in_trtable = as.numeric(tmptable_rows$nodenum_at_top_of_branch)
		#print(tmptable_rows)
		num_newrows = nrow(tmptable_rows)
		tmptable = rbind(tmptable, tmptable_rows)
		} # END for (i in 1:length(events_txt_list))
	events_table = BioGeoBEARS::dfnums_to_numeric(BioGeoBEARS::adf2(tmptable))
	names(events_table) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")	
	
	return(events_table)
	}

events_txt_into_events_table<-function(branch_events_txt)
	{
	words = strsplit(branch_events_txt, split=";")[[1]]
	
	events_table_for_branch = t(sapply(X=words, FUN=event_txt_into_events_row))
	row.names(events_table_for_branch) = NULL
	events_table_for_branch
	
	events_table_for_branch = BioGeoBEARS::adf2(events_table_for_branch)
	events_table_for_branch
	names(events_table_for_branch) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")
	
	return(events_table_for_branch)
	}

event_txt_into_events_row<-function(word)
	{

	split_key_item<-function(word2)
		{
		output_pair = c("", "")
		words3 = strsplit(word2, split=":")[[1]]
		numwords = length(words3)
		
		output_pair[1:numwords] = words3[1:numwords]
		
		return(output_pair)
		}


	words2 = strsplit(word, split=",")[[1]]
	output = sapply(X=words2, FUN=split_key_item)
	tmprow = matrix(data=output[2,], nrow=1)
	return(tmprow)
	}
	