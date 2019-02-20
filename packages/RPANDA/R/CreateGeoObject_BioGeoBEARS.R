
CreateGeoObject_BioGeoBEARS<-function(full.phylo,trimmed.phylo=NULL,ana.events,clado.events,stratified=FALSE){

	if(stratified){
		if(is.null(trimmed.phylo)){
			clado_events_tables<-list()
			clado_events_tables[[1]]<-clado.events
			smap<-.stratified_BGB_to_tables(full.phylo,clado_events_tables,1)
			x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=full.phylo,ana.events=smap$ana.int,clado.events=smap$clado.int,nat.only=FALSE)
			return(x)
			}else{
			clado_events_tables<-list()
			clado_events_tables[[1]]<-clado.events
			smap<-.stratified_BGB_to_tables(full.phylo,clado_events_tables,1)
			x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=trimmed.phylo,ana.events=smap$ana.int,clado.events=smap$clado.int,nat.only=FALSE)
			return(x)		
			}		
	} else{
		if(is.null(trimmed.phylo)){
			x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=full.phylo,ana.events=ana.events,clado.events=clado.events,nat.only=FALSE)
			return(x)
		} else{
			x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=trimmed.phylo,ana.events=ana.events,clado.events=clado.events,nat.only=FALSE)
			return(x)
		}
	}
}