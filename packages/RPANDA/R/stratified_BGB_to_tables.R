.stratified_BGB_to_tables<-function(tree,clado_events_tables,i){
##hack to get ana_events_table for constructing stochastic maps
hold<-.events_txt_list_into_events_table(events_txt_list=clado_events_tables[[i]]$anagenetic_events_txt_below_node, trtable=NULL, recalc_abs_ages=TRUE)
rows<-which(clado_events_tables[[i]][,34]!="none")
if(length(rows)<1){
	ana.int<-NULL
	}else{
row.len<-length(rows)
if(row.len<dim(hold)[1]){
	ho<-clado_events_tables[[i]][,34]
	torep<-vector()
	for(n in 1:length(ho)){
		torep<-c(torep,length(strsplit(ho[n],';')[[1]]))
		}
	torep<-torep[rows]		
	rows2<-vector()
	for(m in 1:row.len){
		rows2<-c(rows2,rep(rows[m],torep[m]))
		}
	rows<-rows2	
	}	

desc<-clado_events_tables[[i]][rows,1]
anc<-clado_events_tables[[i]][rows,5]
anmat<-cbind(desc,anc,hold)
m<-regexpr("->",anmat$event_txt)
regs<-regmatches(as.character(anmat$event_txt),m,invert=T)
new_rangetxt<-unlist(lapply(regs,function(x)x[2]))
nodenum_at_top_of_branch<-anmat[,1]
abs_event_time<-anmat$abs_event_time
ana.int<-data.frame(abs_event_time,nodenum_at_top_of_branch,new_rangetxt)
ana.int$new_rangetxt<-as.character(ana.int$new_rangetxt)
}
    ############################################
    #### create CLADOGENETIC events tables
    #RES_clado_events_tables[[1]][,1:20] #model
    ###impt items needed for *** ;clado.events	$node, $label ,1 ($node),9 ($node_ht),15 (left_desc_nodes),16 (right_desc_nodes),20 (clado_event_txt)
    ############################################ 

clado.int<-clado_events_tables[[i]]
clado.int$clado_event_txt<-as.character(clado.int$clado_event_txt)
if(any(is.na(clado.int$clado_event_txt))){clado.int$clado_event_txt[which(is.na(clado.int$clado_event_txt))]<-''}
clado.int<-clado.int[which(clado.int$clado_event_txt!=''),]
clado.int[,12]<-clado.int$label
node_ht = round(max(nodeHeights(tree))-clado.int$time_bp,5)
clado.int[,9]<-node_ht
lef<-vector()
rig<-vector()
for(i in 1:length(clado.int$node)){
	nin<-clado.int$node[i]
	hol<-tree$edge[tree$edge[,1]==nin,2]
 	lef<-c(lef,hol[1])
	rig<-c(rig,hol[2])
	}
clado.int[,15]<-lef
clado.int[,16]<-rig
clado.int[,20]<-clado.int$clado_event_txt
colnames(clado.int)[c(9,12,15,16,20)]<-c("node_ht","label","left_desc_nodes","right_desc_nodes","clado_event_txt")

int2<-matrix(nrow=length(tree$tip.label),ncol=dim(clado.int)[2]) #for tips
int2<-as.data.frame(int2)
int2[,1]<-1:length(tree$tip.label)
int2[,12]<-tree$tip.label
colnames(int2)<-colnames(clado.int)
int<-rbind(clado.int,int2)
clado.int<-int

return(list(ana.int=ana.int,clado.int=clado.int))
}
