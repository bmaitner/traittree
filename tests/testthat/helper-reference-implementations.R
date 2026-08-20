# Reference implementations: the original loop-based versions of the functions
# that were vectorized, kept so the replacements can be checked against them.
# Taken verbatim from the pre-vectorization source, renamed with a
# `reference_` prefix and with roxygen comments stripped.

reference_scale_branches_by_traits_fastAnc<-function(tree,traits,percent=FALSE){
  
  #Need to have matching species in the phylogeny and the trait data
  #sp_to_remove<-tree$tip.label[which(tree$tip.label%nin%names(traits))]
  
  sp_to_remove <- setdiff(tree$tip.label,names(traits))
  
  
  tree<-ape::drop.tip(tree,sp_to_remove)#
  traits<-traits[tree$tip.label]
  
  
  #Next, estimate ancestral states and character states
  fast_anc_output<-phytools::fastAnc(tree = tree,x = traits)#seems to work well and reasonably fast (seconds for 4351 spp)
  
  #fastAnc does not include extant species, these must be appended to the output
  tip_labels<-as.data.frame(cbind(1:length(tree$tip.label),tree$tip.label))
  names(tip_labels)<-c("tip.number","species")
  trait_df<-as.data.frame(cbind(names(traits),as.numeric(traits)))
  names(trait_df)<-c("species","trait")
  tip.values<-merge(tip_labels,trait_df,by = "species")
  tips_to_append<-as.numeric(as.character(tip.values$trait))
  names(tips_to_append)<-tip.values$tip.number
  fast_anc_output<-c(fast_anc_output,tips_to_append)
  
  #Manual version: for each branch, figure out the corresponding nodes, get the trait values, set the branch length as the difference between those nodes.
  #in the phylo file, phylo$edge is a matrix listing the nodes corresponding to each edge.  the order is identical to that in phylo$edge.length, so changing branch lengths should be easy
  
  if(!percent){
    output_branches<-matrix(data=NA,nrow=length(tree$edge.length),ncol = 1)
    for(i in 1:length(tree$edge.length)){
      node_1<- tree$edge[i,][1]
      node_2<- tree$edge[i,][2]
      
      value_1<-fast_anc_output[as.character(node_1)]
      value_2<-fast_anc_output[as.character(node_2)]
      
      bl<-as.numeric(abs(value_1-value_2))
      output_branches[i]<-bl
      
    }
  }#percent=false
  
  if(percent){
    
    output_branches<-matrix(data=NA,nrow=length(tree$edge.length),ncol = 1)
    for(i in 1:length(tree$edge.length)){
      node_1<- tree$edge[i,][1]
      node_2<- tree$edge[i,][2]
      
      value_1<-fast_anc_output[as.character(node_1)]#beginning of branch
      value_2<-fast_anc_output[as.character(node_2)]#end of branch
      
      bl<-abs(((value_1 - value_2 )/value_1))*100
      
      output_branches[i]<-bl
      
      
    }#i loop 2
    
  }#if percent
  
  
  
  tree$edge.length<-output_branches
  
  return(tree)
  
}

reference_scale_branches_by_traits_rphylopars<-function(tree,traits,percent=FALSE){

invisible("this function should be combined with the other 'scale_branches_rphylopars' functions to eliminate redundancy")  
  
output_trees<-list()    
class(output_trees)<-"multiPhylo"
#First, remove species from trait data that aren't in the phylogeny:
traits<-traits[which(traits$species%in%tree$tip.label),]

#Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
anc_recon<-Rphylopars::phylopars(trait_data = traits,tree = tree,pheno_error = resolve_pheno_error(traits))$anc_recon
row.names(anc_recon)[1:length(tree$tip.label)]<-1:length(tree$tip.label)

if(!percent){
for(x in 1:ncol(anc_recon)){
  trait_x<-colnames(anc_recon)[x]
  output_branches<-matrix(data=NA,nrow=length(tree$edge.length),ncol = 1)
  tree_x<-tree
  for(i in 1:length(tree$edge.length)){
    node_1<- tree$edge[i,][1]
    node_2<- tree$edge[i,][2]
    
    value_1<-anc_recon[,x][which(row.names(anc_recon)==node_1)]
    value_2<-anc_recon[,x][which(row.names(anc_recon)==node_2)]
    
    bl<-as.numeric(abs(value_1-value_2))
    output_branches[i]<-bl
    
  }#i loop
  
  tree_x$edge.length<-output_branches
  output_trees[[x]]<-tree_x
  
  }#x ncol(ie each trait)    

}#if percent=false
  
if(percent){
  for(x in 1:ncol(anc_recon)){
    trait_x<-colnames(anc_recon)[x]
    output_branches<-matrix(data=NA,nrow=length(tree$edge.length),ncol = 1)
    tree_x<-tree
    for(i in 1:length(tree$edge.length)){
      node_1<- tree$edge[i,][1]
      node_2<- tree$edge[i,][2]
      
      value_1<-anc_recon[,x][which(row.names(anc_recon)==node_1)]
      value_2<-anc_recon[,x][which(row.names(anc_recon)==node_2)]
      
      bl<-abs(((value_1 - value_2 )/value_1))*100
      output_branches[i]<-bl
      
    }#i loop
    
    tree_x$edge.length<-output_branches
    output_trees[[x]]<-tree_x
    
  }#x ncol(ie each trait)    
}#if percent=false

  
if(class(output_trees)=="multiPhylo"){
names(output_trees)<-colnames(anc_recon)  
}

return(output_trees)  
  
  
}#function



#################

reference_richness_raster<-function(template.raster,occurrences){
#richness<-NULL
#for( i in 1:ncell(template.raster)){
#  print(i)
#  richness<-c(richness,length(which(occurrencses[,2]==i)))  
#}

#richness<-setValues(x = template.raster,values = richness)

  output_raster<-template.raster
  output_raster<-raster::setValues(x = output_raster,values = NA)
  
  #iterate through all cells with at least one occurrence, record 
  
  output_raster[as.numeric(unique(occurrences[,2]))] <- sapply(X = unique(occurrences[,2]),FUN = function(x){ length(unique(occurrences[which(occurrences[,2]==x),1]))} )

  return(output_raster)

}
