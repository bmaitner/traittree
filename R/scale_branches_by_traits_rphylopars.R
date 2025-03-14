#'Scale the branches of a phylogeny according to phenotypic change.
#'
#'scale_branches_by_traits_rphylopars
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param percent If TRUE, branch lengths returned will reflect percent change, rather than absolute amount of change.
#' @return phylo formatted phylogeny
#' @note This function DOES NOT account for uncertainty in estimated ancestral traits.
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import Rphylopars
scale_branches_by_traits_rphylopars<-function(tree,traits,percent=FALSE){

message("this function should be combined with the other 'scale_branches_rphylopars' functions to eliminate redundancy")  
  
output_trees<-list()    
class(output_trees)<-"multiPhylo"
#First, remove species from trait data that aren't in the phylogeny:
traits<-traits[which(traits$species%in%tree$tip.label),]

#Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
anc_recon<-phylopars(trait_data = traits,tree = tree)$anc_recon
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



