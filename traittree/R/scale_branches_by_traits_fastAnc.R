#'Scale the branches of a phylogeny according to phenotypic change.
#'
#'scale_branches_by_traits_fastAnc
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param percent If TRUE, branch lengths returned will reflect percent change, rather than absolute amount of change.
#' @return phylo formate phylogeny
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @importFrom  "ape" "drop.tip"
#' @import phytools
scale_branches_by_traits_fastAnc<-function(tree,traits,percent=FALSE){
  
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

