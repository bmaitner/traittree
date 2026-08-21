#'Scale the branches of a phylogeny according to phenotypic change.
#'
#'scale_branches_by_traits_rphylopars
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param percent If TRUE, branch lengths returned will reflect percent change, rather than absolute amount of change.
#' @param pheno_error Logical, or NULL (the default). Controls whether ancestral state reconstruction estimates a within-species (phenotypic) error term. When NULL, the setting is taken from the data: TRUE if at least one species has two or more observations of at least one trait, and FALSE otherwise. Supply TRUE or FALSE to override. See \code{scale_branches_multidimensional}.
#' @return phylo formatted phylogeny
#' @note This function DOES NOT account for uncertainty in estimated ancestral traits.
#' @examples
#' # One tree per trait, each scaled by change in that trait alone
#' trees <- scale_branches_by_traits_rphylopars(tree = example_tree,
#'                                              traits = example_traits)
#'
#' names(trees)
#'
#' plot(trees$body_mass_g$edge.length, trees$wing_mm$edge.length,
#'      xlab = "change in body mass", ylab = "change in wing length")
#' @export
#' @import Rphylopars
scale_branches_by_traits_rphylopars<-function(tree,traits,percent=FALSE,pheno_error=NULL){

message("this function should be combined with the other 'scale_branches_rphylopars' functions to eliminate redundancy")  
  
output_trees<-list()    
class(output_trees)<-"multiPhylo"
#First, remove species from trait data that aren't in the phylogeny:
traits<-traits[which(traits$species%in%tree$tip.label),]

#Decide whether a within-species error term should be estimated, after the
#filtering above so that species absent from the tree cannot drive the choice.
pheno_error<-resolve_pheno_error(traits = traits,pheno_error = pheno_error)

#Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
anc_recon<-phylopars(trait_data = traits,tree = tree,pheno_error = pheno_error)$anc_recon
row.names(anc_recon)[1:length(tree$tip.label)]<-1:length(tree$tip.label)

#Look every edge up at once rather than one at a time, which is quadratic in
#tree size, and do it for all traits in one go.
node_ids<-as.integer(row.names(anc_recon))

value_1<-anc_recon[match(tree$edge[,1],node_ids), ,drop=FALSE]#beginning of branch
value_2<-anc_recon[match(tree$edge[,2],node_ids), ,drop=FALSE]#end of branch

if(!percent){
  
  output_branches<-abs(value_1 - value_2)
  
}else{
  
  warn_unstable_percent_change(denominator = value_1)
  
  output_branches<-abs((value_1 - value_2)/value_1)*100
  
}

#One tree per trait, each carrying that trait's column of branch lengths.
for(x in 1:ncol(anc_recon)){
  
  tree_x<-tree
  tree_x$edge.length<-unname(output_branches[,x])
  output_trees[[x]]<-tree_x
  
}#x ncol(ie each trait)    

if(inherits(output_trees,"multiPhylo")){
names(output_trees)<-colnames(anc_recon)  
}

return(output_trees)  
  
  
}#function



#################
