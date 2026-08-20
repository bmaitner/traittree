##########
#'Scale the branches of a phylogeny according to phenotypic change, accounting for uncertainty.
#'
#'scale_branches_multidimensional_with_variation
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param rate If TRUE, branch lengths returned will reflect rates of change, rather than absolute amount of change.
#' @param pheno_error Logical, or NULL (the default). Controls whether ancestral state reconstruction estimates a within-species (phenotypic) error term. When NULL, the setting is taken from the data: TRUE if at least one species has two or more observations of at least one trait, and FALSE otherwise. Supply TRUE or FALSE to override. See \code{scale_branches_multidimensional}.
#' @return phylo formatted phylogeny
#' @note This function accounts for uncertainty in estimated ancestral traits and can be run multiple times to generate a distribution of phylogenies.
#' @note Each node is sampled independently for every edge it touches, so a node
#'   does not carry the same value on all of its branches. The returned tree is
#'   therefore a collection of per-branch draws rather than a single coherent
#'   realization of the ancestral states.
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import Rphylopars
scale_branches_multidimensional_with_variation<-function(tree,traits,rate=FALSE,pheno_error=NULL){
  
  message("Note: currently this function generates one phylogeny at a time, which is dumb.  Brian should modify it to generate a user-supplied number of phylogenies")
    
  #First, remove species from trait data that aren't in the phylogeny:
  traits<-traits[which(traits$species%in%tree$tip.label),]
  
  #Decide whether a within-species error term should be estimated, after the
  #filtering above so that species absent from the tree cannot drive the choice.
  pheno_error<-resolve_pheno_error(traits = traits,pheno_error = pheno_error)
  
  #Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
  phylopars_output <- Rphylopars::phylopars(trait_data = traits,tree = tree,pheno_error = pheno_error)
  anc_recon<- phylopars_output$anc_recon
  anc_var<-phylopars_output$anc_var
  
  row.names(anc_recon)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
  row.names(anc_var)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
  
      tree_x<-tree
      
      #Look every edge up at once rather than one at a time, which is quadratic
      #in tree size.
      node_ids<-as.integer(row.names(anc_recon))
      
      row_1<-match(tree$edge[,1],node_ids)
      row_2<-match(tree$edge[,2],node_ids)
      
      #As before, each end of each edge is drawn independently, so a node
      #touched by several edges is redrawn for each of them.
      value_1<-sample_node_values(means = anc_recon[row_1, ,drop=FALSE],
                                  variances = anc_var[row_1, ,drop=FALSE])
      
      value_2<-sample_node_values(means = anc_recon[row_2, ,drop=FALSE],
                                  variances = anc_var[row_2, ,drop=FALSE])
      
      output_branches<-sqrt(rowSums((value_1 - value_2)^2))
      
      if(rate){
        
        output_branches<-output_branches/tree$edge.length
        
      }
      
      tree_x$edge.length<-unname(output_branches)
      
  return(tree_x)  
  
  
}#function


#' Draw one normal deviate per entry of a matrix of means and variances
#'
#' Internal helper.
#'
#' @param means Matrix of means.
#' @param variances Matrix of variances, matching `means`.
#' @return Matrix of draws, the same shape as `means`.
#' @noRd
sample_node_values <- function(means, variances){

  matrix(data = stats::rnorm(n = length(means),
                             mean = means,
                             sd = sqrt(variances)),
         nrow = nrow(means),
         ncol = ncol(means))

}

########################
