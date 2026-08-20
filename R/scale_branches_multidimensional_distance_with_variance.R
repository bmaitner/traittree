##########
#'Scale the branches of a phylogeny according to phenotypic change, accounting for uncertainty.
#'
#'scale_branches_multidimensional_with_variation
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param rate If TRUE, branch lengths returned will reflect rates of change, rather than absolute amount of change.
#' @param n_trees Number of trees to generate. Defaults to 1.
#' @param pheno_error Logical, or NULL (the default). Controls whether ancestral state reconstruction estimates a within-species (phenotypic) error term. When NULL, the setting is taken from the data: TRUE if at least one species has two or more observations of at least one trait, and FALSE otherwise. Supply TRUE or FALSE to override. See \code{scale_branches_multidimensional}.
#' @return A phylo object when \code{n_trees} is 1, otherwise a multiPhylo object of length \code{n_trees}.
#' @note This function accounts for uncertainty in the estimated ancestral traits. Each tree is a draw from the joint posterior of the ancestral states, so a node takes a single value shared by all of the branches that meet at it, and correlations both between nodes and between traits are respected. Run with \code{n_trees > 1} to obtain a distribution of trees over which downstream analyses can be repeated.
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import Rphylopars
scale_branches_multidimensional_with_variation<-function(tree,
                                                         traits,
                                                         rate=FALSE,
                                                         n_trees=1,
                                                         pheno_error=NULL){
  
  if(!is.numeric(n_trees) || length(n_trees) != 1 || is.na(n_trees) || n_trees < 1){
    
    stop("n_trees must be a single positive number.", call. = FALSE)
    
  }
  
  n_trees <- as.integer(n_trees)
  
  #First, remove species from trait data that aren't in the phylogeny:
  traits<-traits[which(traits$species%in%tree$tip.label),]
  
  #Decide whether a within-species error term should be estimated, after the
  #filtering above so that species absent from the tree cannot drive the choice.
  pheno_error<-resolve_pheno_error(traits = traits,pheno_error = pheno_error)
  
  #Fit the model once.  Every tree is a fresh draw from the posterior it implies,
  #so the reconstruction is not repeated.
  phylopars_output <- Rphylopars::phylopars(trait_data = traits,
                                            tree = tree,
                                            pheno_error = pheno_error)
  
  posterior_mean <- unname(phylopars_output$anc_recon)
  phylocov <- phylopars_output$pars$phylocov
  phenocov <- phylopars_output$pars$phenocov
  
  output_trees <- vector(mode = "list", length = n_trees)
  
  for(i in seq_len(n_trees)){
    
    sampled_states <- sample_ancestral_states(tree = tree,
                                              traits = traits,
                                              posterior_mean = posterior_mean,
                                              phylocov = phylocov,
                                              phenocov = phenocov,
                                              pheno_error = pheno_error)
    
    tree_i <- tree
    
    #A node carries one sampled value, so the branch lengths are the distances
    #between adjacent nodes of a single coherent realization.
    tree_i$edge.length <- edge_distances_from_node_values(tree = tree,
                                                          node_values = sampled_states,
                                                          rate = rate)
    
    output_trees[[i]] <- tree_i
    
  }
  
  if(n_trees == 1){
    
    return(output_trees[[1]])
    
  }
  
  class(output_trees) <- "multiPhylo"
  
  return(output_trees)  
  
  
}#function

########################
