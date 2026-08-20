##########

#' Euclidean distance in trait space across each edge of a tree
#'
#' Internal helper.  Returns, for every edge of `tree`, the Euclidean distance
#' between the reconstructed trait vectors of the two nodes it subtends.
#'
#' This is done in a single vectorized step rather than one edge at a time:
#' looking each node up with `which(row.names(anc_recon) == node)` is a full
#' scan of the row names, so an edge-at-a-time loop is quadratic in tree size.
#'
#' @param tree A phylogeny with branch lengths in units of time.
#' @param anc_recon Matrix of reconstructed trait values, as returned by
#'   `Rphylopars::phylopars()`, with tips in `tree$tip.label` order followed by
#'   internal nodes named by node number.
#' @param rate If TRUE, distances are divided by the original branch lengths,
#'   giving rates of change rather than absolute amounts.
#' @return Numeric vector, one element per edge, in `tree$edge` order.
#' @noRd
edge_trait_distances <- function(tree,
                                 anc_recon,
                                 rate = FALSE){

  row.names(anc_recon)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)

  node_ids <- as.integer(row.names(anc_recon))

  row_1 <- match(tree$edge[,1], node_ids)
  row_2 <- match(tree$edge[,2], node_ids)

  squared_change <- (anc_recon[row_1, , drop = FALSE] -
                       anc_recon[row_2, , drop = FALSE])^2

  branch_lengths <- sqrt(rowSums(squared_change))

  if(rate){

    branch_lengths <- branch_lengths/tree$edge.length

  }

  unname(branch_lengths)

}

#'Scale the branches of a phylogeny according to phenotypic change.
#'
#'scale_branches_multidimensional
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param rate If TRUE, branch lengths returned will reflect rates of change, rather than absolute amount of change.
#' @param return_traits Logical. Defaults = FALSE. IF TRUE, object returned will be a list containing 1) the phylogeny, and 2) traits for the tips.
#' @param pheno_error Logical, or NULL (the default). Controls whether ancestral state reconstruction estimates a within-species (phenotypic) error term. When NULL, the setting is taken from the data: TRUE if at least one species has two or more observations of at least one trait, and FALSE otherwise. Supply TRUE or FALSE to override. Note that this reproduces the results of the Rphylopars defaults rather than changing them: with a single observation per species phylopars() fits no within-species term regardless of the setting. It matters for replicated data, where the two choices give substantially different reconstructions, and a warning is issued if replication is too sparse to be clearly intentional.
#' @return phylo formate phylogeny
#' @note This function DOES NOT account for uncertainty in estimated ancestral traits.
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import Rphylopars
scale_branches_multidimensional <- function(tree,
                                          traits,
                                          rate = FALSE,
                                          return_traits = FALSE,
                                          pheno_error = NULL){
  
  
  #First, remove species from trait data that aren't in the phylogeny:
  
    traits <- traits[which(traits$species%in%tree$tip.label),]
  
  #Decide whether a within-species error term should be estimated.  This is done
  #after the filtering above, so that species absent from the tree cannot drive
  #the choice.
  
    pheno_error <- resolve_pheno_error(traits = traits,
                                       pheno_error = pheno_error)
  
  #Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
  
    anc_recon <- phylopars(trait_data = traits,
                           tree = tree,
                           pheno_error = pheno_error)$anc_recon
    
    if(return_traits){
      
      inferred_traits <- anc_recon[1:length(tree$tip.label),]
      
    }
    
      tree_x <- tree
      
      tree_x$edge.length <- edge_trait_distances(tree = tree,
                                                 anc_recon = anc_recon,
                                                 rate = rate)
      
      # Append trait data if requested
      
      if(return_traits){
        
        output <- list(tree = tree_x,
                       traits = inferred_traits)
        
        return(output)

      }
      
  return(tree_x)  
  
  
}#function
