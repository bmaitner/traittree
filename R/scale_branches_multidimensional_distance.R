##########

#'Scale the branches of a phylogeny according to phenotypic change.
#'
#'scale_branches_multidimensional
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param rate If TRUE, branch lengths returned will reflect rates of change, rather than absolute amount of change.
#' @param return_traits Logical. Defaults = FALSE. IF TRUE, object returned will be a list containing 1) the phylogeny, and 2) traits for the tips.
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
                                          return_traits = FALSE){
  
  
  #First, remove species from trait data that aren't in the phylogeny:
  
    traits <- traits[which(traits$species%in%tree$tip.label),]
  
  #Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
  
    anc_recon <- phylopars(trait_data = traits,tree = tree)$anc_recon
    
    if(return_traits){
      
      inferred_traits <- anc_recon[1:length(tree$tip.label),]
      
    }
    
    row.names(anc_recon)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
  
      tree_x <- tree
      
      # Branch lengths are the Euclidean distance in trait space between the two
      # nodes subtending each edge.  Done in a single vectorized step: looking the
      # nodes up one edge at a time is quadratic in tree size.
      
      node_ids <- as.integer(row.names(anc_recon))
      
      row_1 <- match(tree$edge[,1], node_ids)
      row_2 <- match(tree$edge[,2], node_ids)
      
      squared_change <- (anc_recon[row_1, , drop = FALSE] -
                           anc_recon[row_2, , drop = FALSE])^2
      
      output_branches <- sqrt(rowSums(squared_change))
      
      if(rate){
        
        output_branches <- output_branches/tree$edge.length
        
      }
      
      tree_x$edge.length <- unname(output_branches)
      
      # Append trait data if requested
      
      if(return_traits){
        
        output <- list(tree = tree_x,
                       traits = inferred_traits)
        
        return(output)

      }
      
  return(tree_x)  
  
  
}#function
