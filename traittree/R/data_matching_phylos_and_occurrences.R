#'Prune trait and phylogeny to a common set of species
#'
#'data_matching
#' @param occurrences A set of occurrences in "tidy" format: first column is species name, second is raster cells where the species occurs.
#' @param phylogeny A phylogeny with time-scaled branch lengths
#' @return List the pruned phylogeny and occurrences data
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @importFrom  "ape" "drop.tip"
data_matching<-function(phylogeny,occurrences){
  #removing occurrences that aren't in the phylogeny
  phylo_species<-phylogeny$tip.label
  occurrences[,1]<-gsub(" ","_",occurrences[,1])
  
  occurrences_pruned<-occurrences[which(occurrences[,1]%in%phylo_species),]
  occurrences<-occurrences_pruned
  rm(occurrences_pruned)
  
  #remove species from phylogeny that aren't needed
  species_to_prune<-setdiff(x = phylo_species,y = unique(occurrences[,1]))
  phylogeny_pruned<-ape::drop.tip(phylogeny,species_to_prune)
  rm(species_to_prune)
  phylogeny<-phylogeny_pruned
  rm(phylogeny_pruned)
  rm(phylo_species)
  
  outlist<-list(phylogeny,occurrences)
    
  names(outlist)<-c("phylogeny","occurrences")
  
  return(outlist)
  
  #call<-deparse(sys.call())
  #call<-unlist(strsplit(call,split="[]),(]"))
  #occurrence_file<-unlist(strsplit(grep(pattern="occurrences",x=call,value=TRUE),split="= "))[2]
  #phylogeny_file<-unlist(strsplit(grep(pattern="phylogeny",x=call,value=TRUE),split="= "))[2]
  #assign(x=phylogeny_file,value = phylogeny,envir = parent.env((environment())))
  #assign(x=occurrence_file,value = occurrences,envir = parent.env((environment())))

  
    
}
