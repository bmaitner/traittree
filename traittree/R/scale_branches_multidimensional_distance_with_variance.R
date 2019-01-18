
##########
#'Scale the branches of a phylogeny according to phenotypic change, accounting for uncertainty.
#'
#'scale_branches_multidimensional_with_variation
#' @param tree A phylogeny with branch lengths in units of time.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param rate If TRUE, branch lengths returned will reflect rates of change, rather than absolute amount of change.
#' @return phylo formate phylogeny
#' @note This function accounts for uncertainty in estimated ancestral traits and can be run multiple times to generate a distribution of phylogenies.
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import Rphylopars
#' @importFrom "stats" "dist"
scale_branches_multidimensional_with_variation<-function(tree,traits,rate=F){
  
  message("Note: currently this function generates one phylogeny at a time, which is dumb.  Brian should modify it to generate a user-supplied number of phylogenies")
    
  #First, remove species from trait data that aren't in the phylogeny:
  traits<-traits[which(traits$species%in%tree$tip.label),]
  
  #Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
  phylopars_output <- Rphylopars::phylopars(trait_data = traits,tree = tree)
  anc_recon<- phylopars_output$anc_recon
  anc_var<-phylopars_output$anc_var
  
  row.names(anc_recon)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
  row.names(anc_var)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
  
  
      output_branches<-matrix(data=NA,nrow=length(tree$edge.length),ncol = 1)
      tree_x<-tree
      for(i in 1:length(tree$edge.length)){
        node_1<- tree$edge[i,][1]
        node_2<- tree$edge[i,][2]
        
        #value_1<-anc_recon[which(row.names(anc_recon)==node_1),]
        #value_2<-anc_recon[which(row.names(anc_recon)==node_2),]
        
        value_1<-apply(X = cbind(anc_recon[which(row.names(anc_recon)==node_1),],anc_var[which(row.names(anc_var)==node_1),]),MARGIN = 1,FUN = function(x){stats::rnorm(n = 1,mean = x[1],sd = sqrt(x[2])) })
        value_2<-apply(X = cbind(anc_recon[which(row.names(anc_recon)==node_2),],anc_var[which(row.names(anc_var)==node_2),]),MARGIN = 1,FUN = function(x){stats::rnorm(n = 1,mean = x[1],sd = sqrt(x[2])) })
        
        
        bl<-stats::dist(rbind(value_1,value_2))[1]
        
        if(rate){
        bl<- (bl/tree$edge.length[i])
          
        }
        
        output_branches[i]<-bl
        
      }#i loop
      
      tree_x$edge.length<-output_branches
      
    
  
  return(tree_x)  
  
  
}#function


########################

