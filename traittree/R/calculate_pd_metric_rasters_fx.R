#'Calculate time- and trait-scaled phylogenies and pd/pdi rasters
#'
#'calculate_pd_metric_rasters
#' @param occurrences A set of occurrences in "tidy" format: first column is species name, second is raster cells where the species occurs.
#' @param phylogeny A phylogeny with time-scaled branch lengths
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param template.raster The raster that corresponds to the cell numbers in the second column of the occurrences file
#' @return List containing 4 rasters and 3 phylogenies (original, trait-scaled and rate-scaled)
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import raster
#' @import PhyloMeasures
calculate_pd_metric_rasters<-function(occurrences,phylogeny,traits,template.raster){

  #1) create trait-scaled phylogeny  
  trait_phylo<-scale_branches_multidimensional(tree = phylogeny,traits = traits)
  rate_phylo<-scale_branches_multidimensional(tree=phylogeny,traits = traits,rate = T)
  time_phylo<-phylogeny
  
  #2) Remove species from phylo/occurrences that aren't in both
  
#unlist(base::union(x = trait_phylo$tip.label,y = occurrences$current_species))
#setdiff(y = trait_phylo$tip.label,x = occurrences$current_species)  
  
  trait_phylo_list<-data_matching(phylogeny = trait_phylo,occurrences = occurrences)
  trait_phylo<-trait_phylo_list$phylogeny
  occurrences<-trait_phylo_list$occurrences
  rm(trait_phylo_list)
  
  time_phylo<- data_matching(phylogeny = time_phylo,occurrences = occurrences)$phylogeny
  
  rate_phylo<- data_matching(phylogeny = rate_phylo,occurrences = occurrences)$phylogeny
  
  #3) Iterate through template raster
  
  pd_time<-NULL
  pd_traits<-NULL
  pdi_time<-NULL
  pdi_traits<-NULL
  
  for(i in 1:length(unique(occurrences[,2]))){
  spp_to_include<-occurrences[,1][which(occurrences[,2]==unique(occurrences[,2])[i])]
  print(i)
  if(length(spp_to_include)<2){
    pd_time<-c(pd_time,NA)
    pd_traits<-c(pd_traits,NA)
    pdi_time<-c(pdi_time,NA)
    pdi_traits<-c(pdi_traits,NA)
    
    
    
  }else{
  matrix_i<-matrix(nrow=1,ncol=length(time_phylo$tip.label),data = 0)
  colnames(matrix_i)<-time_phylo$tip.label
  matrix_i[which(colnames(matrix_i)%in%spp_to_include)]<-1
   
  pd_time<-c(pd_time,PhyloMeasures::pd.query(tree = time_phylo,matrix = matrix_i,standardize = F))  
  pdi_time<-c(pdi_time,PhyloMeasures::pd.query(tree = time_phylo,matrix = matrix_i,standardize = T) ) 
  pd_traits<-c(pd_traits,PhyloMeasures::pd.query(tree = trait_phylo,matrix = matrix_i,standardize = F))  
  pdi_traits<-c(pdi_traits,PhyloMeasures::pd.query(tree = trait_phylo,matrix = matrix_i,standardize = T)) 
    
    
  }#else
  }#i loop
  
  pd_time_vals<-rep(x = NA,ncell(template.raster))
  pd_traits_vals<-rep(x = NA,ncell(template.raster))
  pdi_time_vals<-rep(x = NA,ncell(template.raster))
  pdi_traits_vals<-rep(x = NA,ncell(template.raster))
  
  
  pd_time_vals[as.numeric(unique(occurrences[,2]))]<-pd_time
  pd_traits_vals[as.numeric(unique(occurrences[,2]))]<-pd_traits
  pdi_time_vals[as.numeric(unique(occurrences[,2]))]<-pdi_time
  pdi_traits_vals[as.numeric(unique(occurrences[,2]))]<-pdi_traits
  
  pd_time_raster<-setValues(x = template.raster,values = pd_time_vals)
  pdi_time_raster<-setValues(x = template.raster,values = pdi_time_vals)
  pd_traits_raster<-setValues(x = template.raster,values = pd_traits_vals)
  pdi_traits_raster<-setValues(x = template.raster,values = pdi_traits_vals)
    
  pd_stack<-stack(pd_time_raster,pdi_time_raster,pd_traits_raster,pdi_traits_raster)  
  names(pd_stack)<-c("pd_time","pdi_time","pd_traits","pdi_traits")
  
  
  output_list<- list()
  output_list[[1]]<-pd_stack
  output_list[[2]]<-time_phylo
  output_list[[3]]<-trait_phylo
  output_list[[4]]<-rate_phylo
  names(output_list)<-c("pd_stack","time_phylo","trait_phylo","rate_phylo")
  
  return(output_list)
  
}

####################


