# Reference implementations: the original loop-based versions of the functions
# that were vectorized, kept so the replacements can be checked against them.
# Taken verbatim from the pre-vectorization source, renamed with a
# `reference_` prefix and with roxygen comments stripped.

reference_scale_branches_by_traits_fastAnc<-function(tree,traits,percent=FALSE){
  
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

reference_scale_branches_by_traits_rphylopars<-function(tree,traits,percent=FALSE){

invisible("this function should be combined with the other 'scale_branches_rphylopars' functions to eliminate redundancy")  
  
output_trees<-list()    
class(output_trees)<-"multiPhylo"
#First, remove species from trait data that aren't in the phylogeny:
traits<-traits[which(traits$species%in%tree$tip.label),]

#Next, do ancestral state reconstruction on all traits at once using BM with Rphylopars  
anc_recon<-Rphylopars::phylopars(trait_data = traits,tree = tree,pheno_error = resolve_pheno_error(traits))$anc_recon
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

reference_richness_raster<-function(template.raster,occurrences){
#richness<-NULL
#for( i in 1:ncell(template.raster)){
#  print(i)
#  richness<-c(richness,length(which(occurrencses[,2]==i)))  
#}

#richness<-setValues(x = template.raster,values = richness)

  output_raster<-template.raster
  output_raster<-raster::setValues(x = output_raster,values = NA)
  
  #iterate through all cells with at least one occurrence, record 
  
  output_raster[as.numeric(unique(occurrences[,2]))] <- sapply(X = unique(occurrences[,2]),FUN = function(x){ length(unique(occurrences[which(occurrences[,2]==x),1]))} )

  return(output_raster)

}


# Exact joint posterior of the internal node states under Brownian motion with
# an improper prior on the root, built densely from the tree's shared path
# lengths.  Quadratic in tree size, so only usable on small trees, but it is an
# independent check on the linear-time sampler used in the package.
# Univariate: `x_tips` is a named vector, `R` the scalar rate.

exact_joint_posterior <- function(tree, x_tips, R){

  n <- length(tree$tip.label)
  n_nodes <- 2 * n - 1

  # shared root-to-MRCA path length between every pair of nodes
  depth <- ape::node.depth.edgelength(tree)
  C <- matrix(depth[ape::mrca(tree, full = TRUE)], nrow = n_nodes, ncol = n_nodes)

  tips <- 1:n
  internal <- (n + 1):n_nodes

  C_tips_inv <- solve(C[tips, tips])
  C_internal_tips <- C[internal, tips, drop = FALSE]

  ones_tips <- rep(1, n)
  ones_internal <- rep(1, length(internal))

  # improper root prior: root estimated by GLS, with its uncertainty carried
  # into the conditional covariance
  denominator <- as.numeric(t(ones_tips) %*% C_tips_inv %*% ones_tips)
  root <- as.numeric(t(ones_tips) %*% C_tips_inv %*% x_tips[tree$tip.label])/denominator

  weights <- C_internal_tips %*% C_tips_inv
  residual <- ones_internal - as.numeric(weights %*% ones_tips)

  list(mean = as.numeric(root + weights %*% (x_tips[tree$tip.label] - ones_tips * root)),
       cov = (C[internal, internal] -
                weights %*% t(C_internal_tips) +
                outer(residual, residual)/denominator) * R)

}


# Original loop-based implementation of calculate_pd_metric_rasters(), with the
# progress print removed, kept to check the batched rewrite against.

reference_calculate_pd_metric_rasters<-function(occurrences,phylogeny,traits,template.raster){

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
