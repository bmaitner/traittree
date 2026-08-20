#'Calculate time- and trait-scaled phylogenies and pd/pdi rasters
#'
#'calculate_pd_metric_rasters
#' @param occurrences A set of occurrences in "tidy" format: first column is species name, second is raster cells where the species occurs.
#' @param phylogeny A phylogeny with time-scaled branch lengths
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param template.raster The raster that corresponds to the cell numbers in the second column of the occurrences file
#' @param verbose Logical. If TRUE, report progress while querying cells. Defaults to FALSE.
#' @return List containing 4 rasters and 3 phylogenies (original, trait-scaled and rate-scaled)
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import raster
#' @import PhyloMeasures
calculate_pd_metric_rasters<-function(occurrences,phylogeny,traits,template.raster,verbose=FALSE){

  #1) create trait-scaled phylogeny
  trait_phylo<-scale_branches_multidimensional(tree = phylogeny,traits = traits)
  rate_phylo<-scale_branches_multidimensional(tree=phylogeny,traits = traits,rate = TRUE)
  time_phylo<-phylogeny

  #2) Remove species from phylo/occurrences that aren't in both

  trait_phylo_list<-data_matching(phylogeny = trait_phylo,occurrences = occurrences)
  trait_phylo<-trait_phylo_list$phylogeny
  occurrences<-trait_phylo_list$occurrences
  rm(trait_phylo_list)

  time_phylo<- data_matching(phylogeny = time_phylo,occurrences = occurrences)$phylogeny

  rate_phylo<- data_matching(phylogeny = rate_phylo,occurrences = occurrences)$phylogeny

  #3) Query every occupied cell

  #Group the occurrences by cell in a single pass.  Scanning the whole
  #occurrence table once per cell is quadratic in the size of the problem.
  cell_values<-unique(occurrences[,2])

  cell_key<-factor(x = as.character(occurrences[,2]),
                   levels = as.character(cell_values))

  tip_labels<-time_phylo$tip.label

  #Columns of the presence/absence matrix that each occurrence corresponds to.
  species_columns<-match(as.character(occurrences[,1]),tip_labels)

  columns_by_cell<-split(x = species_columns,f = cell_key)

  #As before, cells with fewer than two occurrence records are not queried.
  records_per_cell<-as.numeric(table(cell_key))

  cells_to_query<-which(records_per_cell >= 2)

  pd_time<-rep(NA_real_,length(cell_values))
  pdi_time<-pd_time
  pd_traits<-pd_time
  pdi_traits<-pd_time

  #PhyloMeasures accepts many communities at once, so query in blocks rather
  #than one cell at a time.  The block size keeps the presence/absence matrix to
  #a few million entries regardless of how many species there are.
  cells_per_block<-max(1,floor(2e6/length(tip_labels)))

  if(length(cells_to_query) > 0){

    for(block_start in seq(from = 1,to = length(cells_to_query),by = cells_per_block)){

      block<-cells_to_query[block_start:min(block_start + cells_per_block - 1,
                                            length(cells_to_query))]

      community_matrix<-matrix(data = 0,nrow = length(block),ncol = length(tip_labels))
      colnames(community_matrix)<-tip_labels

      for(row in seq_along(block)){

        #Occurrences of species absent from the tree contribute no column, which
        #is how the original name matching behaved.
        columns<-columns_by_cell[[block[row]]]
        columns<-columns[!is.na(columns)]

        if(length(columns) > 0){

          community_matrix[row,columns]<-1

        }

      }

      pd_time[block]<-PhyloMeasures::pd.query(tree = time_phylo,matrix = community_matrix,standardize = FALSE)
      pdi_time[block]<-PhyloMeasures::pd.query(tree = time_phylo,matrix = community_matrix,standardize = TRUE)
      pd_traits[block]<-PhyloMeasures::pd.query(tree = trait_phylo,matrix = community_matrix,standardize = FALSE)
      pdi_traits[block]<-PhyloMeasures::pd.query(tree = trait_phylo,matrix = community_matrix,standardize = TRUE)

      if(verbose){

        message("Queried ",min(block_start + cells_per_block - 1,length(cells_to_query)),
                " of ",length(cells_to_query)," cells")

      }

    }#block loop

  }

  pd_time_vals<-rep(x = NA,ncell(template.raster))
  pd_traits_vals<-rep(x = NA,ncell(template.raster))
  pdi_time_vals<-rep(x = NA,ncell(template.raster))
  pdi_traits_vals<-rep(x = NA,ncell(template.raster))


  pd_time_vals[as.numeric(cell_values)]<-pd_time
  pd_traits_vals[as.numeric(cell_values)]<-pd_traits
  pdi_time_vals[as.numeric(cell_values)]<-pdi_time
  pdi_traits_vals[as.numeric(cell_values)]<-pdi_traits

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
