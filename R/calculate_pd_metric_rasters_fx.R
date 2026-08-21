#'Calculate time- and trait-scaled phylogenies and pd/pdi rasters
#'
#'calculate_pd_metric_rasters
#' @param occurrences A set of occurrences in "tidy" format: first column is species name, second is raster cells where the species occurs.
#' @param phylogeny A phylogeny with time-scaled branch lengths
#' @param traits a set of trait data where the first column is species name and additional columns are trait data. May be NULL if \code{trait_phylo} is supplied, in which case no rate-scaled phylogeny is returned.
#' @param template.raster The raster that corresponds to the cell numbers in the second column of the occurrences file
#' @param trait_phylo Optional trait-scaled phylogeny, or a multiPhylo of them, to use instead of scaling the phylogeny internally. Defaults to NULL, meaning scale it here.
#' @param n_trees Number of trait-scaled phylogenies to draw from the joint posterior of the ancestral states, so that uncertainty in the reconstruction is carried into the metrics. Defaults to 1, which uses the point estimate and reproduces the previous behaviour. Ignored when \code{trait_phylo} is supplied.
#' @param probs The two quantiles summarising the trait-based metrics across trees. Defaults to the central 95 percent.
#' @param verbose Logical. If TRUE, report progress while querying cells. Defaults to FALSE.
#' @return A list with \code{pd_stack}, \code{time_phylo}, \code{trait_phylo} and \code{rate_phylo}. With a single trait-scaled phylogeny, \code{pd_stack} holds \code{pd_time}, \code{pdi_time}, \code{pd_traits} and \code{pdi_traits}. With several, the trait-based layers are replaced by their mean, standard deviation and the two requested quantiles across trees, and \code{trait_phylo} is the multiPhylo that produced them.
#' @note The time-based metrics do not depend on the traits, so they are computed once and are unaffected by \code{n_trees}.
#' @examples \dontrun{
#' Write example text
#' }
#' @export
#' @import raster
#' @import PhyloMeasures
calculate_pd_metric_rasters<-function(occurrences,
                                      phylogeny,
                                      traits=NULL,
                                      template.raster,
                                      trait_phylo=NULL,
                                      n_trees=1,
                                      probs=c(0.025,0.975),
                                      verbose=FALSE){

  if(!is.numeric(n_trees) || length(n_trees) != 1 || is.na(n_trees) || n_trees < 1){

    stop("n_trees must be a single positive number.", call. = FALSE)

  }

  n_trees <- as.integer(n_trees)

  if(!is.numeric(probs) || length(probs) != 2 || anyNA(probs) ||
     any(probs < 0) || any(probs > 1)){

    stop("probs must be two probabilities between 0 and 1.", call. = FALSE)

  }

  if(is.null(trait_phylo) && is.null(traits)){

    stop("Supply either traits, to scale the phylogeny here, or trait_phylo, ",
         "an already scaled phylogeny.", call. = FALSE)

  }

  #1) obtain the trait-scaled phylogenies

  if(is.null(trait_phylo)){

    if(n_trees == 1){

      trait_trees<-list(scale_branches_multidimensional(tree = phylogeny,traits = traits))

    }else{

      #Each tree is a draw from the joint posterior of the ancestral states, so
      #the spread across them is the reconstruction uncertainty.
      trait_trees<-scale_branches_multidimensional_with_variation(tree = phylogeny,
                                                                  traits = traits,
                                                                  n_trees = n_trees)

    }

  }else{

    if(n_trees != 1){

      stop("n_trees applies only when the phylogeny is scaled here; ",
           "trait_phylo already fixes the trees to use.", call. = FALSE)

    }

    trait_trees<-trait_phylo

  }

  if(inherits(trait_trees,"phylo")){

    trait_trees<-list(trait_trees)

  }

  trait_trees<-unclass(trait_trees)

  if(length(trait_trees) == 0 || !all(vapply(trait_trees,inherits,logical(1),"phylo"))){

    stop("trait_phylo must be a phylo object or a multiPhylo of them.", call. = FALSE)

  }

  #The rate-scaled phylogeny is returned for the caller's own use rather than
  #rasterized here, and needs the traits to build.
  if(!is.null(traits)){

    rate_phylo<-scale_branches_multidimensional(tree = phylogeny,traits = traits,rate = TRUE)

  }else{

    rate_phylo<-NULL

  }

  time_phylo<-phylogeny

  #2) Remove species from phylo/occurrences that aren't in both

  matched<-data_matching(phylogeny = trait_trees[[1]],occurrences = occurrences)

  occurrences<-matched$occurrences

  retained_species<-matched$phylogeny$tip.label

  #Every trait tree carries the same tips, so prune them all to that set.
  trait_trees<-lapply(X = trait_trees,
                      FUN = function(tree_i){

                        to_drop<-setdiff(tree_i$tip.label,retained_species)

                        if(length(to_drop) == 0){

                          return(tree_i)

                        }

                        ape::drop.tip(tree_i,to_drop)

                      })

  time_phylo<- data_matching(phylogeny = time_phylo,occurrences = occurrences)$phylogeny

  if(!is.null(rate_phylo)){

    rate_phylo<- data_matching(phylogeny = rate_phylo,occurrences = occurrences)$phylogeny

  }

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

  #One row per trait tree, one column per cell.
  pd_traits<-matrix(data = NA_real_,nrow = length(trait_trees),ncol = length(cell_values))
  pdi_traits<-pd_traits

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

      #The community matrix does not depend on the tree, so build it once per
      #block and reuse it for every tree.
      pd_time[block]<-PhyloMeasures::pd.query(tree = time_phylo,matrix = community_matrix,standardize = FALSE)
      pdi_time[block]<-PhyloMeasures::pd.query(tree = time_phylo,matrix = community_matrix,standardize = TRUE)

      for(tree_index in seq_along(trait_trees)){

        pd_traits[tree_index,block]<-PhyloMeasures::pd.query(tree = trait_trees[[tree_index]],matrix = community_matrix,standardize = FALSE)
        pdi_traits[tree_index,block]<-PhyloMeasures::pd.query(tree = trait_trees[[tree_index]],matrix = community_matrix,standardize = TRUE)

      }

      if(verbose){

        message("Queried ",min(block_start + cells_per_block - 1,length(cells_to_query)),
                " of ",length(cells_to_query)," cells across ",length(trait_trees),
                if(length(trait_trees) == 1) " tree" else " trees")

      }

    }#block loop

  }

  #4) Turn the per-cell values into rasters

  to_raster<-function(values){

    all_cells<-rep(x = NA_real_,ncell(template.raster))
    all_cells[as.numeric(cell_values)]<-values

    setValues(x = template.raster,values = all_cells)

  }

  layers<-list(pd_time = to_raster(pd_time),
               pdi_time = to_raster(pdi_time))

  if(length(trait_trees) == 1){

    layers$pd_traits<-to_raster(pd_traits[1,])
    layers$pdi_traits<-to_raster(pdi_traits[1,])

  }else{

    #Summarise the trait-based metrics across the posterior sample.  Only the
    #queried cells hold draws; the rest stay NA rather than being summarised.
    quantile_names<-paste0("q",format(100 * probs,trim = TRUE))

    summarise_cells<-function(draws,summary_function,...){

      values<-rep(NA_real_,ncol(draws))

      if(length(cells_to_query) > 0){

        values[cells_to_query]<-apply(X = draws[,cells_to_query,drop = FALSE],
                                      MARGIN = 2,
                                      FUN = summary_function,
                                      ...)

      }

      values

    }

    for(metric in c("pd_traits","pdi_traits")){

      draws<-get(metric)

      layers[[paste0(metric,"_mean")]]<-to_raster(summarise_cells(draws,mean))
      layers[[paste0(metric,"_sd")]]<-to_raster(summarise_cells(draws,stats::sd))

      for(bound in seq_along(probs)){

        layers[[paste0(metric,"_",quantile_names[bound])]]<-
          to_raster(summarise_cells(draws,
                                    function(column) stats::quantile(column,
                                                                     probs = probs[bound],
                                                                     names = FALSE)))

      }

    }

  }

  pd_stack<-stack(layers)
  names(pd_stack)<-names(layers)

  if(length(trait_trees) == 1){

    returned_trait_phylo<-trait_trees[[1]]

  }else{

    returned_trait_phylo<-trait_trees
    class(returned_trait_phylo)<-"multiPhylo"

  }

  output_list<- list(pd_stack = pd_stack,
                     time_phylo = time_phylo,
                     trait_phylo = returned_trait_phylo,
                     rate_phylo = rate_phylo)

  return(output_list)

}

####################
