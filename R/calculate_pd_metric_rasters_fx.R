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
#' @param verbose Logical. If TRUE, report what is being computed. Defaults to FALSE.
#' @param rarefy_to Optional richness to rarefy every cell to before standardizing, which removes the dependence of the index on species richness (see Details). Either a whole number of at least 2, or \code{"min"} for the smallest richness among the reported cells with two or more species. Defaults to NULL, meaning no rarefied index is computed.
#' @return A list with \code{pd_stack}, \code{time_phylo}, \code{trait_phylo} and \code{rate_phylo}. With a single trait-scaled phylogeny, \code{pd_stack} holds \code{pd_time}, \code{pdi_time}, \code{pd_traits} and \code{pdi_traits}. With several, the trait-based layers are replaced by their mean, standard deviation and the two requested quantiles across trees, and \code{trait_phylo} is the multiPhylo that produced them. With \code{rarefy_to}, the rarefied index is appended as \code{pdi_time_rarefied} and \code{pdi_traits_rarefied} (summarised across trees in the same way), and the richness used is returned as \code{rarefied_to}.
#' @note The time-based metrics do not depend on the traits, so they are computed once and are unaffected by \code{n_trees}.
#' @details Phylogenetic diversity (PD) is the total branch length of the minimal subtree connecting the species in a cell; the path from that subtree to the root is not included. The phylogenetic diversity index (PDI) standardizes it: the mean PD of every set of the same number of species drawn from the phylogeny is subtracted, and the result divided by the standard deviation across those sets (Sandel 2018). Both moments are computed exactly rather than by randomization (Nipperess and Matsen 2013), and the species pool is every species left in the phylogeny once those without occurrences are removed. Negative values indicate species more closely related than expected by chance. A cell holding a single species has PD and PDI of zero.
#'
#' The standardization guarantees only that randomly assembled cells score zero whatever their richness. It does not make other values comparable between cells of different richness: under the same degree of phylogenetic clustering, richer cells receive more extreme values (Sandel 2018). Maps of PDI across cells whose richness varies should be read with that in mind. \code{rarefy_to} adds an index without this dependence. Each cell is rarefied to a common richness and compared with the null distribution for that richness; rather than drawing random subsets, the average over every possible rarefaction is computed exactly, as the expected PD of a random subset of the species in the cell, less the null mean, over the null standard deviation. Cells with fewer species than the target are NA. Rarefaction discards information, so the rarefied index is noisier than the plain one, increasingly so as the target falls; a very small target, which \code{"min"} can produce, leaves little power to tell cells apart.
#' @seealso \code{\link{phylo_diversity}}, which does the same calculation for any set of communities, without rasters.
#' @references Nipperess, D. A. and Matsen, F. A. (2013). The mean and variance of phylogenetic diversity under rarefaction. Methods in Ecology and Evolution 4: 566-572. \doi{10.1111/2041-210X.12042}
#'
#' Sandel, B. (2018). Richness-dependence of phylogenetic diversity indices. Ecography 41: 837-844. \doi{10.1111/ecog.02967}
#' @examples
#' template <- raster::raster(nrows = 10, ncols = 10)
#'
#' # Point estimate: one trait-scaled phylogeny
#' result <- calculate_pd_metric_rasters(occurrences = example_occurrences,
#'                                       phylogeny = example_tree,
#'                                       traits = example_traits,
#'                                       template.raster = template)
#'
#' names(result$pd_stack)
#'
#' # Carry the reconstruction uncertainty into the metrics instead
#' sampled <- calculate_pd_metric_rasters(occurrences = example_occurrences,
#'                                        phylogeny = example_tree,
#'                                        traits = example_traits,
#'                                        template.raster = template,
#'                                        n_trees = 10)
#'
#' names(sampled$pd_stack)
#'
#' # The spread across trees is the uncertainty in trait-based PD
#' summary(raster::values(sampled$pd_stack[["pd_traits_sd"]]))
#' @export
#' @import raster
calculate_pd_metric_rasters<-function(occurrences,
                                      phylogeny,
                                      traits=NULL,
                                      template.raster,
                                      trait_phylo=NULL,
                                      n_trees=1,
                                      probs=c(0.025,0.975),
                                      verbose=FALSE,
                                      rarefy_to=NULL){

  if(!is.numeric(n_trees) || length(n_trees) != 1 || is.na(n_trees) || n_trees < 1){

    stop("n_trees must be a single positive number.", call. = FALSE)

  }

  n_trees <- as.integer(n_trees)

  if(!is.numeric(probs) || length(probs) != 2 || anyNA(probs) ||
     any(probs < 0) || any(probs > 1)){

    stop("probs must be two probabilities between 0 and 1.", call. = FALSE)

  }

  if(!is.null(rarefy_to) && !identical(rarefy_to,"min") &&
     !(is.numeric(rarefy_to) && length(rarefy_to) == 1 && !is.na(rarefy_to) &&
       rarefy_to >= 2 && rarefy_to == round(rarefy_to))){

    stop("rarefy_to must be NULL, \"min\", or a single whole number of at least 2.",
         call. = FALSE)

  }

  if(is.null(trait_phylo) && is.null(traits)){

    stop("Supply either traits, to scale the phylogeny here, or trait_phylo, ",
         "an already scaled phylogeny.", call. = FALSE)

  }

  #A supplied trait_phylo is validated here, with the other arguments, rather
  #than where it is used, so that a bad call fails before any of the
  #reconstruction work is done.
  if(!is.null(trait_phylo)){

    if(n_trees != 1){

      stop("n_trees applies only when the phylogeny is scaled here; ",
           "trait_phylo already fixes the trees to use.", call. = FALSE)

    }

    supplied_trees<-if(inherits(trait_phylo,"phylo")){

      list(trait_phylo)

    }else if(inherits(trait_phylo,"multiPhylo")){

      #Shared tip labels may live in an attribute rather than in each tree.
      unclass(ape::.uncompressTipLabel(trait_phylo))

    }else{

      trait_phylo

    }

    if(!is.list(supplied_trees) || length(supplied_trees) == 0 ||
       !all(vapply(supplied_trees,inherits,logical(1),"phylo"))){

      stop("trait_phylo must be a phylo object or a multiPhylo of them.", call. = FALSE)

    }

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

      trait_trees<-unclass(trait_trees)

    }

  }else{

    trait_trees<-supplied_trees

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

  #The tips are the species pool, and the pool defines the null, so every
  #tree has to have the same ones once pruned. A supplied trait tree missing
  #a species would otherwise be measured against a different null, silently.
  pd_check_trees(c(list(time_phylo),trait_trees))

  if(!is.null(rate_phylo)){

    rate_phylo<- data_matching(phylogeny = rate_phylo,occurrences = occurrences)$phylogeny

  }

  #3) Phylogenetic diversity of every occupied cell

  cell_values<-unique(occurrences[,2])

  cell_index<-match(as.character(occurrences[,2]),as.character(cell_values))

  #As before, cells with fewer than two occurrence records are not reported.
  records_per_cell<-tabulate(cell_index,nbins = length(cell_values))

  cells_to_query<-which(records_per_cell >= 2)

  if(verbose){

    message("Computing phylogenetic diversity for ",length(cells_to_query),
            " cells across ",length(trait_trees),
            if(length(trait_trees) == 1) " tree" else " trees")

  }

  #Settle the richness to rarefy to, if any. "min" is the smallest richness
  #among the cells that will be reported and that have a subtree to measure.
  rarefied_to<-NULL

  if(!is.null(rarefy_to)){

    if(identical(rarefy_to,"min")){

      species<-as.character(occurrences[,1])

      counted<-species %in% time_phylo$tip.label &
        !duplicated(data.frame(cell_index,species))

      species_per_cell<-tabulate(cell_index[counted],nbins = length(cell_values))

      eligible<-species_per_cell[records_per_cell >= 2 & species_per_cell >= 2]

      rarefied_to<-if(length(eligible) > 0) min(eligible) else 2

    }else{

      rarefied_to<-as.integer(rarefy_to)

    }

    if(rarefied_to > length(time_phylo$tip.label)){

      stop("rarefy_to (",rarefied_to,") is larger than the number of species in the ",
           "phylogeny with occurrences (",length(time_phylo$tip.label),").",
           call. = FALSE)

    }

  }

  #The time-scaled tree and the trait-scaled trees differ only in their branch
  #lengths, so they are measured together: everything that depends on topology
  #is worked out once and shared.
  diversity<-pd_by_community(trees = c(list(time_phylo),trait_trees),
                             community = cell_index,
                             species = as.character(occurrences[,1]),
                             n_communities = length(cell_values),
                             rarefy_to = rarefied_to)

  not_reported<-records_per_cell < 2

  diversity$pd[not_reported,]<-NA
  diversity$pdi[not_reported,]<-NA

  #4) Turn the per-cell values into rasters

  to_raster<-function(values){

    all_cells<-rep(x = NA_real_,ncell(template.raster))
    all_cells[as.numeric(cell_values)]<-values

    setValues(x = template.raster,values = all_cells)

  }

  #The trait-based layers for one metric: the values themselves for a single
  #tree, or their summary across a posterior sample of trees. Cells without
  #values stay NA rather than being summarised.
  quantile_names<-paste0("q",format(100 * probs,trim = TRUE))

  trait_layers<-function(values,name){

    if(ncol(values) == 1){

      return(stats::setNames(list(to_raster(values[,1])),name))

    }

    complete<-which(!is.na(rowSums(values)))

    summarise_cells<-function(summary_function){

      summary<-rep(NA_real_,nrow(values))

      if(length(complete) > 0){

        summary[complete]<-apply(X = values[complete,,drop = FALSE],
                                 MARGIN = 1,
                                 FUN = summary_function)

      }

      to_raster(summary)

    }

    out<-list(summarise_cells(mean),
              summarise_cells(stats::sd),
              summarise_cells(function(draws) stats::quantile(draws,probs = probs[1],names = FALSE)),
              summarise_cells(function(draws) stats::quantile(draws,probs = probs[2],names = FALSE)))

    stats::setNames(out,paste0(name,"_",c("mean","sd",quantile_names)))

  }

  layers<-c(list(pd_time = to_raster(diversity$pd[,1]),
                 pdi_time = to_raster(diversity$pdi[,1])),
            trait_layers(diversity$pd[,-1,drop = FALSE],"pd_traits"),
            trait_layers(diversity$pdi[,-1,drop = FALSE],"pdi_traits"))

  if(!is.null(rarefied_to)){

    diversity$rarefied_pdi[not_reported,]<-NA

    if(all(is.na(diversity$rarefied_pdi))){

      warning("No cell has ",rarefied_to," or more species, so the rarefied index is NA ",
              "everywhere. Choose a smaller rarefy_to.",call. = FALSE)

    }

    layers<-c(layers,
              list(pdi_time_rarefied = to_raster(diversity$rarefied_pdi[,1])),
              trait_layers(diversity$rarefied_pdi[,-1,drop = FALSE],"pdi_traits_rarefied"))

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

  if(!is.null(rarefied_to)){

    output_list$rarefied_to<-rarefied_to

  }

  return(output_list)

}

####################
