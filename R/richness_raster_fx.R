#'Generate a species richness raster
#'
#'richness_raster
#' @param occurrences A set of occurrences in "tidy" format: first column is species name, second is raster cells where the species occurs.
#' @param template.raster The raster that corresponds to the cell numbers in the second column of the occurrences file
#' @return List containing 4 rasters and 3 phylogenies (original, trait-scaled and rate-scaled)
#' @examples \dontrun{
#' Write example text
#' }
#' @export
richness_raster<-function(template.raster,occurrences){

  output_raster<-template.raster
  output_raster<-raster::setValues(x = output_raster,values = NA)
  
  #Count the distinct species in each occupied cell.  Scanning the occurrences
  #once per cell is quadratic in the size of the problem, so instead drop
  #repeated species-cell pairs and tally the remaining rows in a single pass.
  
  species<-as.character(occurrences[,1])
  cells<-as.character(occurrences[,2])
  
  distinct_pairs<-!duplicated(data.frame(cells,species))
  
  richness<-table(cells[distinct_pairs])
  
  output_raster[as.numeric(names(richness))]<-as.numeric(richness)

  return(output_raster)

}
