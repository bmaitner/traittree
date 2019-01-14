#'Generate a species richness raster
#'
#'richness_raster
#' @param occurrences A set of occurrences in "tidy" format: first column is species name, second is raster cells where the species occurs.
#' @param traits a set of trait data where the first column is species name and additional columns are trait data
#' @param template.raster The raster that corresponds to the cell numbers in the second column of the occurrences file
#' @return List containing 4 rasters and 3 phylogenies (original, trait-scaled and rate-scaled)
#' @examples \dontrun{
#' Write example text
#' }
#' @export
richness_raster<-function(template.raster,occurrences){
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


