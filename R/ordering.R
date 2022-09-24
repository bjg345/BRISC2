BRISC2_order <- function(coords, order = "Sum_coords", verbose = TRUE){
  n <- nrow(coords)
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(order == "AMMD" && nrow(coords) < 65){stop("error: Number of data points must be atleast 65 to use AMMD")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or the number of rows is different than
         data used in the model formula")
  }

  coords <- round(coords, 12)

  if(order != "AMMD" && order != "Sum_coords"){
    stop("error: Please insert a valid ordering scheme choice given by 1 or 2.")
  }
  if(verbose == TRUE){
    cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tOrdering Coordinates"), collapse="   "), "\n")
  }
  if(order == "AMMD"){ord <- orderMaxMinLocal(coords)}
  if(order == "Sum_coords"){ord <- order(coords[,1] + coords[,2])}

  return(ord)
}
