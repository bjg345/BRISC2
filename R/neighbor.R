BRISC2_neighbor <- function(coords, n.neighbors = 15, n_omp = 1,
                             order = "Sum_coords", search.type = "tree", verbose = TRUE, ordering = NULL, tol = 12
){
  n <- nrow(coords)
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(order == "AMMD" && nrow(coords) < 65){stop("error: Number of data points must be atleast 65 to use AMMD")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or the number of rows is different than
         data used in the model formula")
  }

  coords <- round(coords, tol)


  if(is.null(ordering)){
    if(order != "AMMD" && order != "Sum_coords"){
      stop("error: Please insert a valid ordering scheme choice given by 1 or 2.")
    }
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tOrdering Coordinates"), collapse="   "), "\n")
    }
    if(order == "AMMD"){ord <- orderMaxMinLocal(coords)}
    if(order == "Sum_coords"){ord <- order(coords[,1] + coords[,2])}
  }
  if(!is.null(ordering)){
    ord <- ordering
  }
  coords <- coords[ord,]

  ##Search type
  search.type.names <- c("brute", "tree", "cb")
  if(!search.type %in% search.type.names){
    stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
  }
  search.type.indx <- which(search.type == search.type.names)-1
  storage.mode(search.type.indx) <- "integer"

  ##Option for Multithreading if compiled with OpenMp support
  n.omp.threads <- as.integer(n_omp)
  storage.mode(n.omp.threads) <- "integer"


  ##type conversion
  storage.mode(n) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(verbose) <- "integer"



  p1<- proc.time()

  ##estimtion
  result <- .Call("BRISC2_neighborcpp", n, n.neighbors, coords, search.type.indx, n.omp.threads, verbose, PACKAGE = "BRISC2")

  p2 <- proc.time()

  result
}
