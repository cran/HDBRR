#' @rdname predict
#' @export

predict.HDBRR <- function(object,...){
  if(!inherits(object, "HDBRR")) stop("This function only works for objects of class 'HDBRR'\n");
  pred <- object$yhat
  predic <- as.vector(round(pred,5))
  return(predic)
}
