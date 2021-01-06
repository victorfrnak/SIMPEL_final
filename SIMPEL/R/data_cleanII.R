#'Pull out the second field which correpsonds to categroy of sample in the non-stationary labeling experiment
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' data_cleanII()
data_cleanII <- function(x) sapply (strsplit(x , '[_]' ), `[` , 2)
