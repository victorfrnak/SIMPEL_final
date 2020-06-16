#'This function is calculating the number of each of the elements present in the formula
#'
#' @param comp_formula is the formula from the annotation file
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' get_element_count()


get_element_count <- function(comp_formula){
  comp_formula %>%
    makeup() %>%
    as.list() %>%
    return()
}
