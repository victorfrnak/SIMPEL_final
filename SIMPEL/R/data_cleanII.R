#'Pull out the second field which correpsonds to categroy of sample in the non-stationary labeling experiment
data_cleanII <- function(x) sapply (strsplit(x , '[_]' ), `[` , 2)
