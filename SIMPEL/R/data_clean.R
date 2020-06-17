#'Pull out the first field which correpsonds to timepoint of non-stationary labeling experiment
data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)
