#pull out the ending field which correpsonds to replicate
data_cleanIII <- function(x) sapply (strsplit(x , '[_]' ), `[` , 3)

