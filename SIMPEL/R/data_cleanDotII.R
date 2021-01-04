data_cleanDotII <- function(x) sapply (strsplit(x , '[.]' ), `[` , 2)
