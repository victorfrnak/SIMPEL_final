data_cleanDot <- function(x) sapply (strsplit(x , '[.]' ), `[` , 1)
