get_vector_with_proportion <- function(n, p) {
  zeros <- rep(0, n * (1 - p)) ### false positives
  ones <- rep(1, n * p) ### true positives
  redt <- c(zeros, ones)[sample(1:n)] ### randomize them
  redt[is.na(redt)] <- 0 ## substitute occassional NAs by 0s (necessary because of rounding errors)
  redt ### return vector
}