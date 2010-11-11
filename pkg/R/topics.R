topics <- function(K, corpus) {
  if (length(K) != 1 || !is.integer(K)) {
    stop("K must be a length 1 integer.")
  }
  # Check corpus class.
  corpus.class <- get(".cppclass", envir = as.environment(corpus))
  if (.Call("Class__name", corpus.class) != "Corpus") {
    stop("corpus must be of class 'Corpus'.")
  }
  tt <- new(.module$Topics)
  corpus.pointer <- get(".pointer", envir = as.environment(corpus))
  tt$load(K, corpus.pointer)
  return(tt);
}
