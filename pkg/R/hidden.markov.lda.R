hidden.markov.lda <- function(corpus,
                              K,
                              alpha = 0.1,
                              eta = 0.1,
                              initial = NULL) {
  if (length(K) != 1 || !is.integer(K)) {
    stop("K must be a length 1 integer.")
  }

  ## Check corpus class.
  corpus.class <- get(".cppclass", envir = as.environment(corpus))
  if (.Call("Class__name", corpus.class) != "Corpus") {
    stop("corpus must be of class 'Corpus'.")
  }
  corpus.pointer <- get(".pointer", envir = as.environment(corpus))

  tt <- new(.module$HiddenTopicTM)
  tt$load(corpus.pointer, alpha, eta, K)
  if (is.null(initial)) {
    tt$initializeRandom()
  } else {
    if (length(initial) != corpus$getTotalCount()) {
      stop("Incorrect length for 'initial'.")
    }
    
    tt$initialize(initial)
  }
  return(tt)
}
