rtm <- function(corpus,
                links,
                K,
                alpha = 0.1,
                eta = 0.1,
                beta = 3,
                initial = NULL) {
  beta <- rep(beta, length.out = K)
  if (length(K) != 1 || !is.integer(K)) {
    stop("K must be a length 1 integer.")
  }

  ## Check corpus class.
  corpus.class <- get(".cppclass", envir = as.environment(corpus))
  if (.Call("Class__name", corpus.class) != "Corpus") {
    stop("corpus must be of class 'Corpus'.")
  }
  corpus.pointer <- get(".pointer", envir = as.environment(corpus))

  ## Check links class.
  links.class <- get(".cppclass", envir = as.environment(links))
  if (.Call("Class__name", links.class) != "Links") {
    stop("links must be of class 'Links'.")
  }
  links.pointer <- get(".pointer", envir = as.environment(links))
  
  tt <- new(.module$SparseRTM)
  tt$load(corpus.pointer, links.pointer, alpha, eta, beta, K)
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
