lda <- function(corpus,
                K,
                alpha = 0.1,
                eta = 0.1,
                initial = NULL) {
  beta <- rep(0, length.out = K)
  ll <- links(integer(0), rep(0L, corpus$getDocumentCount()))
  return(rtm(corpus, ll, K, alpha, eta, beta, initial))
}
