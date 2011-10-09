expand.documents <- function(documents) {
  lapply(documents, function(x) rep(x[1,], x[2,]))
}

lda.corpus <- function(documents) {
  documents <- lapply(documents, function(x) {
    new(Document, list(new(DiscreteDocumentData, x)))
  })
  new(Corpus, documents)
}

slda.corpus <- function(documents, annotations) {
}
