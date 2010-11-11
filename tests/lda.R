require(rtm)
require(lda)
require(jjplot)
require(reshape)

data(cora.documents)
data(cora.vocab)
set.seed(8675309)

my.corpus <- corpus(documents = cora.documents)
model <- lda(my.corpus, 10L)

system.time(model$iterateCorpus(25L))

top <- model$getTopics()
ass <- model$getAssignments()
doc <- model$getDocumentSums()
top.sum <- model$getTopicSums()
rownames(top) <- cora.vocab

jjplot.stat.head <- function(state, n, sort,
                             decreasing = TRUE) {
  if (!missing(sort)) {
    state$data <- head(state$data[order(state$data[[sort]],
                                        decreasing = decreasing),],
                       n = n)
  } else {
    state$data <- head(state$data, n = n)
  }  
}

top.words <- top.topic.words(t(top), by.score = T)

jjplot(X1 ~ text(label=value) + X1, melt(top.words),
       facet.x = factor(X2), facet.nrow = 2)

## system.time(lda.collapsed.gibbs.sampler(cora.documents,
##                                         100L,
##                                         cora.vocab,
##                                         25L,
##                                         0.1,
##                                         0.1,
##                                         trace = 2L))
