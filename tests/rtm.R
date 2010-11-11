require(rtm)
require(lda)

data(cora.documents)
data(cora.cites)
data(cora.vocab)
set.seed(8675309)

my.corpus <- corpus(documents = cora.documents)
my.links <- links(list = cora.cites)
model <- rtm(my.corpus, my.links, 10L)
system.time(model$iterateCorpus(25L))

top <- model$getTopics()
ass <- model$getAssignments()
doc <- model$getDocumentSums()
top.sum <- model$getTopicSums()
rownames(top) <- cora.vocab

top.words <- top.topic.words(t(top), by.score = T)

## system.time(lda.collapsed.gibbs.sampler(cora.documents,
##                                         100L,
##                                         cora.vocab,
##                                         25L,
##                                         0.1,
##                                         0.1,
##                                         trace = 2L))
