require(rtm)
require(lda)

data(cora.documents)
data(cora.cites)
data(cora.vocab)
set.seed(8675309)

my.corpus <- corpus(documents = cora.documents)
model <- hidden.markov.lda(my.corpus, 10L)
system.time(model$iterateCorpus(25L))

top <- model$getTopics()
ass <- model$getAssignments()
transitions <- model$getTransitionCounts()
top.sum <- model$getTopicSums()
rownames(top) <- cora.vocab

print(transitions)
top.words <- top.topic.words(t(top), by.score = T)
print(top.words)

## system.time(lda.collapsed.gibbs.sampler(cora.documents,
##                                         100L,
##                                         cora.vocab,
##                                         25L,
##                                         0.1,
##                                         0.1,
##                                         trace = 2L))
