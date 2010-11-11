require(rtm)
require(lda)

data(cora.documents)
data(cora.cites)

my.corpus <- corpus(documents=cora.documents)

stopifnot(my.corpus$getWordCount(1L) == 879)
stopifnot(my.corpus$getDocumentCount() == length(cora.documents))
stopifnot(my.corpus$getWord(my.corpus$getIndex(3, 1)) == cora.documents[[4]][1,2])

my.links <- links(list = cora.cites)
stopifnot(my.links$getDocumentCount() == length(cora.cites))
stopifnot(my.links$getNumLinks(0L) == length(cora.cites[[1]]))
stopifnot(my.links$getTarget(3L, 1L) == cora.cites[[4]][2])
