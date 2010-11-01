require(rtm)
require(lda)

data(cora.documents)

my.corpus <- corpus(documents=cora.documents)

stopifnot(my.corpus$getWordCount(1L) == 879)
stopifnot(my.corpus$getDocumentCount() == length(cora.documents))
stopifnot(my.corpus$getWord(my.corpus$getIndex(3, 1)) == cora.documents[[4]][1,2])
