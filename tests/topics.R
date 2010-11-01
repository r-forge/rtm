require(rtm)
require(lda)

data(cora.documents)

my.corpus <- corpus(words = c(0L, 1L, 1L, 1L, 1L, 1L, 2L),
                    lengths = 7L)
my.topics <- topics(4L, my.corpus)
stopifnot(my.topics$getLength(1L) == 4L)
stopifnot(my.topics$getLength(0L) == 1L)

my.topics$updateCount(1L, 1L, 1L)
stopifnot(isTRUE(all.equal(my.topics$getTopics()[2,], c(0, 1, 0, 0))))

my.topics$updateCount(1L, 2L, 1L)
stopifnot(isTRUE(all.equal(my.topics$getTopics()[2,], c(0, 1, 1, 0))))

my.topics$updateCount(1L, 2L, 1L)
stopifnot(isTRUE(all.equal(my.topics$getTopics()[2,], c(0, 1, 2, 0))))

my.topics$updateCount(1L, 2L, -1L)
stopifnot(isTRUE(all.equal(my.topics$getTopics()[2,], c(0, 1, 1, 0))))

my.topics$updateCount(1L, 1L, -1L)
stopifnot(isTRUE(all.equal(my.topics$getTopics()[2,], c(0, 0, 1, 0))))
