#include <vector>
#include <Rcpp.h>

#include "Corpus.h"
#include "Posterior.h"

RCPP_MODULE(rtm) {

  using namespace Rcpp;

  Rcpp::class_<DiscreteDocumentData>("DiscreteDocumentData")
    .constructor<std::vector<int> >()
    .const_method("getWord", &DiscreteDocumentData::getWord)
    .const_method("numWords", &DiscreteDocumentData::numWords);

  Rcpp::class_<ContinuousDocumentData>("ContinuousDocumentData")
    .constructor<double>()
    .const_method("getValue", &ContinuousDocumentData::getValue);

  Rcpp::class_<Document>("Document")
    .constructor<Rcpp::List>()
    .const_method("getData", &Document::getRData);

  Rcpp::class_<Corpus>("Corpus")
    .constructor<Rcpp::List>()
    .method("numDocuments", &Corpus::numDocuments)
    .method("getDocument", &Corpus::getRDocument);
}
