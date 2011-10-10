#include <vector>
#include <Rcpp.h>

#include "Corpus.h"
#include "Posterior.h"
#include "Sampler.h"

RCPP_MODULE(rtm) {

  using namespace Rcpp;

  Rcpp::class_<DiscreteDocumentData>("DiscreteDocumentData")
    .constructor<std::vector<int> >()
    .const_method("getWord", &DiscreteDocumentData::getWord)
    .const_method("numWords", &DiscreteDocumentData::numWords);

  Rcpp::class_<ContinuousAnnotationDocumentData>("ContinuousAnnotationDocumentData")
    .constructor<double>()
    .const_method("getValue", &ContinuousAnnotationDocumentData::getValue);

  Rcpp::class_<Document>("Document")
    .constructor<Rcpp::List>()
    .const_method("getData", &Document::getRData);

  Rcpp::class_<Corpus>("Corpus")
    .constructor<Rcpp::List>()
    .method("numDocuments", &Corpus::numDocuments)
    .method("getDocument", &Corpus::getRDocument);

  Rcpp::class_<DirichletPosterior>("DirichletPosterior")
    .constructor<int, int, double>();

  Rcpp::class_<DiscretePosterior>("DiscretePosterior")
    .constructor<int, int, double>();

  Rcpp::class_<ContinuousAnnotationPosterior>("ContinuousAnnotationPosterior")
    .constructor<std::vector<double>, double>();

//  Rcpp::class_<Sampler>("Sampler")
//    .constructor();
}
