#include <vector>
#include <Rcpp.h>

/**
 * Abstract base class for data types.
 */
class DocumentData {
public:
  DocumentData() {
  }
};

class DiscreteDocumentData : public DocumentData {
public:
  DiscreteDocumentData(std::vector<int> words) :
    DocumentData(), words_(words) {
  }

  const size_t getWord(const size_t index) const {
    return words_[index];
  }

  const size_t numWords() const {
    return words_.size();
  }

private:
  std::vector<int> words_;
};

class ContinuousDocumentData : public DocumentData {
public:
  ContinuousDocumentData(double value) :
    DocumentData(), value_(value) { }

  const double getValue() const {
    return value_;
  }

private:
  double value_;
};

template<typename T>
class WrappedList {
public:
  WrappedList(const Rcpp::List& data) : Rdata_(data), data_(data.size()) {
    for (int ii = 0; ii < data.size(); ++ii) {
      data_.push_back(
        reinterpret_cast<T*>(
          Rcpp::Environment((SEXP)data[ii]).get(".pointer")));
    }
  }

  const T& getData(const size_t index) const {
    return *data_[index];
  }

  SEXP getRData(const size_t index) const {
    return Rdata_[index];
  }

  const size_t size() const {
    return data_.size();
  }

private:
  Rcpp::List Rdata_;
  std::vector<T*> data_;
};

class Document {
public:
  Document(Rcpp::List data) : data_(data) {
  }

  const DocumentData& getData(const size_t index) const {
    return data_.getData(index);
  }

  SEXP getRData(const size_t index) const {
    return data_.getRData(index);
  }

private:
  WrappedList<DocumentData> data_;
};

class Corpus {
public:
  Corpus(Rcpp::List documents) : documents_(documents) {
  }

  const Document& getDocument(const size_t index) const {
    return documents_.getData(index);
  }

  SEXP getRDocument(const size_t index) const {
    return documents_.getRData(index);
  }

  size_t numDocuments() const {
    return documents_.size();
  }

private:
  WrappedList<Document> documents_;
};

