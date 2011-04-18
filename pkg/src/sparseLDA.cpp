#include <Rcpp.h>
#include <stdint.h>
#include <cmath>
#include <assert.h>
#include <vector>
#include <iostream>

// TODO:
//  * Finish/test RTM
//  * Support train/test.
//  * sLDA

#define min(x, y) (((x) < (y)) ? (x) : (y)) 

class Links {
	std::vector<unsigned int> links_;
	std::vector<unsigned int> lengths_;
	std::vector<unsigned int> indices_;
	unsigned int D_;
	
public:
	void load(const std::vector<unsigned int>& links,
						const std::vector<unsigned int>& lengths) {
		links_ = links;
		lengths_ = lengths;

		D_ = lengths.size();
		indices_.resize(D_);
		indices_[0] = 0;
		for (int ii = 1; ii < D_; ++ii) {
			indices_[ii] = indices_[ii - 1] + lengths_[ii - 1];
		}
	}

	unsigned int getNumLinks(int document) const {
		return lengths_[document];
	}

	unsigned int getTarget(int document, int offset) const {
		return links_[indices_[document] + offset];
	}	

	unsigned int getDocumentCount() const {
		return D_;
	}	
};

class Corpus {
	std::vector<unsigned int> words_;
	std::vector<unsigned int> lengths_;
	unsigned int D_;
	unsigned int V_;
	std::vector<unsigned int> counts_;
	std::vector<unsigned int> indices_;
	unsigned int total_count_;

	Rcpp::Matrix<REALSXP> annotations_;
	unsigned int num_annotations_;

 public:
  Corpus() { } 

	void load(const std::vector<unsigned int>& words,
						const std::vector<unsigned int>& lengths,
						unsigned int V) {
		num_annotations_ = 0;
		words_ = words;
		lengths_ = lengths;
		D_ = lengths.size();
		V_ = V; 
		//		words_(words), lengths_(lengths), D_(words.size()), V_(V) {
		indices_.resize(D_);
		indices_[0] = 0;
		total_count_ = lengths[0];
		for (unsigned int ii = 1; ii < D_; ++ii) {
			indices_[ii] = indices_[ii - 1] + lengths[ii - 1];
			total_count_ += lengths[ii];
		}

		counts_.resize(V);
		for (unsigned int ii = 0; ii < V; ++ii) {
			counts_[ii] = 0;
		}
		for (unsigned int ii = 0; ii < total_count_; ++ii) {
			counts_[words[ii]]++;
		}
	}

	void setAnnotations(Rcpp::Matrix<REALSXP> annotations) {
		annotations_ = annotations;
		num_annotations_ = annotations.nrow();
		assert(annotations.ncol() == D_);
	}

	unsigned int getNumAnnotations() const {
		return num_annotations_;		
	}

	std::vector<double> getAnnotations(unsigned int document) const {
		std::vector<double> result;
		for (int ii = 0; ii < num_annotations_; ++ii) {
			result.push_back(annotations_[ii + document * num_annotations_]);
		}
		return result;
	}

	unsigned int getIndex(int document, int word) const {
		assert(word < V_);
		assert(document < D_);
		return indices_[document] + word;
	}

	unsigned int getDocumentCount() const {
		return D_;
	}

	unsigned int getTotalCount() const {
		return total_count_;
	}

	unsigned int getWordCount(int word) const {
		assert(word < V_);
		return counts_[word];
	}

	int getLength(int document) const {
		assert(document < D_);
		return lengths_[document];
	}

	unsigned int getWord(int index) const {
		return words_[index];
	}

	unsigned int getV() const {
		return V_;
	}	
};

class Topics {
	std::vector<uint32_t> data_;
	unsigned int K_;
	unsigned int V_;
	std::vector<uint32_t> indices_;
	std::vector<uint32_t> lengths_;
	unsigned int M_;
	uint32_t mask_;

 public:
  void loadR(unsigned int K, SEXP corpus) {
		Rcpp::XPtr<Corpus> c(corpus);
		load(K, *c);
	}

  void load(unsigned int K, const Corpus& corpus) {
		//	K_(K), V_(corpus.getV()) {
		K_ = K;
		V_ = corpus.getV();

		indices_.resize(V_);
		lengths_.resize(V_);

		indices_[0] = 0;
		lengths_[0] = min(corpus.getWordCount(0), K);
		int total = lengths_[0];

		for (unsigned int ii = 1; ii < V_; ++ii) {
			lengths_[ii] = min(corpus.getWordCount(ii), K);
			indices_[ii] = indices_[ii - 1] + lengths_[ii - 1];
			total += lengths_[ii];
		}

		data_.resize(total);
		for (unsigned int ii = 0; ii < total; ++ii) {
			data_[ii] = 0;
		}
		// Set up M_ and mask_
		M_ = ceil(log2(K));
		mask_ = (1L << (M_)) - 1;
	}

	unsigned int getLength(int word) const {
		return lengths_[word];
	}
	
	class WordIterator {
		unsigned int offset_;
		unsigned int length_;
		unsigned int index_;
		Topics *topics_;
		
		unsigned int getIndexCount(int index) const {
			return topics_->data_[offset_ + index] >> topics_->M_;
		}

	public:
   	WordIterator(int word, Topics* topics) : index_(0), topics_(topics) {
			offset_ = topics->indices_[word];
			length_ = topics->lengths_[word];
		}
		
		void next() {
			++index_;
		}

		bool end() const {
			return index_ >= length_ || getCount() == 0;
		}

		unsigned int getCount() const {
			return getIndexCount(index_);
		}

		unsigned int getTopic() const {
			return topics_->data_[offset_ + index_] & topics_->mask_;
		}		

		unsigned int getIndex() const {
			return index_;
		}

		void swap(int index1, int index2) {
			uint32_t temp = topics_->data_[offset_ + index1];
			topics_->data_[offset_ + index1] = topics_->data_[offset_ + index2];
			topics_->data_[offset_ + index2] = temp;
		}

		void setCountTopic(int count, int topic) {
			assert(index_ < length_);
			topics_->data_[offset_ + index_] = (count << topics_->M_) + topic;
		}

		void dumpState() {
			for (int ii = 0; ii < length_; ++ii) { 
				std::cout << ii << ":" << topics_->data_[offset_ + ii] << std::endl;
			}
		}

		void bubbleSort(int position, int amount) {
			if (position > 0 && amount == 1) {
				if (topics_->data_[offset_ + position] > 
						topics_->data_[offset_ + position - 1]) {
					swap(position, position - 1);
					bubbleSort(position - 1, amount);
					return;
				}				
			} 
			if (position < length_ - 1 && amount == -1) {
				if (topics_->data_[offset_ + position] < 
						topics_->data_[offset_ + position + 1]) {
					swap(position, position + 1);
					bubbleSort(position + 1, amount);
					return;
				}				
			}
		}

		void modifyCount(int amount) {			
			topics_->data_[offset_ + index_] += amount << topics_->M_;
			bubbleSort(index_, amount);
		}
	};

	void updateCount(int word, unsigned int topic, int amount) {
		assert(amount == 1 || amount == -1);
		WordIterator ii(word, this);

		/*
			if (word == 633) {
			std::cout << word << " " << topic << std::endl;
			ii.dumpState();
			}
		*/

		for (; !ii.end(); ii.next()) {
			if (ii.getTopic() == topic) {
				ii.modifyCount(amount);
				return;
			}
		}

		assert(amount == 1);
		ii.setCountTopic(amount, topic);
	}

	Rcpp::Matrix<INTSXP> getTopics() const {		
		Rcpp::Matrix<INTSXP> result(V_, K_);

		for (int vv = 0; vv < V_; ++vv) {
			for (int kk = 0; kk < K_; ++kk) {
				result(vv, kk) = 0;
			}

			for (WordIterator ii(vv, const_cast<Topics*>(this)); 
					 !ii.end(); ii.next()) {
				assert(result(vv, ii.getTopic()) == 0);
				//				ii.dumpState();
				result(vv, ii.getTopic()) = ii.getCount();
			}
		}
		return result;
	}
};

class HiddenMarkovTM {
  std::vector<unsigned int> z_;
  Rcpp::Matrix<INTSXP> transition_counts_;
	std::vector<unsigned int> topic_sums_;

  double alpha_;
	double eta_;
	unsigned int K_;
  unsigned int V_;

  Topics topics_;
	const Corpus* corpus_;

	void updateCounts(int topic, int amount,
                    int prev_topic, int next_topic) {
    if (prev_topic != -1) {
      transition_counts_(prev_topic, topic) += amount;
    }
    if (next_topic != -1) {
      transition_counts_(topic, next_topic) += amount;
    }
		topic_sums_[topic] += amount;
	}

  std::vector<double> tmp_;

  int sampleWord(int word, int prev_z, int next_z) {
    for (int ii = 0; ii < K_; ++ii) {
      tmp_[ii] = eta_;
    }

		for (Topics::WordIterator ii(word, &topics_);
				 !ii.end(); ii.next()) {
			int count = ii.getCount();
			int topic = ii.getTopic();
      tmp_[topic] += count;
		}

    double tmp_sum = 0.0;
    for (int ii = 0; ii < K_; ++ii) {
      if (prev_z != -1) {
        tmp_[ii] *= transition_counts_(prev_z, ii) + alpha_; 
      }
      if (next_z != -1) {
        tmp_[ii] *= (transition_counts_(ii, next_z) + alpha_) / (topic_sums_[ii] + alpha_ * K_);
      }
      tmp_[ii] /= topic_sums_[ii] + V_ * eta_;
      tmp_sum += tmp_[ii];
    }

		double U = unif_rand() * tmp_sum;

		for (unsigned int ii = 0; ii < K_; ++ii) {
      if (U <= tmp_[ii]) {
        return ii;
      }
      U -= tmp_[ii];
    }			
		assert(false);
	}

	void load(const Corpus* corpus, 
						double alpha, double eta, 
						unsigned int K) { 
	  alpha_ = alpha;
	  eta_ = eta;

    K_ = K;
	  V_ = corpus->getV();
	  corpus_ = corpus;

    transition_counts_ = Rcpp::Matrix<INTSXP>(K_, K_);
	  topic_sums_.resize(K);
	  topics_.load(K, *corpus);
		 
	  z_.resize(corpus->getTotalCount());
    tmp_.resize(K_);
	}

public:
  HiddenMarkovTM() { }

	void initializeRandom() {
	  GetRNGstate();
		for (unsigned int dd = 0; dd < corpus_->getDocumentCount(); ++dd) { 
			for (int ww = 0; ww < corpus_->getLength(dd); ++ww) {
				int index = corpus_->getIndex(dd, ww);
				z_[index] = unif_rand() * K_;
				topics_.updateCount(corpus_->getWord(index), z_[index], 1);
			}

			for (int ww = 0; ww < corpus_->getLength(dd); ++ww) {
				int index = corpus_->getIndex(dd, ww);
        int prev_topic = -1;
        if (ww > 0) {
          prev_topic = z_[index - 1];
        }
        int next_topic = -1;
        if (ww < corpus_->getLength(dd) - 1) {
          next_topic= z_[index + 1];
        }
        updateCounts(z_[index], 1, prev_topic, next_topic);
			}
		}
		PutRNGstate();
	}

  void iterateDocument(int document) {
		for (int ww = 0; ww < corpus_->getLength(document); ++ww) {
			int index = corpus_->getIndex(document, ww);
			int old_topic = z_[index];
     
      int prev_topic = -1;
      if (ww > 0) {
        prev_topic = z_[index - 1];
      }
      int next_topic = -1;
      if (ww < corpus_->getLength(document) - 1) {
        next_topic= z_[index + 1];
      }
      
      updateCounts(old_topic, -1, prev_topic, next_topic);
			topics_.updateCount(corpus_->getWord(index), old_topic, -1);

			int new_topic = sampleWord(corpus_->getWord(index), prev_topic, next_topic);
			z_[index] = new_topic;
      updateCounts(new_topic, 1, prev_topic, next_topic);
			topics_.updateCount(corpus_->getWord(index), new_topic, 1);
		}
  }

	void loadR(SEXP corpus, 
						 double alpha, double eta, 
						 unsigned int K) { 
		Rcpp::XPtr<Corpus> c(corpus);
		load(c, alpha, eta, K);
	}

	const std::vector<unsigned int>& getAssignments() const {
		return z_;
	}
	
	const std::vector<unsigned int>& getTopicSums() const {
		return topic_sums_;
	}
	
	Rcpp::Matrix<INTSXP> getTopics() const {
		return topics_.getTopics();
	}
	
	const Rcpp::Matrix<INTSXP>& getTransitionCounts() const {
		return transition_counts_;
	}

	void iterateCorpus(int num_iterations) {
		GetRNGstate();
		for (int ii = 0; ii < num_iterations; ++ii) {
			std::cout << "Iteration " << ii << std::endl;
			for (unsigned int dd = 0; dd < corpus_->getDocumentCount(); ++dd) {
				// 
				//				if (dd % 100 == 0) {
				//					std::cout << ".";
				//				}
				R_CheckUserInterrupt();
				iterateDocument(dd);
			}
		}
		PutRNGstate();
	}
};


class SparseRTM {
	std::vector<double> s_;
	double s_sum_;
	std::vector<double> r_;
	double r_sum_;

	std::vector<unsigned int> z_;

	double alpha_;
	double eta_;
	std::vector<double> beta_;
	
	unsigned int K_;
	unsigned int V_;
	unsigned int D_;
	
	std::vector<unsigned int> topic_sums_;
	Rcpp::Matrix<INTSXP> document_sums_;
	std::vector<double> link_probs_;

  Topics topics_;
	const Corpus* corpus_;
	const Links* links_;

	unsigned int q_assignments;
	unsigned int s_assignments;
	unsigned int r_assignments;

	std::vector<double> annotation_weight;

	void load(const Corpus* corpus, 
						const Links* links,
						double alpha, double eta, 
						std::vector<double> beta, 
						unsigned int K) { 
		 //  	alpha_(alpha), eta_(eta), beta_(beta), K_(K), V_(corpus.getV()), corpus_(corpus) {
	  alpha_ = alpha;
	  eta_ = eta;
	  beta_ = beta;
	  K_ = K;
	  V_ = corpus->getV();
	  corpus_ = corpus;

		annotation_weight.resize(corpus->getNumAnnotations());

		link_probs_.resize(K_);
		links_ = links;
		assert(links->getDocumentCount() == corpus_->getDocumentCount());

		// Init document sums
		document_sums_ = Rcpp::Matrix<INTSXP>(K_, corpus_->getDocumentCount());
	 
	  topic_sums_.resize(K);
	  topics_.load(K, *corpus);
		 
	  z_.resize(corpus->getTotalCount());

		// Initialize s.
		s_.resize(K_);
		s_sum_ = 0.0;
		
		for (unsigned int ii = 0; ii < K_; ++ii) {
			s_[ii] = alpha_ * eta_ / (V_ * eta_ + topic_sums_[ii]);
			s_sum_ += s_[ii];
		}

		// Allocate r.
		r_.resize(K_);
		r_sum_ = 0.0;		
	}

  void initializeR(int document) {
		r_sum_ = 0.0;
		for (unsigned int ii = 0; ii < K_; ++ii) {
			r_[ii] = document_sums_(ii, document) * eta_ / 
				(V_ * eta_ + topic_sums_[ii]);
			r_sum_ += r_[ii];
		}
	}

  int sampleWord(int word) {
		double q_sum = 0.0;
		for (Topics::WordIterator ii(word, &topics_);
				 !ii.end(); ii.next()) {
			int count = ii.getCount();
			int topic = ii.getTopic();
			q_sum += (r_[topic] + s_[topic]) * count;
		}

		double U = unif_rand() * (r_sum_ + s_sum_ + q_sum);

		if (U < s_sum_) {
			s_assignments++;
			for (unsigned int ii = 0; ii < K_; ++ii) {
				if (U <= s_[ii]) {
					return ii;
				}
				U -= s_[ii];
			}
		} else if (U < s_sum_ + r_sum_) {
			r_assignments++;
			U -= s_sum_;
			for (unsigned int ii = 0; ii < K_; ++ii) {
				if (U <= r_[ii]) {
					return ii;
				}
				U -= r_[ii];
			}			
		} else {
			q_assignments++;
			U -= s_sum_ + r_sum_;			
			for (Topics::WordIterator ii(word, &topics_);
					 !ii.end(); ii.next()) {
				int count = ii.getCount();
				int topic = ii.getTopic();
				double q_i = (r_[topic] + s_[topic]) * count;
				if (U <= q_i) {
					return topic;
				}				
				U -= q_i;
			}
		}

		assert(false);
	}

	void initializeLinks(int document) {		
		for (int kk = 0; kk < K_; ++kk) {
			link_probs_[kk] = 1.0;
		}

		for (int ii = 0; ii < links_->getNumLinks(document); +ii) {
			int doc2 = links_->getTarget(document, ii);
			assert(doc2 >= 0 && doc2 < D_);
			for (int kk = 0; kk < K_; ++kk) {
				link_probs_[kk] *= 
					exp(beta_[kk] * document_sums_(kk, doc2) / 
							corpus_->getLength(document) / 
							corpus_->getLength(doc2));
			}
		}
	}

	// Coefficients are gamma.
	// We memoize the sum:
	//   A = (1/N) \sum_w gamma[z_w]

	// We can rewrite the impact of the probability by:
	// exp(-(x - (1 / N) \sum_w gamma[z_w])^2 / (2 * sigma^2))
  //  = exp(-(x - A)^2 / (2 * sigma^2))
	// x^2 term in constant /wrt choice of z_w, so
	// exp(x * A / sigma^2 - A^2 / (2 * sigma^2))

	// Now denote A = A' + gamma[z_w] / N <- word to be sampled.
	// A' is constant.
	// exp(x * (A' + gamma[z_w] / N) / sigma^2 - (A' + gamma[z_w] / N)^2 / (2 * sigma^2))
	// exp(x * gamma[z_w] / N / sigma^2) * exp(-A' * gamma[z_w] / N / sigma^2) * exp(-gamma[z_w]^2 / N^2 / 2 / sigma^2)
	// Let G = gamma[z_w] / N
	// exp((x - A') * G / sigma^2) * exp(-G^2 / 2 /sigma^2)
	// = exp((-G / 2 + (x - A')) * G / sigma^2)

	// annotation_weight = (x - A')
#if 0
	void initializeAnnotations(int document) {
		std::vector<double> annotations = corpus_->getAnnotations(document);
		assert(annotations.size() == corpus_->getNumAnnotations());
		unsigned int N = corpus_->getWordCount(document);
		for (int ii = 0; ii < corpus_->getNumAnnotations(); ++ii) {
			annotation_weight[ii] = annotations[ii];
			for (int ww = 0; ww < N; ++ww) {
				int z_w = z_[getIndex(document, ww)];
				annotation_weight[ii] -= gamma_[z_w] / N;
			}
		}

		exp(G * (annotation_weight[ii] - G / 2) / variance);
	}
#endif

	void updateCounts(int document, int topic, int amount) {
		document_sums_(topic, document) += amount;
		topic_sums_[topic] += amount;
	}

	void updateSR(int document, int topic) {
		s_sum_ -= s_[topic];
		r_sum_ -= r_[topic];
		
		s_[topic] = alpha_ * eta_ / (V_ * eta_ + topic_sums_[topic]);
		r_[topic] = document_sums_(topic,  document) * eta_ / 
			(V_ * eta_ + topic_sums_[topic]);
		
		s_sum_ += s_[topic];
		r_sum_ += r_[topic];
	}

public:
	void loadR(SEXP corpus, 
						 SEXP links,
						 double alpha, double eta, 
						 std::vector<double> beta, unsigned int K) { 
		Rcpp::XPtr<Corpus> c(corpus);
		Rcpp::XPtr<Links> d(links);
		load(c, d, alpha, eta, beta, K);
	}

	void initializeRandom() {
	  GetRNGstate();
		for (unsigned int dd = 0; dd < corpus_->getDocumentCount(); ++dd) { 
			for (int ww = 0; ww < corpus_->getLength(dd); ++ww) {
				int index = corpus_->getIndex(dd, ww);
				z_[index] = unif_rand() * K_;
				topics_.updateCount(corpus_->getWord(index), z_[index], 1);
				updateCounts(dd, z_[index], 1);
			}
		}
		PutRNGstate();
	}

	void initialize(std::vector<unsigned int> assignments) {
		for (unsigned int dd = 0; dd < corpus_->getDocumentCount(); ++dd) { 
			for (int ww = 0; ww < corpus_->getLength(dd); ++ww) {
				int index = corpus_->getIndex(dd, ww);
				z_[index] = assignments[index];
				topics_.updateCount(corpus_->getWord(index), z_[index], 1);
				updateCounts(dd, z_[index], 1);
			}
		}
	}

	void iterateCorpus(int num_iterations) {
		GetRNGstate();
		q_assignments = 0;
		r_assignments = 0;
		s_assignments = 0;
		for (int ii = 0; ii < num_iterations; ++ii) {
			std::cout << "Iteration " << ii;
			for (unsigned int dd = 0; dd < corpus_->getDocumentCount(); ++dd) {
				// 
				//				if (dd % 100 == 0) {
				//					std::cout << ".";
				//				}
				R_CheckUserInterrupt();
				iterateDocument(dd);
			}
			std::cout << "Q: " << q_assignments << "(" << (double)q_assignments / (q_assignments + r_assignments + s_assignments) << "%) ";
			std::cout << "R: " << r_assignments << "(" << (double)r_assignments / (q_assignments + r_assignments + s_assignments) << "%) ";
			std::cout << "S: " << s_assignments << "(" << (double)s_assignments / (q_assignments + r_assignments + s_assignments) << "%) " << std::endl;
		}
		PutRNGstate();
	}

	void iterateDocument(int document) {
		initializeLinks(document);
		initializeR(document);		
		for (int ww = 0; ww < corpus_->getLength(document); ++ww) {
			int index = corpus_->getIndex(document, ww);
			int old_topic = z_[index];
			topics_.updateCount(corpus_->getWord(index), old_topic, -1);
			updateCounts(document, old_topic, -1);
			updateSR(document, old_topic);

			int new_topic = sampleWord(corpus_->getWord(index));
			z_[index] = new_topic;
			topics_.updateCount(corpus_->getWord(index), new_topic, 1);
			updateCounts(document, new_topic, 1);
			updateSR(document, new_topic);
		}
	}

	const std::vector<unsigned int>& getAssignments() const {
		return z_;
	}
	
	const std::vector<unsigned int>& getTopicSums() const {
		return topic_sums_;
	}
	
	Rcpp::Matrix<INTSXP> getTopics() const {
		return topics_.getTopics();
	}
	
	const Rcpp::Matrix<INTSXP>& getDocumentSums() const {
		return document_sums_;
	}
};


RCPP_MODULE(rtm) {
	using namespace Rcpp;

	class_<Links>("Links")		
		.method("load", &Links::load)
		.const_method("getNumLinks", &Links::getNumLinks)
		.const_method("getTarget", &Links::getTarget)
		.const_method("getDocumentCount", &Links::getDocumentCount);

	class_<Corpus>("Corpus")		
    .constructor()
		.method("load", &Corpus::load)
		.method("setAnnotations", &Corpus::setAnnotations)
		.const_method("getNumAnnotations", &Corpus::getNumAnnotations)
		.const_method("getAnnotations", &Corpus::getAnnotations)
		.const_method("getIndex", &Corpus::getIndex)
		.const_method("getDocumentCount", &Corpus::getDocumentCount)
		.const_method("getTotalCount", &Corpus::getTotalCount)
		.const_method("getWordCount", &Corpus::getWordCount)
		.const_method("getLength", &Corpus::getLength)
		.const_method("getV", &Corpus::getV)
		.const_method("getWord", &Corpus::getWord);

	class_<Topics>("Topics")		
		.method("load", &Topics::loadR)
		.const_method("getTopics", &Topics::getTopics)
		.const_method("getLength", &Topics::getLength)
		.method("updateCount", &Topics::updateCount);

	class_<SparseRTM>("SparseRTM")
		.method("load", &SparseRTM::loadR)
		.const_method("getAssignments", &SparseRTM::getAssignments)
		.const_method("getTopics", &SparseRTM::getTopics)
		.const_method("getDocumentSums", &SparseRTM::getDocumentSums)
		.const_method("getTopicSums", &SparseRTM::getTopicSums)
		.method("iterateCorpus", &SparseRTM::iterateCorpus)
		.method("iterateDocument", &SparseRTM::iterateDocument)
		.method("initialize", &SparseRTM::initialize)
		.method("initializeRandom", &SparseRTM::initializeRandom);

  class_<HiddenMarkovTM>("HiddenMarkovTM")
    .constructor()
		.method("load", &HiddenMarkovTM::loadR)
		.const_method("getAssignments", &HiddenMarkovTM::getAssignments)
		.const_method("getTopics", &HiddenMarkovTM::getTopics)
		.const_method("getTransitionCounts", &HiddenMarkovTM::getTransitionCounts)
		.const_method("getTopicSums", &HiddenMarkovTM::getTopicSums)
		.method("iterateCorpus", &HiddenMarkovTM::iterateCorpus)
		.method("iterateDocument", &HiddenMarkovTM::iterateDocument)
		.method("initializeRandom", &HiddenMarkovTM::initializeRandom);
}



