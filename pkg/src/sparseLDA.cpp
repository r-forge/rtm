#include <Rcpp.h>
#include <stdint.h>
#include <cmath>
#include <assert.h>
#include <vector>

#define min(x, y) (((x) < (y)) ? (x) : (y)) 

class Corpus {
	std::vector<unsigned int> words_;
	std::vector<unsigned int> lengths_;
	unsigned int D_;
	unsigned int V_;
	std::vector<unsigned int> counts_;
	std::vector<unsigned int> indices_;
	unsigned int total_count_;

 public:
	void load(const std::vector<unsigned int>& words,
						const std::vector<unsigned int>& lengths,
						unsigned int V) {
		words_ = words;
		lengths_ = lengths;
		D_ = words.size();
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

		for (unsigned int ii = 0; ii < total_count_; ++ii) {
			counts_[words[ii]]++;
		}
	}

	unsigned int getIndex(int document, int word) const {
		return indices_[document] + word;
	}

	unsigned int getDocumentCount() const {
		return D_;
	}

	unsigned int getTotalCount() const {
		return total_count_;
	}

	unsigned int getWordCount(int word) const {
		return counts_[word];
	}

	int getLength(int document) const {
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
	uint32_t* data_;
	unsigned int K_;
	unsigned int V_;
	uint32_t* indices_;
	uint32_t* lengths_;
	unsigned int M_;
	uint32_t mask_;

 public:
  Topics(unsigned int K,
				 const Corpus& corpus) :
	K_(K), V_(corpus.getV()) {
		indices_ = new uint32_t[V_];
		lengths_ = new uint32_t[V_];

		indices_[0] = 0;
		lengths_[0] = min(corpus.getWordCount(0), K);
		int total = lengths_[0];

		for (unsigned int ii = 1; ii < V_; ++ii) {
			lengths_[ii] = min(corpus.getWordCount(ii), K);
			indices_[ii] = indices_[ii - 1] + lengths_[ii - 1];
			total += lengths_[ii];
		}

		data_ = new uint32_t[total];

		// Set up M_ and mask_
		M_ = ceil(log2(K));
		mask_ = (1L << (M_)) - 1;
	}
	
  ~Topics() {
		delete indices_;
		delete data_;
	}

	class WordIterator {
		uint32_t *word_data_;
		unsigned int length_;
		unsigned int index_;
		Topics *topics_;
		
		unsigned int getIndexCount(int index) const {
			return word_data_[index] >> topics_->M_;
		}

	public:
   	WordIterator(int word, Topics* topics) : index_(0), topics_(topics) {
			word_data_ = &topics->data_[topics->indices_[word]];
			length_ = topics->lengths_[word];
		}
		
		void next() {
			++index_;
		}

		bool end() const {
			return index_ >= length_ && getCount() > 0;
		}

		unsigned int getCount() const {
			return getIndexCount(index_);
		}

		unsigned int getTopic() const {
			return word_data_[index_] & topics_->mask_;
		}		

		unsigned int getIndex() const {
			return index_;
		}

		void swap(int index2) {
			uint32_t temp = word_data_[index_];
			word_data_[index_] = word_data_[index2];
			word_data_[index2] = temp;
		}

		void setCountTopic(int count, int topic) {
			assert(index_ < length_);
			word_data_[index_] = (count << topics_->M_) + topic;
		}

		void modifyCount(int amount) {			
			word_data_[index_] += amount << topics_->M_;
			if (index_ > 0 && amount == 1) {
				if (word_data_[index_] > word_data_[index_ - 1]) {
					swap(index_ - 1);
					return;
				}				
			} 
			if (index_ < length_ - 1 && amount == -1) {
				if (word_data_[index_] < word_data_[index_ + 1]) {
					swap(index_ + 1);
					return;
				}				
			}
		}
	};

	void updateCount(int word, unsigned int topic, int amount) {
		assert(amount == 1 || amount == -1);
		WordIterator ii(word, this);
		for (; !ii.end(); ii.next()) {
			if (ii.getTopic() == topic) {
				ii.modifyCount(amount);
			}
		}

		assert(amount == 1);
		ii.setCountTopic(amount, topic);
	}
};


class SparseRTM {
	double *s_;
	double s_sum_;
	double *r_;
	double r_sum_;

	unsigned int *z_;

	double alpha_;
	double eta_;
	double *beta_;
	
	unsigned int K_;
	unsigned int V_;
	unsigned int D_;
	
	uint32_t *topic_sums_;
	uint32_t *document_sums_;

  Topics* topics_;
	const Corpus& corpus_;

   SparseRTM(const Corpus& corpus, double alpha, double eta, double *beta, unsigned int K) :
  	alpha_(alpha), eta_(eta), beta_(beta), K_(K), V_(corpus.getV()), corpus_(corpus) {

	  GetRNGstate();
		topic_sums_ = new uint32_t[K];
		document_sums_ = new uint32_t[D_];
		topics_ = new Topics(K, corpus);

		z_ = new unsigned int[corpus.getTotalCount()];
		initializeZ();

		// Initialize s.
		s_ = new double[K_];
		s_sum_ = 0.0;
		
		for (unsigned int ii = 0; ii < K_; ++ii) {
			s_[ii] = alpha_ * eta_ / (V_ * eta_ + topic_sums_[ii]);
			s_sum_ += s_[ii];
		}

		// Allocate r.
		r_ = new double[K_];
		r_sum_ = 0.0;		

		PutRNGstate();
	}

	~SparseRTM() {
		delete r_;
		delete s_;
	}

	void initializeZ() {
		for (unsigned int dd = 0; dd < corpus_.getDocumentCount(); ++dd) { 
			for (int ww = 0; ww < corpus_.getLength(dd); ++ww) {
				int index = corpus_.getIndex(dd, ww);
				z_[index] = unif_rand() * K_;
				topics_->updateCount(corpus_.getWord(index), z_[index], 1);
				updateCounts(dd, z_[index], 1);
			}
		}
	}
				
  void initializeR(int document) {
		for (unsigned int ii = 0; ii < K_; ++ii) {
			r_[ii] = document_sums_[ii + K_ * document] * eta_ / 
				(V_ * eta_ + topic_sums_[ii]);
			r_sum_ += r_[ii];
		}
	}

	void iterateCorpus(int num_iterations) {
		GetRNGstate();
		for (int ii = 0; ii < num_iterations; ++ii) {
			for (unsigned int dd = 0; dd < corpus_.getDocumentCount(); ++dd) {
				iterateDocument(dd);
			}
		}
		PutRNGstate();
	}

	void iterateDocument(int document) {
		initializeR(document);
		for (int ww = 0; ww < corpus_.getLength(document); ++ww) {
			int index = corpus_.getIndex(document, ww);
			int old_topic = z_[index];
			topics_->updateCount(corpus_.getWord(index), old_topic, -1);
			updateCounts(document, old_topic, -1);
			updateSR(document, old_topic);

			int new_topic = sampleWord(corpus_.getWord(index));
			z_[index] = new_topic;
			topics_->updateCount(corpus_.getWord(index), new_topic, 1);
			updateCounts(document, new_topic, 1);
			updateSR(document, new_topic);
		}
	}

  int sampleWord(int word) {
		double q_sum = 0.0;
		for (Topics::WordIterator ii(word, topics_);
				 !ii.end(); ii.next()) {
			int count = ii.getCount();
			int topic = ii.getTopic();
			q_sum += (r_[topic] + s_[topic]) * count;
		}

		double U = unif_rand() * (r_sum_ + s_sum_ + q_sum);

		if (U < s_sum_) {
			for (unsigned int ii = 0; ii < K_; ++ii) {
				if (U <= s_[ii]) {
					return ii;
				}
				U -= s_[ii];
			}
		} else if (U < s_sum_ + r_sum_) {
			U -= s_sum_;
			for (unsigned int ii = 0; ii < K_; ++ii) {
				if (U <= r_[ii]) {
					return ii;
				}
				U -= r_[ii];
			}			
		} else {
			U -= s_sum_ + r_sum_;			
			for (Topics::WordIterator ii(word, topics_);
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

	void updateCounts(int document, int topic, int amount) {
		document_sums_[topic + K_ * document] += amount;
		topic_sums_[topic] += amount;
	}

	void updateSR(int document, int topic) {
		s_sum_ -= s_[topic];
		r_sum_ -= r_[topic];
		
		s_[topic] = alpha_ * eta_ / (V_ * eta_ + topic_sums_[topic]);
		r_[topic] = document_sums_[topic + K_ * document] * eta_ / 
			(V_ * eta_ + topic_sums_[topic]);
		
		s_sum_ += s_[topic];
		r_sum_ += r_[topic];
	}
};

RCPP_MODULE(rtm) {
	using namespace Rcpp;

	class_<Corpus>("Corpus")
		.const_method("getIndex", &Corpus::getIndex)
		.const_method("getDocumentCount", &Corpus::getDocumentCount)
		.const_method("getTotalCount", &Corpus::getTotalCount)
		.const_method("getWordCount", &Corpus::getWordCount)
		.const_method("getLength", &Corpus::getLength)
		.const_method("getWord", &Corpus::getWord);
}
