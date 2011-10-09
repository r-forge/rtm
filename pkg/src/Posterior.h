#include <cassert>
#include <vector>

#include <Rcpp.h>


typedef std::vector<int> State;

class Position {
  int document;
  int word;
  int index;
public:
  int getDocument() const {
    return document;
  }

  int getWord() const {
    return word;
  }

  int getIndex() const {
    return index;
  }
};

/**
 * Abstract class defining a posterior component.  The full posterior for a
 * collection of assignments is proportional to the product of these posterior
 * components.
 */
class Posterior {
public:
  virtual double topicPosterior(const Position& position,
                                const State& state) = 0;
  virtual void updateState(const Position& position,
                           const State& state, int amount) = 0;
};


/**
 * Abstract class defining a posterior component.  The full posterior for a
 * collection of assignments is proportional to the product of these posterior
 * components.
 */
class BasicPosterior {
public:
  virtual double topicPosterior(const Position& position,
                                const int topic) = 0;
  virtual void updateState(const Position& position,
                           const int topic, const int amount) = 0;

  double topicPosterior(const Position& position,
                        const State& state) {
    assert(state.size() == 1);
    return topicPosterior(position, state[0]);
  }

  void updateState(const Position& position,
                   const State& state, int amount) {
    assert(state.size() == 1);
    return updateState(position, state[0], amount);
  }
};

/**
 * The following are the pieces for LDA.
 */
class DirichletPosterior : public BasicPosterior {
  Rcpp::Matrix<INTSXP> document_sums_;
  double alpha_;

public:
  DirichletPosterior(const int K, const int D, const double alpha) :
    document_sums_(K, D), alpha_(alpha) {
  }

  double topicPosterior(const Position& position, const int topic) {
    return document_sums_(position.getDocument(), topic) + alpha_;
  }

  void updateState(const Position& position,
                   const int topic, const int amount) {
    document_sums_(position.getDocument(), topic) += amount;
  }
};

class DiscretePosterior : public BasicPosterior {
  Rcpp::Matrix<INTSXP> topics_;
  std::vector<int> topic_sums_;
  double eta_;
  int V_;

  std::vector<int> words_;

public:
  DiscretePosterior(int K, int V, double eta, const std::vector<int>& words) :
    topics_(K, V), topic_sums_(K), eta_(eta), V_(V), words_(words) {
  }

  double topicPosterior(const Position& position, const int topic) {
    return (topics_(topic, words_[position.getIndex()]) + eta_) /
           (topic_sums_[topic] + V_ * eta_);
  }

  void updateState(const Position& position,
                   const int topic, const int amount) {
    topic_sums_[topic] += amount;
    topics_(topic, words_[position.getIndex()]) += amount;
  }
};

#if 0
class SwitchingPosterior : public Posterior {
  Posterior* switch_posterior_;
  std::vector<Posterior*> posteriors_;

public:
  SwitchingPosterior(Posterior* switch_posterior,
                     const std::vector<Posterior*>& posteriors) :
    switch_posterior_(switch_posterior), posteriors_(posteriors) {
  }

  double topicPosterior(const Position& position,
                        const State& state) {
    State substate = state;
    std::vector<int> switch_topic;
    switch_topic.push_back(substate[0]);
    substate.pop_front();
    assert(switch_topic[0] >= 0);
    assert(switch_topic[0] < posteriors_.size());
    return switch_posterior_->topicPosterior(position, switch_topic) *
           posteriors_[switch_topic[0]]->topicPosterior(position, substate);
  }

  void updateState(const Position& position,
                   const State& state, const int amount) {
    State substate = state;
    std::vector<int> switch_topic;
    switch_topic.push_back(substate[0]);
    substate.pop_front();
    assert(switch_topic[0] >= 0);
    assert(switch_topic[0] < posteriors_.size());
    switch_posterior_->udpateState(position, switch_topic, amount);
    posteriors_[switch_topic[0]]->udpateState(position, substate, amount);
  }
};

/**
 * Additional types of covariates.
*/
class ContinuousPosterior : public BasicPosterior {
  std::vector<double> words_;
  std::vector<int> n_;
  std::vector<int> x_;
  std::vector<int> x2_;

  double mean0_;
  double var0_;

  std::vector<double> values_;

public:
  ContinuousPosterior(int K, double mean0, double var0) :
    mean0_(mean0), var0_(var0), n_(K), x_(K), x2_(K) {
  }

  double topicPosterior(const Position& position, const int topic) {
    double mu = x_[topic] / n_[topic];
    double sigma2 = x2_[topic] / n_[topic] - mu * mu;

    double mean = ((sigma2 / n_[topic]) * mean0_ + var0_ * mu) /
                  (sigma2 / n_[topic] + var0_);
    double var = 1.0 / (n_[topic] / sigma2 + 1 / var0_);

    double val = values_[position.getIndex()];

    return exp(-(val - mean) * (val - mean) / 2 / var) / sqrt(var);
  }

  void updateState(const Position& position,
                   const int topic, const int amount) {
    double val = values_[position.getIndex()];

    x_[topic] += amount * val;
    x2_[topic] += amount * val * val;
    n_[topic] += amount;
  }
};

class DocumentAnnotationPosterior : public BasicPosterior {
  std::vector<double> words_;
  std::vector<double> beta_;

  std::vector<double> values_;

public:
  DocumentAnnotationPosterior(int K, const std::vector<double> beta) :
    beta_(beta) {
  }

  double topicPosterior(const Position& position, const int topic) {

    return exp(-(val - mean) * (val - mean) / 2 / var) / sqrt(var);
  }

  void updateState(const Position& position,
                   const int topic, const int amount) {
  }
};

class BinaryPosterior : public BasicPosterior {
  Rcpp::Matrix<INTSXP> pos_counts_;
  Rcpp::Matrix<INTSXP> neg_counts_;
  double eta1_;
  double eta0_;
  int V_;

  Rcpp::Matrix<INTSXP> words_;

public:
  BinaryPosterior(int K, int V, double eta1, double eta0,
                  const Rcpp::Matrix<INTSXP>& words) :
    pos_counts_(K, V), neg_counts_(K, V),
    eta1_(eta1), eta0_(eta0), V_(V), words_(words) {
  }

  double topicPosterior(const Position& position, const int topic) {
    double pos_prob = (pos_counts_(topic, position.getWord()) + eta1_) /
      (pos_counts_(topic, position.getWord()) +
       neg_counts_(topic, position.getWord()) +
       eta1_ + eta0_);
    if (words_(position.getDocument(), position.getWord())) {
      return pos_prob;
    } else {
      return 1.0 - pos_prob;
    }
  }

  void updateState(const Position& position,
                   const int topic, const int amount) {
    if (words_(position.getDocument(), position.getWord())) {
      pos_counts_(topic, position.getWord()) += amount;
    } else {
      neg_counts_(topic, position.getWord()) += amount;
    }
  }
};

/**
 * Additional priors.
 */
class LabeledDirichletPosterior : public DirichletPosterior {

protected:
  virtual bool canAssign(const Position& position, int topic) = 0;

public:
  LabeledDirichletPosterior(int K, int D, double alpha) :
    DirichletPosterior(K, D, alpha) {
  }

  double topicPosterior(const Position& position, int topic) {
    if (canAssign(position, topic)) {
      return DirichletPosterior::topicPosterior(position, topic);
    } else {
      return 0;
    }
  }
};

class DocumentLabeledDirichletPosterior : public LabeledDirichletPosterior {
  Rcpp::Matrix<INTSXP> allowed_topics_;

public:
  DocumentLabeledDirichletPosterior(int K, int D, double alpha,
                                    const Rcpp::Matrix<INTSXP>& allowed_topics) :
    LabeledDirichletPosterior(K, D, alpha),
    allowed_topics_(allowed_topics) {
  }

protected:
  bool canAssign(const Position& position, int topic) {
    return allowed_topics_(position.getDocument(), topic);
  }
};

class WordLabeledDirichletPosterior : public LabeledDirichletPosterior {
  Rcpp::Matrix<INTSXP> allowed_topics_;

protected:
  WordLabeledDirichletPosterior(int K, int D, double alpha,
                                const Rcpp::Matrix<INTSXP>& allowed_topics) :
    LabeledDirichletPosterior(K, D, alpha),
    allowed_topics_(allowed_topics) {
  }

  bool canAssign(const Position& position, int topic) {
    return allowed_topics_(position.getIndex(), topic);
  }
};

#endif
