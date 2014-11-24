#ifndef HAM_TRELLIS_H
#define HAM_TRELLIS_H

#include <vector>
#include <stdint.h>
#include <iomanip>

#include "sequences.h"
#include "model.h"
#include "tracebackpath.h"

using namespace std;
namespace ham {

typedef vector<vector<int16_t> > int_2D;
typedef vector<vector<float> > float_2D;
typedef vector<vector<double> > double_2D;

// ----------------------------------------------------------------------------------------
class trellis {
public:
  trellis(Model *hmm, Sequence seq, trellis *cached_trellis=nullptr);
  trellis(Model *hmm, Sequences seqs, trellis *cached_trellis=nullptr);
  void Init();
  ~trellis();

  Model *model() { return hmm_; }
  // double ending_viterbi_log_prob() { return viterbi_log_probs_->at(seqs_->GetSequenceLength() - 1); }
  double ending_viterbi_log_prob() { return ending_viterbi_log_prob_; }
  double ending_viterbi_log_prob(size_t length);  // return the log prob of the most probable path of length <length> NOTE this tacks the ending transition log prob onto whatever was in <viterbi_log_prob_>
  Sequences seqs() { return seqs_; }
  float_2D* forward_table() { return forward_table_; }
  double forward_log_prob() { return ending_forward_log_prob_; }
  int_2D *traceback_table() const { return traceback_table_; }
  size_t viterbi_pointer(size_t length) {  // i.e. the zeroth entry of viterbi_pointers_ corresponds to stopping with sequence of length 1
    assert(length <= viterbi_pointers_->size());
    return (*viterbi_pointers_)[length-1];
  }
  vector<double> *viterbi_log_probs() { return viterbi_log_probs_; }  // TODO this is confusing having the functions above subtract one from the length, you need to chage it
  vector<int> *viterbi_pointers() { return viterbi_pointers_; }

  void Viterbi();
  void Forward();
  void Traceback(TracebackPath &path);

  void Dump();
private:
  Model *hmm_;
  Sequences seqs_;
  int_2D *traceback_table_;

  trellis *cached_trellis_;  // pointer to another trellis that already has its dp table(s) filled in, the idea being this trellis only needs a subset of that table, so we don't need to calculate anything new for this one

  int16_t ending_viterbi_pointer_;
  // float_2D *viterbi_table_;
  float_2D *forward_table_;
  double  ending_viterbi_log_prob_;
  double  ending_forward_log_prob_;

  // internal loop variables TODO wouldn't it make more sense to put these somewhere else?
  vector<double> *viterbi_log_probs_;  // log prob of best path up to and including each position NOTE does *not* include log prob of transition to end
  vector<int> *viterbi_pointers_;  // pointer to the state at which this best log prob occurred
  vector<double> *scoring_current_;
  vector<double> *scoring_previous_;
  vector<double> *swap_ptr_;
};

}
#endif
