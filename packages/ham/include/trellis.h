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
  trellis(Model *hmm, Sequence seq, trellis *cached_trellis = nullptr);
  trellis(Model *hmm, Sequences seqs, trellis *cached_trellis = nullptr);
  void Init();
  ~trellis();

  Model *model() { return hmm_; }
  double ending_viterbi_log_prob() { return ending_viterbi_log_prob_; }
  double ending_viterbi_log_prob(size_t length);  // return the log prob of the most probable path of length <length> NOTE this tacks the ending transition log prob onto whatever was in <viterbi_log_prob_>
  Sequences seqs() { return seqs_; }
  double ending_forward_log_prob() { return ending_forward_log_prob_; }
  double ending_forward_log_prob(size_t length);
  int_2D *traceback_table() const { return traceback_table_; }

  // NOTE (and beware) this is confusing to subtract one from the length. BUT it is totally on purpose: I want the calling code to be able to just worry about how long its sequence is.
  // In other words, I'm pretty sure we'll have to subtract (or add) 1 *somewhere*, and I've chosen to compartmentalize it into trellis.{h,cc}.
  // NOTE also that <viterbi_pointers_> *includes* the ending transition probability at each point.
  size_t viterbi_pointer(size_t length) {  // i.e. the zeroth entry of viterbi_pointers_ corresponds to stopping with sequence of length 1
    assert(length <= viterbi_pointers_->size());
    return (*viterbi_pointers_)[length - 1];
  }
  vector<double> *viterbi_log_probs() { return viterbi_log_probs_; }
  vector<double> *forward_log_probs() { return forward_log_probs_; }
  vector<int> *viterbi_pointers() { return viterbi_pointers_; }

  void SwapColumns(vector<double> *&scoring_previous, vector<double> *&scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states);
  void CacheViterbiVals(size_t position, double dpval, size_t i_st_current);  // if state <i_st_current> is the most likely at <position (i.e. <dpval> is the largest so far), then cache for future retrieval
  void CacheForwardLogProb(size_t position, double dpval, size_t i_st_current);  // add the logprob corresponding to state <i_st_current> at <position> (<dpval>) to the running cached total
  void MiddleVals(string algorithm, vector<double> *scoring_previous, vector<double> *scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states, size_t position);
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
  double  ending_viterbi_log_prob_;
  double  ending_forward_log_prob_;

  vector<double> *viterbi_log_probs_;  // log prob of best path up to and including each position NOTE includes log prob of transition to end
  vector<double> *forward_log_probs_;  // total log prob of all paths up to and including each position NOTE includes log prob of transition to end
  vector<int> *viterbi_pointers_;  // pointer to the state at which the best log prob occurred
  vector<double> *swap_ptr_;
};

}
#endif
