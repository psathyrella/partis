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

// ----------------------------------------------------------------------------------------
class Trellis {
public:
  Trellis(Model *hmm, Sequence seq, Trellis *cached_trellis = nullptr);
  Trellis(Model *hmm, Sequences seqs, Trellis *cached_trellis = nullptr);
  void Init();
  Trellis();
  ~Trellis();

  Model *model() { return hmm_; }
  Sequences seqs() { return seqs_; }
  double ending_viterbi_log_prob() { return ending_viterbi_log_prob_; }  // for full sequence length
  double ending_forward_log_prob() { return ending_forward_log_prob_; }  // for full sequence length
  // NOTE (and beware) this is confusing to subtract one from the length. BUT it is totally on purpose: I want the calling code to be able to just worry about how long its sequence is.
  // In other words, I'm pretty sure we'll have to subtract (or add) 1 *somewhere*, and I've chosen to compartmentalize it into trellis.{h,cc}.
  // NOTE also that <viterbi_indices_> *includes* the ending transition probability at each point.
  double ending_viterbi_log_prob(size_t length) { return viterbi_log_probs_pointer_->at(length - 1); } // the most probable path of length <length> NOTE do *not* use <viterbi_log_probs_>
  double ending_forward_log_prob(size_t length) { return forward_log_probs_pointer_->at(length - 1); } // NOTE do *not* use <forward_log_probs_>
  size_t viterbi_pointer(size_t length) { return viterbi_indices_pointer_->at(length - 1); } // i.e. the zeroth entry of viterbi_indices_ corresponds to stopping with sequence of length 1 NOTE do *not* use <viterbi_indices_>  

  int_2D *traceback_table_pointer() const { return traceback_table_pointer_; }
  vector<double> *viterbi_log_probs_pointer() { return viterbi_log_probs_pointer_; }
  vector<double> *forward_log_probs_pointer() { return forward_log_probs_pointer_; }
  vector<int> *viterbi_indices_pointer() { return viterbi_indices_pointer_; }

  void SwapColumns(vector<double> *&scoring_previous, vector<double> *&scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states);
  void MiddleViterbiVals(vector<double> *scoring_previous, vector<double> *scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states, size_t position);
  void MiddleForwardVals(vector<double> *scoring_previous, vector<double> *scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states, size_t position);
  void CacheViterbiVals(size_t position, double dpval, size_t i_st_current);
  void CacheForwardVals(size_t position, double dpval, size_t i_st_current);
  void Viterbi();
  void Forward();
  void Traceback(TracebackPath &path);

  string SizeString();
  double ApproxBytesUsed();

  void Dump();
private:
  Model *hmm_;
  Sequences seqs_;
  int_2D *traceback_table_pointer_;  // if we have a cached trellis, this points to the cached trellis's table
  int_2D traceback_table_;  // if we have a cached trellis, this isn't initialized

  Trellis *cached_trellis_;  // pointer to another trellis that already has its dp table(s) filled in, the idea being this trellis only needs a subset of that table, so we don't need to calculate anything new for this one

  int16_t ending_viterbi_pointer_;
  double  ending_viterbi_log_prob_;
  double  ending_forward_log_prob_;

  // chunk caching stuff
  vector<double> *viterbi_log_probs_pointer_;  // see notes for traceback_table_
  vector<double> *forward_log_probs_pointer_;  // see notes for traceback_table_
  vector<int> *viterbi_indices_pointer_;  // see notes for traceback_table_
  vector<double> viterbi_log_probs_;  // log prob of best path up to and including each position NOTE includes log prob of transition to end
  vector<double> forward_log_probs_;  // total log prob of all paths up to and including each position NOTE includes log prob of transition to end
  vector<int> viterbi_indices_;  // pointer to the state at which the best log prob occurred

  vector<double> *swap_ptr_;
  vector<double> scoring_current_, scoring_previous_;
};

}
#endif
