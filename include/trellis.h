#ifndef HAM_TRELLIS_H
#define HAM_TRELLIS_H

#include <vector>
#include <stdint.h>
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
  trellis(Model* hmm, Sequence *seq);
  trellis(Model* hmm, Sequences *seqs);
  void Init();
  ~trellis();

  Model *model() { return hmm_; }
  double ending_viterbi_log_prob() { return ending_viterbi_log_prob_; }
  Sequences *seqs() { return seqs_; }
  float_2D* forward_table() { return forward_table_; }
  double forward_log_prob() { return ending_forward_log_prob_; }

  void Viterbi();
  void Forward();
  void Traceback(TracebackPath &path);
private:
  Model *hmm_;
  Sequences *seqs_;
  int_2D *traceback_table_;

  double ending_viterbi_log_prob_;
  int16_t ending_viterbi_pointer_;
  float_2D* forward_table_;
  double  ending_forward_log_prob_;

  // internal loop variables TODO wouldn't it make more sense to put these somewhere else?
  vector<double> *scoring_current_;
  vector<double> *scoring_previous_;
  vector<double> *swap_ptr_;
};

}
#endif
