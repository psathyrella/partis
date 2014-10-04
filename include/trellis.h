#ifndef HAM_TRELLIS_H
#define HAM_TRELLIS_H

#include <vector>
#include <stdint.h>
#include "sequences.h"
#include "hmm.h"
#include "tracebackpath.h"

using namespace std;
namespace ham {
      
typedef vector<vector<int16_t> > int_2D;
typedef vector<vector<float> > float_2D;
typedef vector<vector<double> > double_2D;

// ----------------------------------------------------------------------------------------
class trellis{
public:
  double  ending_viterbi_score;
  trellis(model* h , Sequences* sqs);
  ~trellis();

  inline model* getModel(){return hmm;}
  inline Sequences* getSeq(){return seqs;}
  inline float_2D* getForwardTable(){return forward_score;}
  inline double getForwardProbability(){return ending_forward_prob;}
              
  void viterbi();
  void forward();
  void traceback(TracebackPath& path);
private:
  model* hmm;
  Sequences* seqs; //Digitized Sequences
  int_2D*         traceback_table;        //Simple traceback table
              
  //Score Tables
  // float_2D*       viterbi_score;      //Storing Viterbi scores
  float_2D*       forward_score;      //Storing Forward scores
  float_2D*       backward_score;     //Storing Backward scores
  double_2D*      posterior_score;        //Storing Posterior scores
              
  //Ending Cells
  int16_t ending_viterbi_tb;
  double  ending_forward_prob;
              
  //Cells used for scoring each cell
  vector<double>* scoring_current;
  vector<double>* scoring_previous;
              
  vector<double>* swap_ptr;
};
}
#endif
