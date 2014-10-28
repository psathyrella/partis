#ifndef HAM_TRANSITIONS_H
#define HAM_TRANSITIONS_H

#include <string>
#include <math.h>

using namespace std;

namespace ham {
class State;

// ----------------------------------------------------------------------------------------
class Transition {
  // friend class model;
  // friend class State;
public:
  Transition(string to_state, double prob);

  inline void set_name(string& txt) { to_state_name_ = txt; }
  inline void set_state(State* st) { to_state_ = st; }
  inline string &to_state_name() {return to_state_name_;};
  inline State *to_state() {return to_state_;};

  inline double log_prob() { return log_prob_; }

  void Print();
private:
  string to_state_name_;  //What state we are transitioning to (Fill out when parsing)
  State* to_state_;    // pointer to the to-state (Filled during HMM finalization)
  double log_prob_; //Log of standard transition probability
};

}
#endif
