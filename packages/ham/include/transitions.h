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
  inline void set_to_state(State* st) { to_state_ = st; }
  inline string &to_state_name() {return to_state_name_;};
  inline State *to_state() {return to_state_;};

  inline double log_prob() { return log_prob_; }

  void Print();
private:
  string to_state_name_;  // the state to which we are transitioning (set during parsement)
  State *to_state_;  // pointer to the to-state (set in model::FinalizeState)
  double log_prob_;  // log transition probability
};

}
#endif
