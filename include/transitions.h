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

  void set_to_state(State* st) { to_state_ = st; }
  string &to_state_name() { return to_state_name_; }
  State *to_state() { return to_state_; }
  double log_prob() { return log_prob_; }

  void Print();
private:
  string to_state_name_;  // the state to which we are transitioning (set during parsement)
  State *to_state_;  // pointer to the to-state (set in model::FinalizeState)
  double log_prob_;  // log transition probability
};

}
#endif
