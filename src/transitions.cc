#include "transitions.h"
namespace stochhmm {
  
// ----------------------------------------------------------------------------------------
transition::transition(string to_state, double prob) : stateName(to_state), log_trans(log(prob)) {
  toState = NULL;
}
  
// ----------------------------------------------------------------------------------------
string transition::stringify() {
  return "      " + stateName + " " + double_to_string(exp(log_trans));
}
}
