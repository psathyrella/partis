#include "transitions.h"
namespace stochhmm {
  
// ----------------------------------------------------------------------------------------
transition::transition(string to_state, double prob) : stateName(to_state), log_trans(log(prob)) {
  toState = NULL;
}
  
// ----------------------------------------------------------------------------------------
string transition::stringify() {
  string transString;
  transString+="\t" + stateName + ":\t" + double_to_string(log_trans);
  return transString;
}
}
