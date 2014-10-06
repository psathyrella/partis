#include "transitions.h"
namespace ham {
  
// ----------------------------------------------------------------------------------------
transition::transition(string to_state, double prob) : stateName(to_state), log_trans(log(prob)) {
  toState = NULL;
}
  
// ----------------------------------------------------------------------------------------
void transition::print() {
  cout << "      " << stateName << "     " << double_to_string(exp(log_trans)) << endl;;
}
}
