#include "transitions.h"
namespace ham {
  
// ----------------------------------------------------------------------------------------
transition::transition(string to_state, double prob) : stateName(to_state), log_trans(log(prob)) {
  toState = NULL;
}
  
// ----------------------------------------------------------------------------------------
void transition::print() {
  printf("%30s%14.3e\n", stateName.c_str(), exp(log_trans));
}
}
