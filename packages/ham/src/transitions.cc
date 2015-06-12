#include "transitions.h"
namespace ham {

// ----------------------------------------------------------------------------------------
Transition::Transition(string to_state, double prob) :
  to_state_name_(to_state),
  log_prob_(log(prob))
{
  to_state_ = nullptr;
}

// ----------------------------------------------------------------------------------------
void Transition::Print() {
  printf("%30s%14.3e\n", to_state_name_.c_str(), exp(log_prob_));
}

}
