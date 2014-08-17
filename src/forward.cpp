#include "trellis.h"

namespace StochHMM {
// ----------------------------------------------------------------------------------------
void trellis::forward(model* h, sequences* sqs) {
  hmm = h;
  seqs = sqs;
  seq_size = seqs->GetSequenceLength();
  state_size = hmm->state_size();
  // exDef_defined	= seqs->exDefDefined();
  forward();
}
// ----------------------------------------------------------------------------------------
void trellis::forward() {
  forward_score = new float_2D(seq_size, std::vector<float>(state_size, -INFINITY));
  scoring_current = new std::vector<double> (state_size, -INFINITY);
  scoring_previous= new std::vector<double> (state_size, -INFINITY);
		
  std::bitset<STATE_MAX> next_states;
  std::bitset<STATE_MAX> current_states;
  double  forward_temp(-INFINITY);
  double  emission(-INFINITY);
  state* init = hmm->getInitial();
  std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
  std::bitset<STATE_MAX>* from_trans(NULL);
		
  // calculate forward scores from INIT state, and initialize next_states
  for(size_t st=0; st<state_size; ++st) {
    if ((*initial_to)[st]) {  // if the bitset is set (meaning there is a transition to this state), calculate the viterbi
      double emscore = (*hmm)[st]->emission_score(*seqs, 0);
      // forward_temp = emscore + init->transition_score(st);  //init->getTrans(st)->getTransition(); //getTransition(init, st, 0);
      forward_temp = emscore + init->getTrans(st)->score(); //getTransition(init, st, 0);
      if (forward_temp > -INFINITY) {
	(*forward_score)[0][st] = forward_temp;
	(*scoring_current)[st] = forward_temp;
	next_states |= (*(*hmm)[st]->getTo());
      }
    }
  }
      
  // calculate the rest of the forward scores
  for (size_t position=1; position<seq_size; ++position) {
    // swap current and previous viterbi scores
    scoring_previous->assign(state_size,-INFINITY);
    swap_ptr = scoring_previous;
    scoring_previous = scoring_current;
    scoring_current = swap_ptr;

    //Swap current_states and next states sets
    current_states.reset();
    current_states |= next_states;
    next_states.reset();
          
    for (size_t st_current=0; st_current<state_size; ++st_current) { //i is current state that emits value
      if (!current_states[st_current])
	continue;
              
      emission = (*hmm)[st_current]->emission_score(*seqs, position);
      from_trans = (*hmm)[st_current]->getFrom();
      for (size_t previous=0; previous<state_size; ++previous) {  //j is previous state
	if (!(*from_trans)[previous])
	  continue;
					
	if ((*scoring_previous)[previous] != -INFINITY) {
	  // forward_temp = (*scoring_previous)[previous] + emission + (*hmm)[previous]->transition_score(st_current);  //(*hmm)[previous]->getTrans(st_current)->getTransition(); //getTransition((*hmm)[previous], st_current , position);
	  forward_temp = (*scoring_previous)[previous] + emission + (*hmm)[previous]->getTrans(st_current)->score(); //getTransition((*hmm)[previous], st_current , position);
	  if ((*scoring_current)[st_current] == -INFINITY) {
	    (*scoring_current)[st_current] = forward_temp;
	    (*forward_score)[position][st_current] = forward_temp;
	  } else {
	    (*scoring_current)[st_current] = AddInLogSpace(forward_temp, (*scoring_current)[st_current]);
	    (*forward_score)[position][st_current] = (*scoring_current)[st_current];
	  }
	  next_states |= (*(*hmm)[st_current]->getTo());
	}
      }
    }
  }
		
  // swap current and previous scores
  scoring_previous->assign(state_size, -INFINITY);
  swap_ptr = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr;
		
  ending_forward_prob = -INFINITY;
  for(size_t st_previous=0; st_previous<state_size; ++st_previous) {
    if ((*scoring_previous)[st_previous] != -INFINITY) {
      forward_temp = (*scoring_previous)[st_previous] + (*hmm)[st_previous]->getEndTrans();
      if (forward_temp > -INFINITY) {
	if (ending_forward_prob == -INFINITY){
	  ending_forward_prob = forward_temp;
	} else {
	  ending_forward_prob = AddInLogSpace(ending_forward_prob,forward_temp);
	}
      }
    }
  }
		
  delete scoring_previous;
  delete scoring_current;
  scoring_previous = NULL;
  scoring_current = NULL;
}
}

