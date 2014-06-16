#include "trellis.h"

namespace StochHMM {
      
// ----------------------------------------------------------------------------------------
void trellis::viterbi(model* h, sequences* sqs) {
  hmm = h;
  seqs = sqs;
  seq_size = seqs->getLength();
  state_size = hmm->state_size();
  viterbi();
}
  
// ----------------------------------------------------------------------------------------
void trellis::viterbi() {
  assert(hmm->isBasic());

  //Initialize the traceback table
  if (traceback_table != NULL) {
    delete traceback_table;
  }
              
  traceback_table = new int_2D(seq_size,std::vector<int16_t> (state_size,-1));
  scoring_previous = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
  scoring_current  = new (std::nothrow) std::vector<double> (state_size,-INFINITY);
  assert(scoring_previous && scoring_current && traceback_table);
              
  std::bitset<STATE_MAX> next_states;
  std::bitset<STATE_MAX> current_states;
              
  double  viterbi_temp(-INFINITY);
  double  emission(-INFINITY);
  ending_viterbi_tb = -1;
  ending_viterbi_score = -INFINITY;
              
  state* init = hmm->getInitial();
              
  std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
  std::bitset<STATE_MAX>* from_trans(NULL);
              
  //Calculate Viterbi from transitions from INIT (initial) state
  for(size_t st = 0; st < state_size; ++st) {
    if ((*initial_to)[st]) {  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
      viterbi_temp = (*hmm)[st]->get_emission_prob(*seqs,0) + getTransition(init, st, 0);
      if (viterbi_temp > -INFINITY) {
	if ((*scoring_current)[st] < viterbi_temp) {  // NOTE this is always true since all we've done to scoring_current so far is initialize it to -INFINITY
	  (*scoring_current)[st] = viterbi_temp;
	}
	next_states |= (*(*hmm)[st]->getTo());  // no effect right here, but leaves next_states set to the OR of all transitions out of all states
      }
    }
  }

  // loop over the sequence
  for(size_t position=1; position<seq_size; ++position) {
    //Swap current and previous viterbi scores
    scoring_previous->assign(state_size,-INFINITY); // NOTE I think this can be replaced with
    swap_ptr = scoring_previous;		    // scoring_previous = scoring_current;		
    scoring_previous = scoring_current;	            // scoring_current->assign(state_size,-INFINITY); 
    scoring_current = swap_ptr;
                      
    //Swap current_states and next states sets. ie set current_states to the states which can be transitioned to from *any* of the previous states.
    current_states.reset();
    current_states |= next_states;
    next_states.reset();
                      
    //Current states
    for (size_t st_current = 0; st_current < state_size; ++st_current) { //Current state that emits value
      // Check to see if current state is valid. ie if transition to this state is allowed from *any* state which we passed through at the previous position
      if (!current_states[st_current])
	continue;
                              
      // Get emission of current state
      emission = (*hmm)[st_current]->get_emission_prob(*seqs, position);
                              
      if (emission == -INFINITY)  // zero probabiility, may as will stop
	continue;
                              
      // get list of states that are valid previous states
      from_trans = (*hmm)[st_current]->getFrom();
                              
      for (size_t st_previous=0; st_previous<state_size ; ++st_previous) {  //for previous states
	if (!(*from_trans)[st_previous])
	  continue;
                                      
	//Check that previous state has transition to current state
	//and that the previous viterbi score is not -INFINITY
	if ((*scoring_previous)[st_previous] != -INFINITY) {
	  viterbi_temp = getTransition((*hmm)[st_previous], st_current , position) + emission + (*scoring_previous)[st_previous];
                                              
	  if (viterbi_temp > (*scoring_current)[st_current]) {
	    (*scoring_current)[st_current] = viterbi_temp;
	    (*traceback_table)[position][st_current] = st_previous;
	  }
                                              
	  next_states |= (*(*hmm)[st_current]->getTo());
	}
      }
    }
  }
              
  //TODO:  Calculate ending and set the final viterbi and traceback pointer
  //Swap current and previous viterbi scores
  scoring_previous->assign(state_size,-INFINITY);
  swap_ptr = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr;
              
  //Calculate ending viterbi score and traceback from END state
  for(size_t st_previous = 0; st_previous < state_size ;++st_previous) {
    if ((*scoring_previous)[st_previous] > -INFINITY) {
      viterbi_temp = (*scoring_previous)[st_previous] + (*hmm)[st_previous]->getEndTrans();
                              
      if (viterbi_temp > ending_viterbi_score) {
	ending_viterbi_score = viterbi_temp;
	ending_viterbi_tb = st_previous;
      }
    }
  }
              
  delete scoring_previous;
  delete scoring_current;
  scoring_previous = NULL;
  scoring_current  = NULL;
}
}
