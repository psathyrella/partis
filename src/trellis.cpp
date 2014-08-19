#include "trellis.h"

namespace StochHMM {
      
// ----------------------------------------------------------------------------------------
trellis::trellis(model* h, sequences* sqs){
  hmm = h;
  seqs = sqs;
              
  traceback_table = NULL;
  forward_score = NULL;
  ending_viterbi_tb = -1;
  ending_forward_prob = -INFINITY;

  scoring_current = NULL;
  scoring_previous = NULL;
}

// ----------------------------------------------------------------------------------------
trellis::~trellis() {
  delete traceback_table;
  delete forward_score;
  delete scoring_previous;
  delete scoring_current;
}
      
// ----------------------------------------------------------------------------------------
void trellis::viterbi() {
  //Initialize the traceback table
  if (traceback_table != NULL) {
    delete traceback_table;
  }
              
  traceback_table = new int_2D(seqs->GetSequenceLength(),std::vector<int16_t> (hmm->n_states(),-1));
  scoring_previous = new (std::nothrow) std::vector<double> (hmm->n_states(),-INFINITY);
  scoring_current  = new (std::nothrow) std::vector<double> (hmm->n_states(),-INFINITY);
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
  for(size_t st = 0; st < hmm->n_states(); ++st) {
    if ((*initial_to)[st]) {  //if the bitset is set (meaning there is a transition to this state), calculate the viterbi
      viterbi_temp = (*hmm)[st]->emission_score(*seqs,0) + init->transition_score(st);
      if (viterbi_temp > -INFINITY) {
	if ((*scoring_current)[st] < viterbi_temp) {  // NOTE this is always true since all we've done to scoring_current so far is initialize it to -INFINITY
	  (*scoring_current)[st] = viterbi_temp;
	}
	next_states |= (*(*hmm)[st]->getTo());  // no effect right here, but leaves next_states set to the OR of all transitions out of all states
      }
    }
  }

  // loop over the sequence
  for(size_t position=1; position<seqs->GetSequenceLength(); ++position) {
    //Swap current and previous viterbi scores
    scoring_previous->assign(hmm->n_states(),-INFINITY); // NOTE I think this can be replaced with
    swap_ptr = scoring_previous;		    // scoring_previous = scoring_current;		
    scoring_previous = scoring_current;	            // scoring_current->assign(hmm->n_states(),-INFINITY); EDIT hmm, nope, doesn't seem to work
    scoring_current = swap_ptr;
                      
    //Swap current_states and next states sets. ie set current_states to the states which can be transitioned to from *any* of the previous states.
    current_states.reset();
    current_states |= next_states;
    next_states.reset();
                      
    //Current states
    for (size_t st_current = 0; st_current < hmm->n_states(); ++st_current) { //Current state that emits value
      // Check to see if current state is valid. ie if transition to this state is allowed from *any* state which we passed through at the previous position
      if (!current_states[st_current])
	continue;
                              
      // Get emission of current state
      emission = (*hmm)[st_current]->emission_score(*seqs, position);
                              
      if (emission == -INFINITY)  // zero probabiility, may as will stop
	continue;
                              
      // get list of states that are valid previous states
      from_trans = (*hmm)[st_current]->getFrom();
                              
      for (size_t st_previous=0; st_previous<hmm->n_states() ; ++st_previous) {  //for previous states
	if (!(*from_trans)[st_previous])
	  continue;
                                      
	//Check that previous state has transition to current state
	//and that the previous viterbi score is not -INFINITY
	if ((*scoring_previous)[st_previous] != -INFINITY) {
	  viterbi_temp = (*hmm)[st_previous]->transition_score(st_current) + emission + (*scoring_previous)[st_previous];
                                              
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
  scoring_previous->assign(hmm->n_states(),-INFINITY);
  swap_ptr = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr;
              
  //Calculate ending viterbi score and traceback from END state
  for(size_t st_previous = 0; st_previous < hmm->n_states() ;++st_previous) {
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

// ----------------------------------------------------------------------------------------
void trellis::forward() {
  forward_score = new float_2D(seqs->GetSequenceLength(), std::vector<float>(hmm->n_states(), -INFINITY));
  scoring_current = new std::vector<double> (hmm->n_states(), -INFINITY);
  scoring_previous= new std::vector<double> (hmm->n_states(), -INFINITY);
		
  std::bitset<STATE_MAX> next_states;
  std::bitset<STATE_MAX> current_states;
  double  forward_temp(-INFINITY);
  double  emission(-INFINITY);
  state* init = hmm->getInitial();
  std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
  std::bitset<STATE_MAX>* from_trans(NULL);
		
  // calculate forward scores from INIT state, and initialize next_states
  for(size_t st=0; st<hmm->n_states(); ++st) {
    if ((*initial_to)[st]) {  // if the bitset is set (meaning there is a transition to this state), calculate the viterbi
      double emscore = (*hmm)[st]->emission_score(*seqs, 0);
      forward_temp = emscore + init->getTrans(st)->score();
      if (forward_temp > -INFINITY) {
	(*forward_score)[0][st] = forward_temp;
	(*scoring_current)[st] = forward_temp;
	next_states |= (*(*hmm)[st]->getTo());
      }
    }
  }
      
  // calculate the rest of the forward scores
  for (size_t position=1; position<seqs->GetSequenceLength(); ++position) {
    // swap current and previous viterbi scores
    scoring_previous->assign(hmm->n_states(),-INFINITY);
    swap_ptr = scoring_previous;
    scoring_previous = scoring_current;
    scoring_current = swap_ptr;

    //Swap current_states and next states sets
    current_states.reset();
    current_states |= next_states;
    next_states.reset();
          
    for (size_t st_current=0; st_current<hmm->n_states(); ++st_current) { //i is current state that emits value
      if (!current_states[st_current])
	continue;
              
      emission = (*hmm)[st_current]->emission_score(*seqs, position);
      from_trans = (*hmm)[st_current]->getFrom();
      for (size_t previous=0; previous<hmm->n_states(); ++previous) {  //j is previous state
	if (!(*from_trans)[previous])
	  continue;
					
	if ((*scoring_previous)[previous] != -INFINITY) {
	  forward_temp = (*scoring_previous)[previous] + emission + (*hmm)[previous]->transition_score(st_current);
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
  scoring_previous->assign(hmm->n_states(), -INFINITY);
  swap_ptr = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr;
		
  ending_forward_prob = -INFINITY;
  for(size_t st_previous=0; st_previous<hmm->n_states(); ++st_previous) {
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

// ----------------------------------------------------------------------------------------
//!Perform traceback through trellis
//!\return path trackback_path
void trellis::traceback(traceback_path& path) {
  assert(seqs->GetSequenceLength() != 0);
  assert(traceback_table);
  if (path.getModel() == NULL)
    path.setModel(hmm);
  if(ending_viterbi_score == -INFINITY) return;  // no valid path through this hmm
  path.setScore(ending_viterbi_score);
  path.push_back(ending_viterbi_tb);  // push back the state that led to END state
          
  int16_t pointer = ending_viterbi_tb;
  for(size_t position=seqs->GetSequenceLength()-1; position>0; position--) {
    pointer = (*traceback_table)[position][pointer];
    if (pointer == -1){
	cerr << "No valid path at Position: " << position << endl;
	return;
    }
    path.push_back(pointer);
  }
  assert(path.size() > 0);
}
}
