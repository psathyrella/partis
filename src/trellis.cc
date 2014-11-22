#include "trellis.h"

namespace ham {

// ----------------------------------------------------------------------------------------
trellis::trellis(Model* hmm, Sequence *seq) :
  hmm_(hmm),
  seqs_(NULL) {
  seqs_ = new Sequences;
  seqs_->AddSeq(seq);
  Init();
}

// ----------------------------------------------------------------------------------------
trellis::trellis(Model* hmm, Sequences *seqs) :
  hmm_(hmm),
  seqs_(seqs) {
  Init();
}

// ----------------------------------------------------------------------------------------
void trellis::Init() {
  traceback_table_ = NULL;
  ending_viterbi_log_prob_ = -INFINITY;
  ending_viterbi_pointer_ = -1;
  forward_table_ = NULL;
  ending_forward_log_prob_ = -INFINITY;
  scoring_current_ = NULL;
  scoring_previous_ = NULL;
  swap_ptr_ = NULL;
}

// ----------------------------------------------------------------------------------------
trellis::~trellis() {
  delete traceback_table_;
  delete forward_table_;
  delete scoring_previous_;
  delete scoring_current_;
}

// ----------------------------------------------------------------------------------------
void trellis::Viterbi() {
  viterbi_log_probs_ = new vector<double> (seqs_->GetSequenceLength(), -INFINITY);  // log prob of best path up to (and including) each position in the sequence

  if(traceback_table_ != nullptr)
    delete traceback_table_;
  traceback_table_ = new int_2D(seqs_->GetSequenceLength(), vector<int16_t>(hmm_->n_states(), -1));
  scoring_current_  = new vector<double> (hmm_->n_states(), -INFINITY);  // viterbi values in the current column (i.e. at the current position in the query sequence)
  scoring_previous_ = new vector<double> (hmm_->n_states(), -INFINITY);  // same, but for the previous position

  bitset<STATE_MAX> next_states;
  bitset<STATE_MAX> current_states;

  // calculate Viterbi from transitions from INIT (initial) state
  State *init = hmm_->init_state();
  bitset<STATE_MAX>* initial_to_states = hmm_->initial_to_states();
  for(size_t st = 0; st < hmm_->n_states(); ++st) {
    if(!(*initial_to_states)[st])  // skip <st> if there's no transition to it from <init>
      continue;
    double emission_val = hmm_->state(st)->emission_logprob(*seqs_, 0);  // zeroth position in sequence
    double viterbi_val = emission_val + init->transition_logprob(st);
    if(viterbi_val > -INFINITY) {
      (*scoring_current_)[st] = viterbi_val;
      (*viterbi_log_probs_)[0] = viterbi_val;
      next_states |= (*hmm_->state(st)->to_states());  // add <st>'s outbound transitions to the list of states to check when we get to the next column
                                                       // this leaves <next_states> set to the OR of all states to which we can transition from if start from a state to which we can transition from <init>
    }
  }

  // loop over the rest of the sequence
  bitset<STATE_MAX>* from_trans(NULL);
  for(size_t position = 1; position < seqs_->GetSequenceLength(); ++position) {
    // swap <scoring_current_> and <scoring_previous_>
    scoring_previous_->assign(hmm_->n_states(), -INFINITY); // NOTE I think this can be replaced with
    swap_ptr_ = scoring_previous_;		            // scoring_previous_ = scoring_current_;
    scoring_previous_ = scoring_current_;	            // scoring_current_->assign(hmm_->n_states(),-INFINITY); EDIT nope, doesn't seem to work. arg. dunno why.
    scoring_current_ = swap_ptr_;

    // swap <current_states> and <next_states> sets. ie set current_states to the states to which we can transition from *any* of the previous states.
    current_states.reset();
    current_states |= next_states;
    next_states.reset();

    // current states
    for(size_t st_current = 0; st_current < hmm_->n_states(); ++st_current) {
      if(!current_states[st_current])  // check if transition to this state is allowed from any state through which we passed at the previous position
        continue;

      double emission_val = hmm_->state(st_current)->emission_logprob(*seqs_, position);
      if(emission_val == -INFINITY)
        continue;

      
      from_trans = hmm_->state(st_current)->from_states();  // list of states from which we could've arrive at <st_current>
      for(size_t st_previous = 0; st_previous < hmm_->n_states() ; ++st_previous) { //for previous states
        if(!(*from_trans)[st_previous])
          continue;
        if((*scoring_previous_)[st_previous] == -INFINITY)  // skip if <st_previous> was a dead end, i.e. that row in the previous column had zero probability
	  continue;
	double viterbi_val = (*scoring_previous_)[st_previous] + emission_val + hmm_->state(st_previous)->transition_logprob(st_current);
	if(viterbi_val > (*scoring_current_)[st_current]) {
	  (*scoring_current_)[st_current] = viterbi_val;  // save this value as the best value we've so far come across
	  (*viterbi_log_probs_)[position] = viterbi_val;
	  (*traceback_table_)[position][st_current] = st_previous;  // and mark which state it came from for later traceback
	}
	next_states |= (*hmm_->state(st_current)->to_states());
      }
    }
  }

  // swap <scoring_current_> and <scoring_previous_>
  scoring_previous_->assign(hmm_->n_states(), -INFINITY);
  swap_ptr_ = scoring_previous_;
  scoring_previous_ = scoring_current_;
  scoring_current_ = swap_ptr_;

  // calculate ending viterbi score and traceback from END state
  ending_viterbi_pointer_ = -1;
  ending_viterbi_log_prob_ = -INFINITY;
  for(size_t st_previous = 0; st_previous < hmm_->n_states() ; ++st_previous) {
    if((*scoring_previous_)[st_previous] == -INFINITY)
      continue;
    double viterbi_val = (*scoring_previous_)[st_previous] + hmm_->state(st_previous)->end_transition_logprob();
    if(viterbi_val > ending_viterbi_log_prob_) {
      ending_viterbi_log_prob_ = viterbi_val;
      (*viterbi_log_probs_)[seqs_->GetSequenceLength() - 1] = viterbi_val;
      ending_viterbi_pointer_ = st_previous;
    }
  }

  delete scoring_previous_;
  delete scoring_current_;
  scoring_previous_ = NULL;
  scoring_current_  = NULL;
}

// ----------------------------------------------------------------------------------------
void trellis::Forward() {
  forward_table_ = new float_2D(seqs_->GetSequenceLength(), vector<float>(hmm_->n_states(), -INFINITY));
  scoring_current_ = new vector<double> (hmm_->n_states(), -INFINITY);
  scoring_previous_ = new vector<double> (hmm_->n_states(), -INFINITY);

  bitset<STATE_MAX> next_states;
  bitset<STATE_MAX> current_states;
  double  forward_temp(-INFINITY);
  double  emission(-INFINITY);
  State* init = hmm_->init_state();
  bitset<STATE_MAX>* initial_to = hmm_->initial_to_states();
  bitset<STATE_MAX>* from_trans(NULL);

  // calculate forward scores from INIT state, and initialize next_states
  for(size_t st = 0; st < hmm_->n_states(); ++st) {
    if((*initial_to)[st]) {   // if the bitset is set (meaning there is a transition to this state), calculate the viterbi
      double emscore = hmm_->state(st)->emission_logprob(*seqs_, 0);
      forward_temp = emscore + init->transition(st)->log_prob();
      if(forward_temp > -INFINITY) {
        (*forward_table_)[0][st] = forward_temp;
        (*scoring_current_)[st] = forward_temp;
        next_states |= (*hmm_->state(st)->to_states());
      }
    }
  }

  // calculate the rest of the forward scores
  for(size_t position = 1; position < seqs_->GetSequenceLength(); ++position) {
    // swap current and previous viterbi scores
    scoring_previous_->assign(hmm_->n_states(), -INFINITY);
    swap_ptr_ = scoring_previous_;
    scoring_previous_ = scoring_current_;
    scoring_current_ = swap_ptr_;

    //Swap current_states and next states sets
    current_states.reset();
    current_states |= next_states;
    next_states.reset();

    for(size_t st_current = 0; st_current < hmm_->n_states(); ++st_current) { //i is current state that emits value
      if(!current_states[st_current])
        continue;

      emission = hmm_->state(st_current)->emission_logprob(*seqs_, position);
      from_trans = hmm_->state(st_current)->from_states();
      for(size_t previous = 0; previous < hmm_->n_states(); ++previous) { //j is previous state
        if(!(*from_trans)[previous])
          continue;

        if((*scoring_previous_)[previous] != -INFINITY) {
          forward_temp = (*scoring_previous_)[previous] + emission + hmm_->state(previous)->transition_logprob(st_current);
          if((*scoring_current_)[st_current] == -INFINITY) {
            (*scoring_current_)[st_current] = forward_temp;
            (*forward_table_)[position][st_current] = forward_temp;
          } else {
            (*scoring_current_)[st_current] = AddInLogSpace(forward_temp, (*scoring_current_)[st_current]);
            (*forward_table_)[position][st_current] = (*scoring_current_)[st_current];
          }
          next_states |= (*hmm_->state(st_current)->to_states());
        }
      }
    }
  }

  // swap current and previous scores
  scoring_previous_->assign(hmm_->n_states(), -INFINITY);
  swap_ptr_ = scoring_previous_;
  scoring_previous_ = scoring_current_;
  scoring_current_ = swap_ptr_;

  ending_forward_log_prob_ = -INFINITY;
  for(size_t st_previous = 0; st_previous < hmm_->n_states(); ++st_previous) {
    if((*scoring_previous_)[st_previous] != -INFINITY) {
      forward_temp = (*scoring_previous_)[st_previous] + hmm_->state(st_previous)->end_transition_logprob();
      if(forward_temp > -INFINITY) {
        if(ending_forward_log_prob_ == -INFINITY) {
          ending_forward_log_prob_ = forward_temp;
        } else {
          ending_forward_log_prob_ = AddInLogSpace(ending_forward_log_prob_, forward_temp);
        }
      }
    }
  }

  delete scoring_previous_;
  delete scoring_current_;
  scoring_previous_ = NULL;
  scoring_current_ = NULL;
}

// ----------------------------------------------------------------------------------------
//!Perform traceback through trellis
//!\return path trackback_path
void trellis::Traceback(TracebackPath& path) {
  assert(seqs_->GetSequenceLength() != 0);
  assert(traceback_table_);
  assert(path.model());
  path.set_model(hmm_);
  if(ending_viterbi_log_prob_ == -INFINITY) return;  // no valid path through this hmm
  path.set_score(ending_viterbi_log_prob_);
  path.push_back(ending_viterbi_pointer_);  // push back the state that led to END state

  int16_t pointer(ending_viterbi_pointer_);
  for(size_t position = seqs_->GetSequenceLength() - 1; position > 0; position--) {
    pointer = (*traceback_table_)[position][pointer];
    if(pointer == -1) {
      cerr << "No valid path at Position: " << position << endl;
      return;
    }
    path.push_back(pointer);
  }
  assert(path.size() > 0);
}
}
