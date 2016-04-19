#include "trellis.h"

namespace ham {

// ----------------------------------------------------------------------------------------
trellis::trellis(Model* hmm, Sequence seq, trellis *cached_trellis) :
  hmm_(hmm),
  cached_trellis_(cached_trellis) {
  seqs_.AddSeq(seq);
  Init();
}

// ----------------------------------------------------------------------------------------
trellis::trellis(Model* hmm, Sequences seqs, trellis *cached_trellis) :
  hmm_(hmm),
  seqs_(seqs),
  cached_trellis_(cached_trellis) {
  Init();
}

// ----------------------------------------------------------------------------------------
void trellis::Init() {
  if(cached_trellis_) {
    if(seqs_.GetSequenceLength() > cached_trellis_->seqs().GetSequenceLength())
      throw runtime_error("ERROR cached trellis sequence length " + to_string(cached_trellis_->seqs().GetSequenceLength()) + " smaller than mine " + to_string(seqs_.GetSequenceLength()));
    if(hmm_ != cached_trellis_->model())
      throw runtime_error("ERROR model in cached trellis " + cached_trellis_->model()->name() + " not the same as mine " + hmm_->name());
  }

  traceback_table_ = nullptr;
  viterbi_log_probs_ = nullptr;
  forward_log_probs_ = nullptr;
  viterbi_pointers_ = nullptr;
  swap_ptr_ = nullptr;

  ending_viterbi_log_prob_ = -INFINITY;
  ending_viterbi_pointer_ = -1;
  ending_forward_log_prob_ = -INFINITY;
}

// ----------------------------------------------------------------------------------------
trellis::~trellis() {
  if(viterbi_log_probs_)
    delete viterbi_log_probs_;
  if(forward_log_probs_)
    delete forward_log_probs_;
  if(viterbi_pointers_)
    delete viterbi_pointers_;
  if(traceback_table_ && !cached_trellis_)   // if we have a cached trellis, we used its traceback_table_, so we don't want to delete it
    delete traceback_table_;
}

// ----------------------------------------------------------------------------------------
void trellis::Dump() {
  for(size_t ipos = 0; ipos < seqs_.GetSequenceLength(); ++ipos) {
    cout
        << setw(12) << hmm_->state((*viterbi_pointers_)[ipos])->name()[0]
        << setw(12) << (*viterbi_log_probs_)[ipos];
    cout << endl;
  }
}
// ----------------------------------------------------------------------------------------
double trellis::ending_viterbi_log_prob(size_t length) {  // NOTE this is the length of the sequence, so we return position length - 1
  assert(length <= viterbi_log_probs_->size());
  return viterbi_log_probs_->at(length - 1);
}

// ----------------------------------------------------------------------------------------
double trellis::ending_forward_log_prob(size_t length) {  // NOTE this is the length of the sequence, so we return position length - 1
  assert(length <= forward_log_probs_->size());
  return forward_log_probs_->at(length - 1);
}

// ----------------------------------------------------------------------------------------
void trellis::Viterbi() {
  // initialize stored values for chunk caching
  assert(viterbi_log_probs_  == nullptr);
  assert(viterbi_pointers_ == nullptr);
  viterbi_log_probs_ = new vector<double> (seqs_.GetSequenceLength(), -INFINITY);
  viterbi_pointers_ = new vector<int> (seqs_.GetSequenceLength(), -1);

  if(cached_trellis_) {   // ok, rad, we have another trellis with the dp table already filled in, so we can just poach the values we need from there
    traceback_table_ = cached_trellis_->traceback_table();  // note that the table from the cached trellis is larger than we need right now (that's the whole point, after all)
    ending_viterbi_pointer_ = cached_trellis_->viterbi_pointer(seqs_.GetSequenceLength());
    ending_viterbi_log_prob_ = cached_trellis_->ending_viterbi_log_prob(seqs_.GetSequenceLength());
    // and also set things to allow this trellis to be passed as a cached trellis
    for(size_t ip = 0; ip < seqs_.GetSequenceLength(); ++ip) {
      (*viterbi_log_probs_)[ip] = cached_trellis_->viterbi_log_probs()->at(ip);
      (*viterbi_pointers_)[ip] = cached_trellis_->viterbi_pointers()->at(ip);
    }
    return;
  }

  assert(!traceback_table_);
  traceback_table_ = new int_2D(seqs_.GetSequenceLength(), vector<int16_t>(hmm_->n_states(), -1));

  vector<double> *scoring_current  = new vector<double> (hmm_->n_states(), -INFINITY);  // dp table values in the current column (i.e. at the current position in the query sequence)
  vector<double> *scoring_previous = new vector<double> (hmm_->n_states(), -INFINITY);  // same, but for the previous position
  bitset<STATE_MAX> next_states;  // bitset of states which we need to check at the next position (column)
  bitset<STATE_MAX> current_states;  // 

  // first calculate log probs for first position in sequence
  size_t position(0);
  State *init = hmm_->init_state();
  bitset<STATE_MAX> *initial_to_states = hmm_->initial_to_states();
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!(*initial_to_states)[i_st_current])  // skip <i_st_current> if there's no transition to it from <init>
      continue;
    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    double dpval = emission_val + init->transition_logprob(i_st_current);
    if(dpval == -INFINITY)
      continue;
    (*scoring_current)[i_st_current] = dpval;
    CacheViterbiVals(position, dpval, i_st_current);
    next_states |= (*hmm_->state(i_st_current)->to_states());  // add <i_st_current>'s outbound transitions to the list of states to check when we get to the next position (column)
  }


  // then loop over the rest of the sequence
  for(size_t position = 1; position < seqs_.GetSequenceLength(); ++position) {
    // swap <scoring_current> and <scoring_previous>, and set <scoring_current> values to -INFINITY
    swap_ptr_ = scoring_previous;
    scoring_previous = scoring_current;
    scoring_current = swap_ptr_;
    scoring_current->assign(hmm_->n_states(), -INFINITY);

    // swap the <current_states> and <next_states> bitsets (ie set current_states to the states to which we can transition from *any* of the previous states)
    current_states.reset();
    current_states |= next_states;
    next_states.reset();

    MiddleViterbi(scoring_previous, scoring_current, current_states, next_states, position);
  }

  // swap <scoring_current> and <scoring_previous>
  scoring_previous->assign(hmm_->n_states(), -INFINITY);
  swap_ptr_ = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr_;

  // NOTE now that I've got the chunk caching info, it may be possible to remove this
  // calculate ending probability and get final traceback pointer
  ending_viterbi_pointer_ = -1;
  ending_viterbi_log_prob_ = -INFINITY;
  for(size_t st_previous = 0; st_previous < hmm_->n_states(); ++st_previous) {
    if((*scoring_previous)[st_previous] == -INFINITY)
      continue;
    double dpval = (*scoring_previous)[st_previous] + hmm_->state(st_previous)->end_transition_logprob();
    if(dpval > ending_viterbi_log_prob_) {
      ending_viterbi_log_prob_ = dpval;  // NOTE should *not* be replaced by last entry in viterbi_log_probs_, since that does not include the ending transition
      ending_viterbi_pointer_ = st_previous;
    }
  }

  delete scoring_previous;
  delete scoring_current;
}

// ----------------------------------------------------------------------------------------
void trellis::MiddleViterbi(vector<double> *scoring_previous, vector<double> *scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states, size_t position) {
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!current_states[i_st_current])  // check if transition to this state is allowed from any state through which we passed at the previous position
      continue;

    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    if(emission_val == -INFINITY)
      continue;
    bitset<STATE_MAX> *from_trans = hmm_->state(i_st_current)->from_states();  // list of states from which we could've arrive at <i_st_current>

    for(size_t i_st_previous = 0; i_st_previous < hmm_->n_states(); ++i_st_previous) {
      if(!(*from_trans)[i_st_previous])
	continue;
      if((*scoring_previous)[i_st_previous] == -INFINITY)  // skip if <i_st_previous> was a dead end, i.e. that row in the previous column had zero probability
	continue;
      double dpval = (*scoring_previous)[i_st_previous] + emission_val + hmm_->state(i_st_previous)->transition_logprob(i_st_current);
      if(dpval > (*scoring_current)[i_st_current]) {
	(*scoring_current)[i_st_current] = dpval;  // save this value as the best value we've so far come across
	(*traceback_table_)[position][i_st_current] = i_st_previous;  // and mark which state it came from for later traceback
      }
      CacheViterbiVals(position, dpval, i_st_current);
      next_states |= (*hmm_->state(i_st_current)->to_states());  // TODO wait, doesn't this want to be outside the i_st_previous loop?
    }
  }
}

// ----------------------------------------------------------------------------------------
void trellis::MiddleLogProbs(vector<double> *scoring_previous, vector<double> *scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states, size_t position) {
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!current_states[i_st_current])  // check if transition to this state is allowed from any state through which we passed at the previous position
      continue;

    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    if(emission_val == -INFINITY)
      continue;
    bitset<STATE_MAX> *from_trans = hmm_->state(i_st_current)->from_states();  // list of states from which we could've arrived at <i_st_current>

    for(size_t i_st_previous = 0; i_st_previous < hmm_->n_states(); ++i_st_previous) {
      if(!(*from_trans)[i_st_previous])
	continue;
      if((*scoring_previous)[i_st_previous] == -INFINITY)  // skip if <i_st_previous> was a dead end, i.e. that row in the previous column had zero probability
	continue;
      double dpval = (*scoring_previous)[i_st_previous] + emission_val + hmm_->state(i_st_previous)->transition_logprob(i_st_current);
      (*scoring_current)[i_st_current] = AddInLogSpace(dpval, (*scoring_current)[i_st_current]);
      CacheForwardLogProb(position, dpval, i_st_current);
      next_states |= (*hmm_->state(i_st_current)->to_states());
    }
  }
}

// ----------------------------------------------------------------------------------------
void trellis::CacheViterbiVals(size_t position, double dpval, size_t i_st_current) {
    double end_trans_val = hmm_->state(i_st_current)->end_transition_logprob();
    if(dpval + end_trans_val > (*viterbi_log_probs_)[position]) {
      (*viterbi_log_probs_)[position] = dpval + end_trans_val;  // since this is the log prob of *ending* at this point, we have to add on the prob of going to the end state from this state
      (*viterbi_pointers_)[position] = i_st_current;
    }
}

// ----------------------------------------------------------------------------------------
void trellis::CacheForwardLogProb(size_t position, double dpval, size_t i_st_current) {
  double end_trans_val = hmm_->state(i_st_current)->end_transition_logprob();
  double logprob = dpval + end_trans_val;
  (*forward_log_probs_)[position] = AddInLogSpace(logprob, (*forward_log_probs_)[position]);
}

// ----------------------------------------------------------------------------------------
void trellis::Forward() {
  assert(forward_log_probs_ == nullptr);  // you can't be too careful
  forward_log_probs_ = new vector<double> (seqs_.GetSequenceLength(), -INFINITY);

  if(cached_trellis_) {
    if(!cached_trellis_->forward_log_probs())
      throw runtime_error("ERROR I got a trellis to cache that didn't have any forward information");
    ending_forward_log_prob_ = cached_trellis_->ending_forward_log_prob(seqs_.GetSequenceLength());
    for(size_t ip = 0; ip < seqs_.GetSequenceLength(); ++ip)  // and also set things to allow this trellis to be passed as a cached trellis
      (*forward_log_probs_)[ip] = cached_trellis_->forward_log_probs()->at(ip);
    return;
  }

  vector<double> *scoring_current  = new vector<double> (hmm_->n_states(), -INFINITY);  // dp table values in the current column (i.e. at the current position in the query sequence)
  vector<double> *scoring_previous = new vector<double> (hmm_->n_states(), -INFINITY);  // same, but for the previous position
  bitset<STATE_MAX> next_states;
  bitset<STATE_MAX> current_states;

  // first calculate log probs for first position in sequence
  size_t position(0);
  State *init = hmm_->init_state();
  bitset<STATE_MAX> *initial_to_states = hmm_->initial_to_states();  // states to which we can transition from the initial state
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!(*initial_to_states)[i_st_current])  // skip <i_st_current> if there's no transition to it from <init>
      continue;
    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    double dpval = emission_val + init->transition_logprob(i_st_current);
    if(dpval == -INFINITY)
      continue;
    (*scoring_current)[i_st_current] = dpval;
    next_states |= (*hmm_->state(i_st_current)->to_states());  // add <i_st_current>'s outbound transitions to the list of states to check when we get to the next column. This leaves <next_states> set to the OR of all states to which we can transition from if start from a state to which we can transition from <init>
    CacheForwardLogProb(position, dpval, i_st_current);
  }

  // then loop over the rest of the sequence
  for(position = 1; position < seqs_.GetSequenceLength(); ++position) {
    // swap <scoring_current> and <scoring_previous>, and set <scoring_current> values to -INFINITY
    swap_ptr_ = scoring_previous;
    scoring_previous = scoring_current;
    scoring_current = swap_ptr_;
    scoring_current->assign(hmm_->n_states(), -INFINITY);

    // swap the <current_states> and <next_states> bitsets (ie set current_states to the states to which we can transition from *any* of the previous states)
    current_states.reset();
    current_states |= next_states;
    next_states.reset();

    MiddleLogProbs(scoring_previous, scoring_current, current_states, next_states, position);
  }

  // swap <scoring_current> and <scoring_previous>
  scoring_previous->assign(hmm_->n_states(), -INFINITY);
  swap_ptr_ = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr_;

  ending_forward_log_prob_ = -INFINITY;
  for(size_t st_previous = 0; st_previous < hmm_->n_states(); ++st_previous) {
    if((*scoring_previous)[st_previous] == -INFINITY)
      continue;
    double dpval = (*scoring_previous)[st_previous] + hmm_->state(st_previous)->end_transition_logprob();
    if(dpval == -INFINITY)
      continue;
    ending_forward_log_prob_ = AddInLogSpace(ending_forward_log_prob_, dpval);
  }

  delete scoring_previous;
  delete scoring_current;
}

// ----------------------------------------------------------------------------------------
void trellis::Traceback(TracebackPath& path) {
  assert(seqs_.GetSequenceLength() != 0);
  assert(traceback_table_);
  assert(path.model());
  path.set_model(hmm_);
  if(ending_viterbi_log_prob_ == -INFINITY) return;  // no valid path through this hmm
  path.set_score(ending_viterbi_log_prob_);
  path.push_back(ending_viterbi_pointer_);  // push back the state that led to END state

  int16_t pointer(ending_viterbi_pointer_);
  for(size_t position = seqs_.GetSequenceLength() - 1; position > 0; position--) {
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
