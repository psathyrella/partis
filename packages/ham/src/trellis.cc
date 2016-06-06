#include "trellis.h"

namespace ham {

// ----------------------------------------------------------------------------------------
double Trellis::ApproxBytesUsed() {
  double bytes(0.);
  // NOTE doesn't include traceback table!
  bytes += sizeof(double) * viterbi_log_probs_pointer_->size();
  bytes += sizeof(double) * forward_log_probs_pointer_->size();
  bytes += sizeof(int) * viterbi_indices_.size();
  return bytes;
}

// ----------------------------------------------------------------------------------------
string Trellis::SizeString() {
  char buffer[2000];
  sprintf(buffer, "%8zu  %8zu  %8zu  %8zu",
	  viterbi_log_probs_pointer_->size(),
	  forward_log_probs_pointer_->size(),
	  viterbi_indices_.size(),
	  swap_ptr_ ? swap_ptr_->size() : 0);
  return string(buffer);
}

// ----------------------------------------------------------------------------------------
Trellis::Trellis(Model* hmm, Sequence seq, Trellis *cached_trellis) :
  hmm_(hmm),
  cached_trellis_(cached_trellis),
  scoring_current_(hmm_->n_states(), -INFINITY),
  scoring_previous_(hmm_->n_states(), -INFINITY)
{
  seqs_.AddSeq(seq);
  Init();
}

// ----------------------------------------------------------------------------------------
Trellis::Trellis(Model* hmm, Sequences seqs, Trellis *cached_trellis) :
  hmm_(hmm),
  seqs_(seqs),
  cached_trellis_(cached_trellis),
  scoring_current_(hmm_->n_states(), -INFINITY),
  scoring_previous_(hmm_->n_states(), -INFINITY)
{
  Init();
}

// ----------------------------------------------------------------------------------------
Trellis::Trellis() : hmm_(nullptr), cached_trellis_(nullptr)
{
  Init();
}

// ----------------------------------------------------------------------------------------
void Trellis::Init() {
  if(cached_trellis_) {
    if(seqs_.GetSequenceLength() > cached_trellis_->seqs().GetSequenceLength())
      throw runtime_error("ERROR cached trellis sequence length " + to_string(cached_trellis_->seqs().GetSequenceLength()) + " smaller than mine " + to_string(seqs_.GetSequenceLength()));
    if(hmm_ != cached_trellis_->model())
      throw runtime_error("ERROR model in cached trellis " + cached_trellis_->model()->name() + " not the same as mine " + hmm_->name());
  }

  traceback_table_pointer_ = nullptr;
  viterbi_log_probs_pointer_ = nullptr;
  forward_log_probs_pointer_ = nullptr;
  viterbi_indices_pointer_ = nullptr;
  swap_ptr_ = nullptr;

  ending_viterbi_log_prob_ = -INFINITY;
  ending_viterbi_pointer_ = -1;
  ending_forward_log_prob_ = -INFINITY;
}

// ----------------------------------------------------------------------------------------
Trellis::~Trellis() {
}

// ----------------------------------------------------------------------------------------
void Trellis::Dump() {
  for(size_t ipos = 0; ipos < seqs_.GetSequenceLength(); ++ipos) {
    cout
        << setw(12) << hmm_->state(viterbi_indices_[ipos])->name()[0]
        << setw(12) << viterbi_log_probs_pointer_->at(ipos);
    cout << endl;
  }
}

// ----------------------------------------------------------------------------------------
void Trellis::MiddleViterbiVals(vector<double> *scoring_previous, vector<double> *scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states, size_t position) {
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!current_states[i_st_current])  // check if transition to this state is allowed from any state through which we passed at the previous position
      continue;

    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    if(emission_val == -INFINITY)
      continue;

    for(auto &i_st_previous : *hmm_->state(i_st_current)->from_state_indices()) {  // list of states from which we could've arrived at <i_st_current>
      if((*scoring_previous)[i_st_previous] == -INFINITY)  // skip if <i_st_previous> was a dead end, i.e. that row in the previous column had zero probability
	continue;
      double dpval = (*scoring_previous)[i_st_previous] + emission_val + hmm_->state(i_st_previous)->transition_logprob(i_st_current);
      if(dpval > (*scoring_current)[i_st_current]) {
	(*scoring_current)[i_st_current] = dpval;  // save this value as the best value we've so far come across
	(*traceback_table_pointer_)[position][i_st_current] = i_st_previous;  // and mark which state it came from for later traceback NOTE do *not* use <traceback_table_>, since we want the cached trellis's table if we have a cached trellis)
      }
      CacheViterbiVals(position, dpval, i_st_current);
      next_states |= (*hmm_->state(i_st_current)->to_states());  // NOTE we want this *inside* the <i_st_previous> loop because we only want to include previous states that are really needed
    }
  }
}

// ----------------------------------------------------------------------------------------
void Trellis::MiddleForwardVals(vector<double> *scoring_previous, vector<double> *scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states, size_t position) {
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!current_states[i_st_current])  // check if transition to this state is allowed from any state through which we passed at the previous position
      continue;

    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    if(emission_val == -INFINITY)
      continue;

    for(auto &i_st_previous : *hmm_->state(i_st_current)->from_state_indices()) {  // list of states from which we could've arrived at <i_st_current>
      if((*scoring_previous)[i_st_previous] == -INFINITY)  // skip if <i_st_previous> was a dead end, i.e. that row in the previous column had zero probability
	continue;
      double dpval = (*scoring_previous)[i_st_previous] + emission_val + hmm_->state(i_st_previous)->transition_logprob(i_st_current);
      (*scoring_current)[i_st_current] = AddInLogSpace(dpval, (*scoring_current)[i_st_current]);
      CacheForwardVals(position, dpval, i_st_current);
      next_states |= (*hmm_->state(i_st_current)->to_states());  // NOTE we want this *inside* the <i_st_previous> loop because we only want to include previous states that are really needed
    }
  }
}

// ----------------------------------------------------------------------------------------
void Trellis::SwapColumns(vector<double> *&scoring_previous, vector<double> *&scoring_current, bitset<STATE_MAX> &current_states, bitset<STATE_MAX> &next_states) {
  // swap <scoring_current> and <scoring_previous>, and set <scoring_current> values to -INFINITY
  swap_ptr_ = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr_;
  scoring_current->assign(hmm_->n_states(), -INFINITY);
  swap_ptr_ = nullptr;

  // swap the <current_states> and <next_states> bitsets (ie set current_states to the states to which we can transition from *any* of the previous states)
  current_states.reset();
  current_states |= next_states;
  next_states.reset();
}

// ----------------------------------------------------------------------------------------
void Trellis::CacheViterbiVals(size_t position, double dpval, size_t i_st_current) {
  double end_trans_val = hmm_->state(i_st_current)->end_transition_logprob();
  double logprob = dpval + end_trans_val;
  if(logprob > viterbi_log_probs_[position]) {
    viterbi_log_probs_[position] = logprob;  // since this is the log prob of *ending* at this point, we have to add on the prob of going to the end state from this state
    viterbi_indices_[position] = i_st_current;
  }
}

// ----------------------------------------------------------------------------------------
void Trellis::CacheForwardVals(size_t position, double dpval, size_t i_st_current) {
  double end_trans_val = hmm_->state(i_st_current)->end_transition_logprob();
  double logprob = dpval + end_trans_val;
  forward_log_probs_[position] = AddInLogSpace(logprob, forward_log_probs_[position]);
}

// ----------------------------------------------------------------------------------------
void Trellis::Viterbi() {
  if(cached_trellis_) {   // ok, rad, we have another trellis with the dp table already filled in, so we can just poach the values we need from there
    traceback_table_pointer_ = cached_trellis_->traceback_table_pointer();  // note that the table from the cached trellis is larger than we need right now (that's the whole point, after all)
    ending_viterbi_pointer_ = cached_trellis_->viterbi_pointer(seqs_.GetSequenceLength());
    ending_viterbi_log_prob_ = cached_trellis_->ending_viterbi_log_prob(seqs_.GetSequenceLength());
    viterbi_log_probs_pointer_ = cached_trellis_->viterbi_log_probs_pointer();
    viterbi_indices_pointer_ = cached_trellis_->viterbi_indices_pointer();
    return;
  }

  // initialize stored values for chunk caching
  viterbi_log_probs_.resize(seqs_.GetSequenceLength(), -INFINITY);
  viterbi_indices_.resize(seqs_.GetSequenceLength(), -1);
  viterbi_log_probs_pointer_ = &viterbi_log_probs_;
  viterbi_indices_pointer_ = &viterbi_indices_;

  traceback_table_ = int_2D(seqs_.GetSequenceLength(), vector<int16_t>(hmm_->n_states(), -1));
  traceback_table_pointer_ = &traceback_table_;

  vector<double> *scoring_current = &scoring_current_;  // dp table values in the current column (i.e. at the current position in the query sequence)
  vector<double> *scoring_previous = &scoring_previous_;  // same, but for the previous position
  scoring_current->assign(scoring_current->size(), -INFINITY);
  scoring_previous->assign(scoring_previous->size(), -INFINITY);
  bitset<STATE_MAX> next_states, current_states;  // bitset of states which we need to check at the next/current position

  // first calculate log probs for first position in sequence
  size_t position(0);
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!(*hmm_->initial_to_states())[i_st_current])  // skip <i_st_current> if there's no transition to it from <init>
      continue;
    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    double dpval = emission_val + hmm_->init_state()->transition_logprob(i_st_current);
    if(dpval == -INFINITY)
      continue;
    (*scoring_current)[i_st_current] = dpval;
    CacheViterbiVals(position, dpval, i_st_current);
    next_states |= (*hmm_->state(i_st_current)->to_states());  // add <i_st_current>'s outbound transitions to the list of states to check when we get to the next position (column)
  }


  // then loop over the rest of the sequence
  for(size_t position = 1; position < seqs_.GetSequenceLength(); ++position) {
    SwapColumns(scoring_previous, scoring_current, current_states, next_states);
    MiddleViterbiVals(scoring_previous, scoring_current, current_states, next_states, position);
  }

  SwapColumns(scoring_previous, scoring_current, current_states, next_states);

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
}

// ----------------------------------------------------------------------------------------
void Trellis::Forward() {
  if(cached_trellis_) {
    ending_forward_log_prob_ = cached_trellis_->ending_forward_log_prob(seqs_.GetSequenceLength());
    forward_log_probs_pointer_ = cached_trellis_->forward_log_probs_pointer();
    return;
  }

  // initialize stored values for chunk caching
  forward_log_probs_.resize(seqs_.GetSequenceLength(), -INFINITY);
  forward_log_probs_pointer_ = &forward_log_probs_;

  vector<double> *scoring_current = &scoring_current_;  // dp table values in the current column (i.e. at the current position in the query sequence)
  vector<double> *scoring_previous = &scoring_previous_;  // same, but for the previous position
  scoring_current->assign(scoring_current->size(), -INFINITY);
  scoring_previous->assign(scoring_previous->size(), -INFINITY);
  bitset<STATE_MAX> next_states, current_states;  // bitset of states which we need to check at the next/current position

  // first calculate log probs for first position in sequence
  size_t position(0);
  for(size_t i_st_current = 0; i_st_current < hmm_->n_states(); ++i_st_current) {
    if(!(*hmm_->initial_to_states())[i_st_current])  // skip <i_st_current> if there's no transition to it from <init>
      continue;
    double emission_val = hmm_->state(i_st_current)->EmissionLogprob(&seqs_, position);
    double dpval = emission_val + hmm_->init_state()->transition_logprob(i_st_current);
    if(dpval == -INFINITY)
      continue;
    (*scoring_current)[i_st_current] = dpval;
    next_states |= (*hmm_->state(i_st_current)->to_states());  // add <i_st_current>'s outbound transitions to the list of states to check when we get to the next column. This leaves <next_states> set to the OR of all states to which we can transition from if start from a state to which we can transition from <init>
    CacheForwardVals(position, dpval, i_st_current);
  }

  // then loop over the rest of the sequence
  for(position = 1; position < seqs_.GetSequenceLength(); ++position) {
    SwapColumns(scoring_previous, scoring_current, current_states, next_states);
    MiddleForwardVals(scoring_previous, scoring_current, current_states, next_states, position);
  }

  SwapColumns(scoring_previous, scoring_current, current_states, next_states);

  ending_forward_log_prob_ = -INFINITY;
  for(size_t st_previous = 0; st_previous < hmm_->n_states(); ++st_previous) {
    if((*scoring_previous)[st_previous] == -INFINITY)
      continue;
    double dpval = (*scoring_previous)[st_previous] + hmm_->state(st_previous)->end_transition_logprob();
    if(dpval == -INFINITY)
      continue;
    ending_forward_log_prob_ = AddInLogSpace(ending_forward_log_prob_, dpval);
  }
}

// ----------------------------------------------------------------------------------------
void Trellis::Traceback(TracebackPath& path) {
  assert(seqs_.GetSequenceLength() != 0);
  assert(path.model());
  path.set_model(hmm_);
  if(ending_viterbi_log_prob_ == -INFINITY) return;  // no valid path through this hmm
  path.set_score(ending_viterbi_log_prob_);
  path.push_back(ending_viterbi_pointer_);  // push back the state that led to END state

  int16_t pointer(ending_viterbi_pointer_);
  for(size_t position = seqs_.GetSequenceLength() - 1; position > 0; position--) {
    pointer = (*traceback_table_pointer_)[position][pointer];  // NOTE do *not* use <traceback_table_>, since we want the cached trellis's table if we have a cached trellis)
    if(pointer == -1) {
      cerr << "No valid path at Position: " << position << endl;
      return;
    }
    path.push_back(pointer);
  }
  assert(path.size() > 0);  // NOTE don't remove this! dphandler assumes paths are invalid/not set if path size is zero
}
}
