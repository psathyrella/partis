#include "state.h"

namespace ham {

// ----------------------------------------------------------------------------------------
State::State() : name_(""), germline_nuc_(""), ambiguous_emission_logprob_(-INFINITY), ambiguous_char_(""), trans_to_end_(nullptr), index_(SIZE_MAX) {
  transitions_ = new vector<Transition*>;
}

// ----------------------------------------------------------------------------------------
State::~State() {
  delete transitions_;
  transitions_ = nullptr;
}

// ----------------------------------------------------------------------------------------
void State::Parse(YAML::Node node, vector<string> state_names, Tracks trks) {
  name_ = node["name"].as<string>();
  assert(name_.size() > 0);
  if(node["extras"]["germline"])
    germline_nuc_ = node["extras"]["germline"].as<string>();
  if(node["extras"]["ambiguous_emission_prob"])
    ambiguous_emission_logprob_ = log(node["extras"]["ambiguous_emission_prob"].as<double>());
  if(node["extras"]["ambiguous_char"])
    ambiguous_char_ = node["extras"]["ambiguous_char"].as<string>();

  double total(0.0); // make sure things add to 1.0
  for(YAML::const_iterator it = node["transitions"].begin(); it != node["transitions"].end(); ++it) {
    string to_state(it->first.as<string>());
    if(to_state != "end" && find(state_names.begin(), state_names.end(), to_state) == state_names.end()) {   // make sure transition is either to "end", or to a state that we know about
      cout << "ERROR attempted to add transition to unknown state \"" << to_state << "\"" << endl;
      throw runtime_error("configuration");
    }
    double prob(it->second.as<double>());
    total += prob;
    Transition *trans = new Transition(to_state, prob);
    if(trans->to_state_name() == "end")
      trans_to_end_ = trans;
    else
      transitions_->push_back(trans);
  }
  // NOTE it would be better to use something cleverer than a hard coded EPS that I just pulled ooma
  if(fabs(total - 1.0) >= EPS) { // make sure transition probs sum to 1.0
    cerr << "ERROR normalization failed on transitions in state \"" << name_ << "\"" << endl;
    cerr << node << endl;
    throw runtime_error("configuration");
  }

  // emissions
  if(name_ == "init")
    return;

  // make sure at least one emission was specified
  if(node["emissions"].IsNull()) {
    stringstream node_ss;
    node_ss << node;
    throw runtime_error("ERROR no emissions found in " + node_ss.str());
  }
  emission_.Parse(node["emissions"], trks);
}

// ----------------------------------------------------------------------------------------
void State::RescaleOverallMuteFreq(double factor) {
  if(germline_nuc_ == ambiguous_char_ || germline_nuc_ == "")  // if the germline state is N, or if this state has no germline (most likely fv or jf insertion)
    return;

  if(factor <= 0.0 || factor > 10.0)  // ten is a hack... but boy, you probably don't really want to multiply by more than 10
    throw runtime_error("ERROR State::RescaleOverallMuteFreq got a bad factor: " + to_string(factor) + "\n");

  vector<vector<double> > new_log_probs(emission_.log_probs());
  assert(new_log_probs.size() == 1);  // only support one column a.t.m.
  size_t icol(0);
  assert(new_log_probs[0].size() == emission_.track()->alphabet_size());
  for(size_t ip=0; ip<new_log_probs[icol].size(); ++ip) {
    // cout << emission_.track()->symbol(ip) << " " << exp(emission_.score(ip)) << endl;
    // TODO this is wasteful to go out of and back into log space
    double old_emit_prob = exp(new_log_probs[icol][ip]);
    double old_mute_freq;
    bool is_germline;
    // if(germline_nuc_ == ambiguous_char_ || germline_nuc_ == "") {
    //   is_germline = true;  // just arbitrarily say everything is germline if the germline is ambiguous
    // } else {
      if(emission_.track()->symbol_index(germline_nuc_) >= emission_.track()->alphabet_size())  // this'll throw an exception on the symbol_index call if the germline nuc is bad
	throw runtime_error("bad symbol");  // ...so this should never be reached
      is_germline = emission_.track()->symbol(ip) == germline_nuc_;
    // }
    if(is_germline)
      old_mute_freq = 1. - old_emit_prob;
    else
      old_mute_freq = 3. * old_emit_prob;
    double new_mute_freq = min(0.95, factor*old_mute_freq);  // .95 is kind of arbitrary, but from looking at lots of plots, the only cases where the extrapolation flies above 1.0 is where we have little information, so .95 is probably a good compromise
    if(new_mute_freq <= 0.0 || new_mute_freq >= 1.0)
      throw runtime_error("ERROR new_mute_freq not in (0,1) (" + to_string(new_mute_freq) + ") in State::RescaleOverallMuteFreq old: " + to_string(old_mute_freq) + " factor: " + to_string(factor) + " is_germline: " + to_string(is_germline));
    if(is_germline)
      new_log_probs[icol][ip] = log(1.0 - new_mute_freq);
    else
      new_log_probs[icol][ip] = log(new_mute_freq / 3.);
  }

  emission_.ReplaceLogProbs(new_log_probs);
}
// ----------------------------------------------------------------------------------------
double State::emission_logprob(Sequences *seqs, size_t pos) {
  if(seqs->n_seqs() == 1) {  // NOTE they're log probs, not scores, but I haven't yet managed to eliminate all the old 'score' names
    if(ambiguous_char_ != "" && (*seqs)[0][pos] == emission_.track()->ambiguous_index())
      return ambiguous_emission_logprob_;
    else
      return emission_.score(seqs->get_ptr(0), pos);  // get_ptr shenaniganery is to avoid performance hit from pass-by-value. Yes, I will at some point make it more elegant!
  } else {
    assert(seqs->n_seqs() > 1);

    // initialize <log_prob> for the emission from the first sequence
    double log_prob(-INFINITY);
    if(ambiguous_char_ != "" && (*seqs)[0][pos] == emission_.track()->ambiguous_index())
      log_prob = ambiguous_emission_logprob_;
    else
      log_prob = emission_.score(seqs->get_ptr(0), pos);

    // then loop over the rest of the sequences
    for(size_t iseq = 1; iseq < seqs->n_seqs(); ++iseq) {
      // add to <log_prob> the emission log prob for the <iseq>th sequence, i.e. prob1 *and* prob2
      double this_log_prob(-INFINITY);
      if(ambiguous_char_ != "" && (*seqs)[iseq][pos] == emission_.track()->ambiguous_index())
	this_log_prob = ambiguous_emission_logprob_;
      else
	this_log_prob = emission_.score(seqs->get_ptr(iseq), pos);
      log_prob = AddWithMinusInfinities(log_prob, this_log_prob);
    }

    return log_prob;
  }
}

// ----------------------------------------------------------------------------------------
void State::Print() {
  cout << "state: " << name_;
  if(germline_nuc_ != "")
    cout << " (" << germline_nuc_ << ")";
  cout << endl;

  cout << "  transitions:" << endl;;
  for(size_t i = 0; i < transitions_->size(); ++i) {
    if((*transitions_)[i] == nullptr)  // reminder: this is a bitset over all states, with the bit set (well, with a non-null-pointer set) if we can transition to the corresponding state; *not* a vector of allowed transitions
      continue;
    else
      (*transitions_)[i]->Print();
  }

  if(trans_to_end_)
    trans_to_end_->Print();

  if(name_ == "init")
    return;

  cout << "  emissions:" << endl;;
  emission_.Print();
}

// ----------------------------------------------------------------------------------------
//! Get the log probability transitioning to end from the state.
double State::end_transition_logprob() {
  if(trans_to_end_ == nullptr)
    return -INFINITY;
  return trans_to_end_->log_prob();
}

// ----------------------------------------------------------------------------------------
// On initial import of the states the allowed transitions are pushed onto <transitions_> in
// the order written in the model file. But later on we need them to be in the order specified by <index_>, in a vector of length <n_states>.
// So here we make a new vector <fixed_transitions> with the proper ordering and length <n_states> (each set to nullptr by default), and replace <transitions_> with this
// new vector.
void State::ReorderTransitions(map<string, State*> &state_indices) {
  size_t n_states(state_indices.size());
  vector<Transition*> *fixed_transitions = new vector<Transition*>(n_states - 1, nullptr); // subtract 1 because initial state is kept separate

  // find the proper place for the transition and put it in the correct position
  for(size_t i = 0; i < transitions_->size(); ++i) {  // reminder: transitions_ and fixed_transitions are not the same length
    Transition* temp = (*transitions_)[i];
    string to_state_name(temp->to_state_name());
    assert(state_indices.count(to_state_name));
    State *st(state_indices[to_state_name]);
    (*fixed_transitions)[st->index()] = temp;
  }

  delete transitions_;
  transitions_ = fixed_transitions;
}

}
