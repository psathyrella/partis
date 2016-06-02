#include "state.h"

namespace ham {

// ----------------------------------------------------------------------------------------
State::State() :
  name_(""),
  germline_nuc_(""),
  ambiguous_emission_logprob_(-INFINITY),
  ambiguous_char_(""),
  trans_to_end_(nullptr),
  index_(SIZE_MAX)
{
  transitions_ = new vector<Transition*>;
}

// ----------------------------------------------------------------------------------------
State::~State() {
  for(size_t it=0; it<transitions_->size(); ++it)
    delete (*transitions_)[it];
  delete transitions_;
  transitions_ = nullptr;
  if(trans_to_end_ != nullptr)
    delete trans_to_end_;
}

// ----------------------------------------------------------------------------------------
void State::Parse(YAML::Node node, vector<string> state_names, Track *track) {
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
    throw runtime_error("no emissions found in " + node_ss.str());
  }
  emission_.Parse(node["emissions"], track);
}

// ----------------------------------------------------------------------------------------
void State::RescaleOverallMuteFreq(double factor) {
  if(germline_nuc_ == ambiguous_char_ || germline_nuc_ == "")  // if the germline state is N, or if this state has no germline (most likely fv or jf insertion)
    return;

  if(factor <= 0.0 || factor > 15.0)  // ten is a hack... but boy, you probably don't really want to multiply by more than 10
    cout << "very large factor in State::RescaleOverallMuteFreq: " << to_string(factor) << endl;

  assert(emission_.track()->symbol_index(germline_nuc_) < emission_.track()->alphabet_size());  // this'll throw an exception on the symbol_index call if the germline nuc is bad

  vector<double> new_log_probs(emission_.log_probs());
  assert(new_log_probs.size() == emission_.track()->alphabet_size());
  for(size_t ip=0; ip<new_log_probs.size(); ++ip) {
    // NOTE this is wasteful to go out of and back into log space (but doesn't matter at all in actual practice)
    double old_emit_prob = exp(new_log_probs[ip]);
    bool is_germline(emission_.track()->symbol(ip) == germline_nuc_);
    double old_mute_freq;
    if(is_germline)
      old_mute_freq = 1. - old_emit_prob;
    else
      old_mute_freq = 3. * old_emit_prob;
    double new_mute_freq = min(0.95, factor*old_mute_freq);  // .95 is kind of arbitrary, but from looking at lots of plots, the only cases where the extrapolation flies above 1.0 is where we have little information, so .95 is probably a good compromise
    if(new_mute_freq <= 0.0 || new_mute_freq >= 1.0)
      throw runtime_error("ERROR new_mute_freq not in (0,1) (" + to_string(new_mute_freq) + ") in State::RescaleOverallMuteFreq old: "
			  + to_string(old_mute_freq) + " factor: " + to_string(factor) + " is_germline: " + to_string(is_germline));
    if(is_germline)
      new_log_probs[ip] = log(1.0 - new_mute_freq);
    else
      new_log_probs[ip] = log(new_mute_freq / 3.);
  }

  emission_.ReplaceLogProbs(new_log_probs);
}

// ----------------------------------------------------------------------------------------
void State::UnRescaleOverallMuteFreq() {
  if(germline_nuc_ == ambiguous_char_ || germline_nuc_ == "")  // if the germline state is N, or if this state has no germline (most likely fv or jf insertion)
    return;
  emission_.UnReplaceLogProbs();
}

// ----------------------------------------------------------------------------------------
double State::EmissionLogprob(uint8_t ch) {
  if(ambiguous_char_ != "" && ch == emission_.track()->ambiguous_index())
    return ambiguous_emission_logprob_;
  else
    return emission_.score(ch);
}

// ----------------------------------------------------------------------------------------
double State::EmissionLogprob(Sequences *seqs, size_t pos) {
  double logprob(0.);  // multiplying probabilities, so initial prob value should be 1.
  for(size_t iseq=0; iseq<seqs->n_seqs(); ++iseq)
    logprob = AddWithMinusInfinities(logprob, EmissionLogprob((*seqs->get_ptr(iseq))[pos]));
  return logprob;
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
  vector<Transition*> *fixed_transitions = new vector<Transition*>(n_states - 1, nullptr);  // subtract 1 because initial state is kept separate

  // find the proper place for the transition and put it in the correct position
  for(size_t i = 0; i < transitions_->size(); ++i) {  // reminder: transitions_ and fixed_transitions are not the same length
    Transition* tmp_trans = (*transitions_)[i];
    string to_state_name(tmp_trans->to_state_name());
    assert(state_indices.count(to_state_name));
    State *to_state(state_indices[to_state_name]);
    (*fixed_transitions)[to_state->index()] = tmp_trans;
  }

  delete transitions_;
  transitions_ = fixed_transitions;
}

// ----------------------------------------------------------------------------------------
void State::SetFromStateIndices() {
  for(size_t istate=0; istate<from_states_.size(); ++istate)
    if(from_states_[istate])
      from_state_indices_.push_back(istate);
}

}
