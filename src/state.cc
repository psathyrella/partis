#include "state.h"

namespace ham {

// ----------------------------------------------------------------------------------------
State::State() : trans_to_end_(nullptr), index_(SIZE_MAX) {
  transitions_ = new(nothrow) vector<Transition*>;
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
double State::emission_logprob(Sequences *seqs, size_t pos) {
  if(seqs->n_seqs() == 1) {
    // NOTE they're log probs, not scores, but I haven't yet managed to eliminate all the old 'score' names
    return emission_.score(seqs->get_ptr(0), pos);  // get_ptr shenaniganery is to avoid performance hit from pass-by-value. Yes, I will at some point make it more elegant!
  } else {
    assert(seqs->n_seqs() > 1);
    double log_prob = emission_.score(seqs->get_ptr(0), pos);  // initialize <log_prob> for the emission from the first sequence
    for(size_t iseq = 1; iseq < seqs->n_seqs(); ++iseq)   // then loop over the rest of the sequences
      // add to <log_prob> the emission log prob for the <iseq>th sequence, i.e. prob1 *and* prob2
      log_prob = AddWithMinusInfinities(log_prob, emission_.score(seqs->get_ptr(iseq), pos));
    return log_prob;
  }
}

// ----------------------------------------------------------------------------------------
void State::Print() {
  cout << "state: " << name_ << endl;

  cout << "  transitions:" << endl;;
  for(size_t i = 0; i < transitions_->size(); ++i) {
    if((*transitions_)[i] == nullptr) { assert(0); continue;} // wait wtf would this happen?
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
// On initial import of the states they are pushed onto <transitions_> in
// the order written in the model file. But later on we need them to be in the order specified by <index_>.
// So here we make a new vector <fixed_transitions> with the proper ordering and replace <transitions_> with this
// new vector.
void State::ReorderTransitions(map<string, State*> &state_indices) {
  size_t n_states(state_indices.size());
  vector<Transition*> *fixed_transitions = new vector<Transition*>(n_states - 1, nullptr); // subtract 1 because initial state is kept separate

  // find the proper place for the transition and put it in the correct position
  for(size_t i = 0; i < transitions_->size(); ++i) {
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
