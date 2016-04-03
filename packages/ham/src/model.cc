#include "model.h"

namespace ham {
// ----------------------------------------------------------------------------------------
Model::Model() :
  overall_prob_(0.0),
  original_overall_mute_freq_(0.0),
  rescale_ratio_(-INFINITY),
  ambiguous_char_(""),
  track_(nullptr),
  initial_(nullptr),
  finalized_(false)
{
  ending_ = new State;
}

// ----------------------------------------------------------------------------------------
Model::~Model() {
  delete ending_;
  ending_ = nullptr;
  for(auto &kv : states_by_name_)
    delete kv.second;
  delete track_;
}

// ----------------------------------------------------------------------------------------
void Model::Parse(string infname) {
  if(!ifstream(infname))
    throw runtime_error("input file " + infname + " does not exist.");

  // load yaml
  YAML::Node config = YAML::LoadFile(infname);
  // first get model-wide information
  try {
    name_ = config["name"].as<string>();
    if(config["extras"] && config["extras"]["gene_prob"])
      overall_prob_ = config["extras"]["gene_prob"].as<double>();
    if(config["extras"] && config["extras"]["overall_mute_freq"])
      original_overall_mute_freq_ = config["extras"]["overall_mute_freq"].as<double>();
    if(config["extras"]["ambiguous_char"])
      ambiguous_char_ = config["extras"]["ambiguous_char"].as<string>();
  } catch(...) {
    cerr << "ERROR invalid model header info in " << infname << endl;
    throw;
  }

  try {
    // and the tracks
    YAML::Node tracks(config["tracks"]);
    assert(tracks.size() == 1);  // don't at the moment support multiple tracks
    for(YAML::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
      assert(track_ == nullptr);  // shouldn't already be initialized
      track_ = new Track;
      track_->set_name(it->first.as<string>());
      if(ambiguous_char_ != "")
	track_->SetAmbiguous(ambiguous_char_);
      for(size_t ic = 0; ic < it->second.size(); ++ic)
        track_->AddSymbol(it->second[ic].as<string>());
    }
  } catch(...) {
    cerr << "ERROR invalid track specifications in " << infname << endl;
    throw;
  }

  // then push back each state name
  vector<string> state_names;
  for(size_t is = 0; is < config["states"].size(); ++is) {
    string name;
    try {
      name = config["states"][is]["name"].as<string>();
    } catch(...) {
      cerr << "ERROR invalid state name in " << infname << endl;
      throw;
    }
    for(auto chkname : state_names) {
      if(name == chkname) {
        cerr << "ERROR added two states with name '" << name << "'" << endl;
        throw;
      }
    }
    state_names.push_back(name);
  }

  // then actually parse the info for each state
  for(size_t ist = 0; ist < state_names.size(); ++ist) {
    State *state(new State);
    try {
      state->Parse(config["states"][ist], state_names, track_);
    } catch(...) {
      cerr << "ERROR invalid specification for state '" << state_names[ist] << "' in " << infname << endl;
      throw;
    }

    if(state->name() == "init") {
      initial_ = state;
    } else {
      assert(states_.size() < STATE_MAX);
      states_.push_back(state);
    }
    states_by_name_[state->name()] = state;
  }

  Finalize(); // post process states and/to create an end state with only transitions-from
}

// ----------------------------------------------------------------------------------------
void Model::AddState(State* state) {
  throw runtime_error("do I ever get here?");
  assert(states_.size() < STATE_MAX);
  states_.push_back(state);
  states_by_name_[state->name()] = state;
  return;
}

// ----------------------------------------------------------------------------------------
void Model::RescaleOverallMuteFreq(double overall_mute_freq) {
  assert(overall_mute_freq != -INFINITY);
  // cout << "rescaling " << name_ << " from " << original_overall_mute_freq_ << " to " << overall_mute_freq << endl;
  for(auto &state : states_) {
    // NOTE it is arguable that the denominator here should be the original mute freq only over the sequences that had
    //  *this* germline gene (rather than over all sequence in the data set). However, it'd be a bunch more work to do
    //  it that way, and even if it's more correcter, I don't think it'd make much difference
    double factor = max(0.01, overall_mute_freq) / original_overall_mute_freq_;  // NOTE the 1% is kind of a hack (to protect against zero) -- but it's roughly equal to the uncertainty on our mute freq estimates, so it's reasonable
    state->RescaleOverallMuteFreq(factor);  // REMINDER still not in log space
  }
}

// ----------------------------------------------------------------------------------------
void Model::UnRescaleOverallMuteFreq() {
  // cout << "  unrescaling" << endl;
  for(auto &state : states_)
    state->UnRescaleOverallMuteFreq();
}

// ----------------------------------------------------------------------------------------
// set transitions and perform some other checks
void Model::Finalize() {
  assert(!finalized_);  // well it wouldn't *hurt* to call this twice, but you still *oughtn't* to

  // set each state's index within this model
  for(size_t i = 0; i < states_.size(); ++i)
    states_[i]->SetIndex(i);
  // set transition information (NOTE does not only modify states_[i])
  for(size_t i = 0; i < states_.size(); ++i)
    FinalizeState(states_[i]);
  if(!initial_)
    throw runtime_error("ERROR no 'init' state was specified");
  FinalizeState(initial_);

  // reorder state::transitions_ so it agrees with state:index_
  for(size_t i = 0; i < states_.size(); ++i)
    states_[i]->ReorderTransitions(states_by_name_);
  initial_->ReorderTransitions(states_by_name_);

  CheckTopology();

  finalized_ = true;
}


// ----------------------------------------------------------------------------------------
void Model::FinalizeState(State *st) {
  // Modiry to_state_ and from_state_ bitsets in <st> and its transition partners
  vector<Transition*>* transitions(st->transitions());
  for(size_t it = 0; it < transitions->size(); ++it) { // loops over the transitions out of <st>
    string to_state_name(transitions->at(it)->to_state_name());
    assert(states_by_name_.count(to_state_name));
    State *to_state(states_by_name_[to_state_name]);

    st->AddToState(to_state); // add <to_state> to the list of states which can be reached from <st>
    transitions->at(it)->set_to_state(to_state);  // set to_state pointer for the it'th transition
    if(st != initial_)
      to_state->AddFromState(st);  // add <st> to the list of states from which you can reach <to_state>
  }

  if(st->trans_to_end())
    ending_->AddFromState(st);
}

// ----------------------------------------------------------------------------------------
void Model::CheckTopology() {
  // check for states with
  //   - zero outbound transitions
  //   - only self-transitions
  //   - zero inbound transitions
  // make sure there's an end state
  // make sure all states are reachable from init

  vector<uint16_t> states_to_check;  // Dynamic vector of states to which we've managed to get (starting from init).
  // i.e. we push onto <states_to_check> when we first encounter a state,
  // and pop that state back off when we've verified the state has a non-self transition

  AddToStateIndices(initial_, states_to_check);  // push init's transitions onto <states_to_check>

  vector<bool> checked_states(states_.size(), false);  // states that 1) are reachable from init and 2) have a non-self transition
  while(states_to_check.size() > 0) {
    uint16_t icheck(states_to_check.back());  // index of the state we're now checking
    states_to_check.pop_back();  // we're checking it now, so it no longer needs to be in <states_to_check>. NOTE states can in general appear in <states_to_check> more than once

    if(checked_states[icheck])  // skip it if we're already sure this state is ok
      continue;

    vector<uint16_t> tmp_visited;  // vector of the states to which we can transition from the <icheck>th state
    AddToStateIndices(states_[icheck], tmp_visited);
    size_t num_visited(tmp_visited.size());  // number of states to which we can transition from <icheck>. NOTE recall that transitions are not included in the vector of transitions (I don't know if there's a real reason for this)

    // make sure <icheck>th state has at least one non-self transition
    if(num_visited == 0) {  // if we didn't visit any, <icheck> better have a transition to 'end'
      if(states_[icheck]->trans_to_end() == nullptr)
        throw runtime_error("ERROR state '" + states_[icheck]->name() + "' in '" + name_ + "' has no transitions");
    } else if(num_visited == 1 && tmp_visited[0] == icheck) {  // if <icheck> only visited itself
      if(states_[icheck]->trans_to_end() == nullptr)
        throw runtime_error("ERROR state '"  + states_[icheck]->name() + "' in '" + name_ + "' has only a transition to itself");
    }

    checked_states[icheck] = true;

    // push onto <states_to_check> the states to which we can transition from <icheck>
    for(size_t i = 0; i < tmp_visited.size(); ++i) {
      if(!checked_states[tmp_visited[i]])  // skip 'em if we've already verified that they're reachable
        states_to_check.push_back(tmp_visited[i]);
    }
  }

  // make sure there's at least one transition to end
  bool found_end_trans(false);
  for(size_t i = 0; i < states_.size(); ++i) {
    if(states_[i]->trans_to_end()) {
      found_end_trans = true;
      break;
    }
  }
  if(!found_end_trans)
    throw runtime_error("ERROR no transtion to 'end' in '" + name_ + "'");

  // and finally, make sure we actually reached every state when we began from init and traversed through the entire transition network
  for(size_t i = 0; i < checked_states.size(); ++i) {
    if(!checked_states[i])
      throw runtime_error("ERROR state '" + states_[i]->name() + "' is not reachable from init.");
  }
}

// ----------------------------------------------------------------------------------------
void Model::AddToStateIndices(State *st, vector<uint16_t> &visited) {
  // push onto <visited> the index of each state to which we can transition from <st>
  for(size_t i = 0; i < st->transitions()->size(); ++i) {
    if(st->transitions()->at(i))
      visited.push_back(st->transitions()->at(i)->to_state()->index());
  }
}
}
