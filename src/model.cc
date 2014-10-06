#include "model.h"

namespace ham {
// ----------------------------------------------------------------------------------------
model::model() : overall_prob_(0.0), initial_(NULL), finalized_(false) {
  ending_ = new State;
}

// ----------------------------------------------------------------------------------------
void model::parse(string infname) {
  YAML::Node config = YAML::LoadFile(infname);
  name_ = config["name"].as<string>();
  overall_prob_ = config["extras"]["gene_prob"].as<double>();
  
  YAML::Node tracks(config["tracks"]);
  for (YAML::const_iterator it=tracks.begin(); it!=tracks.end(); ++it) {
    Track *trk = new Track;
    trk->setName(it->first.as<string>());
    for (size_t ic=0; ic<it->second.size(); ++ic)
      trk->addAlphabetChar(it->second[ic].as<string>());
    tracks_.push_back(trk);
  }

  vector<string> state_names;
  for (size_t is=0; is<config["states"].size(); ++is) {
    string name = config["states"][is]["name"].as<string>();
    for (auto chkname: state_names) {
      if (name == chkname) {
	cerr << "ERROR added two states with name \"" << name << "\"" << endl;
	assert(0);
      }
    }
    state_names.push_back(name);
  }

  for (size_t ist=0; ist<state_names.size(); ++ist) {
    State *st(new State);
    st->parse(config["states"][ist], state_names, tracks_);
    // st->print();

    if (st->name() == "init") {
      initial_ = st;
      states_by_name_[st->name()] = st;
    } else {
      assert(states_.size() < STATE_MAX);
      states_.push_back(st);
      states_by_name_[st->name()] = st;
    }
  }
      
  finalize(); // post process states and/to create an end state with only transitions-from
}
      
// ----------------------------------------------------------------------------------------
void model::add_state(State* state) {
  assert(states_.size() < STATE_MAX);
  states_.push_back(state);
  states_by_name_[state->name()] = state;
  return;
}
      
// ----------------------------------------------------------------------------------------
// set transitions, check labels, and perform other checks
void model::finalize() {
  assert(!finalized_);  // well it wouldn't *hurt* to call this twice, but you still *oughtn't* to

  // set each state's index within this model
  for (size_t i=0; i<states_.size(); ++i)
    states_[i]->set_index(i);
  // set each state's various transition pointers
  for(size_t i=0; i<states_.size(); ++i)
    finalize_state(states_[i]);
  finalize_state(initial_);
                      
  // now we need to fix state::transitions_ so that it agrees with state:index_
  for(size_t i=0; i<states_.size(); ++i)
    states_[i]->reorder_transitions(states_by_name_);
  initial_->reorder_transitions(states_by_name_);
          
  check_topology();
                      
  finalized_ = true;
}
      

// ----------------------------------------------------------------------------------------
void model::finalize_state(State *st) {
  // Set the bitsets in <st> that say which states we can go to from <st>, and from which we can arrive at <st>.
  // Also set to-state pointers in <st>'s transitions
  vector<Transition*>* transitions(st->transitions());
  for(size_t it=0; it<transitions->size(); ++it) {  // loops over the transitions out of <st>
    string to_state_name(transitions->at(it)->to_state_name());
    assert(states_by_name_.count(to_state_name));
    State* to_state(states_by_name_[to_state_name]);
    st->add_to_state(to_state); // add <to_state> to the list of states that you can go to from <st>
    transitions->at(it)->set_state(to_state);  // set <to_state> as the state corresponding to to_state_name_ in <transitions>
    if (st != initial_)
      to_state->add_from_state(st);  // add <st> to the list of states from which you can reach <to_state>
  }
      
  if (st->end_trans())
    ending_->add_from_state(st);
}
  
// ----------------------------------------------------------------------------------------
void model::check_topology(){
  //!Check model topology
  //!Iterates through all states to check to see if there are any:
  //! 1. Orphaned States
  //! 2. Dead end States
  //! 3. Uncompleted States
              
  vector<uint16_t> visited;
  add_to_state_indices(initial_, visited);  // add <initial_>'s transition to-states to <visited>
              
  vector<bool> states_visited(states_.size(), false);
  while (visited.size()>0) {
    uint16_t st_iter(visited.back());
    visited.pop_back();
    if (!states_visited[st_iter]) {
      vector<uint16_t> tmp_visited;
      add_to_state_indices(states_[st_iter], tmp_visited);
      size_t num_visited(tmp_visited.size());
                              
      // check orphaned
      if (num_visited == 0 ){
	// we get here if the state only has a transition to the end state. TODO reinstate this check
        // cerr << "Warning: State: "  << states_[st_iter]->name() << " has no transitions defined\n";
      } else if (num_visited == 1 && tmp_visited[0] == st_iter) {
        if (states_[st_iter]->end_trans() == NULL) {
          cerr << "State: "  << states_[st_iter]->name() << " is an orphaned state that has only transition to itself" << endl;
	} else {
	  cerr << "State: "  << states_[st_iter]->name() << " may be an orphaned state that only has transitions to itself and END state." << endl;
	}
      }
                              
      for(size_t i=0; i<tmp_visited.size(); ++i) {
        if (!states_visited[tmp_visited[i]])
          visited.push_back(tmp_visited[i]);
      }

      states_visited[st_iter] = true;
    }
  }
              
  // make sure that ending is defined
  bool ending_defined(false);
  for(size_t i=0; i<states_.size(); ++i) {
    if (states_[i]->end_trans()) {
      ending_defined = true;
      break;
    }
  }
  if (!ending_defined) {
    cerr << "No END state defined in the model\n";
    assert(0);
  }
              
  for(size_t i=0; i<states_visited.size(); ++i) {
    if (!states_visited[i]) {
      cerr << "ERROR state \""  << states_[i]->name() << "\" has bad topology. Check its transitions.\n";
      assert(0);
    }
  }
}
      
// ----------------------------------------------------------------------------------------
void model::add_to_state_indices(State *st, vector<uint16_t> &visited) {
  // for each transition out of <st>, add the index of its to-state to <visited>
  for(size_t i=0; i<st->transitions()->size(); ++i) {
    if (st->transitions()->at(i))
      visited.push_back(st->transitions()->at(i)->to_state()->index());
  }
}
}
