#ifndef HAM_MODEL_H
#define HAM_MODEL_H

#include <fstream>
#include "state.h"
#include "yaml-cpp/yaml.h"

using namespace std;
namespace ham {

class Model {
public:
  Model();
  void Parse(string);
  void AddState(State*);
  void RescaleOverallMuteFreq(double overall_mute_freq);  // Rescale emissions to reflect <overall_mute_freq>, unless <overall_mute_freq> is -INFINITY, in which case we *re*-rescale them to what they were originally
  void UnRescaleOverallMuteFreq();  // Undo the above
  void Finalize();

  inline string &name() { return name_; }
  inline Track *track(size_t it) { return tracks_[it]; }
  inline size_t n_states() { return states_.size(); }
  inline State *state(string name) { assert(states_by_name_.count(name)); return states_by_name_[name]; }
  inline State *state(size_t ist) {
    if(ist >= states_.size())
      throw runtime_error("ERROR state index too large in model::state(): " + to_string(ist));
    return states_[ist];
  }
  inline bitset<STATE_MAX> *initial_to_states() { return initial_->to_states(); }  //!Get vector of states that the initial state transitions to
  inline State *init_state() { return initial_; }  //!Get pointer to the initial state
  inline double overall_prob() { return overall_prob_; }
  inline double original_overall_mute_freq() { return original_overall_mute_freq_; }

private:
  string name_;
  double overall_prob_;  // overall probability of this hmm/gene (not the same 'overall' as <overall_mute_freq_>)
  double original_overall_mute_freq_;  // mean mutation frequency, over v, d and j (not insertions), for the sequences in the data set
                                       // from which this hmm was derived. Reiterating: mean over all genes and all regions, *not* just this gene.
                                       // Note, this is the *original* one, i.e. we don't reset it when we reset the mute freqs
  double rescale_ratio_;  // ratio by which we have rescaled the emission probabilities (-INFINITY if we haven't rescaled them, i.e. if they correspond to <original_overall_mute_freq_>)
  string ambiguous_char_;
  Tracks tracks_;
  vector<State*> states_; //!  All the states contained in the model
  map<string, State*> states_by_name_; //Ptr to state stored by State name;
  State* initial_;
  State* ending_;
  bool finalized_;

  void FinalizeState(State *st);
  void CheckTopology();
  void AddToStateIndices(State* st, vector<uint16_t>& visited); // that's 'to-state', as in, 'here we push back the to-state indices onto <visited>'
};

}
#endif
