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
  ~Model();
  void Parse(string);
  void AddState(State*);
  void RescaleOverallMuteFreq(double overall_mute_freq);  // Rescale emissions to reflect <overall_mute_freq>, unless <overall_mute_freq> is -INFINITY, in which case we *re*-rescale them to what they were originally
  void UnRescaleOverallMuteFreq();  // Undo the above
  void Finalize();
  void AddMaybeFasterFromStateStuff();

  string &name() { return name_; }
  Track *track() { return track_; }
  size_t n_states() { return states_.size(); }
  State *state(string name) { assert(states_by_name_.count(name)); return states_by_name_[name]; }
  State *state(size_t ist) { assert(ist < states_.size()); return states_[ist]; }
  bitset<STATE_MAX> *initial_to_states() { return initial_->to_states(); }  // get vector of states to which the initial state may transition
  State *init_state() { return initial_; }
  double overall_prob() { return overall_prob_; }
  double original_overall_mute_freq() { return original_overall_mute_freq_; }

private:
  void FinalizeState(State *st);
  void CheckTopology();
  void AddToStateIndices(State* st, vector<uint16_t>& visited); // that's 'to-state', as in, 'here we push back the to-state indices onto <visited>'

  string name_;
  double overall_prob_;  // overall probability of this hmm/gene (not the same 'overall' as <overall_mute_freq_>)
  double original_overall_mute_freq_;  // mean mutation frequency, over v, d and j (not insertions), for the sequences in the data set
                                       // from which this hmm was derived. Reiterating: mean over all genes and all regions, *not* just this gene.
                                       // Note, this is the *original* one, i.e. we don't reset it when we reset the mute freqs
  double rescale_ratio_;  // ratio by which we have rescaled the emission probabilities (-INFINITY if we haven't rescaled them, i.e. if they correspond to <original_overall_mute_freq_>)
  string ambiguous_char_;
  Track *track_;
  vector<State*> states_; //!  All the states contained in the model
  map<string, State*> states_by_name_; //Ptr to state stored by State name;
  State *initial_;
  State *ending_;
  bool finalized_;
};

}
#endif
