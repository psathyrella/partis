#ifndef STOCHHMM_EMM_H
#define STOCHHMM_EMM_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "yaml-cpp/yaml.h"
#include "lexicaltable.h"

using namespace std;
namespace stochhmm {
class emm {
  friend class state;
  friend class model;
public:
  emm();
  bool parse(YAML::Node config, string is_pair, tracks model_tracks);
  ~emm();
              
  // accessors
  inline bool isSimple() { return true; }
  inline bool isComplex() { return false; }
  bool pair() { return pair_; }
  size_t get_n_tracks() { return tracks_->size(); }
  inline double score(Sequence& seq, size_t pos) { return scores.getValue(seq, pos); }
  inline double score(Sequences& seqs, size_t pos) { return scores.getValue(seqs, pos); }
  
  inline LexicalTable* getTables(){ return &scores; }

  inline void print(){ cout << stringify() << endl; }
  string stringify();
private:
  bool pair_;
  LexicalTable scores;
  vector<track*>* tracks_;             //Tracks used
  vector<size_t>* track_indices; //Indices of tracks used
};
}
#endif
