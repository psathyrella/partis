#ifndef HAM_EMISSION_H
#define HAM_EMISSION_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

#include "yaml-cpp/yaml.h"
#include "lexicaltable.h"

using namespace std;

namespace ham {

class Emission {
public:
  Emission();
  void Parse(YAML::Node config, Tracks model_tracks);
  void ReplaceLogProbs(vector<vector<double> > new_log_probs) { scores_.ReplaceLogProbs(new_log_probs); }
  void UnReplaceLogProbs() { scores_.UnReplaceLogProbs(); }
  size_t lps() { return scores_.log_probs().size(); }  // TODO remove
  ~Emission();

  inline double score(Sequence *seq, size_t pos) { return scores_.LogProb(seq, pos); }
  // inline double score(Sequence *seq, size_t pos) {     cout << "0---" << endl; cout << scores_.LogProb(seq, pos) << endl; cout << "1---" << endl; return scores_.LogProb(seq, pos); }
  // double score(size_t letter) { return scores_.LogProb(letter); }
  Track *track() { assert(tracks_->size() == 1); return tracks_->at(0); }  // I should really eliminate the possiblity of more than one track, since it's not supported any more
  vector<vector<double> > log_probs() { return scores_.log_probs(); }

  void Print();
private:
  double total_;
  LexicalTable scores_;
  vector<Track*>* tracks_; // tracks on to which this emission is wont to spit
};

}
#endif
