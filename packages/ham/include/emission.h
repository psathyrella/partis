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
  void Parse(YAML::Node config, Track *track);
  void ReplaceLogProbs(vector<double> new_log_probs) { scores_.ReplaceLogProbs(new_log_probs); }
  void UnReplaceLogProbs() { scores_.UnReplaceLogProbs(); }
  ~Emission();

  double score(Sequence *seq, size_t pos) { return scores_.LogProb(seq, pos); }
  inline double score(uint8_t index) { return scores_.LogProb(index); }
  Track *track() { return track_; }
  vector<double> log_probs() { return scores_.log_probs(); }

  void Print();
private:
  double total_;
  LexicalTable scores_;
  Track *track_;
};

}
#endif
