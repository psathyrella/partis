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
  ~Emission();

  inline double score(Sequence *seq, size_t pos) { return scores_.LogProb(seq, pos); }

  void Print();
private:
  double total_;
  LexicalTable scores_;
  vector<Track*>* tracks_; // tracks on to which this emission is wont to spit
};

}
#endif
