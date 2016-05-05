#ifndef HAM_CLUSTERPATH_H
#define HAM_CLUSTERPATH_H

#include <set>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "text.h"

using namespace std;
namespace ham {

typedef set<string> Partition;

// ----------------------------------------------------------------------------------------
class ClusterPath {  // sequence of gradually coalescing partitions, with associated info
public:
  ClusterPath() {}
  ClusterPath(Partition initial_partition, double initial_logprob=-INFINITY);
  void AddPartition(Partition partition, double logprob);  // , double max_drop);
  int PotentialNumberOfParents(Partition &partition, bool debug=false);  // number of partitions from which we could have arrived at this partition (i.e. number of ways to split it)
  Partition &CurrentPartition() { return partitions_.back(); }  // return current (most recent) partition
  double CurrentLogProb() { return logprobs_.back(); }  // return logprob of current (most recent partition)
  vector<Partition> &partitions() { return partitions_; }
  void set_logprob(size_t il, double logprob) { logprobs_[il] = logprob; }
  vector<double> &logprobs() { return logprobs_; }
  bool finished_;
  int initial_path_index_;  // index (in the batch of last glomeration steps) of the path which gave rise to this path [if you have to ask, you really don't want to know]
private:
  vector<Partition> partitions_;
  vector<double> logprobs_;

  double max_log_prob_of_partition_;
  Partition best_partition_;
};
}
#endif
