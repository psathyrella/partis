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
  ClusterPath(Partition initial_partition, double initial_logprob, double initial_logweight);
  void AddPartition(Partition partition, double logprob);
  int PotentialNumberOfParents(Partition &partition, bool debug=false);  // number of partitions from which we could have arrived at this partition (i.e. number of ways to split it)
  Partition &CurrentPartition() { return partitions_.back(); }  // return current (most recent) partition
  double CurrentLogProb() { return logprobs_.back(); }  // return logprob of current (most recent partition)
  vector<Partition> &partitions() { return partitions_; }
  vector<double> &logprobs() { return logprobs_; }
  vector<double> &logweights() { return logweights_; }
  double CurrentLogWeight() { return logweights_.back(); }
  bool finished_;
  int tmpi;
private:
  vector<Partition> partitions_;
  vector<double> logprobs_;
  vector<double> logweights_;  // (log of the) product of the inverse of the number of potential parents up to each step

  double max_log_prob_of_partition_;
  Partition best_partition_;
};
}
#endif
