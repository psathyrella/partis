#ifndef HAM_CLUSTERPATH_H
#define HAM_CLUSTERPATH_H

#include <set>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

using namespace std;
namespace ham {

typedef set<string> Partition;

// ----------------------------------------------------------------------------------------
class ClusterPath {  // sequence of gradually coalescing partitions, with associated info
public:
  ClusterPath() {}
  ClusterPath(Partition initial_partition, double initial_logprob);
  void AddPartition(Partition partition, double logprob);
  Partition &CurrentPartition() { return partitions_.back(); }  // return current (most recent) partition
  double CurrentLogProb() { return logprobs_.back(); }  // return logprob of current (most recent partition)
  vector<Partition> &partitions() { return partitions_; }
  vector<double> &logprobs() { return logprobs_; }
  bool finished_;
private:
  vector<Partition> partitions_;
  vector<double> logprobs_;

  double max_log_prob_of_partition_;
  Partition best_partition_;

};
}
#endif
