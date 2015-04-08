#ifndef HAM_CLUSTERPATH_H
#define HAM_CLUSTERPATH_H
#include <set>
#include <vector>

using namespace std;
namespace ham {

typedef set<string> Partition;

// ----------------------------------------------------------------------------------------
class ClusterPath {  // sequence of gradually coalescing partitions, with associated info
public:
  ClusterPath(Partition initial_partition, double initial_logprob);
  void AddPartition(Partition partition, double logprob);
  Partition &CurrentPartition() { return partitions_[partitions_.size()-1]; }  // return current (most recent) partition
  bool finished_;
private:
  vector<Partition> partitions_;
  vector<double> logprobs_;

  double max_log_prob_of_partition_;
  Partition best_partition_;
};
}
#endif
