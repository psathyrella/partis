#include "clusterpath.h"
namespace ham {
// ----------------------------------------------------------------------------------------
ClusterPath::ClusterPath(Partition initial_partition, double initial_logprob):
  max_log_prob_of_partition_(-INFINITY),
  finished_(false)
{
  partitions_.push_back(initial_partition);
  logprobs_.push_back(initial_logprob);
}

// ----------------------------------------------------------------------------------------
void ClusterPath::AddPartition(Partition partition, double logprob) {
  if(args_->debug())
    PrintPartition(partition, "current");

  partitions_.push_back(partition);
  logprobs_.push_back(logprob);

  if(logprob > max_log_prob_of_partition_) {
    max_log_prob_of_partition_ = logprob;
    best_partition_ = partition;
  }

  if(max_log_prob_of_partition_ - logprob > 1000.0) {  // stop if we've moved too far past the maximum
    cout << "    stopping after drop " << max_log_prob_of_partition_ << " --> " << logprob << endl;
    finished_ = true;  // NOTE this will not play well with multiple maxima, but I'm pretty sure we shouldn't be getting those
  }
}
}
