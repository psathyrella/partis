#include "clusterpath.h"
namespace ham {
// ----------------------------------------------------------------------------------------
ClusterPath::ClusterPath(Partition initial_partition, double initial_logprob):
  finished_(false),
  initial_path_index_(0),  // NOTE changing default from -1 to 0... which I think is ok, but it may screw something up
  max_log_prob_of_partition_(-INFINITY)
{
  partitions_.push_back(initial_partition);
  logprobs_.push_back(initial_logprob);
}

// ----------------------------------------------------------------------------------------
void ClusterPath::set_logprob(size_t il, double logprob) {
  logprobs_[il] = logprob;

  // NOTE see also AddPartition()
  if(max_log_prob_of_partition_ == -INFINITY or logprob > max_log_prob_of_partition_) {  // partitions usually get added with -INFINITY, in which case <max_log_prob_of_partition_> will stay -INFINITY, so <i_best_> will be the last one
    max_log_prob_of_partition_ = logprob;
    i_best_ = il;
    // best_partition_ = partition;
  }
  
}

// ----------------------------------------------------------------------------------------
void ClusterPath::AddPartition(Partition partition, double logprob, size_t n_max_partitions)  { // , double max_drop) {
  partitions_.push_back(partition);
  logprobs_.push_back(logprob);

  // int pot_parents = PotentialNumberOfParents(partition);
  // double combifactor(1.);
  // if(pot_parents > 0)  // if all partitions consist of one sequence, the above formula gives zero, but we want one
  //   combifactor = double(1.) / pot_parents;
  // logweights_.push_back(logweights_.back() + log(combifactor));


  // NOTE see also set_logprob()
  if(max_log_prob_of_partition_ == -INFINITY or logprob > max_log_prob_of_partition_) {  // partitions usually get added with -INFINITY, in which case <max_log_prob_of_partition_> will stay -INFINITY, so <i_best_> will be the last one
    max_log_prob_of_partition_ = logprob;
    i_best_ = logprobs_.size() - 1;
    // best_partition_ = partition;
  }

  // if(max_log_prob_of_partition_ - logprob > max_drop) {  // stop if we've moved too far past the maximum
  //   cout << "        stopping after drop " << max_log_prob_of_partition_ << " --> " << logprob << endl;
  //   finished_ = true;  // NOTE this will not play well with multiple maxima, but I'm pretty sure we shouldn't be getting those
  // }

  if(n_max_partitions > 0 && partitions_.size() > n_max_partitions) {  // NOTE we don't check here that we're not removing the best partition
    partitions_.pop_front();
    logprobs_.pop_front();
  }

}

// // ----------------------------------------------------------------------------------------
// int ClusterPath::PotentialNumberOfParents(Partition &partition, bool debug) {
//   int combifactor(0);
//   for(auto &cluster : partition) {
//     int n_k(SplitString(cluster, ":").size());
//     combifactor += pow(2, n_k - 1) - 1;
//   }

//   // if(debug) {
//   //   for(auto &key : partition)
//   //     cout << "          x " << key << endl;
//   //   cout << "    combi " << combifactor << endl;
//   // }

//   return combifactor;
// }
}
