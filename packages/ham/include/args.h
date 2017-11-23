#ifndef HAM_ARGS_H
#define HAM_ARGS_H
#include <algorithm>
#include <map>
#include <set>
#include <fstream>
#include <cassert>
#include <ctime>
#include <cmath>

#include <text.h>
#include "tclap/CmdLine.h"
using namespace TCLAP;
using namespace std;
namespace ham {
// ----------------------------------------------------------------------------------------
// input processing class
// NOTE some input is passed on the command line (global configuration), while some is passed in a csv file (stuff that depends on each (pair of) sequence(s)).
class Args {
public:
  Args(int argc, const char * argv[]);
  // void Check();  // make sure everything's the same length (i.e. the input file had all the expected columns)

  string hmmdir() { return hmmdir_arg_.getValue(); }
  string datadir() { return datadir_arg_.getValue(); }
  string infile() { return infile_arg_.getValue(); }
  string outfile() { return outfile_arg_.getValue(); }
  string annotationfile() { return annotationfile_arg_.getValue(); }
  string input_cachefname() { return input_cachefname_arg_.getValue(); }
  string output_cachefname() { return output_cachefname_arg_.getValue(); }
  string locus() { return locus_arg_.getValue(); }
  float hamming_fraction_bound_lo() { return hamming_fraction_bound_lo_arg_.getValue(); }
  float hamming_fraction_bound_hi() { return hamming_fraction_bound_hi_arg_.getValue(); }
  float logprob_ratio_threshold() { return logprob_ratio_threshold_arg_.getValue(); }
  float max_logprob_drop() { return max_logprob_drop_arg_.getValue(); }
  string algorithm() { return algorithm_arg_.getValue(); }
  string ambig_base() { return ambig_base_arg_.getValue(); }
  string seed_unique_id() { return seed_unique_id_arg_.getValue(); }
  int debug() { return debug_arg_.getValue(); }
  int naive_hamming_cluster() { return naive_hamming_cluster_arg_.getValue(); }
  int biggest_naive_seq_cluster_to_calculate() { return biggest_naive_seq_cluster_to_calculate_arg_.getValue(); }
  int biggest_logprob_cluster_to_calculate() { return biggest_logprob_cluster_to_calculate_arg_.getValue(); }
  int n_partitions_to_write() { return n_partitions_to_write_arg_.getValue(); }
  unsigned n_final_clusters() { return n_final_clusters_arg_.getValue(); }
  unsigned min_largest_cluster_size() { return min_largest_cluster_size_arg_.getValue(); }
  unsigned max_cluster_size() { return max_cluster_size_arg_.getValue(); }
  unsigned random_seed() { return random_seed_arg_.getValue(); }
  bool no_chunk_cache() { return no_chunk_cache_arg_.getValue(); }
  bool partition() { return partition_arg_.getValue(); }
  bool dont_rescale_emissions() { return dont_rescale_emissions_arg_.getValue(); }
  bool cache_naive_seqs() { return cache_naive_seqs_arg_.getValue(); }
  bool cache_naive_hfracs() { return cache_naive_hfracs_arg_.getValue(); }
  bool only_cache_new_vals() { return only_cache_new_vals_arg_.getValue(); }
  bool write_logprob_for_each_partition() { return write_logprob_for_each_partition_arg_.getValue(); }
 
  // command line arguments
  vector<string> algo_strings_;
  vector<int> debug_ints_;
  ValuesConstraint<string> algo_vals_;
  ValuesConstraint<int> debug_vals_;
  ValueArg<string> hmmdir_arg_, datadir_arg_, infile_arg_, outfile_arg_, annotationfile_arg_, input_cachefname_arg_, output_cachefname_arg_, locus_arg_, algorithm_arg_, ambig_base_arg_, seed_unique_id_arg_;
  ValueArg<float> hamming_fraction_bound_lo_arg_, hamming_fraction_bound_hi_arg_, logprob_ratio_threshold_arg_, max_logprob_drop_arg_;
  ValueArg<int> debug_arg_, naive_hamming_cluster_arg_, biggest_naive_seq_cluster_to_calculate_arg_, biggest_logprob_cluster_to_calculate_arg_, n_partitions_to_write_arg_;
  ValueArg<unsigned> n_final_clusters_arg_, min_largest_cluster_size_arg_, max_cluster_size_arg_, random_seed_arg_;
  SwitchArg no_chunk_cache_arg_, partition_arg_, dont_rescale_emissions_arg_, cache_naive_seqs_arg_, cache_naive_hfracs_arg_, only_cache_new_vals_arg_, write_logprob_for_each_partition_arg_;

  // arguments read from csv input file
  map<string, vector<string> > strings_;
  map<string, vector<int> > integers_;
  map<string, vector<double> > floats_;
  map<string, vector<vector<string> > > str_lists_;
  map<string, vector<vector<int> > > int_lists_;
  map<string, vector<vector<double> > > float_lists_;
  set<string> str_headers_, int_headers_, float_headers_, str_list_headers_, int_list_headers_, float_list_headers_;
};
}
#endif
