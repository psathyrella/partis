#include "args.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Args::Args(int argc, const char * argv[]):
  algo_strings_ {"viterbi", "forward"},
  debug_ints_ {0, 1, 2},
  algo_vals_(algo_strings_),
  debug_vals_(debug_ints_),
  hmmdir_arg_("", "hmmdir", "directory in which to look for hmm model files", true, "", "string"),
  datadir_arg_("", "datadir", "directory in which to look for non-sample-specific data (eg human germline seqs)", true, "", "string"),
  infile_arg_("", "infile", "input (whitespace-separated) file", true, "", "string"),
  outfile_arg_("", "outfile", "output csv file", true, "", "string"),
  annotationfile_arg_("", "annotationfile", "if specified, write annotations for each cluster to here", false, "", "string"),
  input_cachefname_arg_("", "input-cachefname", "input cached log prob/naive seq csv file", false, "", "string"),
  output_cachefname_arg_("", "output-cachefname", "output cached log prob/naive seq csv file", false, "", "string"),
  locus_arg_("", "locus", "ig{h,k,l} or tr{a,b,g,d}", true, "", "string"),
  algorithm_arg_("", "algorithm", "algorithm to run", true, "", &algo_vals_),
  ambig_base_arg_("", "ambig-base", "ambiguous base", false, "", "string"),
  seed_unique_id_arg_("", "seed-unique-id", "seed unique id", false, "", "string"),
  hamming_fraction_bound_lo_arg_("", "hamming-fraction-bound-lo", "if hamming fraction for a pair is smaller than this, merge them without calculating lratio", false, 0.0, "float"),
  hamming_fraction_bound_hi_arg_("", "hamming-fraction-bound-hi", "if hamming fraction for a pair is larger than this, skip without calculating lratio", false, 1.0, "float"),
  logprob_ratio_threshold_arg_("", "logprob-ratio-threshold", "", false, -INFINITY, "float"),
  max_logprob_drop_arg_("", "max-logprob-drop", "stop glomerating when the total logprob has dropped by this much", false, -1.0, "float"),
  debug_arg_("", "debug", "debug level", false, 0, &debug_vals_),
  naive_hamming_cluster_arg_("", "naive-hamming-cluster", "cluster sequences using naive hamming distance", false, 0, "int"),
  biggest_naive_seq_cluster_to_calculate_arg_("", "biggest-naive-seq-cluster-to-calculate", "", false, 99999, "int"),
  biggest_logprob_cluster_to_calculate_arg_("", "biggest-logprob-cluster-to-calculate", "", false, 99999, "int"),
  n_partitions_to_write_arg_("", "n-partitions-to-write", "how many partitions, before the best one, should we write to the output file", false, 99999, "int"),
  n_final_clusters_arg_("", "n-final-clusters", "instead of stopping at the most likely partition, stop when you have this many clusters", false, 0, "unsigned"),
  min_largest_cluster_size_arg_("", "min-largest-cluster-size", "instead of stopping at the most likely partition, stop when your largest cluster is this big", false, 0, "unsigned"),
  max_cluster_size_arg_("", "max-cluster-size", "if any cluster gets bigger than this, stop clustering", false, 0, "unsigned"),
  random_seed_arg_("", "random-seed", "", false, time(NULL), "unsigned"),
  no_chunk_cache_arg_("", "no-chunk-cache", "don't perform chunk caching?", false),
  partition_arg_("", "partition", "", false),
  dont_rescale_emissions_arg_("", "dont-rescale-emissions", "", false),
  cache_naive_seqs_arg_("", "cache-naive-seqs", "cache all naive sequences", false),
  cache_naive_hfracs_arg_("", "cache-naive-hfracs", "cache naive hamming fraction between sequence sets (in addition to log probs and naive seqs)", false),
  only_cache_new_vals_arg_("", "only-cache-new-vals", "only write sequence sets with newly-calculated values to cache file", false),
  write_logprob_for_each_partition_arg_("", "write-logprob-for-each-partition", "By default, we don't know the total logprob of each partition (since many merges are by naive hfrac). This argument tells us that this is the last time through (with one process) and we want to know the total probability of each partition.", false),
  str_headers_ {},
  int_headers_ {"k_v_min", "k_v_max", "k_d_min", "k_d_max", "cdr3_length"},
  float_headers_ {"mut_freq"},
  str_list_headers_ {"names", "seqs", "only_genes"},  // passed as colon-separated lists of strings
  int_list_headers_ {},  // passed as colon-separated lists of ints
  float_list_headers_ {}  // passed as colon-separated lists of floats
{
  try {
    CmdLine cmd("bcrham -- the fantabulous HMM compiler goes to B-Cellville", ' ', "");
    cmd.add(hmmdir_arg_);
    cmd.add(datadir_arg_);
    cmd.add(infile_arg_);
    cmd.add(outfile_arg_);
    cmd.add(annotationfile_arg_);
    cmd.add(input_cachefname_arg_);
    cmd.add(output_cachefname_arg_);
    cmd.add(locus_arg_);
    cmd.add(hamming_fraction_bound_lo_arg_);
    cmd.add(hamming_fraction_bound_hi_arg_);
    cmd.add(logprob_ratio_threshold_arg_);
    cmd.add(max_logprob_drop_arg_);
    cmd.add(algorithm_arg_);
    cmd.add(ambig_base_arg_);
    cmd.add(seed_unique_id_arg_);
    cmd.add(debug_arg_);
    cmd.add(naive_hamming_cluster_arg_);
    cmd.add(biggest_naive_seq_cluster_to_calculate_arg_);
    cmd.add(biggest_logprob_cluster_to_calculate_arg_);
    cmd.add(n_partitions_to_write_arg_);
    cmd.add(n_final_clusters_arg_);
    cmd.add(min_largest_cluster_size_arg_);
    cmd.add(max_cluster_size_arg_);
    cmd.add(random_seed_arg_);
    cmd.add(no_chunk_cache_arg_);
    cmd.add(cache_naive_seqs_arg_);
    cmd.add(cache_naive_hfracs_arg_);
    cmd.add(only_cache_new_vals_arg_);
    cmd.add(write_logprob_for_each_partition_arg_);
    cmd.add(partition_arg_);
    cmd.add(dont_rescale_emissions_arg_);

    cmd.parse(argc, argv);

  } catch(ArgException &e) {
    cerr << "ERROR: " << e.error() << " for argument " << e.argId() << endl;
    throw;
  }

  for(auto & head : str_headers_)
    strings_[head] = vector<string>();
  for(auto & head : int_headers_)
    integers_[head] = vector<int>();
  for(auto & head : float_headers_)
    floats_[head] = vector<double>();

  vector<string> loci{"igh", "igk", "igl", "tra", "trb", "trg", "trd"};  // this is ugly... but oh, well
  if(find(loci.begin(), loci.end(), locus()) == loci.end())
    throw runtime_error("--locus argument '" + locus() + "' not among ig{h,k,l} or tr{a,b,g,d}");

  ifstream ifs(infile());
  if(!ifs.is_open())
    throw runtime_error("args.cc: bcrham input file '" + infile() + "' d.n.e.\n");
  string line;
  // get header line
  getline(ifs, line);
  stringstream ss(line);
  vector<string> headers;  // keep track of the file's column order
  while(!ss.eof()) {
    string head;
    ss >> head;
    if(head != "")   // make sure not to add a blank header at the end of the file
      headers.push_back(head);
  }
  while(getline(ifs, line)) {
    if(line.size() < 10)  // 10 is kinda arbitrary, but we just want to skip blank lines
      continue;
    stringstream ss(line);
    string tmpstr;
    int tmpint;
    double tmpfloat;
    for(auto & head : headers) {
      if(str_headers_.count(head)) {
        ss >> tmpstr;
        strings_[head].push_back(tmpstr);
      } else if(str_list_headers_.count(head)) {
        ss >> tmpstr;
        str_lists_[head].push_back(SplitString(tmpstr, ":"));
      } else if(int_list_headers_.count(head)) {
        ss >> tmpstr;
        int_lists_[head].push_back(Intify(SplitString(tmpstr, ":")));
      } else if(float_list_headers_.count(head)) {
        ss >> tmpstr;
	float_lists_[head].push_back(Floatify(SplitString(tmpstr, ":")));
      } else if(int_headers_.count(head)) {
        ss >> tmpint;
        integers_[head].push_back(tmpint);
      } else if(float_headers_.count(head)) {
        ss >> tmpfloat;
        floats_[head].push_back(tmpfloat);
      } else {
        throw runtime_error("args.cc: found unexpected header " + head + "' in input file " + infile());
      }
    }
  }

  // Check();
}

// // ----------------------------------------------------------------------------------------
//   // oh, wait, this isn't right
//   // need to rearrange a couple things
// void Args::Check() {
//   assert(strings_.size() == integers_.size());
//   assert(strings_.size() == floats_.size());
//   assert(strings_.size() == str_lists_.size());
//   assert(strings_.size() == int_lists_.size());
//   assert(strings_.size() == float_lists_.size());
//   for(auto &kv : str_lists_) {  // loop over column headers
//     string column(kv.first);
//     vector<vector<string> > &sub_str_lists(kv.second);
//     for(size_t il=0; il<sub_str_lists.size(); ++il) {  // loop over queries
//       vector<string> &str_list(sub_str_lists[il]);
//       vector<int> &int_list(int_lists_[column][il]);
//       vector<double> &float_list(float_lists_[column][il]);
//       // make sure each list for this query has the same length
//       assert(str_list.size() == int_list.size());
//       assert(str_list.size() == float_list.size());
//     }
//   }
// }
}
