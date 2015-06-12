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
  cachefile_arg_("", "cachefile", "input (and output) cache log prob csv file", false, "", "string"),
  algorithm_arg_("", "algorithm", "algorithm to run", true, "", &algo_vals_),
  hamming_fraction_cutoff_arg_("", "hamming-fraction-cutoff", "hamming fraction cutoff for clustering", false, -1.0, "float"),
  debug_arg_("", "debug", "debug level", false, 0, &debug_vals_),
  n_best_events_arg_("", "n_best_events", "number of candidate recombination events to write to file", true, -1, "int"),
  smc_particles_arg_("", "smc-particles", "number of particles (paths) to run in sequential monte carlo (do not run smc if < 2)", false, 1, "int"),
  chunk_cache_arg_("", "chunk-cache", "perform chunk caching?", false),
  partition_arg_("", "partition", "", false),
  truncate_seqs_arg_("", "truncate-seqs", "truncate sequences to the same length on either side of the conserved cysteine. NOTE this disables caching, so do *not* turn it on unless you really have different-length sequences", false),
  rescale_emissions_arg_("", "rescale-emissions", "", false),
  unphysical_insertions_arg_("", "unphysical-insertions", "", false),
  str_headers_ {},
  int_headers_ {"path_index", "k_v_min", "k_v_max", "k_d_min", "k_d_max"},
  float_headers_ {"logweight"},
  str_list_headers_ {"names", "seqs", "only_genes"},  // passed as colon-separated lists of strings
  int_list_headers_ {"cyst_positions"},  // passed as colon-separated lists of ints
  float_list_headers_ {"mute_freqs"}  // passed as colon-separated lists of floats
{
  try {
    CmdLine cmd("bcrham -- the fantabulous HMM compiler goes to B-Cellville", ' ', "");
    cmd.add(hmmdir_arg_);
    cmd.add(datadir_arg_);
    cmd.add(infile_arg_);
    cmd.add(outfile_arg_);
    cmd.add(cachefile_arg_);
    cmd.add(hamming_fraction_cutoff_arg_);
    cmd.add(algorithm_arg_);
    cmd.add(debug_arg_);
    cmd.add(n_best_events_arg_);
    cmd.add(smc_particles_arg_);
    cmd.add(chunk_cache_arg_);
    cmd.add(partition_arg_);
    cmd.add(truncate_seqs_arg_);
    cmd.add(rescale_emissions_arg_);
    cmd.add(unphysical_insertions_arg_);

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

  ifstream ifs(infile());
  if(!ifs.is_open())
    throw runtime_error("ERROR infile (" + infile() + ") d.n.e.\n");
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
        throw runtime_error("ERROR header " + head + "' not found");
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
