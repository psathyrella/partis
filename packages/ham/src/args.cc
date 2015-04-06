#include "args.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Args::Args(int argc, const char * argv[]):
  algo_strings_ {"viterbi", "forward"},
  debug_ints_ {0, 1, 2},
  algo_vals_(algo_strings_),
  debug_vals_(debug_ints_),
  hmmdir_arg_("m", "hmmdir", "directory in which to look for hmm model files", true, "", "string"),
  datadir_arg_("d", "datadir", "directory in which to look for non-sample-specific data (eg human germline seqs)", true, "", "string"),
  infile_arg_("i", "infile", "input (whitespace-separated) file", true, "", "string"),
  outfile_arg_("o", "outfile", "output csv file", true, "", "string"),
  cachefile_arg_("u", "cachefile", "input (and output) cache log prob csv file", false, "", "string"),
  algorithm_arg_("a", "algorithm", "algorithm to run", true, "", &algo_vals_),
  hamming_fraction_cutoff_arg_("f", "hamming-fraction-cutoff", "hamming fraction cutoff for clustering", false, -1.0, "float"),
  debug_arg_("g", "debug", "debug level", false, 0, &debug_vals_),
  n_best_events_arg_("n", "n_best_events", "number of candidate recombination events to write to file", true, -1, "int"),
  chunk_cache_arg_("c", "chunk-cache", "perform chunk caching?", false),
  partition_arg_("z", "partition", "", false),
  rescale_emissions_arg_("q", "rescale-emissions", "", false),
  str_headers_ {},
  int_headers_ {"k_v_min", "k_v_max", "k_d_min", "k_d_max"},
  str_list_headers_ {"names", "seqs", "only_genes"},  // passed as colon-separated lists of strings
  float_list_headers_ {"mute_freqs"}  // passed as colon-separated lists of floats
{
  try {
    CmdLine cmd("ham -- the fantastic HMM compiler", ' ', "");
    cmd.add(hmmdir_arg_);
    cmd.add(datadir_arg_);
    cmd.add(infile_arg_);
    cmd.add(outfile_arg_);
    cmd.add(cachefile_arg_);
    cmd.add(hamming_fraction_cutoff_arg_);
    cmd.add(algorithm_arg_);
    cmd.add(debug_arg_);
    cmd.add(n_best_events_arg_);
    cmd.add(chunk_cache_arg_);
    cmd.add(partition_arg_);
    cmd.add(rescale_emissions_arg_);

    cmd.parse(argc, argv);

    // algorithm_ = algorithm_arg_.getValue();
    // debug_ = debug_arg_.getValue();
  } catch(ArgException &e) {
    cerr << "ERROR: " << e.error() << " for argument " << e.argId() << endl;
    throw;
  }

  for(auto & head : str_headers_)
    strings_[head] = vector<string>();
  for(auto & head : int_headers_)
    integers_[head] = vector<int>();

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
    stringstream ss(line);
    string tmpstr;
    int tmpint;
    for(auto & head : headers) {
      if(str_headers_.find(head) != str_headers_.end()) {
        ss >> tmpstr;
        strings_[head].push_back(tmpstr);
      } else if(str_list_headers_.find(head) != str_list_headers_.end()) {
        ss >> tmpstr;
        str_lists_[head].push_back(SplitString(tmpstr, ":"));
      } else if(float_list_headers_.count(head)) {
        ss >> tmpstr;
	float_lists_[head].push_back(Floatify(SplitString(tmpstr, ":")));
      } else if(int_headers_.find(head) != int_headers_.end()) {
        ss >> tmpint;
        integers_[head].push_back(tmpint);
      } else {
        throw runtime_error("ERROR header " + head + "' not found");
      }
    }
  }
}

}
