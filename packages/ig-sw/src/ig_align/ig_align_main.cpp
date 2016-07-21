#include "tclap/CmdLine.h"
#include <ValuesConstraint.h>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include "ig_align.h"
}

// Gets file name from locus and directory.
std::vector<std::string> GetFileName(std::string vdj_dir, std::string locus) {
  std::vector<char> genes;
  std::vector<std::string> paths;
  if (locus == "IGH") {
    genes.push_back('v');
    genes.push_back('d');
    genes.push_back('j');
  }
  if (locus == "IGL" || locus == "IGK") {
    genes.push_back('v');
    genes.push_back('j');
  }
  // DJ rearrangements use the heavy chain genes,
  // but are incompletely rearranged
  if (locus == "DJ") {
    genes.push_back('d');
    genes.push_back('j');
    locus = "IGH";
  }
  // Converts to lowercase.
  std::transform(locus.begin(), locus.end(), locus.begin(), ::tolower);
  // Assembles file path.
  for (int i = 0; i < genes.size(); i++) {
    std::string path = vdj_dir + locus + genes[i] + ".fasta";
    paths.push_back(path);
  }
  return paths;
}

int main(int argc, const char *argv[]) {
  try {
    TCLAP::CmdLine cmd("Aligns reads from fasta files", ' ', "1");

    // tclap command line arguments
    TCLAP::UnlabeledValueArg<std::string> qry_path_opt(
        "query-path", "The path to the query file", true, "", "string");
    cmd.add(qry_path_opt);

    TCLAP::UnlabeledValueArg<std::string> output_path_opt(
        "output-path", "The path to the output file", true, "output.sam",
        "string");
    cmd.add(output_path_opt);

    TCLAP::ValueArg<int> match_opt("m", "match", "Match score: default 2",
                                   false, 2, "int");
    cmd.add(match_opt);

    TCLAP::ValueArg<int> mismatch_opt(
        "u", "mismatch", "Mismatch score: default 2", false, 2, "int");
    cmd.add(mismatch_opt);

    TCLAP::ValueArg<int> gap_o_opt(
        "o", "gap-open", "Gap opening penalty: default 3", false, 3, "int");
    cmd.add(gap_o_opt);

    TCLAP::ValueArg<int> gap_e_opt(
        "e", "gap-extend", "Gap extension penalty: default 1", false, 1, "int");
    cmd.add(gap_e_opt);

    TCLAP::ValueArg<unsigned> max_drop_opt("d", "max-drop",
                                           "Max drop: default 1000", false,
                                           1000, "int (unsigned)");
    cmd.add(max_drop_opt);

    TCLAP::ValueArg<int> min_score_opt("s", "min-score", "Min score: default 0",
                                       false, 0, "int");
    cmd.add(min_score_opt);

    TCLAP::ValueArg<unsigned> bandwidth_opt("b", "bandwidth",
                                            "Bandwidth: default 150", false,
                                            150, "int (unsigned)");
    cmd.add(bandwidth_opt);

    TCLAP::ValueArg<uint8_t> n_threads_opt(
        "j", "threads", "Number of threads: default 1", false, 1, "int");
    cmd.add(n_threads_opt);

    std::vector<std::string> options;
    options.push_back("IGH");
    options.push_back("IGL");
    options.push_back("IGK");
    options.push_back("DJ");
    TCLAP::ValuesConstraint<std::string> allowedVals(options);

    TCLAP::ValueArg<std::string> locus_opt(
        "l", "locus", "Locus to align reads against: default IGH", false, "IGH",
        &allowedVals);
    cmd.add(locus_opt);

    TCLAP::ValueArg<std::string> vdj_dir_opt(
        "p", "vdj-dir", "Directory from which to read germline genes", true, "",
        "string");
    cmd.add(vdj_dir_opt);

    cmd.parse(argc, argv);

    // qry_path
    std::string str_qry_path = qry_path_opt.getValue();
    // char *qry_path = ToCharArray(str_qry_path);
    char *qry_path = &str_qry_path[0];
    // output_path
    std::string str_out_path = output_path_opt.getValue();
    char *output_path = &str_out_path[0];
    // match
    int match = match_opt.getValue();
    // mismatch
    int mismatch = mismatch_opt.getValue();
    // gap o
    int gap_o = gap_o_opt.getValue();
    // gap e
    int gap_e = gap_e_opt.getValue();
    // max_drop
    unsigned max_drop = max_drop_opt.getValue();
    // min_score
    int min_score = min_score_opt.getValue();
    // bandwidth
    int bandwidth = bandwidth_opt.getValue();
    // n_threads
    uint8_t n_threads = n_threads_opt.getValue();
    // locus
    std::string locus = locus_opt.getValue();
    // vdj_dir
    std::string vdj_dir = vdj_dir_opt.getValue();

    // assemble paths
    std::vector<std::string> paths = GetFileName(vdj_dir, locus);
    std::string str_ref_path = paths[0];
    // convert to char*
    char *ref_path = &str_ref_path[0];

    int n_extra_refs = paths.size() - 1;
    char *extra_paths[n_extra_refs];

    for (int i = 1; i < paths.size(); i++) {
      extra_paths[i - 1] = &(paths[i])[0];
    }
    const char **extra_ref_paths = (const char **)extra_paths;

    ig_align_reads(ref_path, n_extra_refs, extra_ref_paths, qry_path,
                   output_path, match, mismatch, gap_o, gap_e, max_drop,
                   min_score, bandwidth, n_threads, NULL, NULL);

  } catch (TCLAP::ArgException &e) // catch any exception
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId()
              << std::endl;
  }
  return 0;
}
