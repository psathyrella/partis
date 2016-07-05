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
// Converts string to char array
char *ToCharArray(std::string str) {
  char *c_array = new char[str.length() + 1];
  std::strcpy(c_array, str.c_str());
  return c_array;
}

// Converts string with spaces as delimiters to array of char arrays
// void PathArray(char **c_array, std::string str) {
//   char *message = ToCharArray(str);
//   int i = 0;
//   c_array[i] = strtok(message, " ");
//   i++;
//   while ((c_array[i] = strtok(NULL, " ")) != NULL) {
//     i++;
//   }
// }

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
  if (locus == "DJ") {
    // files probably don't start with dj
    genes.push_back('d');
    genes.push_back('j');
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
    // TCLAP::ValueArg<std::string> ref_path_opt(
    //     "r", "ref-path", "The path to the reference file", true, "",
    //     "string");
    // cmd.add(ref_path_opt);
    //
    // TCLAP::ValueArg<int> n_extra_refs_opt("n", "n-extra-refs",
    //                                       "Number of extra reference files",
    //                                       false, 0, "int");
    // cmd.add(n_extra_refs_opt);
    //
    // TCLAP::ValueArg<std::string> extra_ref_paths_opt(
    //     "x", "extra-ref-paths", "The paths to the extra reference files",
    //     false,
    //     "", "string (space delimiter)");
    // cmd.add(extra_ref_paths_opt);

    TCLAP::ValueArg<std::string> qry_path_opt(
        "q", "query-path", "The path to the query file", true, "", "string");
    cmd.add(qry_path_opt);

    TCLAP::ValueArg<std::string> output_path_opt("O", "output-path",
                                                 "The path to the output file",
                                                 true, "output.sam", "string");
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
    // the files are probably not going to be named djj.fasta and djd.fasta, so
    // maybe something else would be better here
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

    // TCLAP::ValueArg<std::string> locus_opt(
    //     "l", "locus", "Locus to align reads against: default IGH", false,
    //     "IGH",
    //     "string");
    // cmd.add(locus_opt);

    cmd.parse(argc, argv);

    // // ref_path
    // std::string str_ref_path = ref_path_opt.getValue();
    // char *ref_path = ToCharArray(str_ref_path);
    // // n_extra_refs
    // int n_extra_refs = n_extra_refs_opt.getValue();
    // // extra_ref_paths
    // std::string paths = extra_ref_paths_opt.getValue();
    // char *extra_paths[n_extra_refs + 1];
    // PathArray(extra_paths, extra_ref_paths_opt.getValue());
    // const char **extra_ref_paths = (const char **)extra_paths;

    // qry_path
    std::string str_qry_path = qry_path_opt.getValue();
    char *qry_path = ToCharArray(str_qry_path);
    // output_path
    std::string str_out_path = output_path_opt.getValue();
    char *output_path = ToCharArray(str_out_path);
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

    std::vector<std::string> paths = GetFileName(vdj_dir, locus);
    std::string str_ref_path = paths[0];
    char *ref_path = ToCharArray(str_ref_path);

    int n_extra_refs = paths.size() - 1;
    char *extra_paths[n_extra_refs];

    for (int i = 1; i < paths.size(); i++) {
      extra_paths[i - 1] = ToCharArray(paths[i]);
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
