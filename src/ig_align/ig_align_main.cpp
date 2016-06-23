#include "tclap/CmdLine.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>

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
void PathArray(char **c_array, std::string str) {
  char *message = ToCharArray(str);
  int i = 0;
  c_array[i] = strtok(message, " ");
  i = i + 1;
  while ((c_array[i] = strtok(NULL, " ")) != NULL) {
    i = i + 1;
  }
}

int main(int argc, const char *argv[]) {
  try {
    TCLAP::CmdLine cmd("Aligns reads from fasta files", ' ', "1");

    // tclap command line arguments
    TCLAP::ValueArg<std::string> ref_path(
        "r", "ref_path", "The path to the reference file", true,
        "test_data/igkv.fasta", "string");
    cmd.add(ref_path);

    TCLAP::ValueArg<int> n_extra_refs("n", "n_extra_refs",
                                      "Number of extra reference files", false,
                                      0, "int");
    cmd.add(n_extra_refs);

    TCLAP::ValueArg<std::string> extra_ref_paths(
        "x", "extra_ref_paths", "The paths to the extra reference files", false,
        "", "string (space delimiter)");
    cmd.add(extra_ref_paths);

    TCLAP::ValueArg<std::string> qry_path("q", "qry_path",
                                          "The path to the query file", true,
                                          "test_data/iglv.fasta", "string");
    cmd.add(qry_path);

    TCLAP::ValueArg<std::string> output_path("o", "output_path",
                                             "The path to the output file",
                                             true, "output.txt", "string");
    cmd.add(output_path);

    TCLAP::ValueArg<int> match("m", "match", "Match: default 2", false, 2,
                               "int");
    cmd.add(match);

    TCLAP::ValueArg<int> mismatch("i", "mismatch", "Mismatch: default 2", false,
                                  2, "int");
    cmd.add(mismatch);

    TCLAP::ValueArg<int> gap_o("g", "gap_o", "Gap o: default 3", false, 3,
                               "int");
    cmd.add(gap_o);

    TCLAP::ValueArg<int> gap_e("e", "gap_e", "Gap e: default 1", false, 1,
                               "int");
    cmd.add(gap_e);

    TCLAP::ValueArg<unsigned> max_drop("d", "max_drop",
                                       "Max drop: default 1000", false, 1000,
                                       "int (unsigned)");
    cmd.add(max_drop);

    TCLAP::ValueArg<int> min_score("s", "min_score", "Min score: default 0",
                                   false, 0, "int");
    cmd.add(min_score);

    TCLAP::ValueArg<unsigned> bandwidth("b", "bandwidth",
                                        "Bandwidth: default 150", false, 150,
                                        "int (unsigned)");
    cmd.add(bandwidth);

    TCLAP::ValueArg<uint8_t> n_threads(
        "t", "n_threads", "Number of threads: default 1", false, 1, "int");
    cmd.add(n_threads);

    cmd.parse(argc, argv);

    // ref_path
    std::string str_ref_path = ref_path.getValue();
    char *r_path = ToCharArray(str_ref_path);
    // n_extra_refs
    int n_ref = n_extra_refs.getValue();
    // extra_ref_paths
    std::string paths = extra_ref_paths.getValue();
    char *extra_paths[n_ref + 1];
    PathArray(extra_paths, extra_ref_paths.getValue());
    const char **e_paths = (const char **)extra_paths;
    // qry_path
    std::string str_qry_path = qry_path.getValue();
    char *q_path = ToCharArray(str_qry_path);
    // output_path
    std::string str_out_path = output_path.getValue();
    char *o_path = ToCharArray(str_out_path);
    // match
    int mtch = match.getValue();
    // mismatch
    int mismtch = mismatch.getValue();
    // gap o
    int go = gap_o.getValue();
    // gap e
    int ge = gap_e.getValue();
    // max_drop
    unsigned m_drop = max_drop.getValue();
    // min_score
    int m_score = min_score.getValue();
    // bandwidth
    int bwidth = bandwidth.getValue();
    // n_threads
    uint8_t n_thrds = n_threads.getValue();

    ig_align_reads(r_path, n_ref, e_paths, q_path, o_path, mtch, mismtch, go,
                   ge, m_drop, m_score, bwidth, n_thrds, NULL, NULL);

  } catch (TCLAP::ArgException &e) // catch any exception
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId()
              << std::endl;
  }
  return 0;
}
