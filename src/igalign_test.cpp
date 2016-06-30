// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "ig_align/ig_align.h"

// The return of open file: opens file and writes it to string vector
void OpenFile(std::string filename, std::vector<std::string> &original_file) {
  std::ifstream test_file;
  test_file.open(filename.c_str());
  std::string read_file;
  if (test_file.is_open()) {
    while (test_file.good()) {
      test_file >> read_file;
      if (read_file.find("VN:") != 0) // Ignores version numbers
        original_file.push_back(read_file);
      if (test_file.eof())
        break;
      read_file = "";
    }
  }
  test_file.close();
}

// Compares results from ig_align_main.cpp to those from ighutil
TEST_CASE("Testing align_reads") {
  system("./\"ig_align/ig_align\" -q test_data/short_iglv.fasta -o "
         "test_data/testoutput.sam -r "
         "test_data/ighv.fasta -n 2 -x "
         "\"test_data/ighd.fasta test_data/ighj.fasta\" -m 1 "
         "-i 1 -g 7 -e 1 -d 0"); // Paths and parameters
  std::vector<std::string> original_file;
  OpenFile("test_data/testoutput.sam", original_file);
  std::vector<std::string> ighutil_file;
  OpenFile("test_data/vdjalign_output.sam", ighutil_file);
  REQUIRE(original_file.size() == ighutil_file.size());
  for (int i = 0; i < original_file.size(); i++) {
    REQUIRE(original_file[i] == ighutil_file[i]);
  }
}
