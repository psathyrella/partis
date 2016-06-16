#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "kseq.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <zlib.h>

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

extern "C" void ReadFile() {
  gzFile fp;
  kseq_t *seq;
  int l;
  FILE *writefp;
  writefp = fopen("readtest.txt", "w");
  fp = gzopen("test.fastq", "r");     // STEP 2: open the file handler
  seq = kseq_init(fp);                // STEP 3: initialize seq
  while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
    fprintf(writefp, "name: @%s\n", seq->name.s);
    if (seq->comment.l)
      fprintf(writefp, "comment: %s\n", seq->comment.s);
    fprintf(writefp, "seq: %s\n", seq->seq.s);
    if (seq->qual.l)
      fprintf(writefp, "qual: %s\n", seq->qual.s);
  }
  // printf("return value: %d\n", l);
  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp);       // STEP 6: close the file handler
}

std::vector<std::string> OpenFile(std::string filename, bool label) {
  std::ifstream test_file;
  // std::ifstream readtest;
  // std::ofstream write;
  std::vector<std::string> original_file;
  // write.open ("opentest.txt");
  // readtest.open ("readtest.txt");
  test_file.open(filename.c_str());
  std::string read_file;
  if (test_file.is_open()) {
    while (test_file.good()) {
      test_file >> read_file;
      if (label) {
        // strip labels
        if ("name:" != read_file && "comment:" != read_file &&
            "seq:" != read_file && "qual:" != read_file)
          original_file.push_back(read_file);

        // std::string n = "name:";
        // std::string::size_type i = read_file.find(n);
        // if (i != std::string::npos)
        //   read_file.erase(i, n.length());
        //
        // original_file.push_back(read_file);

      }

      else {
        if ("+" != read_file)
          original_file.push_back(read_file);
      }

      // write << std::endl;
      if (test_file.eof())
        break;
      read_file = "";
    }
    // Print vector
    for (int i = 0; i < original_file.size(); i++) {
      // std::cout << original_file[i] + "\n";
    }
  }
  return original_file;
}

// void Compare(std::vector<std::string> original_vector,
// std::vector<std::string> parsed_vector)
// {
//   if (original_vector.size() != parsed_vector.size())
//     throw "Files are not the same length.";
//   for (int i = 0; i < original_vector.size(); i++)
//   {
//
//   }
// }

TEST_CASE("Comparison") {
  // ReadFile();
  // sleep(5);
  std::vector<std::string> parsed_vector = OpenFile("readtest.txt", true);
  std::vector<std::string> original_vector = OpenFile("test.fastq", false);
  REQUIRE(original_vector.size() == parsed_vector.size());
  for (int i = 0; i < original_vector.size(); i++) {
    REQUIRE(original_vector[i] == parsed_vector[i]);
  }
}

// int main(int argc, char const *argv[])
// {
//   //ReadFile();
//   std::vector<std::string> parsed_vector = OpenFile("readtest.txt", true);
//   std::cout << "\n";
//   std::vector<std::string> original_vector = OpenFile("test.fastq", false);
//
//   // std::cout << "max size o " << original_vector.max_size() << "\n";
//   // std::cout << "capacity o " << original_vector.capacity() << "\n";
//   // std::cout << "max size p " << parsed_vector.max_size() << "\n";
//   // std::cout << "capacity p " << parsed_vector.capacity() << "\n";
//   return 0;
// }
