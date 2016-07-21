/*
Tests:
Issue 2: Tests if Catch is included properly with a trivial test case
Issue 3: Tests if klib parsing of fastq files is implemented correctly
Issues 5 and 6: Tests if results from running new ig_align are the same as when
running old files from ighutil
Issue 16: Tests if results from running ig_align with new flags are the same as
the results from ighutil
*/
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "ig_align/ig_align.h"
#include "ig_align/kseq.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <zlib.h>

// Helper functions

// ISSUE 3
KSEQ_INIT(gzFile, gzread)
// Reads fastq file in and parses it.
void ReadFile(char file_name[], std::vector<std::string> &output_seq) {
  gzFile fp;
  kseq_t *seq;
  int l;
  FILE *writefp;
  fp = gzopen(file_name, "r"); // open the file handler
  seq = kseq_init(fp);         // initialize seq
  int i = 0;
  int length = 0;
  while ((l = kseq_read(seq)) >= 0) { // read sequence
    output_seq.push_back(std::string("name: ") + seq->name.s + "\n");
    if (seq->comment.l) {
      output_seq.push_back(std::string("comment: ") + seq->comment.s + "\n");
    }
    output_seq.push_back(std::string("seq: ") + seq->seq.s + "\n");
    if (seq->qual.l) {
      output_seq.push_back(std::string("qual: ") + seq->qual.s + "\n");
    }
  }
  kseq_destroy(seq); // destroy seq
  gzclose(fp);       // close the file handler
}

// Strips labels from parsed vector.
void StripLabels(std::vector<std::string> &parsed_vector) {
  for (int i = 0; i < parsed_vector.size(); i++) {
    std::string temp = parsed_vector[i];
    char temp_char[1024];
    strcpy(temp_char, temp.c_str());
    std::string label = std::strtok(temp_char, " ");
    if ("name:" == label || "comment:" == label || "seq:" == label ||
        "qual:" == label) {
      parsed_vector[i] = strtok(NULL, " ");
    }
    parsed_vector[i].erase(
        std::remove(parsed_vector[i].begin(), parsed_vector[i].end(), '\n'),
        parsed_vector[i].end());
  }
}

// ISSUES 5 AND 6

// The return of OpenFile! Opens a file and writes it to a vector of strings
// Also removes version numbers
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

// Compares two files to see if they are identcal using TCLAP. Calls OpenFile.
void CompareFiles(std::string igsw_file, std::string ighutil_file) {
  std::vector<std::string> original_file;
  OpenFile(igsw_file, original_file);
  std::vector<std::string> new_file;
  OpenFile(ighutil_file, new_file);
  REQUIRE(original_file.size() == new_file.size());
  for (int i = 0; i < original_file.size(); i++) {
    REQUIRE(original_file[i] == new_file[i]);
  }
}

// Test cases

TEST_CASE("Issue 2: Trivial pass", "[trivial]") { REQUIRE(1 == 1); }

TEST_CASE("Issue 3: Comparing parsed and original files") {
  /*Test for checking if data from a fastq file is the same after being parsed.
  Place file "shorttest.fastq" in current working directory or change the path
  in
  the code. To test using a different fastq file, file path, length of file, and
  copied and pasted lines will need to be changed.*/
  char filename[] = "test_data/shorttest.fastq";
  std::vector<std::string> parsed_vector;
  ReadFile(filename, parsed_vector);
  StripLabels(parsed_vector);
  // Comparing expected values (from shorttest.fastq) to the parsed vector
  REQUIRE("GTACTCTGGTTTGTC|PRCONS=20080924-IGHM|CONSCOUNT=17|DUPCOUNT=9" ==
          parsed_vector[0]);
  REQUIRE("NNNNNNNNNNNNNNNNNNNNNNNTCCTACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCAC"
          "GTTCTCTGGGCTCTCACTCAGTACTAGTGAAGTGGCTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGG"
          "CCCGGGAGTGGCTTGCACTCCTTTATTGGAATGATGATAAGTACTATAGTCCATCTCTGAAGAGCAGG"
          "CTCACCATTACTAAGGACACCTCCGAAAATCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGA"
          "CACAGGCACATATTATTGTGCACACGACGTGACAAGAAGTCGGTATGGGATGGACGTCTGGGGCCAAG"
          "GGACCACGGTCACCGTCTCATCAGGGA" == parsed_vector[1]);
  REQUIRE("!!!!!!!!!!!!!!!!!!!!!!!HFHHFFHCD:AFFFGHHHHHHEEEHHHHHFHHHF@"
          "EECFFHFFHHHHFFHHFFGHHCFFGHHBFHDFHHHFHFHHFHHFFH;;5@=C=,=+4?DEDED@@"
          "DEEEE@AC*;B,;CE:'.2'8888CCAAA0.*A**?EECCC:?AEFEECA::?EA:*0ACE?:C?C?"
          "1*CFFEFCEFFDFFE=:?EEA*EFFFFFCC?A?B8'?<CFFFECFFCECAA?A::CAFFCCFEFEEE?"
          "A*EEAEEFFFFFFFE?AFEFEFFECEEFFEEE?ECAA>;"
          "EFEEBECAEFFEEEEEFDDDBFFEDEEEFFFFHHHHHHHHHHHHFFHHHHIHFC<"
          "HEHGDGHGEGFFFDH" == parsed_vector[2]);
  REQUIRE("GGTTTGCAATGGTTT|PRCONS=20080924-IGHG|CONSCOUNT=3|DUPCOUNT=2" ==
          parsed_vector[3]);
  REQUIRE("NNNNNNNNNNNNNNNNNNNNNNNGGGAGGCGTGGTCCAGCCTGGGACGTCCCTGAGACTCTCCTGTGC"
          "AGCGTCTGGATTCACCTTCAGTAGTTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGG"
          "AGTGGGTGGCACTCGCATGGTATGATGGAAGTACTGCATACTATACAAACTCCGTGAAGGGCCGATTC"
          "ACCATCTCCAGAGACAATTCCAAGAACACGCTGTCTTTGCAAATGAACAGCCTGACAGCCGAGGACAC"
          "GGCTGTGTATTACTGTGCGAGAGGCCACATCCCCTATGCCTACAATTACCTTTTTGACTACTGGGGCC"
          "AGGGAACCCTGGTCACCGTCTCCTCAGCCC" == parsed_vector[4]);
  REQUIRE("!!!!!!!!!!!!!!!!!!!!!!!"
          "oioooolooiihnmmoooooooomoaloomnhmiklnilmmmlnnnomlkkkomimomkmmmmogmkk"
          "kkkkkkkkiihiiihgiiiiiWfiZdiiicegigigNeii]^biiVV_eeiLi\\eai_"
          "eeKNeiiigmmkgikmkgiikggcJc^kkii^"
          "kcgejammmmmmmmkekmmmkki\\mmkmkmmmmmiimmmmkimmmmmkii^"
          "jkmkkmmmmmmmmmmmmmimmmmmmmkmmmmmmmmmkmkmmmmmmmmmmmmmmmmkmmmmmmmmmmkk"
          "mkmmmmmmmmooomoomokfoooooooooooooooqqponpqqpnloqqqqqqpoooonklpqlooon"
          "oo" == parsed_vector[5]);
  REQUIRE("GCGTAGCACCCGTGG|PRCONS=20080924-IGHG|CONSCOUNT=2|DUPCOUNT=1" ==
          parsed_vector[6]);
  REQUIRE("NNNNNNNNNNNNNNNNNNNNNNNNCCTACGCTGGTGAAACCCACACAGACCCTCACGCTGACCTGCAC"
          "CTTCTCTGGGTTCTCACTCAGCACTGGTGGAGTGGGTGTGGGCTGGATCCGTCAGCCCCCAGGAAAGG"
          "CCCTGGAGTGGCTTGCACTCATTTATTGGAATGATGATAAGCGCTACAGCCCATCTCTGAAGAGCAGG"
          "CTCACCATCACCAAGGACACCTCCAAAAACCAGGTGGTCCTTACAATGACCAACATGGACCCTGTGGA"
          "CACAGCCACATATTACTGTGCACACAGGGGATACGGTCACCTTGCAGCAACTGGTACTGAGTGGTTCG"
          "ACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGCCT" == parsed_vector[7]);
  REQUIRE("!!!!!!!!!!!!!!!!!!!!!!!!HHooloomoeckkkV`"
          "hlloooomVjooiooommhiikooomjnoooooooikmolmjl`"
          "jhmlhklkkmmfihkokkmollihhikii^_hiighii\\gic^f_f_Ycfgfcgiiiiiiieageg["
          "NXLCPJZag3IZXZINH]gcigZeZXM\\aaNkkkimki`iekegiigciecikgSN_"
          "amiicUeikmigegaikggmikmg\\kkkmmmmmkmmmkmkmmikcmki`"
          "ikmmkkkkkgkkmmmkgmkkmmmmmimmkkiimmmmkgmkikmmkiiPmkg]["
          "Y\\\\kkjkmkmmmmmmmmmmmlmmoooolooooooooopqpqpqqpqppoqqqpoookjddhhemml"
          "npolmi" == parsed_vector[8]);
}

TEST_CASE("Issue 5 and 6: Testing align_reads") {
  // Compares results from ig_align_main.cpp to those from ighutil.
  // Run scons in the src/ig_align/ directory before running this command. It
  // will create the ig_align executable.
  system("./\"ig_align/ig-sw\"  test_data/short_iglv.fasta "
         "test_data/testoutput.sam -p test_data/ -m 1 -u 1 -o 7 -e 1 -d 0");
  CompareFiles("test_data/testoutput.sam", "test_data/vdjalign_output.sam");
}

TEST_CASE("Issue 16: Testing with new flags") {
  // This requires scons as well.
  system("./\"ig_align/ig-sw\" test_data/16test_100.fastq "
         "test_data/16test_100output.sam -p test_data/ -m 1 -u 1 -o 7 -e 1 -d "
         "0");
  CompareFiles("test_data/16test_100output.sam",
               "test_data/vdjalign_16test_100output.sam");
}
