/*Test for checking if data from a fastq file is the same after being parsed.
Place file "shorttest.fastq" in current working directory or change the path in
the code. To test using a different fastq file, file path, length of file, and
copied and pasted lines will need to be changed.*/

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "kseq.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)
// Populates array with char arrays.
extern "C" void WriteToArray(char *temp, char **output_seq, int i) {
  int length = strlen(temp);
  output_seq[i] = (char *)malloc(sizeof(char) * (length + 1));
  strcpy(output_seq[i], temp);
}

// Reads fastq file in and parses it. Calls WriteToArray.
extern "C" void ReadFile(char file_name[], char **output_seq) {
  gzFile fp;
  kseq_t *seq;
  int l;
  FILE *writefp;
  fp = gzopen(file_name, "r"); // open the file handler
  seq = kseq_init(fp);         // initialize seq
  int i = 0;
  int length = 0;
  while ((l = kseq_read(seq)) >= 0) { // read sequence
    char temp[500]; // longer than the longest line in the fastq file
    sprintf(temp, "name: @%s\n", seq->name.s);
    WriteToArray(temp, output_seq, i);
    i++;
    if (seq->comment.l) {
      sprintf(temp, "comment: %s\n", seq->comment.s);
      WriteToArray(temp, output_seq, i);
      i++;
    }
    sprintf(temp, "seq: %s\n", seq->seq.s);
    WriteToArray(temp, output_seq, i);
    i++;
    if (seq->qual.l) {
      sprintf(temp, "qual: %s\n", seq->qual.s);
      WriteToArray(temp, output_seq, i);
      i++;
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

TEST_CASE("Issue 3: Comparing parsed and original files") {
  char filename[] = "test_data/shorttest.fastq";
  int output_length = 9; // Number of lines after parsing
  char *output_seq[output_length];
  ReadFile(filename, output_seq);
  std::vector<std::string> parsed_vector;
  for (int i = 0; i < output_length; i++) {
    const char *temp = output_seq[i];
    std::string str(temp);
    parsed_vector.push_back(str);
  }
  StripLabels(parsed_vector);
  // Comparing expected values (from shorttest.fastq) to the parsed vector
  REQUIRE("@GTACTCTGGTTTGTC|PRCONS=20080924-IGHM|CONSCOUNT=17|DUPCOUNT=9" ==
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
  REQUIRE("@GGTTTGCAATGGTTT|PRCONS=20080924-IGHG|CONSCOUNT=3|DUPCOUNT=2" ==
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
  REQUIRE("@GCGTAGCACCCGTGG|PRCONS=20080924-IGHG|CONSCOUNT=2|DUPCOUNT=1" ==
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

TEST_CASE("Trivial pass", "[trivial]") { REQUIRE(1 == 1); }
