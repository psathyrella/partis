#include "bcrutils.h"

namespace ham {

// ----------------------------------------------------------------------------------------
TermColors::TermColors() {
  codes_["bold"] = "\033[1m";
  codes_["reverse"] = "\033[7m";
  codes_["purple"] = "\033[95m";
  codes_["blue"] = "\033[94m";
  codes_["green"] = "\033[92m";
  codes_["yellow"] = "\033[93m";
  codes_["red"] = "\033[91m";
  codes_["end"] = "\033[0m";
}

// ----------------------------------------------------------------------------------------
string TermColors::Color(string col, string seq) {
  assert(codes_.find(col) != codes_.end());
  return codes_[col] + seq + codes_["end"];
}

// ----------------------------------------------------------------------------------------
string TermColors::RedifyIfMuted(char germline_nuc, char nuc) {
  if(nuc == germline_nuc)
    return string(1, nuc);
  else
    return Color("red", string(1, nuc));
}

// ----------------------------------------------------------------------------------------
// get region from gene name, e.g. IGHD1-1*01 --> d
string TermColors::GetRegion(string gene) {
  assert(gene.find("IGH") != string::npos);
  char region_char(tolower(gene[3]));
  string region(1, region_char);
  assert(region == "v" || region == "d" || region == "j");
  return region;
}

// ----------------------------------------------------------------------------------------
string TermColors::ColorMutants(string color, string seq, string ref_1, vector<string> other_refs, string ambiguous_char) {
  // Return <seq> with mutant bases w.r.t. <ref_1> escaped to appear red (and 'i', inserts, yellow) in bash terminal.
  // If <refs> are specified, use bold text and reverse video to show if <seq> is muted w/respect to more than one ref
  if(ref_1 != "")
    other_refs.push_back(ref_1);  // only doing it this way so we can call it without specifying <other_refs>
  for(auto & ref : other_refs) {
    if(ref.size() != seq.size()) {
      throw runtime_error("ERROR seqs not same length in color_mutants: " + ref + "\n                                              " + seq);
    }
  }

  string return_str;
  for(size_t inuc = 0; inuc < seq.size(); ++inuc) {
    if(seq[inuc] == 'i') {
      return_str += Color("yellow", seq.substr(inuc, 1));
    } else {
      int ndiff(0);  // number of reference sequences that differ at base inuc
      for(auto & ref : other_refs) {
	if(ambiguous_char != "" and ref[inuc] == ambiguous_char[0])  // don't count ambiguous characters as mutated
	  continue;
        if(seq[inuc] != ref[inuc])
          ndiff += 1;
      }
      if(ndiff == 0)
        return_str += seq[inuc];
      else if(ndiff == 1)
        return_str += Color(color, seq.substr(inuc, 1));
      else
        return_str += Color("reverse", Color(color, seq.substr(inuc, 1)));
    }
  }
  return return_str;
}

// ----------------------------------------------------------------------------------------
string TermColors::ColorGene(string gene) {
  string return_str = Color("red", GetRegion(gene));  //gene.substr(0, 3) + Color("purple", gene.substr(3, 1));
  string n_version = gene.substr(4, gene.find("-") - 4);
  string n_subversion = gene.substr(gene.find("-") + 1, gene.find("*") - gene.find("-") - 1);
  if(GetRegion(gene) == "j") {
    n_version = gene.substr(4, gene.find("*") - 4);
    n_subversion = "";
    return_str += "" + Color("purple", n_version);
  } else {
    return_str += "" + Color("purple", n_version) + "-" + Color("purple", n_subversion);
  }
  string allele = gene.substr(gene.find("*") + 1, gene.find("_") - gene.find("*") - 1);
  return_str += "" + Color("yellow", allele);
  if(gene.find("_") != string::npos)   // _F or _P in j gene names
    return_str += gene.substr(gene.find("_"));
  return return_str;
}

// ========================================================================================
GermLines::GermLines(string input_dir):
  regions_({"v", "d", "j"})
{
  for(auto & region : regions_) {
    names_[region] = vector<string>();
    ifstream ifs(input_dir + "/igh" + region + ".fasta");
    assert(ifs.is_open());
    string line, name, seq;
    while(getline(ifs, line)) {
      if(line[0] == '>') {   // read header lines
        assert(line[1] != ' ');   // make my life hard, will you?
        name = line.substr(1, line.find(" ") - 1);  // skip the '>', and run until the first blank. It *should* be the gene name. We'll find out later when we look for the file.
        names_[region].push_back(name);
      } else {
        // line.replace(line.find("\n"), 1, "");
        seq = (seq == "") ? line : seq + line;
        if(ifs.peek() == '>' || ifs.peek() == EOF) { // if we're done with this sequence, add it to the map
          seqs_[name] = seq;
          seq = "";
        }
      }
    }
    ifs.close();
  }

  // get cyst and tryp info
  ifstream ifs;
  string line;
  vector<string> header;

  // get cyst info
  ifs.open(input_dir + "/v-meta.csv");
  assert(ifs.is_open());
  // check header
  getline(ifs, line);
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  header = (SplitString(line, ","));
  assert(header[0] == "gene");
  assert(header[1] == "cyst_start");
  // get info
  while(getline(ifs, line)) {
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    vector<string> info(SplitString(line, ","));
    assert(info[0].find("IGH") == 0);
    cyst_positions_[info[0]] = atoi(info[1].c_str());
  }
  ifs.close();

  // get tryp info
  ifs.open(input_dir + "/j_tryp.csv");
  assert(ifs.is_open());
  // check header
  getline(ifs, line);
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  header = SplitString(line, ",");
  assert(header[0] == "gene");
  assert(header[1] == "tryp_start");
  // get info
  while(getline(ifs, line)) {
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    vector<string> info(SplitString(line, ","));
    assert(info[0].find("IGH") == 0);
    tryp_positions_[info[0]] = atoi(info[1].c_str());
  }
  
}

// ----------------------------------------------------------------------------------------
// replace * with _star_ and /OR15 with _slash_
string GermLines::SanitizeName(string gene_name) {
  size_t istar(gene_name.find("*"));
  if(istar != string::npos)
    gene_name.replace(istar, 1, "_star_");

  size_t islash(gene_name.find("/"));
  if(islash != string::npos)
    gene_name.replace(islash, 1, "_slash_");

  return gene_name;
}

// ----------------------------------------------------------------------------------------
// get region from gene name, e.g. IGHD1-1*01 --> d
string GermLines::GetRegion(string gene) {
  assert(gene.find("IGH") != string::npos);
  char region_char(tolower(gene[3]));
  string region(1, region_char);
  assert(region == "v" || region == "d" || region == "j");
  return region;
}

// ========================================================================================
// ----------------------------------------------------------------------------------------
RecoEvent::RecoEvent() : score_(999)
{
}

// ----------------------------------------------------------------------------------------
void RecoEvent::SetNaiveSeq(GermLines &gl) {
  map<string, string> original_seqs, eroded_seqs;
  map<string, int> lengths;
  for(auto & region : gl.regions_) {
    int del_5p = deletions_[region + "_5p"];
    int del_3p = deletions_[region + "_3p"];
    original_seqs[region] = gl.seqs_[genes_[region]];
    lengths[region] = original_seqs[region].size() - del_5p - del_3p;
    eroded_seqs[region] = original_seqs[region].substr(del_5p, lengths[region]);
  }
  naive_seq_ = insertions_["fv"] + eroded_seqs["v"] + insertions_["vd"] + eroded_seqs["d"] + insertions_["dj"] + eroded_seqs["j"] + insertions_["jf"];

  int eroded_gl_cpos = gl.cyst_positions_[genes_["v"]] - deletions_["v_5p"] + insertions_["fv"].size();
  int eroded_gl_tpos = gl.tryp_positions_[genes_["j"]] - deletions_["j_5p"];
  int tpos_in_joined_seq = eroded_gl_tpos + insertions_["fv"].size() + eroded_seqs["v"].size() + insertions_["vd"].size() + eroded_seqs["d"].size() + insertions_["dj"].size();
  cyst_position_ = eroded_gl_cpos;
  tryp_position_ = tpos_in_joined_seq;
  cdr3_length_ = tpos_in_joined_seq - eroded_gl_cpos + 3;
}

// ----------------------------------------------------------------------------------------
void RecoEvent::AddAuxiliarySeqs(string name, string seq) {
  assert(seq.size() == seq_.size());  // make sure we already have a first seq, and also that the new seq is the same size
  auxiliary_seq_names_.push_back(name);
  auxiliary_seqs_.push_back(seq);
}

// ----------------------------------------------------------------------------------------
void RecoEvent::Print(GermLines &germlines, size_t cyst_position, size_t final_tryp_position, bool one_line, string extra_indent) {
  assert(0);  // this needs to be updated to match the python version
  string print_seq = seq_;  // this variable is a relic of past transgressions, should be redundantized when I finish what I'm doing
  // need to update this to allow printing of auxiliary sequences... or not. Maybe just do it in the python version. I don't really use this anymore, actually
  // if(print_second_seq)
  //   print_seq = second_seq_;

  if(score_ == -INFINITY) {
    cout << "    " << extra_indent << score_ << endl;
    return;
  }
  // make sure we remembered to fill all the necessary information
  vector<string> deletion_names {"v_3p", "d_5p", "d_3p", "j_5p"};
  vector<string> insertion_names {"vd", "dj"};
  for(auto & region : germlines.regions_) assert(genes_.find(region) != genes_.end());
  for(auto & name : deletion_names) assert(deletions_.find(name) != deletions_.end());
  for(auto & name : insertion_names) assert(insertions_.find(name) != insertions_.end());
  assert(print_seq.size() > 0);

  map<string, string> original_seqs;
  for(auto & region : germlines.regions_)
    original_seqs[region] = germlines.seqs_[genes_[region]];

  if(deletions_.find("v_5p") == deletions_.end()) {  // try to infer the left-hand v "deletion"
    deletions_["v_5p"] = original_seqs["v"].size() + original_seqs["d"].size() + original_seqs["j"].size()
                         - deletions_["v_3p"] - deletions_["d_5p"] - deletions_["d_3p"] - deletions_["j_5p"]
                         + insertions_["vd"].size() + insertions_["dj"].size() - print_seq.size();
  }
  original_seqs["v"] = original_seqs["v"].substr(deletions_["v_5p"]);

  size_t v_length = original_seqs["v"].size() - deletions_["v_3p"];
  size_t d_length = original_seqs["d"].size() - deletions_["d_5p"] - deletions_["d_3p"];
  size_t j_length = original_seqs["j"].size() - deletions_["j_5p"];

  map<string, string> eroded_seqs;
  eroded_seqs["v"] = original_seqs["v"].substr(0, original_seqs["v"].size() - deletions_["v_3p"]);
  eroded_seqs["d"] = original_seqs["d"].substr(deletions_["d_5p"], original_seqs["d"].size() - deletions_["d_3p"] - deletions_["d_5p"]);
  eroded_seqs["j"] = original_seqs["j"].substr(deletions_["j_5p"]);

  size_t germline_v_end = original_seqs["v"].size() - 1;
  size_t germline_d_start = original_seqs["v"].size() - deletions_["v_3p"] + insertions_["vd"].size() - deletions_["d_5p"];
  size_t germline_d_end = germline_d_start + original_seqs["d"].size();
  size_t germline_j_start = germline_d_end + 1 - deletions_["d_3p"] + insertions_["dj"].size() - deletions_["j_5p"];

  // see if we've "eroded" left side of v or right side of j, and if so add some more dots
  // if (deletions_.find("v_5p") != deletions_.end()) {  // truncating this part for the moment
  //   for (size_t i=0; i<deletions_["v_5p"]; ++i)
  // 	print_seq = "." + print_seq;
  // }
  if(deletions_.find("j_3p") != deletions_.end()) {
    for(size_t i = 0; i < deletions_["j_3p"]; ++i)
      print_seq += ".";
  }

  string final_seq;
  TermColors tc;
  for(size_t inuc = 0; inuc < print_seq.size(); ++inuc) {
    size_t ilocal = inuc;
    string new_nuc;
    if(ilocal < v_length) {
      new_nuc = tc.RedifyIfMuted(eroded_seqs["v"][ilocal], print_seq[inuc]);
    } else {
      ilocal -= v_length;
      if(ilocal < insertions_["vd"].size()) {
        new_nuc = tc.RedifyIfMuted(insertions_["vd"][ilocal], print_seq[inuc]);
      } else {
        ilocal -= insertions_["vd"].size();
        if(ilocal < d_length) {
          new_nuc = tc.RedifyIfMuted(eroded_seqs["d"][ilocal], print_seq[inuc]);
        } else {
          ilocal -= d_length;
          if(ilocal < insertions_["dj"].size()) {
            new_nuc = tc.RedifyIfMuted(insertions_["dj"][ilocal], print_seq[inuc]);
          } else {
            ilocal -= insertions_["dj"].size();
            new_nuc = tc.RedifyIfMuted(eroded_seqs["j"][ilocal], print_seq[inuc]);
          }
        }
      }
    }
    // reverse video the conserved codons, if we have them
    if(cyst_position > 0 && final_tryp_position > 0) {
      if(inuc == cyst_position || inuc == final_tryp_position) {
        new_nuc = "\033[7m" + new_nuc;
      } else if(inuc == cyst_position + 2 || inuc == final_tryp_position + 2) {
        new_nuc = new_nuc + "\033[m";
      }
    }
    // and finally tack it onto the final sequence
    final_seq += new_nuc;
  }

  // pad with dots
  for(size_t i = 0; i < deletions_["v_3p"]; ++i)
    eroded_seqs["v"] += ".";
  for(size_t i = 0; i < deletions_["d_5p"]; ++i)
    eroded_seqs["d"] = "." + eroded_seqs["d"];
  for(size_t i = 0; i < deletions_["d_3p"]; ++i)
    eroded_seqs["d"] += ".";
  for(size_t i = 0; i < deletions_["j_5p"]; ++i)
    eroded_seqs["j"] = "." + eroded_seqs["j"];

  string insertion_str;
  for(size_t iv = 0; iv < v_length; ++iv)
    insertion_str += " ";
  insertion_str += insertions_["vd"];
  for(size_t id = 0; id < d_length; ++id)
    insertion_str += " ";
  insertion_str += insertions_["dj"];
  for(size_t ij = 0; ij < j_length; ++ij)
    insertion_str += " ";
  string d_str;
  for(size_t ig = 0; ig < germline_d_start; ++ig)
    d_str += " ";
  d_str += eroded_seqs["d"];
  int ij_max = original_seqs["j"].size() - deletions_["j_5p"] + insertions_["dj"].size() - deletions_["d_3p"];
  // assert(ij_max>=0);
  if(ij_max < 0)
    ij_max = 0;
  for(int ij = 0; ij < ij_max; ++ij)
    d_str += " ";
  string vj_str(eroded_seqs["v"]);
  int ivj_max = germline_j_start - germline_v_end - 2;
  // assert(ivj_max >= 0);
  if(ivj_max < 0)
    ivj_max = 0;
  for(int ivj = 0; ivj < ivj_max; ++ivj)
    vj_str += " ";
  vj_str += eroded_seqs["j"];

  if(!one_line) {
    cout << extra_indent << "    " << insertion_str << "   inserts" << endl;
    cout << extra_indent <<  "    " << d_str << "   " << tc.ColorGene(genes_["d"]) << endl;
    cout << extra_indent << "    " << vj_str << "   " << tc.ColorGene(genes_["v"]) << "," << tc.ColorGene(genes_["j"]) << endl;
  }
  cout << extra_indent << "    " << final_seq << "   " << score_ << endl;
}


// ----------------------------------------------------------------------------------------
KBounds KBounds::LogicalOr(KBounds rhs) {
  KBounds kbr(rhs); // return value
  if(vmin < kbr.vmin) kbr.vmin = vmin;
  if(dmin < kbr.dmin) kbr.dmin = dmin;
  if(vmax > kbr.vmax) kbr.vmax = vmax;
  if(dmax > kbr.dmax) kbr.dmax = dmax;
  return kbr;
}

// ----------------------------------------------------------------------------------------
void Result::check_boundaries(KSet best, KBounds kbounds) {
  // if(kbounds.vmax - kbounds.vmin <= 1 || kbounds.dmax - kbounds.dmin <= 2) return; // if k space is very narrow, we expect the max to be on the boundary, so ignore boundary errors

  int delta(2);  // be VERY VERY CAREFUL about subtracting move than the value off of a size_t. Yes, I know I should follow the google standards and not use unsigned integers but it wasn't my choice at the start... I'll switch eventually

  // see if we need to expand
  if(best.v == kbounds.vmin) {
    boundary_error_ = true;
    better_kbounds_.vmin = max((int)1, (int)kbounds.vmin - delta);
  }
  if(best.v == kbounds.vmax - 1) {
    boundary_error_ = true;
    better_kbounds_.vmax = kbounds.vmax + delta;
  }
  if(best.d == kbounds.dmin) {
    boundary_error_ = true;
    better_kbounds_.dmin = max((int)1, (int)kbounds.dmin - delta);
  }
  if(best.d == kbounds.dmax - 1) {
    boundary_error_ = true;
    better_kbounds_.dmax = kbounds.dmax + delta;
  }

  if(boundary_error_ && better_kbounds_.equals(kbounds))
    could_not_expand_ = true;
}

// ----------------------------------------------------------------------------------------
void HMMHolder::CacheAll() {
  for(auto & region : gl_.regions_) {
    for(auto & gene : gl_.names_[region]) {
      string infname(hmm_dir_ + "/" + gl_.SanitizeName(gene) + ".yaml");
      if(ifstream(infname)) {
        cout << "    read " << infname << endl;
        hmms_[gene] = new Model;
        hmms_[gene]->Parse(infname);
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
Model *HMMHolder::Get(string gene, bool debug) {
  if(hmms_.find(gene) == hmms_.end()) {   // if we don't already have it, read it from disk
    hmms_[gene] = new Model;
    string infname(hmm_dir_ + "/" + gl_.SanitizeName(gene) + ".yaml");
    // if (true) cout << "    read " << infname << endl;
    hmms_[gene]->Parse(infname);
  }
  return hmms_[gene];
}

// ----------------------------------------------------------------------------------------
void HMMHolder::RescaleOverallMuteFreqs(map<string, set<string> > &only_genes, double overall_mute_freq) {
  // WOE BETIDE THEE WHO FORGETETH TO RE-RESET THESE
  // Seriously! If you don't re-rescale 'em when you're done with the sequences to which <overall_mute_freq> correspond, the mute freqs in the hmms will be *wrong*

  // then actually do the rescaling for each necessary gene
  for(auto &region : gl_.regions_) {
    for(auto &gene : only_genes[region]) {
      Get(gene, false)->RescaleOverallMuteFreq(overall_mute_freq);
    }
  }
}

// ----------------------------------------------------------------------------------------
void HMMHolder::UnRescaleOverallMuteFreqs(map<string, set<string> > &only_genes) {
  for(auto &region : gl_.regions_) {
    for(auto &gene : only_genes[region]) {
      Get(gene, false)->UnRescaleOverallMuteFreq();
    }
  }
}

// ----------------------------------------------------------------------------------------
HMMHolder::~HMMHolder() {
  for(auto & entry : hmms_)
    delete entry.second;
}

// ----------------------------------------------------------------------------------------
string HMMHolder::NameString(map<string, set<string> > *only_genes, int max_to_print) {
  // NOTE this doesn't check that we actually have xeverybody in <only_genes>
  TermColors tc;
  map<string, string> region_strs;
  map<string, int> n_genes;
  int n_total_genes(0);
  for(auto &region : gl_.regions_) {
    region_strs[region] = "";
    n_genes[region] = 0;
    for(auto &kv : hmms_) {
      string gene(kv.first);
      if(tc.GetRegion(gene) != region)  // skip genes from other regions
	continue;
      if(only_genes && (*only_genes)[region].count(gene) == 0)  // skip genes not in <only_genes>
	continue;
      n_genes[region] += 1;
      if(region_strs[region].size() > 0)
	region_strs[region] += ":";
      region_strs[region] +=  tc.ColorGene(gene);
    }
    n_total_genes += n_genes[region];
  }

  string return_str;
  for(auto &region : gl_.regions_) {
    if(max_to_print < 0 || n_total_genes <= max_to_print)
      return_str += region_strs[region];
    else
      return_str += to_string(n_genes[region]) + region;
    if(region == "v" || region == "d")
      return_str += "  ";
  }

  return return_str;
}

// ----------------------------------------------------------------------------------------
// truncate sequences in <seqs> to the same length (on both ends of the conserved cysteine), and correspondingly modify <kbvector>
void TruncateSeqs(vector<Sequence> &seqs, vector<KBounds> &kbvector, bool debug) {
  assert(seqs.size() == kbvector.size());  // one kbound for each sequence
  if(debug)
    cout << "    truncating      dleft dright   (cpos, len)" << endl;

  // first find min length to left and right of the cysteine position
  int min_left(-1), min_right(-1);
  for(size_t is=0; is<seqs.size(); ++is) {
    Sequence *seq(&seqs[is]);
    int cpos(seq->cyst_position());
    if(cpos < 0 || cpos >= (int)seq->size())
      throw runtime_error("cpos " + to_string(cpos) + " invalid for " + seq->name() + " (" + seq->undigitized() + ")");
    int dleft = cpos;  // NOTE <dright> includes <cpos>, i.e. dleft + dright = len(seq)
    int dright = seq->size() - cpos;
    if(debug)
      printf("                    %4d %4d      (%-4d, %4d)\n", dleft, dright, (int)cpos, (int)seq->size()); //, seq->name().c_str());
    if(min_left == -1 || dleft < min_left) {
      min_left = dleft;
    }
    if(min_right == -1 || dright < min_right) {
      min_right = dright;
    }
  }
  assert(min_left >= 0 && min_right >= 0);
  // printf("  min left %d right %d\n", min_left, min_right);

  // then truncate all the sequences to these lengths
  if(debug)
    cout << "          chops:" << endl;
  for(size_t is=0; is<seqs.size(); ++is) {
    Sequence *seq(&seqs[is]);
    int cpos(seq->cyst_position());
    int istart(cpos - min_left);
    int istop(cpos + min_right);
    int chopleft(istart);
    int chopright(seq->size() - istop);

    // string truncated_str(seq->undigitized().substr(istart, istop - istart));
    // Sequence trunc_seq(seq->track(), seq->name(), truncated_str, cpos - istart);
    // seqs[is] = trunc_seq;  // replace the old sequence with the truncated one
    seqs[is] = Sequence(*seq, istart, istop - istart);  // replace the old sequence with the truncated one (this sets cpos in the new sequence to cpos in the old sequence minus istart)
    assert(chopleft < (int)kbvector[is].vmin);  // kinda nonsensical if we start chopping off the entire v
    kbvector[is].vmin -= chopleft;
    kbvector[is].vmax -= chopleft;
    
    // printf("%s", seq->name());
    if(debug)
      printf("                    %4d %4d\n", chopleft, chopright);  //, seq->name().c_str());
    // printf("  before", self.sw_info[name]['k_v']['min'], self.sw_info[name]['k_v']['max'], self.sw_info[name]['v_5p_del'], self.sw_info[name]['j_3p_del'], self.sw_info[name]['seq']
    // self.sw_info[name]['seq'] = seq[istart : istop]
    // self.sw_info[name]['k_v']['min'] -= chopleft
    // self.sw_info[name]['k_v']['max'] -= chopleft
    // self.sw_info[name]['v_5p_del'] += chopleft
    // self.sw_info[name]['j_3p_del'] += chopright
    // print '   after', self.sw_info[name]['k_v']['min'], self.sw_info[name]['k_v']['max'], self.sw_info[name]['v_5p_del'], self.sw_info[name]['j_3p_del'], self.sw_info[name]['seq']
  }
}

}
