#include "bcrutils.h"

namespace ham {

// ----------------------------------------------------------------------------------------
TermColors::TermColors() {
  codes_["bold"] = "\033[1m";
  codes_["reverse"] = "\033[7m";
  codes_["purple"] = "\033[95m";
  codes_["blue"] = "\033[94m";
  codes_["light_blue"] = "\033[1;34m";
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
// NOTE GermLines class also has its own GetRegion()
string TermColors::GetRegion(string gene) {
  if(gene.find("IG") != 0 && gene.find("TR") != 0)
    throw runtime_error("gene name '" + gene + "' doesn't begin with 'IG' or 'TR'.");
  char region_char(tolower(gene[3]));
  string region(1, region_char);
  assert(region == "v" || region == "d" || region == "j");
  return region;
}

// ----------------------------------------------------------------------------------------
string TermColors::ColorChars(char ch_to_color, string color, string seq) {
  string return_str;
  for(size_t ic=0; ic<seq.size(); ++ic) {
    if(seq[ic] == ch_to_color)
      return_str += Color(color, seq.substr(ic, 1));
    else
      return_str += seq[ic];
  }
  return return_str;
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
	if(ambiguous_char != "" && (ref[inuc] == ambiguous_char[0] || seq[inuc] == ambiguous_char[0]))  // don't count ambiguous characters as mutated
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
GermLines::GermLines(string gldir, string locus):
  locus_(locus),
  regions_({"v", "d", "j"})
{
  if(!HasDGene(locus_)) {
    string upperlocus(locus_);
    transform(upperlocus.begin(), upperlocus.end(), upperlocus.begin(), ::toupper);
    dummy_d_gene = upperlocus + "Dx-x*x";
  }

  for(auto & region : regions_) {
    names_[region] = vector<string>();
    if(!HasDGene(locus_) && region == "d") {
      names_[region].push_back(dummy_d_gene);
      seqs_[dummy_d_gene] = "A";  // NOTE this choice is also set in python/glutils.py
      continue;
    }
    string infname(gldir + "/" + locus_ + "/" + locus_ + region + ".fasta");
    ifstream ifs(infname);
    if(!ifs.is_open())
      throw runtime_error("germline file " + infname + " d.n.e.");
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

  // get cyst/tryp info
  string infname(gldir + "/" + locus_ + "/extras.csv");
  ifs.open(infname);
  if(!ifs.is_open())
    throw runtime_error("germline file " + infname + " d.n.e.");
  // check header
  getline(ifs, line);
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  header = (SplitString(line, ","));
  assert(header[0] == "gene");
  assert(header[1] == "cyst_position");
  assert(header[2] == "tryp_position");
  assert(header[3] == "phen_position");
  // get info
  while(getline(ifs, line)) {
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    vector<string> info(SplitString(line, ","));
    assert(info[0].find("IG") == 0 || info[0].find("TR") == 0);
    if(info[1] != "")
      cyst_positions_[info[0]] = atoi(info[1].c_str());
    else if(info[2] != "")
      tryp_positions_[info[0]] = atoi(info[2].c_str());
    else if(info[3] != "")  // put the phens into tryp_positions_ for now
      tryp_positions_[info[0]] = atoi(info[3].c_str());
  }
  ifs.close();
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
// NOTE TermColors class also has its own GetRegion()
string GermLines::GetRegion(string gene) {
  if(gene.find("IG") != 0 && gene.find("TR") != 0)
    throw runtime_error("gene name '" + gene + "' doesn't begin with 'IG' or 'TR'.");
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
KBounds KBounds::LogicalOr(KBounds rhs) {
  KBounds kbr(rhs); // return value
  if(vmin < kbr.vmin) kbr.vmin = vmin;
  if(dmin < kbr.dmin) kbr.dmin = dmin;
  if(vmax > kbr.vmax) kbr.vmax = vmax;
  if(dmax > kbr.dmax) kbr.dmax = dmax;
  return kbr;
}

// ----------------------------------------------------------------------------------------
void Result::Finalize(GermLines &gl, map<string, double> &unsorted_per_gene_support, KSet best_kset, KBounds kbounds) {
  assert(!finalized_);

  // sort vector of events by score (i.e. find the best path over ksets)
  sort(events_.begin(), events_.end());
  reverse(events_.begin(), events_.end());
  best_event_ = events_[0];

  // set per-gene support (really just rearranging and sorting the values in DPHandler::per_gene_support_) NOTE make sure to do this *after* sorting
  for(auto &region : gl.regions_) {
    vector<SupportPair> support;  // sorted list of (gene, logprob) pairs for this region
    for(auto &pgit : unsorted_per_gene_support) {
      string gene = pgit.first;
      double logprob = pgit.second;
      if(gl.GetRegion(gene) != region)
	continue;
      support.push_back(SupportPair(gene, logprob));
    }
    sort(support.begin(), support.end());
    reverse(support.begin(), support.end());
    // NOTE we *only* want the *best* event to have its per-gene support set -- because the ones in <events_> only correspond to one kset, it doesn't make sense to have their per-gene supports set (well, they'd just be trivial)
    best_event_.per_gene_support_[region] = support;  // NOTE organization is totally different to that of DPHandler::per_gene_support_
  }

  check_boundaries(best_kset, kbounds);

  finalized_ = true;
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
  if(HasDGene(locus_) && best.d == kbounds.dmin) {  // for light chains we expect/require k_d = 1 with no fuzz
    boundary_error_ = true;
    better_kbounds_.dmin = max((int)1, (int)kbounds.dmin - delta);
  }
  if(HasDGene(locus_) && best.d == kbounds.dmax - 1) {
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
Model *HMMHolder::Get(string gene) {
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
      Get(gene)->RescaleOverallMuteFreq(overall_mute_freq);
    }
  }
}

// ----------------------------------------------------------------------------------------
void HMMHolder::UnRescaleOverallMuteFreqs(map<string, set<string> > &only_genes) {
  for(auto &region : gl_.regions_) {
    for(auto &gene : only_genes[region]) {
      Get(gene)->UnRescaleOverallMuteFreq();
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
void StreamHeader(ofstream &ofs, string algorithm) {
  // NOTE make sure to change this in StreamErrorput() and StreamViterbiOutput() (or StreamForwardOutput()) below!
  if(algorithm == "viterbi")
    ofs << "unique_ids,v_gene,d_gene,j_gene,fv_insertion,vd_insertion,dj_insertion,jf_insertion,v_5p_del,v_3p_del,d_5p_del,d_3p_del,j_5p_del,j_3p_del,logprob,seqs,v_per_gene_support,d_per_gene_support,j_per_gene_support,errors" << endl;
  else if(algorithm == "forward")
    ofs << "unique_ids,logprob,errors" << endl;
  else
    throw runtime_error("bad algorithm " + algorithm);
}

// ----------------------------------------------------------------------------------------
void StreamErrorput(ofstream &ofs, string algorithm, vector<Sequence*> &pseqs, string errors) {
  vector<Sequence> seqs(GetSeqVector(pseqs));
  StreamErrorput(ofs, algorithm, seqs, errors);
}

// ----------------------------------------------------------------------------------------
void StreamErrorput(ofstream &ofs, string algorithm, vector<Sequence> &seqs, string errors) {
  if(algorithm == "viterbi") {
    ofs  // be very, very careful to change this *and* the csv header above at the same time
      << SeqNameStr(seqs, ":")
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << ","
      << "," << SeqStr(seqs, ":")
      << ","
      << ","
      << ","
      << "," << errors
      << endl;
  } else {
    ofs
      << SeqNameStr(seqs, ":")
      << ","
      << "," << errors
      << endl;
  }

}

// ----------------------------------------------------------------------------------------
string PerGeneSupportString(vector<SupportPair> &support) {
  string return_str;
  for(size_t is=0; is<support.size(); ++is) {
    if(is > 0)
      return_str += ";";
    return_str += support[is].gene() + ":" + to_string(support[is].logprob());
  }
  return return_str;
}

// ----------------------------------------------------------------------------------------
void StreamViterbiOutput(ofstream &ofs, RecoEvent &event, vector<Sequence*> &pseqs, string errors) {
  vector<Sequence> seqs(GetSeqVector(pseqs));
  StreamViterbiOutput(ofs, event, seqs, errors);
}

// ----------------------------------------------------------------------------------------
void StreamViterbiOutput(ofstream &ofs, RecoEvent &event, vector<Sequence> &seqs, string errors) {
  string second_seq_name, second_seq;
  ofs  // be very, very careful to change this *and* the csv header above at the same time
    << SeqNameStr(seqs, ":")
    << "," << event.genes_["v"]
    << "," << event.genes_["d"]
    << "," << event.genes_["j"]
    << "," << event.insertions_["fv"]
    << "," << event.insertions_["vd"]
    << "," << event.insertions_["dj"]
    << "," << event.insertions_["jf"]
    << "," << event.deletions_["v_5p"]
    << "," << event.deletions_["v_3p"]
    << "," << event.deletions_["d_5p"]
    << "," << event.deletions_["d_3p"]
    << "," << event.deletions_["j_5p"]
    << "," << event.deletions_["j_3p"]
    << "," << event.score_
    << "," << SeqStr(seqs, ":")
    << "," << PerGeneSupportString(event.per_gene_support_["v"])
    << "," << PerGeneSupportString(event.per_gene_support_["d"])
    << "," << PerGeneSupportString(event.per_gene_support_["j"])
    << "," << errors
    << endl;
}

// ----------------------------------------------------------------------------------------
void StreamForwardOutput(ofstream &ofs, vector<Sequence*> &pseqs, double total_score, string errors) {
  vector<Sequence> seqs(GetSeqVector(pseqs));
  StreamForwardOutput(ofs, seqs, total_score, errors);
}

// ----------------------------------------------------------------------------------------
void StreamForwardOutput(ofstream &ofs, vector<Sequence> &seqs, double total_score, string errors) {
  ofs  // be very, very careful to change this *and* the csv header above at the same time
    << SeqNameStr(seqs, ":")
    << "," << total_score
    << "," << errors
    << endl;
}

// ----------------------------------------------------------------------------------------
string SeqStr(vector<Sequence*> &pseqs, string delimiter) {
  vector<Sequence> seqs(GetSeqVector(pseqs));
  return SeqStr(seqs, delimiter);
}

// ----------------------------------------------------------------------------------------
string SeqStr(vector<Sequence> &seqs, string delimiter) {
  string seq_str;
  for(size_t iseq = 0; iseq < seqs.size(); ++iseq) {
    if(iseq > 0) seq_str += delimiter;
    seq_str += seqs[iseq].undigitized();
  }
  return seq_str;
}

// ----------------------------------------------------------------------------------------
string SeqNameStr(vector<Sequence*> pseqs, string delimiter) {
  vector<Sequence> seqs(GetSeqVector(pseqs));
  return SeqNameStr(seqs, delimiter);
}

// ----------------------------------------------------------------------------------------
string SeqNameStr(vector<Sequence> &seqs, string delimiter) {
  string name_str;
  for(size_t iseq = 0; iseq < seqs.size(); ++iseq) {
    if(iseq > 0) name_str += delimiter;
    name_str += seqs[iseq].name();
  }
  return name_str;
}

// ----------------------------------------------------------------------------------------
bool HasDGene(string locus) {
  if(locus == "igh" || locus == "trb" || locus == "trd")
    return true;
  else if(locus == "igk" || locus == "igl" || locus == "tra" || locus == "trg")
    return false;
  else
    throw runtime_error("unkown locus " + locus);
}

// ----------------------------------------------------------------------------------------
bool FishyMultiSeqAnnotation(size_t n_seqs, RecoEvent &event) {
  if(n_seqs < 3)
    return false;
  vector<string> real_deletions{"v_3p", "d_5p", "d_3p", "j_5p"};
  for(auto &delname : real_deletions) {
    if(event.deletions_[delname] > 2)
      return true;
  }
  return false;
}

// ----------------------------------------------------------------------------------------
vector<Sequence> GetSeqVector(vector<Sequence*> pseqvector) {
  vector<Sequence> seqvector(pseqvector.size());
  for(size_t is=0; is<pseqvector.size(); ++is) {
    if(pseqvector[is] == nullptr)
      throw runtime_error("null sequence pointer in vector of length " + to_string(pseqvector.size()));
    seqvector[is] = *pseqvector[is];
  }
  return seqvector;
}

// ----------------------------------------------------------------------------------------
void runps() {  // NOTE this probably isn't worth using any more, the /proc/self/statm call in glomerator.cc is better
  const int MAX_BUFFER = 255;
  char buffer[MAX_BUFFER];
  FILE *fstr = popen("ps -eo rss,pcpu,pmem,user,stime,args --sort pmem | grep bcrham | grep -v grep", "r");
  string outstr;
  while (fgets(buffer, MAX_BUFFER, fstr) != NULL)
    outstr.append(buffer);
  pclose(fstr);
  cout << endl << "ps: " << endl << outstr << endl;
}

// ----------------------------------------------------------------------------------------
int GetMemVal(string path, string name) {
  ifstream ifs("/proc/" + path);
  string tmpline;
  int memval_kb(-1);
  while(getline(ifs, tmpline)) {
    if(tmpline.find(name) == string::npos)
      continue;
    if(memval_kb != -1)
      throw runtime_error("two values for " + name + " in /proc/" + path);
    vector<string> valstrs(PythonSplit(tmpline));
    if(valstrs.size() != 3)
      throw runtime_error("couldn't split line '" + tmpline + "' into 3 words (got " + to_string(valstrs.size()) + ")");
    assert(valstrs[0] == name + ":");
    assert(valstrs[2] == "kB");
    memval_kb = atoi(valstrs[1].c_str());
  }
  return memval_kb;
}

// ----------------------------------------------------------------------------------------
int GetRss() {
  return GetMemVal("self/status", "VmRSS");
}

// ----------------------------------------------------------------------------------------
int GetMemTot() {
  return GetMemVal("meminfo", "MemTotal");
}

}
