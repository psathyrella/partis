#include "dphandler.h"
namespace ham {
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DPHandler::DPHandler(string algorithm, Args *args, GermLines &gl, HMMHolder &hmms):
  algorithm_(algorithm),
  args_(args),
  gl_(gl),
  hmms_(hmms)
{
}

// ----------------------------------------------------------------------------------------
DPHandler::~DPHandler() {
  Clear();
}

// ----------------------------------------------------------------------------------------
void DPHandler::Clear() {
  scratch_cachefo_.clear();
  paths_.clear();
  scores_.clear();
  per_gene_support_.clear();
}

// ----------------------------------------------------------------------------------------
Sequences DPHandler::GetSubSeqs(Sequences &seqs, KSet kset, string region) {
  // get subsequences for one region
  size_t k_v(kset.v), k_d(kset.d);
  if(region == "v")
    return Sequences(seqs, 0, k_v);  // v region (plus vd insert) runs from zero up to k_v
  else if(region == "d")
    return Sequences(seqs, k_v, k_d);  // d region (plus dj insert) runs from k_v up to k_v + k_d
  else if(region == "j")
    return Sequences(seqs, k_v + k_d, seqs.GetSequenceLength() - k_v - k_d);  // j region runs from k_v + k_d to end
  else
    assert(0);
}

// ----------------------------------------------------------------------------------------
map<string, Sequences> DPHandler::GetSubSeqs(Sequences &seqs, KSet kset) {
  // get subsequences for all regions
  map<string, Sequences> subseqs;
  for(auto & region : gl_.regions_)
    subseqs[region] = GetSubSeqs(seqs, kset, region);
  return subseqs;
}


// ----------------------------------------------------------------------------------------
Result DPHandler::Run(Sequence seq, KBounds kbounds, vector<string> only_gene_list, double overall_mute_freq, bool clear_cache) {
  vector<Sequence> seqvector{seq};
  return Run(seqvector, kbounds, only_gene_list, overall_mute_freq, clear_cache);
}

// ----------------------------------------------------------------------------------------
Result DPHandler::Run(vector<Sequence*> pseqvector, KBounds kbounds, vector<string> only_gene_list, double overall_mute_freq, bool clear_cache) {
  vector<Sequence> seqvector(GetSeqVector(pseqvector));
  return Run(seqvector, kbounds, only_gene_list, overall_mute_freq, clear_cache);
}

// ----------------------------------------------------------------------------------------
Result DPHandler::Run(vector<Sequence> seqvector, KBounds kbounds, vector<string> only_gene_list, double overall_mute_freq, bool clear_cache) {
  clock_t run_start(clock());

  Sequences seqs;
  for(auto &seq : seqvector)
    seqs.AddSeq(seq);

  // convert <only_gene_list> to a set for each region
  map<string, set<string> > only_genes;
  if(only_gene_list.size() > 0) {
    for(auto & region : gl_.regions_)  // initialize
      only_genes[region] = set<string>();
    for(auto & gene : only_gene_list)  // insert each gene in the proper region
      only_genes[gl_.GetRegion(gene)].insert(gene);
    for(auto & region : gl_.regions_) // then make sure we have at least one gene for each region
      if(only_genes[region].size() == 0)
        throw runtime_error("ERROR dphandler didn't get any genes for " + region + " region");
  }

  if(kbounds.vmin == 0 || kbounds.dmin == 0 || kbounds.vmax <= kbounds.vmin || kbounds.dmax <= kbounds.dmin) // make sure max values for k_v and k_d are greater than their min values (it at least used to seg fault if you passed in one of them as zero)
    throw runtime_error("k bounds trivial, nonsensical, or include zero (v: " + to_string(kbounds.vmin) + " " + to_string(kbounds.vmax) + "  d: " + to_string(kbounds.dmin) + " " + to_string(kbounds.dmax) + ")");
  if(clear_cache)  // default is true, and be VERY FUCKING CAREFUL if you change that
    Clear();  // delete all existing trellisi, paths, and logprobs NOTE in principal it kinda ought to be faster to keep everything cached between calls to Run()... but in practice there's a fair bit of overhead to keeping all that stuff hanging around, and it's much more efficient to do the caching in Glomerator (which we already do). So, in sum, it's generally faster to Clear() right here. One exception is if you, say, run viterbi on the same sequence fifty times in a row... then you want to keep the cache around. But why would you do that? In practice the only time you're running on the same sequence many times is in Glomerator, and there we're already doing caching more efficiently at a higher level.
  map<KSet, double> best_scores; // best score for each kset (summed over regions)
  map<KSet, double> total_scores; // total score for each kset (summed over regions)
  map<KSet, map<string, string> > best_genes; // map from a kset to its corresponding triplet of best genes
  if(!args_->dont_rescale_emissions()) {  // reset the emission probabilities in the hmms to reflect the frequences in this particular set of sequences
    assert(overall_mute_freq != -INFINITY);  // make sure the caller remembered to set it
    // NOTE it's super important to *un*set them after you're done
    hmms_.RescaleOverallMuteFreqs(only_genes, overall_mute_freq);
  }

  Result result(kbounds, args_->locus());

  // loop over k_v k_d space
  double best_score(-INFINITY);
  KSet best_kset(0, 0);
  double *total_score = &result.total_score_;  // total score for all ksets
  int n_too_long(0), n_run(0), n_total(0);
  for(size_t k_v = kbounds.vmax - 1; k_v >= kbounds.vmin; --k_v) {  // loop in reverse order to facilitate chunk caching: in principle we calculate V once the first time through, and after that can just copy over pieces of the first dp table (roughly the same for D and J)
    for(size_t k_d = kbounds.dmax - 1; k_d >= kbounds.dmin; --k_d) {
      ++n_total;
      if(k_v + k_d >= seqs.GetSequenceLength()) {
        ++n_too_long;
        continue;
      }
      KSet kset(k_v, k_d);
      RunKSet(seqs, kset, only_genes, &best_scores, &total_scores, &best_genes);
      ++n_run;
      *total_score = AddInLogSpace(total_scores[kset], *total_score);  // sum up the probabilities for each kset, log P_tot = log \sum_i P_k_i
      if(args_->debug() == 2 && algorithm_ == "forward") printf("            %9.2f (%.1e)  tot: %7.2f\n", total_scores[kset], exp(total_scores[kset]), *total_score);
      if(best_scores[kset] > best_score) {
        best_score = best_scores[kset];
        best_kset = kset;
      }
      if(algorithm_ == "viterbi" && best_scores[kset] != -INFINITY)  // add event to the vector in <result>
        result.PushBackRecoEvent(FillRecoEvent(seqs, kset, best_genes[kset], best_scores[kset]));
    }
  }
  if(args_->debug() && n_too_long > 0) cout << "      skipped " << n_too_long << " (of " << n_total << ") k sets 'cause they were longer than the sequence (ran " << n_run << ")" << endl;

  // return if no valid path
  if(best_kset.v == 0 && best_kset.d == 0) {
    cout << "    no valid paths for query " << seqs.name_str() << endl;
    result.no_path_ = true;
    return result;
  }

  if(algorithm_ == "viterbi")
    result.Finalize(gl_, per_gene_support_, best_kset, kbounds);

  // print debug info
  if(args_->debug()) {
    double prob;
    string alg_str;
    char kstr[300];
    if(algorithm_ == "viterbi") {
      prob = best_score;
      alg_str = "vtb";
      sprintf(kstr, "%zu [%zu-%zu)  %zu [%zu-%zu)", best_kset.v, kbounds.vmin, kbounds.vmax, best_kset.d, kbounds.dmin, kbounds.dmax);
    } else {
      prob = *total_score;
      alg_str = "fwd";
      sprintf(kstr, "    [%zu-%zu)     [%zu-%zu)", kbounds.vmin, kbounds.vmax, kbounds.dmin, kbounds.dmax);
    }
    double cpu_seconds(((clock() - run_start) / (double)CLOCKS_PER_SEC));
    printf("           %s %12.3f   %-25s  %2zuv %2zud %2zuj  %5.2fs   %s\n", alg_str.c_str(), prob, kstr,
	   only_genes["v"].size(), only_genes["d"].size(), only_genes["j"].size(),  // hmms_.NameString(&only_genes, 30)
	   cpu_seconds, seqs.name_str(":").c_str());

    if(result.boundary_error()) {   // not necessarily a big deal yet -- the bounds get automatical expanded
      // cout << "             max at boundary:"
      // 	   << " " << best_kset.v << " [" << kbounds.vmin << "-" << kbounds.vmax  << ")"
      // 	   << ", " << best_kset.d << " [" << kbounds.dmin << "-" << kbounds.dmax  << ")"
      // 	   << "    better: " << result.better_kbounds().stringify();
      // if(result.could_not_expand())
      // 	cout << " (could not expand)     ";
      // cout << "    " << seqs.name_str()  << endl;
    }
  }

  if(!args_->dont_rescale_emissions())  // if we rescaled them above, re-rescale the overall mean mute freqs
    hmms_.UnRescaleOverallMuteFreqs(only_genes);

  return result;
}

// ----------------------------------------------------------------------------------------
void DPHandler::HandleFishyAnnotations(Result &multi_seq_result, vector<Sequence*> pqry_seqs, KBounds kbounds, vector<string> only_gene_list, double overall_mute_freq) {
  vector<Sequence> qry_seqs(GetSeqVector(pqry_seqs));
  return HandleFishyAnnotations(multi_seq_result, qry_seqs, kbounds, only_gene_list, overall_mute_freq);
}

// ----------------------------------------------------------------------------------------
void DPHandler::HandleFishyAnnotations(Result &multi_seq_result, vector<Sequence> qry_seqs, KBounds kbounds, vector<string> only_gene_list, double overall_mute_freq) {
  Sequence naive_seq(qry_seqs[0].track(), "naive-seq", multi_seq_result.best_event().naive_seq_);
  Result naive_result = Run(naive_seq, kbounds, only_gene_list, overall_mute_freq);
  RecoEvent &naive_event(naive_result.best_event());
  RecoEvent &multi_event(multi_seq_result.best_event());
  for(auto &region : gl_.regions_)
    multi_event.SetGene(region, naive_event.genes_[region]);

  // vector<string> real_deletions{"v_3p", "d_5p", "d_3p", "j_5p"};
  vector<string> all_deletions{"v_5p", "v_3p", "d_5p", "d_3p", "j_5p", "j_3p"};
  for(auto &delname : all_deletions)
    multi_event.SetDeletion(delname, naive_event.deletions_[delname]);  // NOTE they might be the same, but I think the multi-event (real) insertions would be better here, but it'd take some effort to get the right pieces of them

  // vector<string> real_insertions{"vd", "dj"};
  vector<string> all_insertions{"fv", "vd", "dj", "jf"};
  for(auto &insert_name : all_insertions)
    multi_event.SetInsertion(insert_name, naive_event.insertions_[insert_name]);

  multi_event.per_gene_support_ = naive_event.per_gene_support_;
}

// ----------------------------------------------------------------------------------------
void DPHandler::PrintCachedTrellisSize() {
  // double str_bytes(0.), tr_bytes(0.), path_bytes(0.);
  // for(auto &kv_a : cachefo_) {
  //   string gene(kv_a.first);
  //   for(auto &kv_b : cachefo_[gene]) {
  //     vector<string> cached_query_strs(kv_b.first);
  //     Trellis ts(kv_b.second.trellis_);
  //     for(auto &str : cached_query_strs)
  // 	str_bytes += sizeof(string::value_type) * str.size();
  //     tr_bytes += ts.ApproxBytesUsed();
  //     path_bytes += sizeof(int) * kv_b.second.path_.size();  // size is zero if not set (e.g. for fwd)
  //   }
  // }

  // printf("   TOT %.0e   str %.0e   trellis %.0e     path %.0e\n", tr_bytes + str_bytes + path_bytes, tr_bytes, str_bytes, path_bytes);
}

// ----------------------------------------------------------------------------------------
void DPHandler::FillTrellis(KSet kset, Sequences query_seqs, vector<string> query_strs, string gene, string &origin) {

  Trellis *cached_trellis(nullptr);
  if(!args_->no_chunk_cache()) {   // figure out if we've already got a trellis with a dp table which includes the one we're about to calculate (we should, unless this is the first kset)
    // NOTE we're no longer looking through previously chunk cached cachefo here. Which I think is ok, but possible only because we loop over ksets in decreasing order (?)
    for(auto &kv : scratch_cachefo_[gene]) {  // kv: (query string vector, trellis)
      vector<string> cached_query_strs(kv.first);
      if(cached_query_strs.size() != query_strs.size())  // have to have same number of sequences (it'd be much harder for this to happen now that I'm now reusing dphandlers)
	continue;

      // loop over all the query strings for this trellis to see if they all match
      bool found_match(true);
      for(size_t iseq = 0; iseq < cached_query_strs.size(); ++iseq) {  // NOTE this starts to seem like it might be bottlenecking me when I'm applying it for short d sequences
        if(cached_query_strs[iseq].find(query_strs[iseq]) != 0) {  // if <query_strs[iseq]> (the current query) doesn't appear starting at position zero in <cached_query_strs[iseq]> (a previously cached query), we'll need to recalculate
          found_match = false;
          break;
        }
      }

      // if they all match, then use it
      if(found_match) {
	cached_trellis = &scratch_cachefo_[gene][cached_query_strs];  // will copy over the required chunk of the old trellis into a new trellis for the current query
        break;
      }
    }
  }

  Trellis tmptrell(hmms_.Get(gene), query_seqs, cached_trellis);  // NOTE chunk cached trellisi don't get kept around -- we should be able to always just go back to the original one
  Trellis *trell(&tmptrell);  // convenience pointer
  if(cached_trellis == nullptr) {   // if we didn't find a suitable chunk cached trellis
    scratch_cachefo_[gene][query_strs] = Trellis(hmms_.Get(gene), query_seqs);
    trell = &scratch_cachefo_[gene][query_strs];
    origin = "scratch";
  } else {
    origin = "chunk";
  }

  // run the actual dp algorithms
  double uncorrected_score;  // still need to tack on the gene choice prob to this score
  if(algorithm_ == "viterbi") {
    trell->Viterbi();
    uncorrected_score = trell->ending_viterbi_log_prob();
    paths_[gene][kset] = TracebackPath(hmms_.Get(gene));
    if(uncorrected_score != -INFINITY)   // if there's a valid path
      trell->Traceback(paths_[gene][kset]);
  } else if(algorithm_ == "forward") {
    trell->Forward();
    uncorrected_score = trell->ending_forward_log_prob();
  } else {
    assert(0);
  }

  // correct the score for gene choice probs
  double gene_choice_score = log(hmms_.Get(gene)->overall_prob());
  scores_[gene][kset] = AddWithMinusInfinities(uncorrected_score, gene_choice_score);
}

// ----------------------------------------------------------------------------------------
void DPHandler::PrintPath(KSet kset, vector<string> query_strs, string gene, double score, string extra_str) {  // NOTE query_str is seq1xseq2 for pair hmm
  if(score == -INFINITY) {
    // cout << "                    " << gene << " " << score << endl;
    return;
  }
  vector<string> path_names = paths_[gene][kset].name_vector();
  if(path_names.size() == 0) {
    if(args_->debug()) cout << "                     " << gene << " has no valid path" << endl;
    return;
  }
  assert(path_names.size() > 0);  // this will happen if the ending viterbi prob is 0, i.e. if there's no valid path through the hmm (probably the sequence or hmm lengths are screwed up)
  assert(path_names.size() == query_strs[0].size());
  // for(auto &p : path_names)
  //   cout << p;
  // cout << endl;
  string left_insert = GetInsertion("left", path_names);
  string right_insert = GetInsertion("right", path_names);
  size_t left_erosion_length = GetErosionLength("left", path_names, gene);
  size_t right_erosion_length = GetErosionLength("right", path_names, gene);

  TermColors tc;

  // make a string for the germline match
  string germline(gl_.seqs_[gene]);
  string modified_germline = germline.substr(left_erosion_length, germline.size() - right_erosion_length - left_erosion_length);  // remove deletions
  modified_germline = left_insert + modified_germline + right_insert;  // add insertions to either end
  assert(modified_germline.size() == query_strs[0].size());
  assert(germline.size() + left_insert.size() - left_erosion_length - right_erosion_length + right_insert.size() == query_strs[0].size());
  string match_str = tc.ColorMutants("red", modified_germline, "", query_strs, hmms_.track()->ambiguous_char());  // color it relative to the query strings
  match_str = tc.ColorChars(hmms_.track()->ambiguous_char()[0], "light_blue", match_str);
  string l_e_str(to_string(left_erosion_length)), r_e_str(to_string(right_erosion_length));
  if(left_erosion_length > 0)
    match_str = string(max(0, 2 - (int)l_e_str.size()), ' ') + "." + l_e_str + "." + match_str;
  else
    match_str = "    " + match_str;
  if(right_erosion_length > 0)
    match_str += string(max(0, 2 - (int)r_e_str.size()), ' ') + "." + r_e_str + ".";
  else
    match_str += "    ";

  // if there were insertions, we indicate this on a separate line below
  string insert_str(germline.size() - right_erosion_length - left_erosion_length, ' ');
  insert_str = string(left_insert.size(), 'i') + insert_str + string(right_insert.size(), 'i');

  // NOTE this doesn't include the overall gene prob!
  cout << "                    " << match_str << "  " << extra_str << setw(12) << score << setw(25) << gene << endl;
  if(left_insert.size() + right_insert.size() > 0)
    cout << "                      " << "  " << tc.Color("yellow", insert_str) << "  " << endl;
}

// ----------------------------------------------------------------------------------------
RecoEvent DPHandler::FillRecoEvent(Sequences &seqs, KSet kset, map<string, string> &best_genes, double score) {
  RecoEvent event;
  vector<string> seq_strs(seqs.n_seqs(), "");  // build up these strings summing over each regions
  for(auto & region : gl_.regions_) {
    vector<string> query_strs(GetQueryStrs(seqs, kset, region));
    if(best_genes.find(region) == best_genes.end()) {
      seqs.Print();
    }
    assert(best_genes.find(region) != best_genes.end());
    string gene(best_genes[region]);
    vector<string> path_names = paths_[gene][kset].name_vector();
    if(path_names.size() == 0) {
      if(args_->debug()) cout << "                     " << gene << " has no valid path" << endl;
      event.SetScore(-INFINITY);
      return event;
    }
    assert(path_names.size() > 0);
    assert(path_names.size() == query_strs[0].size());
    event.SetGene(region, gene);

    // set right-hand deletions
    event.SetDeletion(region + "_3p", GetErosionLength("right", path_names, gene));
    // and left-hand deletions
    event.SetDeletion(region + "_5p", GetErosionLength("left", path_names, gene));

    SetInsertions(region, path_names, &event);  // NOTE this sets the insertion *only* according to the *first* sequence. Which makes sense at the moment, since the RecoEvent class is only designed to represent a single sequence

    for(size_t iseq = 0; iseq < seq_strs.size(); ++iseq)
      seq_strs[iseq] += query_strs[iseq];
  }

  event.SetScore(score);
  event.SetNaiveSeq(gl_);

  // NOTE we do *not* want to set the event's per_gene_support_, since the event we're filling is only for one kset, while per_gene_support_ should be summed over ksets
  
  return event;
}

// ----------------------------------------------------------------------------------------
vector<string> DPHandler::GetQueryStrs(Sequences &seqs, KSet kset, string region) {
  Sequences query_seqs(GetSubSeqs(seqs, kset, region));
  vector<string> query_strs;
  for(size_t iseq = 0; iseq < seqs.n_seqs(); ++iseq)
    query_strs.push_back(query_seqs[iseq].undigitized());
  return query_strs;
}

// ----------------------------------------------------------------------------------------
KSet DPHandler::FindPartialCacheMatch(string region, string gene, KSet kset) {
  // this is just to avoid having to store all the query string vectors (they get big)
  if(scores_[gene].find(kset) != scores_[gene].end())  // the exact same kset shouldn't actually be in there (except maybe if we're rerunning with expanded boundaries?) but I think we may as well check
    return kset;
  if(region == "v") {
    for(auto &kv : scores_[gene]) {  // kv: (KSet, double)
      if(kv.first.v == kset.v)  // for v, we just need k_v to be the same
	return kv.first;
    }
  } else if(region == "j") {
    for(auto &kv : scores_[gene]) {  // kv: (KSet, double)
      if(kv.first.v + kv.first.d == kset.v + kset.d)  // for j, we need k_v and k_d to sum to the same thing
	return kv.first;
    }
  }

  return KSet(0, 0);
}

// ----------------------------------------------------------------------------------------
void DPHandler::InitCache(string gene) {
  if(scores_.find(gene) == scores_.end()) {
    scratch_cachefo_[gene] = map<vector<string>, Trellis>();
    paths_[gene] = map<KSet, TracebackPath>();
    scores_[gene] = map<KSet, double>();
  }
}

// ----------------------------------------------------------------------------------------
void DPHandler::RunKSet(Sequences &seqs, KSet kset, map<string, set<string> > &only_genes, map<KSet, double> *best_scores, map<KSet, double> *total_scores, map<KSet, map<string, string> > *best_genes) {
  map<string, Sequences> subseqs(GetSubSeqs(seqs, kset));
  (*best_scores)[kset] = -INFINITY;
  (*total_scores)[kset] = -INFINITY;  // total log prob of this kset, i.e. log(P_v * P_d * P_j), where e.g. P_v = \sum_i P(v_i k_v)
  (*best_genes)[kset] = map<string, string>();
  map<string, double> regional_best_scores; // the best score for each region
  map<string, double> regional_total_scores; // the total score for each region, i.e. log P_v
  map<string, double> per_gene_support_this_kset;
  if(args_->debug() == 2) {
    printf("         %3d%3d", (int)kset.v, (int)kset.d);
    if(algorithm_ == "forward")
      printf(" %6s %9s  %7s  %7s", "prob", "logprob", "total", "origin");
    printf(" %s\n", "---------------");
  }
  for(auto & region : gl_.regions_) {
    vector<string> query_strs(GetQueryStrs(seqs, kset, region));

    TermColors tc;
    if(args_->debug() == 2) {
      if(algorithm_ == "viterbi") {
        cout << "                " << region << " query " << tc.ColorChars(hmms_.track()->ambiguous_char()[0], "light_blue", query_strs[0]) << endl;
        for(size_t is = 1; is < query_strs.size(); ++is)
          cout << "                " << region << " query " << tc.ColorChars(hmms_.track()->ambiguous_char()[0], "light_blue", tc.ColorMutants("purple", query_strs[is], "", query_strs, hmms_.track()->ambiguous_char())) << endl;  // use the first query_str as reference sequence... could just as well use any other
      } else {
        cout << "              " << region << endl;
      }
    }

    regional_best_scores[region] = -INFINITY;
    regional_total_scores[region] = -INFINITY;
    for(auto & gene : only_genes[region]) {
      InitCache(gene);
      string origin;
      KSet partial_cache_match(FindPartialCacheMatch(region, gene, kset));  // "partial" in the sense that only this region's query sequence(s) need to be the same
      if(!partial_cache_match.isnull()) {  // first see if we have a match for these exact strings
	paths_[gene][kset] = paths_[gene][partial_cache_match];
	scores_[gene][kset] = scores_[gene][partial_cache_match];
	// NOTE that we don't put anything about this gene/kset combo into the trellis caches. Which is fine now, since later we'll only need the path and score info
	origin = "cached";
      } else {  // no exact cache match, so proceed to check for chunk caching (if that fails it'll actually calculate things)
	FillTrellis(kset, subseqs[region], query_strs, gene, origin);
      }

      double gene_score(scores_[gene][kset]);  // convenience variable
      if(args_->debug() == 2 && algorithm_ == "viterbi")
        PrintPath(kset, query_strs, gene, gene_score, origin);

      // add this score to the regional total score
      regional_total_scores[region] = AddInLogSpace(gene_score, regional_total_scores[region]);  // (log a, log b) --> log a+b, i.e. here we are summing probabilities in log space, i.e. a *or* b
      if(args_->debug() == 2 && algorithm_ == "forward")
        printf("                %6.0e %9.2f  %7.2f  %s  %s\n", exp(gene_score), gene_score, regional_total_scores[region], origin.c_str(), tc.ColorGene(gene).c_str());

      // set best regional scores (and the best gene for this kset)
      if(gene_score > regional_best_scores[region]) {
        regional_best_scores[region] = gene_score;
        (*best_genes)[kset][region] = gene;
      }

      // watch this space for something pithy
      per_gene_support_this_kset[gene] = gene_score;
    }

    // return if we didn't find a valid path for this region
    if((*best_genes)[kset].find(region) == (*best_genes)[kset].end()) {
      if(args_->debug() == 2)
        cout << "                  found no gene for " << region << " so skip" << endl;
      return;
    }
  }

  // store the results
  (*best_scores)[kset] = AddWithMinusInfinities(regional_best_scores["v"], AddWithMinusInfinities(regional_best_scores["d"], regional_best_scores["j"]));  // i.e. best_prob = v_prob * d_prob * j_prob (v *and* d *and* j)
  (*total_scores)[kset] = AddWithMinusInfinities(regional_total_scores["v"], AddWithMinusInfinities(regional_total_scores["d"], regional_total_scores["j"]));

  // work out per-gene support
  for(auto &region : gl_.regions_) {  // we have to do this in a separate loop because we need to know what the regional_best_scores are for the other regions
    for(auto &gene : only_genes[region]) {
      // first multiply the prob for this kset by the *total* for the other two regions
      double score_this_kset(0);  // not -INFINITY, since we're multiplying probabilities
      for(auto &tmpreg : gl_.regions_) {
      	if(tmpreg == region)
      	  score_this_kset = AddWithMinusInfinities(score_this_kset, per_gene_support_this_kset[gene]);
      	else
      	  score_this_kset = AddWithMinusInfinities(score_this_kset, regional_best_scores[tmpreg]);  // i.e. we use the best genes in the other two regions, but single out this gene in its region
      }

      if(per_gene_support_.count(gene) == 0)
      	per_gene_support_[gene] = -INFINITY;
      // per_gene_support_[gene] = AddInLogSpace(per_gene_support_[gene], score_this_kset);  // also, if you do it this way, a large fraction of the events have different viterbi and best-supported d genes
      if(score_this_kset > per_gene_support_[gene])  // NOTE we could also add up the scores for every kset, but what we want to compare to is the viterbi prob for the best annotation, so this is cleaner and clearer, i.e. it doesn't muddle up viterbi and forward probs
      	per_gene_support_[gene] = score_this_kset;
    }
  }
}

// ----------------------------------------------------------------------------------------
void DPHandler::SetInsertions(string region, vector<string> path_names, RecoEvent *event) {
  Insertions ins;
  for(auto & insertion : ins[region]) {  // loop over the boundaries (vd and dj)
    string side(insertion == "jf" ? "right" : "left");
    string inserted_bases = GetInsertion(side, path_names);
    event->SetInsertion(insertion, inserted_bases);
  }
}

// ----------------------------------------------------------------------------------------
size_t DPHandler::GetInsertStart(string side, size_t path_length, size_t insert_length) {
  if(side == "left") {
    return 0;
  } else if(side == "right") {
    return path_length - insert_length;
  } else {
    throw runtime_error("ERROR side must be left or right, not \"" + side + "\"");
  }

}

// ----------------------------------------------------------------------------------------
string DPHandler::GetInsertion(string side, vector<string> names) {
  string inserted_bases;
  if(side == "left") {
    for(auto & name : names) {
      if(name.find("insert") == 0)
        inserted_bases = inserted_bases + name.back();  // last character is "germline-like" base, e.g. insert_left_C
      else
        break;
    }
  } else if(side == "right") {
    for(size_t ip = names.size() - 1; true; --ip)
      if(names[ip].find("insert") == 0)
        inserted_bases = names[ip].back() + inserted_bases;  // last character is "germline-like" base, e.g. insert_left_C
      else
        break;
  } else {
    throw runtime_error("ERROR side must be left or right, not \"" + side + "\"");
  }

  for(size_t ib=0; ib<inserted_bases.size(); ++ib)
    hmms_.track()->symbol_index(inserted_bases.substr(ib, 1));  // throws exception if we've got a character that isn't handled by the track
    
  return inserted_bases;
}

// ----------------------------------------------------------------------------------------
size_t DPHandler::GetErosionLength(string side, vector<string> names, string gene_name) {
  // NOTE this does *not* count a bunch of Ns at the end as an erosion, that interpretation is made in partitiondriver.py

  string germline(gl_.seqs_[gene_name]);

  // first check if we eroded the entire sequence. If so we can't say how much was left and how much was right, so just (integer) divide by two (arbitrarily giving one side the odd base if necessary)
  bool its_inserts_all_the_way_down(true);
  for(auto & name : names) {
    if(name.find("insert") != 0) {
      assert(name.find("IG") == 0 || name.find("TR") == 0);  // Trust but verify, my little ducky, trust but verify.
      its_inserts_all_the_way_down = false;
      break;
    }
  }
  if(its_inserts_all_the_way_down) {   // entire sequence is inserts, so there's no way to tell which part is a left erosion and which is a right erosion
    if(side == "left")
      return floor(float(germline.size()) / 2);
    else if(side == "right")
      return ceil(float(germline.size()) / 2);
    else
      throw runtime_error("ERROR bad side: " + side);
  }

  // find the index in <names> up to which we eroded
  size_t istate(0);  // index (in path) of first non-eroded state
  if(side == "left") { // to get left erosion length we look at the first non-insert state in the path
    for(size_t il = 0; il < names.size(); ++il) { // loop over each state from left to right
      if(names[il].find("insert") == 0) {   // skip any insert states on the left
        continue;
      } else {  // found the leftmost non-insert state -- that's the one we want
        istate = il;
        break;
      }
    }
  } else if(side == "right") { // and for the righthand one we need the last non-insert state
    for(size_t il = names.size() - 1; true; --il) {
      if(names[il].find("insert") == 0) {   // skip any insert states on the left
        continue;
      } else {  // found the leftmost non-insert state -- that's the one we want
        istate = il;
        break;
      }
    }
  } else {
    assert(0);
  }

  // then find the state number (in the hmm's state numbering scheme) of the state found at that index in the viterbi path
  assert(istate < names.size());
  if(names[istate].find("IG") != 0 && names[istate].find("TR") != 0)  // start of state name should be {IG,TR}[HKL][VDJ]
    throw runtime_error("state not of the form {IG,TR}[HKL]<gene>_<position>: " + names[istate]);
  string state_index_str = names[istate].substr(names[istate].find_last_of("_") + 1);
  size_t state_index = atoi(state_index_str.c_str());

  size_t length(0);
  if(side == "left") {
    length = state_index;
  } else if(side == "right") {
    size_t germline_length = gl_.seqs_[gene_name].size();
    length = germline_length - state_index - 1;
  } else {
    assert(0);
  }

  return length;
}

}
