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
  // PrintCachedTrellisSize();
  Clear();
}

// ----------------------------------------------------------------------------------------
void DPHandler::Clear() {
  // PrintCachedTrellisSize();
  chunk_cache_info_.clear();
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

  assert(kbounds.vmax > kbounds.vmin && kbounds.dmax > kbounds.dmin); // make sure max values for k_v and k_d are greater than their min values
  assert(kbounds.vmin > 0 && kbounds.dmin > 0);  // you get the loveliest little seg fault if you accidentally pass in zero for a lower bound
  if(clear_cache)  // default is true, and you pobably want to leave it that way
    Clear();  // delete all existing trellisi, paths, and logprobs NOTE in principal it kinda ought to be faster to keep everything cached between calls to Run()... but in practice there's a fair bit of overhead to keeping all that stuff hanging around, and it's much more efficient to do the caching in Glomerator (which we already do). So, in sum, it's generally faster to Clear() right here. One exception is if you, say, run viterbi on the same sequence fifty times in a row... then you want to keep the cache around. But why would you do that? In practice the only time you're running on the same sequence many times is in Glomerator, and there we're already doing caching more efficiently at a higher level.
  map<KSet, double> best_scores; // best score for each kset (summed over regions)
  map<KSet, double> total_scores; // total score for each kset (summed over regions)
  map<KSet, map<string, string> > best_genes; // map from a kset to its corresponding triplet of best genes
  if(!args_->dont_rescale_emissions()) {  // reset the emission probabilities in the hmms to reflect the frequences in this particular set of sequences
    assert(overall_mute_freq != -INFINITY);  // make sure the caller remembered to set it
    // NOTE it's super important to *un*set them after you're done
    hmms_.RescaleOverallMuteFreqs(only_genes, overall_mute_freq);
  }

  Result result(kbounds);

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
      if(algorithm_ == "viterbi" && best_scores[kset] != -INFINITY)  // add event to the vector in <result> (we get the event from <paths_>, using <kset> and <best_genes_> as indices)
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
      sprintf(kstr, "%zu [%zu-%zu]  %zu [%zu-%zu]", best_kset.v, kbounds.vmin, kbounds.vmax-1, best_kset.d, kbounds.dmin, kbounds.dmax-1);
    } else {
      prob = *total_score;
      alg_str = "fwd";
      sprintf(kstr, "    [%zu-%zu]     [%zu-%zu]", kbounds.vmin, kbounds.vmax-1, kbounds.dmin, kbounds.dmax-1);
    }
    double cpu_seconds(((clock() - run_start) / (double)CLOCKS_PER_SEC));
    printf("           %s %12.3f   %-25s  %2zuv %2zud %2zuj  %5.1fs   %s\n", alg_str.c_str(), prob, kstr,
	   only_genes["v"].size(), only_genes["d"].size(), only_genes["j"].size(),  // hmms_.NameString(&only_genes, 30)
	   cpu_seconds, seqs.name_str(":").c_str());

    if(result.boundary_error()) {   // not necessarily a big deal yet -- the bounds get automatical expanded
      cout << "             max at boundary:"
	   << " " << best_kset.v << " (" << kbounds.vmin << "-" << kbounds.vmax - 1 << ")"
	   << ", " << best_kset.d << " (" << kbounds.dmin << "-" << kbounds.dmax - 1 << ")"
	   << "    better: " << result.better_kbounds().stringify();
      if(result.could_not_expand())
	cout << " (could not expand)     ";
      cout << "    " << seqs.name_str()  << endl;
    }
  }

  if(!args_->dont_rescale_emissions())  // if we rescaled them above, re-rescale the overall mean mute freqs
    hmms_.UnRescaleOverallMuteFreqs(only_genes);

  return result;
}

// ----------------------------------------------------------------------------------------
void DPHandler::PrintCachedTrellisSize() {
  assert(0);
  // double str_bytes(0.), tr_bytes(0.), path_bytes(0.);
  // for(auto &kv_a : trellisi_) {
  //   string gene(kv_a.first);
  //   for(auto &kv_b : trellisi_[gene]) {
  //     vector<string> cached_query_strs(kv_b.first);
  //     trellis ts(kv_b.second);
  //     for(auto &str : cached_query_strs)
  // 	str_bytes += sizeof(string::value_type) * str.size();
  //     tr_bytes += ts.ApproxBytesUsed();
  //     path_bytes += sizeof(int) * paths_[gene][cached_query_strs].size();  // size is zero if not set (e.g. for fwd)
  //   }
  // }

  // printf("   TOT %.0e   str %.0e   trellis %.0e     path %.0e\n", tr_bytes + str_bytes + path_bytes, tr_bytes, str_bytes, path_bytes);
}

// ----------------------------------------------------------------------------------------
void DPHandler::FillTrellis(Sequences query_seqs, vector<string> query_strs, string gene, double *score, string &origin) {
  *score = -INFINITY;
  // initialize trellis and path
  if(chunk_cache_info_.find(gene) == chunk_cache_info_.end())
    chunk_cache_info_[gene] = map<vector<string>, CacheInfoHolder>();
  origin = "scratch";
  if(!args_->no_chunk_cache()) {   // figure out if we've already got a trellis with a dp table which includes the one we're about to calculate (we should, unless this is the first kset)
    for(auto &kv : chunk_cache_info_[gene]) {
      vector<string> cached_query_strs(kv.first);
      if(cached_query_strs.size() != query_strs.size())
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
        for(size_t iseq = 0; iseq < cached_query_strs.size(); ++iseq)  // neurotic double check
          assert(cached_query_strs[iseq].find(query_strs[iseq]) == 0);
        chunk_cache_info_[gene][query_strs].trellis_ = trellis(hmms_.Get(gene, args_->debug()), query_seqs, &chunk_cache_info_[gene][cached_query_strs].trellis_);  // copy over the required chunk of the old trellis into a new trellis for the current query
        origin = "chunk";
        break;
      }
    }
  }

  if(chunk_cache_info_[gene].count(query_strs) == 0)   // if didn't find a suitable chunk cached trellis
    chunk_cache_info_[gene][query_strs].trellis_ = trellis(hmms_.Get(gene, args_->debug()), query_seqs);
  trellis *trell = &chunk_cache_info_[gene][query_strs].trellis_; // this pointer's just to keep the name short

  // run the actual dp algorithms
  if(algorithm_ == "viterbi") {
    trell->Viterbi();
    *score = trell->ending_viterbi_log_prob();  // NOTE still need to add the gene choice prob to this score (it's done in RunKSet)
    if(trell->ending_viterbi_log_prob() == -INFINITY) {   // no valid path through hmm
      chunk_cache_info_[gene][query_strs].path_ = TracebackPath();
      if(args_->debug() == 2) cout << "                    arg " << gene << " " << *score << " " << origin << endl;
    } else {
      chunk_cache_info_[gene][query_strs].path_ = TracebackPath(hmms_.Get(gene, args_->debug()));
      trell->Traceback(chunk_cache_info_[gene][query_strs].path_);
      assert(trell->ending_viterbi_log_prob() == chunk_cache_info_[gene][query_strs].path_.score());  // NOTE it would be better to just not store the darn score in both the places to start with, rather than worry here about them being the same
    }
    assert(fabs(*score) > 1e-200);
    assert(*score == -INFINITY || chunk_cache_info_[gene][query_strs].path_.size() > 0);
  } else if(algorithm_ == "forward") {
    trell->Forward();
    chunk_cache_info_[gene][query_strs].path_ = TracebackPath();  // avoids violating the assumption that paths_ and trellisi_ have the same entries
    *score = trell->ending_forward_log_prob();  // NOTE still need to add the gene choice prob to this score
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------------------------
void DPHandler::PrintPath(vector<string> query_strs, string gene, double score, string extra_str) {  // NOTE query_str is seq1xseq2 for pair hmm
  if(score == -INFINITY) {
    // cout << "                    " << gene << " " << score << endl;
    return;
  }
  vector<string> path_names = chunk_cache_info_[gene][query_strs].path_.name_vector();
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
    vector<string> path_names = chunk_cache_info_[gene][query_strs].path_.name_vector();
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

  event.SetSeq(seqs[0].name(), seq_strs[0]);
  for(size_t iseq = 1; iseq < seq_strs.size(); ++iseq)
    event.AddAuxiliarySeqs(seqs[iseq].name(), seq_strs[iseq]);
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
          cout << "              " << region << " query " << tc.ColorChars(hmms_.track()->ambiguous_char()[0], "light_blue", tc.ColorMutants("purple", query_strs[is], "", query_strs, hmms_.track()->ambiguous_char())) << endl;  // use the first query_str as reference sequence... could just as well use any other
      } else {
        cout << "              " << region << endl;
      }
    }

    regional_best_scores[region] = -INFINITY;
    regional_total_scores[region] = -INFINITY;
    for(auto & gene : only_genes[region]) {
      double *gene_score(&chunk_cache_info_[gene][query_strs].score_);  // pointed-to value is already set if we have this trellis cached, otherwise not
      bool already_cached = chunk_cache_info_.find(gene) != chunk_cache_info_.end() && chunk_cache_info_[gene].find(query_strs) != chunk_cache_info_[gene].end();
      string origin("ARG");
      if(already_cached) {
        origin = "cached";
      } else {
        FillTrellis(subseqs[region], query_strs, gene, gene_score, origin);  // sets *gene_score to uncorrected score
        double gene_choice_score = log(hmms_.Get(gene, args_->debug())->overall_prob());
        *gene_score = AddWithMinusInfinities(*gene_score, gene_choice_score);  // then correct it for gene choice probs
      }
      if(args_->debug() == 2 && algorithm_ == "viterbi")
        PrintPath(query_strs, gene, *gene_score, origin);

      // add this score to the regional total score
      regional_total_scores[region] = AddInLogSpace(*gene_score, regional_total_scores[region]);  // (log a, log b) --> log a+b, i.e. here we are summing probabilities in log space, i.e. a *or* b
      if(args_->debug() == 2 && algorithm_ == "forward")
        printf("                %6.0e %9.2f  %7.2f  %s  %s\n", exp(*gene_score), *gene_score, regional_total_scores[region], origin.c_str(), tc.ColorGene(gene).c_str());

      // set best regional scores (and the best gene for this kset)
      if(*gene_score > regional_best_scores[region]) {
        regional_best_scores[region] = *gene_score;
        (*best_genes)[kset][region] = gene;
      }

      // watch this space for something pithy
      per_gene_support_this_kset[gene] = *gene_score;
    }

    // return if we didn't find a valid path for this region
    if((*best_genes)[kset].find(region) == (*best_genes)[kset].end()) {
      if(args_->debug() == 2)
        cout << "                  found no gene for " << region << " so skip" << endl;
      return;
    }
  }

  (*best_scores)[kset] = AddWithMinusInfinities(regional_best_scores["v"], AddWithMinusInfinities(regional_best_scores["d"], regional_best_scores["j"]));  // i.e. best_prob = v_prob * d_prob * j_prob (v *and* d *and* j)
  (*total_scores)[kset] = AddWithMinusInfinities(regional_total_scores["v"], AddWithMinusInfinities(regional_total_scores["d"], regional_total_scores["j"]));

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
    for(size_t ip = names.size() - 1; ip >= 0; --ip)  // NOTE unsigned comparison is always true, so could really replace with while(true)
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
      assert(name.find("IGH") == 0);  // Trust but verify, my little ducky, trust but verify.
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
    for(size_t il = names.size() - 1; il >= 0; --il) {  // NOTE unsigned comparison is always true, so could really replace with while(true)
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
  assert(istate >= 0 && istate < names.size());
  if(names[istate].find("IGH") != 0)  // start of state name should be IGH[VDJ]
    throw runtime_error("state not of the form IGH<gene>_<position>: " + names[istate]);
  assert(names[istate].find("IGH") == 0);  // start of state name should be IGH[VDJ]
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
