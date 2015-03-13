#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <map>
#include <iomanip>
// #include <ctime>
#include <fstream>
#include <cfenv>
#include "jobholder.h"
#include "germlines.h"
#include "text.h"
#include "args.h"
#include "tclap/CmdLine.h"
using namespace TCLAP;
using namespace ham;
using namespace std;

// class Glomerator {
// public:
//   Glomerator();
// private:
//   map<string, Sequences> current_partition_;
//   map<string, vector<string> > only_genes_;
//   map<string, KBounds> kbinfo_;
//   map<string, double> cached_log_probs_;
//   vector<double> 
// }

// ----------------------------------------------------------------------------------------
// read input sequences from file and return as vector of sequences
vector<Sequences> GetSeqs(Args &args, Track *trk) {
  vector<Sequences> all_seqs;
  set<string> all_names;  // keeps track which sequences we've added to make sure we weren't passed duplicates
  for(size_t iqry = 0; iqry < args.str_lists_["names"].size(); ++iqry) { // loop over queries, where each query can be composed of one, two, or k sequences
    Sequences seqs;
    assert(args.str_lists_["names"][iqry].size() == args.str_lists_["seqs"][iqry].size());
    for(size_t iseq = 0; iseq < args.str_lists_["names"][iqry].size(); ++iseq) { // loop over each sequence in that query
      Sequence sq(trk, args.str_lists_["names"][iqry][iseq], args.str_lists_["seqs"][iqry][iseq]);

      if(all_names.count(sq.name()))
	throw runtime_error("ERROR tried to add sequence with name " + sq.name() + " twice in bcrham::GetSeqs");
      else
	all_names.insert(sq.name());

      seqs.AddSeq(sq);
    }
    all_seqs.push_back(seqs);
  }
  assert(all_seqs.size() == args.str_lists_["names"].size());
  return all_seqs;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score, string errors);
void print_forward_scores(double numerator, vector<double> single_scores, double bayes_factor);
void hierarch_agglom(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args &args, ofstream &ofs);

// // ----------------------------------------------------------------------------------------
// void check_cache(string name, map<string, double> &cache) {
//   if(cached.count(name))
//     cout << "    cached " << cache[name] << " - " << cvals[ic] << " = " << cached_log_probs[cnames[ic]] - cvals[ic] << endl;
//   else
//     cached_log_probs[cnames[ic]] = cvals[ic];
  
// // ----------------------------------------------------------------------------------------
// string sort_name_list(string unsorted_str) {
//   // alphabetically sort the space-separated name list in <unsorted_str>
//   vector<string> unsorted_vector(SplitString(unsorted_str, " "));
//   sort(unsorted_vector.begin(), unsorted_vector.end());
//   string return_str;
//   for(size_t ic=0; ic<unsorted_vector.size(); ++ic) {
//     if(ic > 0)
//       return_str += " ";
//     return_str += unsorted_vector[ic];
//   }
//   return return_str;    
// }

// ----------------------------------------------------------------------------------------
void get_result(Args &args, JobHolder &jh, string name, Sequences &seqs, KBounds &kbounds, map<string, double> &cached_log_probs, string &errors) {
  if(cached_log_probs.count(name))  // already did it
    return;
    
  Result result(kbounds);
  bool stop(false);
  do {
    result = jh.Run(seqs, kbounds);
    kbounds = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args.debug() && !stop)
      cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  cached_log_probs[name] = result.total_score();
  if(result.boundary_error())
    errors = "boundary";
}

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  srand(time(NULL));
  Args args(argc, argv);
  // write csv output headers
  ofstream ofs;
  ofs.open(args.outfile());
  assert(ofs.is_open());
  if(args.algorithm() == "viterbi")
    ofs << "unique_ids,v_gene,d_gene,j_gene,fv_insertion,vd_insertion,dj_insertion,jf_insertion,v_5p_del,v_3p_del,d_5p_del,d_3p_del,j_5p_del,j_3p_del,score,seqs,errors" << endl;
  else if(args.partition())
    ofs << "partition,score,errors" << endl;
  else if(args.algorithm() == "forward")
    ofs << "unique_ids,score,errors" << endl;

  // init some infrastructure
  vector<string> characters {"A", "C", "G", "T"};
  Track trk("NUKES", characters);
  vector<Sequences> qry_seq_list(GetSeqs(args, &trk));
  GermLines gl(args.datadir());
  HMMHolder hmms(args.hmmdir(), gl);
  // hmms.CacheAll();

  if(args.partition()) {
    hierarch_agglom(hmms, gl, qry_seq_list, args, ofs);
    ofs.close();
    return 0;
  }
  
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    if(args.debug()) cout << "  ---------" << endl;
    KSet kmin(args.integers_["k_v_min"][iqry], args.integers_["k_d_min"][iqry]);
    KSet kmax(args.integers_["k_v_max"][iqry], args.integers_["k_d_max"][iqry]);
    KBounds kbounds(kmin, kmax);
    Sequences qry_seqs(qry_seq_list[iqry]);

    JobHolder jh(gl, hmms, args.algorithm(), args.str_lists_["only_genes"][iqry]);
    jh.SetDebug(args.debug());
    jh.SetChunkCache(args.chunk_cache());
    jh.SetNBestEvents(args.n_best_events());

    Result result(kbounds);
    vector<Result> denom_results(qry_seqs.n_seqs(), result);  // only used for forward if n_seqs > 1
    double numerator(-INFINITY);  // numerator in P(A,B,C,...) / (P(A)P(B)P(C)...)
    double bayes_factor(-INFINITY); // final result
    vector<double> single_scores(qry_seqs.n_seqs(), -INFINITY);  // NOTE log probs, not scores, but I haven't managed to finish switching over to the new terminology
    bool stop(false);
    string errors;
    do {
      errors = "";
      // clock_t run_start(clock());
      if(args.debug()) cout << "       ----" << endl;
      result = jh.Run(qry_seqs, kbounds);
      numerator = result.total_score();
      bayes_factor = numerator;
      if(args.algorithm() == "forward" && qry_seqs.n_seqs() > 1) {  // calculate factors for denominator
        for(size_t iseq = 0; iseq < qry_seqs.n_seqs(); ++iseq) {
          denom_results[iseq] = jh.Run(qry_seqs[iseq], kbounds);  // result for a single sequence
          single_scores[iseq] = denom_results[iseq].total_score();
          bayes_factor -= single_scores[iseq];
        }
      }

      kbounds = result.better_kbounds();
      for(auto & res : denom_results)
        kbounds = kbounds.LogicalOr(res.better_kbounds());

      stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
      for(auto & res : denom_results)
        stop &= !res.boundary_error() || res.could_not_expand();
      if(args.debug() && !stop)
        cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
      // cout << "      time " << ((clock() - run_start) / (double)CLOCKS_PER_SEC) << endl;
      if(result.boundary_error())
        errors = "boundary";
      for(auto & res : denom_results) // NOTE <errors> will still just have one "boundary" in it even if multiple results had boundary errors
        if(res.boundary_error())
          errors = "boundary";
    } while(!stop);

    // if(result.could_not_expand())
    //   cout << "WARNING " << qry_seqs.name_str() << " couldn't expand k bounds for " << kbounds.stringify() << endl;
    // if(single_result.boundary_error())
    //   cout << "WARNING boundary errors for " << qry_seqs[iseq].name() << " when together with " << qry_seqs.name_str() << endl;

    if(args.debug() && args.algorithm() == "forward" && qry_seqs.n_seqs() > 1)
      print_forward_scores(numerator, single_scores, bayes_factor);

    if(args.algorithm() == "viterbi" && size_t(args.n_best_events()) > result.events_.size()) {   // if we were asked for more events than we found
      if(result.events_.size() > 0)
        cout << "WARNING asked for " << args.n_best_events() << " events but only found " << result.events_.size() << endl;
      else
        assert(result.no_path_);  // if there's some *other* way we can end up with no events, I want to know about it
    }
    StreamOutput(ofs, args, result.events_, qry_seqs, bayes_factor, errors);
  }

  ofs.close();
  return 0;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score, string errors) {
  if(args.algorithm() == "viterbi") {
    size_t n_max = min(size_t(args.n_best_events()), events.size());
    for(size_t ievt = 0; ievt < n_max; ++ievt) {
      RecoEvent *event = &events[ievt];
      string second_seq_name, second_seq;
      ofs  // be very, very careful to change this *and* the csv header above at the same time
          << seqs.name_str(":")
          << "," << event->genes_["v"]
          << "," << event->genes_["d"]
          << "," << event->genes_["j"]
          << "," << event->insertions_["fv"]
          << "," << event->insertions_["vd"]
          << "," << event->insertions_["dj"]
          << "," << event->insertions_["jf"]
          << "," << event->deletions_["v_5p"]
          << "," << event->deletions_["v_3p"]
          << "," << event->deletions_["d_5p"]
          << "," << event->deletions_["d_3p"]
          << "," << event->deletions_["j_5p"]
          << "," << event->deletions_["j_3p"]
          << "," << event->score_
          << "," << seqs.seq_str(":")
          << "," << errors
          << endl;
    }
  } else {
    ofs
        << seqs.name_str(":")
        << "," << total_score
        << "," << errors
        << endl;
  }
}
// ----------------------------------------------------------------------------------------
void print_forward_scores(double numerator, vector<double> single_scores, double bayes_factor) {
  printf("   %8.3f = ", bayes_factor);
  printf("%2s %8.2f", "", numerator);
  for(auto & score : single_scores)
    printf(" - %8.2f", score);
  printf("\n");
}
// ----------------------------------------------------------------------------------------
int minimal_hamming_distance(Sequences &seqs_a, Sequences &seqs_b) {
  // Minimal hamming distance between any sequence in <seqs_a> and any sequence in <seqs_b>
  // NOTE for now, we require sequences are the same length (otherwise we have to deal with alignming them which is what we would call a CAN OF WORMS.
  assert(seqs_a.n_seqs() > 0 && seqs_b.n_seqs() > 0);
  int min_distance(seqs_a[0].size());
  for(size_t is=0; is<seqs_a.n_seqs(); ++is) {
    for(size_t js=0; js<seqs_b.n_seqs(); ++js) {
      Sequence seq_a = seqs_a[is];
      Sequence seq_b = seqs_b[js];
      assert(seq_a.size() == seq_b.size());
      int distance(0);
      for(size_t ic=0; ic<seq_a.size(); ++ic) {
	if(seq_a[ic] != seq_b[ic])
	  ++distance;
      }
      if(distance < min_distance)
	min_distance = distance;
    }
  }
  return min_distance;
}

// ----------------------------------------------------------------------------------------
vector<string> get_cluster_list(map<string, Sequences> &partinfo) {
  vector<string> clusters;
  for(auto &kv : partinfo)
    clusters.push_back(kv.first);
  return clusters;
}

// ----------------------------------------------------------------------------------------
double log_prob_of_partition(vector<string> &clusters, map<string, double> &log_probs) {
  // get log prob of entire partition given by the keys in <partinfo> using the individual log probs in <log_probs>
  double total_log_prob(0.0);
  for(auto &key : clusters) {
    if(log_probs.count(key) == 0)
      throw runtime_error("ERROR couldn't find key " + key + " in cached log probs\n");
    total_log_prob = AddWithMinusInfinities(total_log_prob, log_probs[key]);
  }
  return total_log_prob;
}

// ----------------------------------------------------------------------------------------
void print_partition(vector<string> &clusters, map<string, double> &log_probs, string extrastr) {
  const char *extra_cstr(extrastr.c_str());
  if(log_probs.size() == 0)
    printf("    %s partition\n", extra_cstr);
  else
    printf("    %-8.2f %s partition\n", log_prob_of_partition(clusters, log_probs), extra_cstr);
  for(auto &key : clusters)
    cout << "          " << key << endl;
}

// ----------------------------------------------------------------------------------------
void glomerate(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args &args, ofstream &ofs, map<string, Sequences> &info,
	       map<string, KBounds> &kbinfo, map<string, vector<string> > &only_genes, map<string, double> &cached_log_probs,
	       vector<double> &list_of_log_probs, vector<vector<string> > &list_of_partitions,
	       double &max_log_prob_of_partition, vector<string> &best_partition, bool &finished) {  // reminder: <qry_seq_list> is a list of lists

  double max_log_prob(-INFINITY);
  pair<string, string> max_pair; // pair of clusters with largest log prob (i.e. the ratio of their prob together to their prob apart is largest)
  // Sequences max_seqs;
  vector<string> max_only_genes;
  KBounds max_kbounds;

  set<string> already_done;  // keeps track of which  a-b pairs we've already done
  for(auto &kv_a : info) {  // note that c++ maps are ordered
    for(auto &kv_b : info) {
      if(kv_a.first == kv_b.first) continue;
      vector<string> names{kv_a.first, kv_b.first};
      sort(names.begin(), names.end());
      string bothnamestr(names[0] + ":" + names[1]);
      if(already_done.count(bothnamestr))  // already did this pair
	continue;
      else
	already_done.insert(bothnamestr);

      Sequences a_seqs(kv_a.second), b_seqs(kv_b.second);
      // TODO cache hamming fraction as well
      double hamming_fraction = float(minimal_hamming_distance(a_seqs, b_seqs)) / a_seqs[0].size();  // minimal_hamming_distance() will fail if the seqs aren't all the same length
      bool TMP_only_get_single_log_probs(false);
      if(hamming_fraction > args.hamming_fraction_cutoff()) {  // if all sequences in a are too far away from all sequences in b
	if(cached_log_probs.count(kv_a.first) == 0 || cached_log_probs.count(kv_b.first) == 0)
	  TMP_only_get_single_log_probs = true;
	else
	  continue;
      }

      // TODO skip all this stuff if we already have all three of 'em cached
      Sequences ab_seqs(a_seqs.Union(b_seqs));
      vector<string> ab_only_genes = only_genes[kv_a.first];
      for(auto &g : only_genes[kv_b.first])  // NOTE this will add duplicates (that's no big deal, though) OPTIMIZATION
	ab_only_genes.push_back(g);
      KBounds ab_kbounds = kbinfo[kv_a.first].LogicalOr(kbinfo[kv_b.first]);

      JobHolder jh(gl, hmms, args.algorithm(), ab_only_genes);  // NOTE it's an ok approximation to compare log probs for sequence sets that were run with different kbounds, but (I'm pretty sure) we do need to run them with the same set of genes. EDIT hm, well, maybe not. Anywa, it's ok for now
      // TODO make sure that using <ab_only_genes> doesn't introduce some bias
      jh.SetDebug(args.debug());
      jh.SetChunkCache(args.chunk_cache());
      jh.SetNBestEvents(args.n_best_events());

      string errors;
      // NOTE error from using the single kbounds rather than the OR seems to be around a part in a thousand or less
      get_result(args, jh, kv_a.first, a_seqs, kbinfo[kv_a.first], cached_log_probs, errors);
      get_result(args, jh, kv_b.first, b_seqs, kbinfo[kv_b.first], cached_log_probs, errors);
      if(TMP_only_get_single_log_probs)
	continue;
      get_result(args, jh, bothnamestr, ab_seqs, ab_kbounds, cached_log_probs, errors);

      // clock_t run_start(clock());
      // cout << "      time " << ((clock() - run_start) / (double)CLOCKS_PER_SEC) << endl;

      double bayes_factor(cached_log_probs[bothnamestr] - cached_log_probs[kv_a.first] - cached_log_probs[kv_b.first]);  // REMINDER a, b not necessarily same order as names[0], names[1]
      if(args.debug())
	print_forward_scores(cached_log_probs[bothnamestr], {cached_log_probs[kv_a.first], cached_log_probs[kv_b.first]}, bayes_factor);

      if(bayes_factor > max_log_prob) {
	max_log_prob = bayes_factor;
	max_pair = pair<string, string>(kv_a.first, kv_b.first);  // REMINDER not always same order as names[0], names[1]
	// max_seqs = ab_seqs;
	max_only_genes = ab_only_genes;
	max_kbounds = ab_kbounds;
      }
    }
  }

  // if <info> only has one cluster, if hamming is too large between all remaining clusters/remaining bayes factors are -INFINITY
  if(max_log_prob == -INFINITY) {
    finished = true;
    return;
  }

  // then merge the two best clusters
  vector<string> max_names{max_pair.first, max_pair.second};
  sort(max_names.begin(), max_names.end());
  Sequences max_seqs(info[max_names[0]].Union(info[max_names[1]]));  // NOTE this will give the ordering {<first seqs>, <second seqs>}, which should be the same as in <max_name_str>. But I don't think it'd hurt anything if the sequences and names were in a different order
  string max_name_str(max_names[0] + ":" + max_names[1]);  // NOTE the names[i] are not sorted *within* themselves, but <names> itself is sorted. This is ok, because we will never again encounter these sequences separately
  info[max_name_str] = max_seqs;
  kbinfo[max_name_str] = max_kbounds;
  only_genes[max_name_str] = max_only_genes;

  info.erase(max_pair.first);
  info.erase(max_pair.second);
  kbinfo.erase(max_pair.first);
  kbinfo.erase(max_pair.second);
  only_genes.erase(max_pair.first);
  only_genes.erase(max_pair.second);

  if(args.debug())
    printf("       merged %-8.2f %s and %s\n", max_log_prob, max_pair.first.c_str(), max_pair.second.c_str());

  vector<string> partition(get_cluster_list(info));
  if(args.debug())
    print_partition(partition, cached_log_probs, "current");
  double total_log_prob = log_prob_of_partition(partition, cached_log_probs);

  list_of_log_probs.push_back(total_log_prob);
  list_of_partitions.push_back(partition);

  if(total_log_prob > max_log_prob_of_partition) {
    best_partition = partition;
    max_log_prob_of_partition = total_log_prob;
  }

  if(max_log_prob_of_partition - total_log_prob > 1000.0) {  // stop if we've moved too far past the maximum
    cout << "    stopping after drop " << max_log_prob_of_partition << " --> " << total_log_prob << endl;
    finished = true;  // NOTE this will not play well with multiple maxima, but I'm pretty sure we shouldn't be getting those
  }
}

// ----------------------------------------------------------------------------------------
void write_cached_log_probs(ofstream &log_prob_ofs, map<string, double> &cached_log_probs) {
  log_prob_ofs << "unique_ids,score" << endl;
  for(auto &kv : cached_log_probs)
    log_prob_ofs << kv.first << "," << kv.second << endl;
}

// ----------------------------------------------------------------------------------------
void write_partition(ofstream &ofs, vector<string> partition, double log_prob) {
  for(size_t ic=0; ic<partition.size(); ++ic) {
    if(ic > 0)
      ofs << ";";
    ofs << partition[ic];
  }
  ofs << ",";
  ofs << log_prob << ",";
  ofs << "n/a\n";
}

// ----------------------------------------------------------------------------------------
map<string, double> read_cached_log_probs(string fname) {
  map<string, double> cached_log_probs;
  if(fname == "")
    return cached_log_probs;
  ifstream ifs(fname);
  assert(ifs.is_open());
  string line;

  // check the header is right TODO should write a general csv reader
  getline(ifs, line);
  vector<string> headstrs(SplitString(line, ","));
  cout << "x" << headstrs[0] << "x" << headstrs[1] << "x" << endl;
  assert(headstrs[0].find("unique_ids") == 0);
  assert(headstrs[1].find("score") == 0);

  while(getline(ifs, line)) {
    vector<string> column_list = SplitString(line, ",");
    assert(column_list.size() == 2);
    string unique_ids(column_list[0]);
    double logprob(stof(column_list[1]));
    cached_log_probs[unique_ids] = logprob;
  }

  return cached_log_probs;
}
// ----------------------------------------------------------------------------------------
void hierarch_agglom(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args &args, ofstream &ofs) {  // reminder: <qry_seq_list> is a list of lists
  // first convert the vector to a map 
  map<string, Sequences> info;
  map<string, vector<string> > only_genes;
  map<string, KBounds> kbinfo;
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    string key(qry_seq_list[iqry].name_str(":"));

    KSet kmin(args.integers_["k_v_min"][iqry], args.integers_["k_d_min"][iqry]);
    KSet kmax(args.integers_["k_v_max"][iqry], args.integers_["k_d_max"][iqry]);

    KBounds kb(kmin, kmax);
    info[key] = qry_seq_list[iqry];  // NOTE this is probably kind of inefficient to remake the Sequences all the time
    only_genes[key] = args.str_lists_["only_genes"][iqry];
    kbinfo[key] = kb;
  }

  vector<string> initial_partition(get_cluster_list(info));
  // then glomerate 'em
  map<string, double> cached_log_probs = read_cached_log_probs(args.incachefile());
  vector<double> list_of_log_probs;  // TODO I think I don't need this any more
  vector<vector<string> > list_of_partitions;  // TODO I think I don't need this any more
  double max_log_prob_of_partition(-INFINITY);  // 
  vector<string> best_partition;
  if(args.debug())
    print_partition(initial_partition, cached_log_probs, "initial");
  bool finished(false);
  do {
    glomerate(hmms, gl, qry_seq_list, args, ofs, info, kbinfo, only_genes, cached_log_probs, list_of_log_probs, list_of_partitions, max_log_prob_of_partition, best_partition, finished);
  } while(!finished);

  // assert(list_of_partitions.size() == max_log_probs.size());
  if(args.debug()) {
    cout << "-----------------" << endl;  
    print_partition(best_partition, cached_log_probs, "best");
  }

  // TODO oh, damn, if the initial partition is better than any subsequent ones this breaks
  // TODO it might be mroe efficient to write the partitions as I go so I don't have to keep them around in memory
  write_partition(ofs, initial_partition, log_prob_of_partition(initial_partition, cached_log_probs));
  for(size_t ip=0; ip<list_of_partitions.size(); ++ip)
    write_partition(ofs, list_of_partitions[ip], list_of_log_probs[ip]);

  // TODO really pass the cachefile as arg? maybe do something cleaner
  ofstream log_prob_ofs;
  log_prob_ofs.open(args.outcachefile());
  write_cached_log_probs(log_prob_ofs, cached_log_probs);
  log_prob_ofs.close();
}
