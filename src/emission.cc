#include "emission.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Emission::Emission() {
  total_ = 0.0;
  tracks_ = NULL;
  pair_ = false;
}

// ----------------------------------------------------------------------------------------
Emission::~Emission() {
  if(tracks_) delete tracks_;
}

// ----------------------------------------------------------------------------------------
void Emission::Parse(YAML::Node config, string is_pair, Tracks model_tracks) {
  scores_.Init();
  // NOTE at this point we only allow one track per emission (in particular, we require that pair emissions be on the same track). kinda TODO This'd be easy to change later, of course
  tracks_ = new vector<Track*>();  // list of the tracks used by *this* emission. Note that this may not be all the tracks used in the model.
  if(is_pair == "single") {
    pair_ = false;
    Track *tk(model_tracks.track(config["track"].as<string>()));
    assert(tk);  // assures we actualy found the track in model_tracks
    tracks_->push_back(tk);
    scores_.AddTrack(tk, 0);

    YAML::Node probs(config["probs"]);
    assert(probs.size() == scores_.alphabet_size(0));  // TODO actually I don't need these either, since I'm looping over the track
    assert(scores_.alphabet_size(0) == tk->alphabet_size()); // TODO arg I shouldn't need this. so complicated...
    vector<double> log_probs;
    total_ = 0.0; // make sure things add to 1.0
    for(size_t ip = 0; ip < scores_.alphabet_size(0); ++ip) {
      double prob(probs[tk->symbol(ip)].as<double>());  // NOTE probs are stored as dicts in the file, so <probs> is unordered
      log_probs.push_back(log(prob));
      total_ += prob;
    }
    if(fabs(total_ - 1.0) >= EPS) { // make sure emissions probs sum to 1.0
      cerr << "ERROR normalization failed for emissions:" << endl;
      cerr << config << endl;
      throw runtime_error("configuration");
    }
    scores_.AddColumn(log_probs);  // NOTE <log_probs> must already be logged
  } else if(is_pair == "pair") {
    pair_ = true;
    assert(config["tracks"].size() == 2);
    for(size_t it = 0; it < config["tracks"].size(); ++it) {
      Track *tk(model_tracks.track(config["tracks"][it].as<string>()));
      assert(tk);  // assures we actualy found the track in model_tracks
      tracks_->push_back(tk);
      scores_.AddTrack(tk, 0);
    }

    YAML::Node probs(config["probs"]);
    assert(probs.size() == scores_.alphabet_size(0));  // TODO actually I don't need these either, since I'm looping over the track
    assert(tracks_->size() == 2);
    assert(scores_.alphabet_size(0) == (*tracks_)[0]->alphabet_size()); // TODO arg I shouldn't need this. so complicated...
    assert(scores_.alphabet_size(0) == (*tracks_)[1]->alphabet_size()); // TODO arg I shouldn't need this. so complicated...
    total_ = 0.0; // make sure things add to 1.0
    for(size_t ip = 0; ip < scores_.alphabet_size(0); ++ip) {
      YAML::Node these_probs(config["probs"][(*tracks_)[0]->symbol(ip)]);
      assert(these_probs.size() == scores_.alphabet_size(0));
      vector<double> log_probs;
      for(size_t ipp = 0; ipp < scores_.alphabet_size(0); ++ipp) {
        double prob(these_probs[(*tracks_)[1]->symbol(ipp)].as<double>());  // NOTE probs are stored as dicts in the file, so <probs> is unordered
        log_probs.push_back(log(prob));
        total_ += prob;
      }
      scores_.AddColumn(log_probs);  // NOTE <log_probs> must already be logged. also NOTE that a column in <scores_> is maybe a row in the yaml file. I didn't choose it!
    }
    // TODO use something cleverer than a random hard coded EPS
    if(fabs(total_ - 1.0) >= EPS) { // make sure emissions probs sum to 1.0
      cerr << "ERROR normalization failed for" << endl;
      cerr << config << endl;
      throw runtime_error("configuration");
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------------------------
void Emission::Print() {
  for(size_t i = 0; i < scores_.n_tracks(); ++i) // TODO dammit I'm confused by having a vector of tracks in the lexical table *and* every freaking where else
    cout << "    " << scores_.track(i)->name();
  cout << "     (normed to within at least " << EPS << ")" << endl;

  printf("%16s", "");
  // print name of each emission
  for(size_t ir = 0; ir < (*tracks_)[0]->alphabet_size(); ++ir)
    printf("%12s", (*tracks_)[0]->symbol(ir).c_str());
  // printf("        (total - 1.0 = %.5e - 1.0 = %.5e\n", total_, total_ - 1.0);
  cout << endl;

  if(!pair_) {
    assert(tracks_->size() == 1);
    printf("%16s", "");
    for(size_t ir = 0; ir < (*tracks_)[0]->alphabet_size(); ++ir) {
      double prob = exp(scores_.LogProb(ir));
      if(prob < 0.01)
        printf("%12.3e", exp(scores_.LogProb(ir)));
      else
        printf("%12.3f", exp(scores_.LogProb(ir)));
    }
    cout << "\n";
  } else {
    assert(tracks_->size() == 2);
    for(size_t ir = 0; ir < (*tracks_)[0]->alphabet_size(); ++ir) {
      printf("%16s", (*tracks_)[0]->symbol(ir).c_str());
      for(size_t ic = 0; ic < (*tracks_)[1]->alphabet_size(); ++ic)
        printf("%12.3e", exp(scores_.LogProb(ir, ic)));
      cout << "\n";
    }
  }
}

}
