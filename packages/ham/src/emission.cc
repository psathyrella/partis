#include "emission.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Emission::Emission() {
  total_ = 0.0;
  tracks_ = NULL;
}

// ----------------------------------------------------------------------------------------
Emission::~Emission() {
  if(tracks_) delete tracks_;
}

// ----------------------------------------------------------------------------------------
void Emission::Parse(YAML::Node config, Tracks model_tracks) {
  scores_.Init();  // TODO doesn't seem to do anything?
  tracks_ = new vector<Track*>();  // list of the tracks used by *this* emission. Note that this may not be all the tracks used in the model.
  Track *tk(model_tracks.track(config["track"].as<string>()));
  assert(tk);  // assures we actualy found the track in model_tracks
  tracks_->push_back(tk);
  scores_.AddTrack(tk, 0);  // TODO what's with the zero yo?

  YAML::Node probs(config["probs"]);
  if(probs.size() != scores_.alphabet_size(0))
    throw runtime_error("ERROR emission probabilities (" + to_string(probs.size()) + ") not the same length as the alphavet passed to lexical table constructor (" + to_string(scores_.alphabet_size(0)) + "). Yes, a v! Hobgoblins and smallness, yo.");
  if(scores_.alphabet_size(0) != tk->alphabet_size())  // TODO what's with the zero?
    throw runtime_error("ERROR alphabet passed to lexical table (" + to_string(scores_.alphabet_size(0)) + ") not the same size as the track's alphabet (" + to_string(tk->alphabet_size()) + ").");
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
}

// ----------------------------------------------------------------------------------------
void Emission::Print() {
  for(size_t i = 0; i < scores_.n_tracks(); ++i)
    cout << "    " << scores_.track(i)->name();
  cout << "     (normed to within at least " << EPS << ")" << endl;

  printf("%16s", "");
  // print name of each emission
  for(size_t ir = 0; ir < (*tracks_)[0]->alphabet_size(); ++ir)
    printf("%12s", (*tracks_)[0]->symbol(ir).c_str());
  // printf("        (total - 1.0 = %.5e - 1.0 = %.5e\n", total_, total_ - 1.0);
  cout << endl;

  assert(tracks_->size() == 1);
  printf("%16s", "");
  for(size_t ir = 0; ir < (*tracks_)[0]->alphabet_size(); ++ir) {
    double prob = exp(scores_.LogProb(ir));
    if(prob < 0.01)
      printf("%22.9e", prob);
    else
      printf("%22.9f", prob);
  }
  cout << "\n";
}

}
