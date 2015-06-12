#include "emission.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Emission::Emission() {
  total_ = 0.0;
  track_ = nullptr;
}

// ----------------------------------------------------------------------------------------
Emission::~Emission() {
  // track_ is owned by the model
}

// ----------------------------------------------------------------------------------------
void Emission::Parse(YAML::Node config, Track *track) {
  track_ = track;
  scores_.Init(track_);

  YAML::Node probs(config["probs"]);
  if(probs.size() != track_->alphabet_size())
    throw runtime_error("ERROR emission probabilities (" + to_string(probs.size()) + ") not the same length as the alphavet passed to lexical table constructor (" + to_string(track_->alphabet_size()) + "). Yes, a v! Hobgoblins and smallness, yo.");
  if(scores_.track()->alphabet_size() != track_->alphabet_size())
    throw runtime_error("ERROR alphabet passed to lexical table (" + to_string(scores_.track()->alphabet_size()) + ") not the same size as the track's alphabet (" + to_string(track_->alphabet_size()) + ").");
  vector<double> log_probs;
  total_ = 0.0; // make sure things add to 1.0
  for(size_t ip = 0; ip < track_->alphabet_size(); ++ip) {
    double prob(probs[track_->symbol(ip)].as<double>());  // NOTE probs are stored as dicts in the file, so <probs> is unordered
    log_probs.push_back(log(prob));
    total_ += prob;
  }
  if(fabs(total_ - 1.0) >= EPS) { // make sure emissions probs sum to 1.0
    cerr << "ERROR normalization failed for emissions:" << endl;
    cerr << config << endl;
    throw runtime_error("configuration");
  }
  scores_.SetLogProbs(log_probs);  // NOTE <log_probs> must already be logged
}

// ----------------------------------------------------------------------------------------
void Emission::Print() {
  cout << "    " << track_->name() << "     (normed to within at least " << EPS << ")" << endl;

  printf("%16s", "");
  // print name of each emission
  for(size_t ir = 0; ir < track_->alphabet_size(); ++ir)
    printf("%12s", track_->symbol(ir).c_str());
  // printf("        (total - 1.0 = %.5e - 1.0 = %.5e\n", total_, total_ - 1.0);
  cout << endl;

  printf("%16s", "");
  for(size_t ir = 0; ir < track_->alphabet_size(); ++ir) {
    double prob = exp(scores_.LogProb(ir));
    if(prob < 0.01)
      printf("%22.9e", prob);
    else
      printf("%22.9f", prob);
  }
  cout << "\n";
}

}
