#include "emm.h"

namespace stochhmm {

// ----------------------------------------------------------------------------------------
emm::emm() {
  tracks_ = NULL;
  track_indices = NULL;
  pair_ = false;
}

// ----------------------------------------------------------------------------------------
emm::~emm(){
  if (tracks_) delete tracks_;
}
  
// ----------------------------------------------------------------------------------------
bool emm::parse(YAML::Node config, string is_pair, tracks model_tracks) {
  pair_ = is_pair;

  scores.init();
  // push back tracks (should only be one a.t.m.)
  // NOTE at this point we only allow one track per emission (in particular, we require that pair emissions be on the same track). kinda TODO This'd be easy to change later, of course
  tracks_ = new vector<track*>();  // list of the tracks used by *this* emission. Note that this may not be all the tracks used in the model.
  if (is_pair=="single") {
    track *tk(model_tracks.getTrack(config["track"]));
    assert(tk);  // assures we actualy found the track in model_tracks
    tracks_->push_back(tk);
    scores.addTrack(tk, 0);
    
    YAML::Node probs(config["probs"]);
    assert(probs.size() == scores.getAlphaSize(0));  // TODO actually I don't need these either, since I'm looping over the track
    assert(scores.getAlphaSize(0) == tk->getAlphaSize()); // TODO arg I shouldn't need this. so complicated...
    vector<double> log_probs;
    for (size_t ip=0; ip<scores.getAlphaSize(0); ++ip) {
      double prob(probs[tk->getAlpha(ip)]);  // NOTE probs are stored as dicts in the file, so <probs> is unordered
      log_probs.push_back(log(prob));
    }
    scores.AddColumn(tmp_vec);  // NOTE <tmp_vec> must already be logged
  } else if (is_pair=="pair") {
    assert(config["tracks"].size() == 2);
    for (size_t it=0; it<config["tracks"].size(); ++it) {
      track *tk(model_tracks.getTrack(config["tracks"][it]));
      assert(tk);  // assures we actualy found the track in model_tracks
      tracks_->push_back(tk);
      scores.addTrack(tk, 0);
    }
    
    YAML::Node probs(config["probs"]);
    assert(probs.size() == scores.getAlphaSize(0));  // TODO actually I don't need these either, since I'm looping over the track
    assert(scores.getAlphaSize(0) == tk->getAlphaSize()); // TODO arg I shouldn't need this. so complicated...
    for (size_t ip=0; ip<scores.getAlphaSize(0); ++ip) {
      YAML::Node these_probs(config["probs"][tk->getAlpha(ip)]);
      assert(these_probs.size() == scores.getAlphaSize(0));
      vector<double> log_probs;
      for (size_t ipp=0; ipp<scores.getAlphaSize(0); ++ipp) {
	double prob(these_probs[tk->getAlpha(ipp)]);  // NOTE probs are stored as dicts in the file, so <probs> is unordered
	log_probs.push_back(log(prob));
      }
      scores.AddColumn(tmp_vec);  // NOTE <tmp_vec> must already be logged. also NOTE that a column in <scores> is maybe a row in the yaml file. I didn't choose it!
    }
  } else {
    assert(0);
  }
}
      
// // ----------------------------------------------------------------------------------------
// double emm::get_emission(sequences& seqs,size_t pos) {
//   return scores.getValue(seqs, pos);
// }
      
// ----------------------------------------------------------------------------------------
string emm::stringify() {
  string emissionString("EMISSION:\t");
  for(size_t i=0;i<scores.getNTracks();i++){
    if (i>0){
      emissionString+=",";
    }
    emissionString+=scores.getTrack(i)->getName();
  }
          
  emissionString+=":\t";
          
  emissionString+="LOG";
          
  emissionString+="\n\tORDER:\t";
          
  for(size_t i=0;i<scores.getNTracks();i++){
    if (i>0){
      emissionString+=",";
    }
    emissionString+=int_to_string(0// scores.getOrder(i)
				  );
  }
          
  emissionString+="\n";

  emissionString+="FOOP";//scores.stringify();
  emissionString+="\n";
      
  return emissionString;
}
}
