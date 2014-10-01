#include "emm.h"

namespace StochHMM{

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
bool emm::parse(string& txt, tracks& model_tracks) {
  stringList lines;  // a list of the lines in <txt>, which specifies the entire emission
  lines.splitString(txt,"\n");
  size_t idx;  // set <idx> to start of emission definition
  if (lines.contains("EMISSION")) {
    idx = lines.indexOf("EMISSION");
  } else {
    cerr << "Missing EMISSION tag from emission. Please check the formatting.   This is what was handed to the emission class:\n " <<  txt << endl;
    assert(0);
  }
      
  // Determine Emission Type and set appropriate flags
  stringList line;
  line.splitString(lines[idx], "\t,: ");  // split lines[idx] (idx probably zero) using any of these delimineters, and store in <line>
  size_t typeBegin(0);
  valueType valtyp(PROBABILITY);
  if (line.contains("P(X)")) {
    typeBegin = line.indexOf("P(X)");
    valtyp = PROBABILITY;
  } else if (line.contains("LOG")) {
    typeBegin = line.indexOf("LOG");
    valtyp = LOG_PROB;
  } else if (line.contains("COUNTS")) {
    typeBegin = line.indexOf("COUNTS");
    valtyp = COUNTS ;
  } else {
    string info = "Couldn't parse Value type in the Emission: " + txt  + " Please check the formatting.   The allowed types are: P(X), LOG, COUNTS, or REAL_NUMBER. \n";
    cerr << info << endl;
  }

  pair_ = line.contains("PAIR");

  scores.init();
  // push back tracks (should only be one a.t.m.)
  tracks_ = new vector<track*>();  // list of the tracks used by *this* emission. Note that this may not be all the tracks used in the model.
  for (size_t i=1; i<typeBegin; i++) {  // loop over the words in <line> from 1 to the start of the <valtype> specification (i.e. over what should be a list of all tracks for this emission)
    track* tk(model_tracks.getTrack(line[i]));
    if (tk) {
      tracks_->push_back(tk);
      scores.addTrack(tk, 0);
    } else {
      cerr << "Emissions tried to add a track named: " << line[i] << " . However, there isn't a matching track in the model.  Please check model formatting.\n";
      assert(0);
    }
  }
      
  size_t expectedColumns(scores.getAlphaSize(0));
  size_t expectedRows(pair_ ? scores.getAlphaSize(0) : 1);

  // check emission labels are written correctly
  stringstream ss(lines[1]);
  string letter;
  ss >> letter;
  if (letter != "@") {
    cerr << "ERROR in emission label line " << lines[1] << endl;
    cerr << "    " << letter << " != @" << endl;
    assert(0);
  }
  for (size_t in=0; in<expectedColumns; ++in) {
    ss >> letter;
    assert(letter == (*tracks_)[0]->getAlpha(in));
  }

  size_t n_lines_to_skip(2);  // NOTE skip header lines
  for (size_t il=n_lines_to_skip; il<lines.size(); il++) {
    line.splitString(lines[il],"\t ");  // reset <line> to a list consisting of lines[il] split by white space

    if (pair_) {
      if (line[0] != (*tracks_)[0]->getAlpha(il - n_lines_to_skip)) {
	cerr << "ERROR beginning of emission line does not match expected letter: " << line[0] << " != " << (*tracks_)[0]->getAlpha(il - n_lines_to_skip) << endl;
	assert(0);
      }
      line.pop_ith(0);
    }

    vector<double> tmp_vec = line.toVecDouble();
    if (tmp_vec.size() != expectedColumns) {
      cerr << "The following line with " << tmp_vec.size() << " columns couldn't be parsed into the required number of columns (" + int_to_string(expectedColumns) + ")\n"
	   << lines[il] << endl;
      assert(0);
    }
    if (valtyp == PROBABILITY) {
      logVector(tmp_vec);
    } else if (valtyp == COUNTS) {
      probVector(tmp_vec);
      logVector(tmp_vec);
    }
    scores.AddColumn(tmp_vec);
  }
          
  if (scores.getLogProbabilityTable()->size() != expectedRows) {
    cerr << " The Emission table doesn't contain enough rows.  Found " << scores.getLogProbabilityTable()->size() << " but expected " << expectedRows << " \n Please check the Emission Table and formatting for " <<  txt << endl;
    assert(0);
  }
                      
  return true;
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
