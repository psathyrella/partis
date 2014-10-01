#include "sequence.h"

namespace stochhmm {

// ----------------------------------------------------------------------------------------
Sequence::Sequence(string& seq_str, track* trk, string name):
  name_(name),
  undigitized_(seq_str),
  track_(trk)
{
  seq_ = new(nothrow) vector<uint8_t>;
  _digitize();
}
  
// ----------------------------------------------------------------------------------------
Sequence::Sequence(const Sequence& rhs) {
  name_ = rhs.name_;
  header_  = rhs.header_;
  track_  = rhs.track_;
  undigitized_ = rhs.undigitized_;
  seq_ = rhs.seq_;
  if (rhs.seq_)
    seq_ = new vector<uint8_t>(*rhs.seq_);
}
      
// ----------------------------------------------------------------------------------------
Sequence::~Sequence() {
  delete seq_;
}
  
// ----------------------------------------------------------------------------------------
void Sequence::clear(){
  header_ = "";
  undigitized_="";
  if (seq_)
    seq_->clear();
  track_ = NULL;
}
  
// ----------------------------------------------------------------------------------------
string* Sequence::getUndigitized() {
  if (!undigitized_.empty() || seq_->empty()) {
    return &undigitized_;
  } else {
    undigitized_ = undigitize();
    return &undigitized_;
  }
}

// ----------------------------------------------------------------------------------------
//! \return a copy of the sequence from <pos> of size <len>
Sequence Sequence::getSubSequence(size_t pos, size_t len) {
  assert(pos < undigitized_.size());
  assert(len < undigitized_.size());
  assert(pos+len <= undigitized_.size());  // arg. if pos+len overflows this still passes
  string subseq_str = undigitized_.substr(pos, len);
  return Sequence(subseq_str, track_, name_);
}

// ----------------------------------------------------------------------------------------
string Sequence::stringify(){
  string output;
  if (!header_.empty())
    output += header_ + "\n";
  output += undigitized_ + "\n";
  return output;
}
      
// ----------------------------------------------------------------------------------------
string Sequence::stringifyWOHeader() {
  return undigitized_ + "\n";
}
  
// ----------------------------------------------------------------------------------------
string Sequence::undigitize() {
  if (!seq_){  //If the sequence is not digitized yet.  Return the undigitized version 
    return undigitized_;
  }
              
  assert(track_);
  string output;
  size_t alphaMax = track_->getAlphaMax();
  for (size_t i=0; i<size(); i++) {
    output+=track_->getAlpha((*seq_)[i]);
    if (alphaMax!=1) {
      output+=" ";
    }
  }
  return output;
}
  
// ----------------------------------------------------------------------------------------
bool Sequence::getFasta(ifstream& file, track* trk) {
  if (seq_) this->clear();
  track_ = trk;
  assert(file.good());

  //Find next header mark 
  while (file.peek() != '>') {
    string temp;
    getline(file, temp, '\n');
    assert(file.good());  // Sequence doesn't contain a header
  }
      
  getline(file, header_, '\n');  // get header line
  // make sure the track name in the fasta header matches that in the track we've been passed
  stringstream ss(header_);
  string name,tmp_track_name;
  ss >> name >> tmp_track_name;
  assert(name != ">");  // make sure there wasn't a space after the '>'. yeah, I know, I should rewrite this
  name_ = name.substr(1);  // first char is the '>' marking this as a header line
  assert(tmp_track_name == trk->getName());
  
  bool success(false);
  //get sequence
  string line;
  while(getline(file, line, '\n')){
    undigitized_ += line;

    char nl_peek = file.peek();
    if (nl_peek=='>') {  // if there's a new sequence header on the next line, digitize what we have then break
      success = _digitize();
      break;
    } else if (nl_peek==EOF) {  // if we're at the end of the file, digitize and break
      getline(file, line, '\n');
      success = _digitize();
      break;
    } else {
      continue;
    }
  }

  assert(success);
  return success;
}
  
// ----------------------------------------------------------------------------------------
bool Sequence::_digitize() {
  assert(track_);
      
  stringList lst;
  clear_whitespace(undigitized_,"\n");
  if (!seq_) {
    seq_ = new vector<uint8_t>(undigitized_.size());
  } else {
    seq_->assign(undigitized_.size(),0);
  }
  for (size_t i=0; i<undigitized_.size(); i++) {
    uint8_t symbl = track_->symbolIndex(undigitized_[i]);
    (*seq_)[i] = symbl;
  }
  return true;
}
  
// ----------------------------------------------------------------------------------------
//! Import one fastq entry from the file
//!FastQ format:
//!Line 1: Start with @
//!Line 2: Sequence  , Can be multiple lines
//!Line 3: Start with +
//!Line 4: Quality Score , Can be multiple lines
//! \param file File stream to file
//! \param trk Track to used to digitize
//! \return true if the sequence was successfully imported
bool Sequence::getFastq(ifstream& file, track* trk){
              
  if (seq_!=NULL){
    this->clear();
  }
              
  track_=trk;
          
  if (!file.good()){
    return false;
  }
      
  //Find first line that starts with "@" and Get header
      
  //Move down until the next line has a "@"
  while(file.peek() != '@' && file.good()){
    string temp;
    getline(file,temp,'\n');
  }
      
  //Get Header (One line)
  if (file.good()){
    getline(file,header_,'\n');
  }
  else{
    return false;
  }
      
      
  string sequence="";
  //Get sequence (Multiple Lines)
  if (file.good()){
    while(getline(file,sequence,'\n')){
      undigitized_+=sequence;

      char nl_peek=file.peek();  // see if we have new sequence header on the next line
      if (nl_peek=='+'){  //If next line begins with + then we are at a 
        break;
      }
      else if (nl_peek==EOF){
        break;
      }
      else{
        continue;
      }
    }
  }
  else{
    return false;
  }
      
  _digitize();
      
  string quality_string;
      
  //Get Quality String (Multiple Lines)
  if (file.good()){
    while(getline(file,sequence,'\n')){
      quality_string+=sequence;
              
      char nl_peek=file.peek();  // see if we have new sequence header on the next line
      if (nl_peek=='@'){  //If next line begins with + then we are at a
        if (quality_string.size() < undigitized_.size()){
          continue;
        }
                  
        break;
      }
      else if (nl_peek==EOF){
        break;
      }
      else{
        continue;
      }
    }
  }
  else{
    return false;
  }
      
  // sequence_length_=seq->size();
  return true;
}
  
//! \return subsequence from (and including) istart to (but excluding) istop
vector<uint8_t> Sequence::getDigitalSubSeq(size_t istart, size_t istop) {
  vector<uint8_t> subseq;
  for (vector<uint8_t>::iterator it = seq_->begin()+istart; it != seq_->begin()+istop; ++it)
    subseq.push_back(*it);
  return subseq;
}
}
