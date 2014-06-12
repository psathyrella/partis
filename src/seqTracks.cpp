//
//  seqTracks.cpp
//Copyright (c) 2007-2012 Paul C Lott 
//University of California, Davis
//Genome and Biomedical Sciences Facility
//UC Davis Genome Center
//Ian Korf Lab
//Website: www.korflab.ucdavis.edu
//Email: lottpaul@gmail.com
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//the Software, and to permit persons to whom the Software is furnished to do so,
//subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



#include "seqTracks.h"
namespace StochHMM{

  //!Create an empty initialized seqTrack
  seqTracks::seqTracks(){
    _init();
  }
    
  //!Create, initialize, and start loading sequence file
  //!\param mdl Single model
  //!\param filename  Sequence File filename
  //!\param format  Sequence file Format (FASTA or FASTQ)
  //!\param trFuncs TrackFunc* functions to create tracks based on imported sequences
  seqTracks::seqTracks(model& mdl , std::string& filename , SeqFileFormat format, TrackFuncs* trFuncs){
    _init();
    loadSeqs(mdl, filename,format, trFuncs);
  }
    
  //!Destroy seqTracks
  seqTracks::~seqTracks(){
    for(size_t i=0;i<filehandles.size();i++){
      if (filehandles[i]!=NULL){
	filehandles[i]->close();
      }
    }
        
    hmms=NULL;
    trackFunctions = NULL;
    attribModelFunc = NULL;
        
    while(!jobQueue.empty()){
      seqJob* element=jobQueue.front();
      jobQueue.pop();
      delete element;
      element = NULL;
    }
    return;
  }
  
  // ----------------------------------------------------------------------------------------
  void seqTracks::_init(){
    numImportJobs=1;
    n_jobs=0;
    hmms    = NULL;
    hmm     = NULL;
    trackFunctions   = NULL;
    attribModelFunc  = NULL;
    seqFormat = FASTA;  //Set default file format to fasta
    fileType = SINGLE_TRACK;
    good=false;
        
    //Make seqTrack thread-safe
#ifdef THREADS
    pthread_mutex_init(&exit_thread_flag_mutex, NULL);
    pthread_cond_init(&exit_thread_flag_cv, NULL);
#endif
    exit_thread=1;
  }
    
  //Reset the Queue, Files, and Import Tracks
  //Does not reset track Functions or Attribute model selection functions 
  void seqTracks::_reset(){
    for(size_t i=0;i<filehandles.size();i++){
      delete filehandles[i];
    }
    filehandles.clear();
    seqFilenames.clear();
    hmms    = NULL;
    hmm     = NULL;
    seqFormat = FASTA;
    fileType  = SINGLE_TRACK;
    good = false;
    exit_thread = 1;
    importTracks.clear();
  }

  // ----------------------------------------------------------------------------------------
  //! Load the fasta sequence file
  //! \param mod  Model to be used
  //! \param seqFile  Sequence filename
  void seqTracks::loadSeqs(model& mod, std::string& seqFile, SeqFileFormat format, TrackFuncs* trFuncs){
    if (filehandles.size()>0) _reset();
        
    hmm = &mod;
    seqFormat = format;
    seqFilenames.push_back(seqFile);
    if (trFuncs) trackFunctions = trFuncs;
        
    //Get State Information and Determine # of tracks to import
    info = mod.getStateInfo();
    _initImportTrackInfo();
    fileType = (importTracks.size()>1) ? MULTI_TRACK : SINGLE_TRACK;
    _open();
    while (good) {  // while we're not at eof
      getNext();
    }
  }
    
  // ----------------------------------------------------------------------------------------
  //!Get the next job from job queue
  seqJob* seqTracks::getJob() {
    seqJob *jb(0);
    if (n_jobs>0) {
      jb = jobQueue.front();
      jobQueue.pop();
      n_jobs--;
    }
    return jb;
  }
    
  // ----------------------------------------------------------------------------------------
  // TODO: fix PT2TRACKFUNC function assignment
  bool seqTracks::_initImportTrackInfo() {
    model* tmp_model(0);
        
    if (hmms!=NULL) {  // if vector of models was initialized
      tmp_model = (*hmms)[0];
      assert(tmp_model);
    } else if (hmm!=NULL) {  // else if single model was initialized
      tmp_model = hmm;
    } else {
      std::cerr << "seqTracks initialization error:  Model is not defined. Therefore, can't initiate seqTrack with necessary model inforamation" << std::endl;
      return false;
    }
        
    modelTracks = tmp_model->getTracks();
        
    //Determine which tracks to import and which to get by using track functions
    track* tmp_track;
    trackCount = tmp_model->track_size();
    for(size_t i=0; i<trackCount; i++) {
      tmp_track = tmp_model->getTrack(i);
      if (tmp_track->isTrackFuncDefined()) {  // get this track using track functions
	ppTrack tmp;
	tmp.func = NULL;
	tmp.trackNumber=i;
	tmp.trackToUse=tmp_model->getTrackIter(tmp_track->getTrackToUse());
	std::string functionTouse=tmp_track->getTrackFunction();
	//tmp.func = funcs->getFunction(functionTouse);
	postprocessTracks.push_back(tmp);
      } else {  // import this track from the model
	importTracks.push_back(std::make_pair(i, tmp_track->getAlphaType()));
      }
    }
    return true;
  }

  // ----------------------------------------------------------------------------------------
  //!Print the seqTracks to stdout
  void seqTracks::print(){
    for (size_t i=0;i<jobQueue.size();i++){
      //jobQueue[i]->print_seq();
    }
    return;
  }

  // ----------------------------------------------------------------------------------------
  //! Read one job from the input file, i.e. add one set of sequences and the model to a new seqJob and push the
  //! job on the queue.
  bool seqTracks::getNext(){
    // create a new seqJob and add a sequence to it for each track defined in the model
    seqJob* temp_job = new(std::nothrow) seqJob(trackCount);
    assert(temp_job); // a lot shorter, innit?
    bool valid=true;  // set to false if sequence import failed for some reason
    sequence* sq;
    for(size_t i=0;i<importTracks.size();i++){  // push back a sequence (from the input file) for each track which was defined in the model
      sq = new(std::nothrow) sequence(false);
      assert(sq);
      bool success;
      size_t ifilehandle(0);  // a.t.m. we're limiting to one file at a time
      if (importTracks[i].second == REAL) {
	success = sq->getReal(*filehandles[ifilehandle], (*modelTracks)[importTracks[i].first]);
      } else if (seqFormat == FASTA) { // AlphaNum and Fasta
	success = sq->getFasta(*filehandles[ifilehandle], (*modelTracks)[importTracks[i].first], info);
      } else {
	success = sq->getFastq(*filehandles[ifilehandle], (*modelTracks)[importTracks[i].first]);
      }
      assert(success && sq);
      good = filehandles[ifilehandle]->good(); 

      if (sq->exDefDefined()) {  //If exDef is defined in sequence put it in sequences
	temp_job->set->setExDef(sq->getExDef());
      }
      temp_job->set->addSeq(sq, importTracks[i].first);
      temp_job->setSeqFilename(seqFilenames[ifilehandle]);
    }
        
    // get sequences defined by external function
    for (size_t i =0; i<postprocessTracks.size(); i++) {
      std::vector<double>* rl = NULL;
      if (postprocessTracks[i].func != NULL ){
	rl = (*postprocessTracks[i].func)(temp_job->set->getUndigitized(postprocessTracks[i].trackToUse));
      } else {
	std::cerr << "Sequence external function not defined for track number: " << postprocessTracks[i].trackNumber << std::endl;
	std::cerr << "Using Sequences from track: " << postprocessTracks[i].trackToUse << std::endl;
	rl = new(std::nothrow) std::vector<double>;
	assert(rl);
      }
      sequence* sq = new(std::nothrow) sequence(rl , (*modelTracks)[postprocessTracks[i].trackNumber]);
      assert(sq);
      temp_job->set->addSeq(sq,postprocessTracks[i].trackNumber);
    }

    temp_job->hmm = hmm;
            
    //Check that all sequences are same length.
    size_t lengthOfAll=SIZE_MAX;
    for (size_t i=0; i<trackCount; i++) {
      size_t length = temp_job->set->getLength(i);
      if (lengthOfAll==SIZE_MAX) {
	lengthOfAll = length;
      } else if (lengthOfAll!=length) {
	std::cerr << "Sequence Lengths not the same" <<std::endl;
	delete temp_job;
	return false;
      } else {
	continue;
      }
    }
    temp_job->set->setLength(lengthOfAll);

    // all is well! actually add the job to the queue
    jobQueue.push(temp_job);
    n_jobs++;

    return true;
  }
    
  // ----------------------------------------------------------------------------------------
  bool seqTracks::_open() {
    assert(seqFilenames.size() > 0);
    assert(importTracks.size() > 0);
    std::ifstream *SEQ = new(std::nothrow) std::ifstream;
    assert(SEQ);
    filehandles.push_back(SEQ);
    filehandles[0]->open(seqFilenames[0].c_str());
    assert(filehandles[0]->is_open());
    assert(filehandles[0]->good());
    good = true;
    return true;
  }
                   
  // ----------------------------------------------------------------------------------------
  bool seqTracks::_close(){
    filehandles[0]->clear();
    delete filehandles[0];
    filehandles.erase(filehandles.begin());
    seqFilenames.erase(seqFilenames.begin());
    return true;
  }
                    
}
