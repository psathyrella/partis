//
//  Lexical.cpp
//  StochHMM
//
//  Created by Paul Lott on 4/2/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include "lexicalTable.h"



namespace StochHMM{
    
  lexicalTable::lexicalTable(){
    max_order=0;
    logProb=NULL;
    counts = NULL;
    prob = NULL;
    x_subarray=NULL;
    y_subarray=NULL;
    return;
  }
    
  lexicalTable::~lexicalTable(){
    delete logProb;
    delete prob;
    delete counts;
    delete x_subarray;
    delete y_subarray;
    logProb=NULL;
    prob=NULL;
    counts=NULL;
    x_subarray=NULL;
    y_subarray=NULL;
  }
    
  // ----------------------------------------------------------------------------------------
  void lexicalTable::initialize_emission_table(){
    if (logProb == NULL){
      std::cerr << "Cannot initialize emission table until after the tables have been assigned";
      exit(2);
    }

    assert(tracks.size() == tracks.size());  // not sure why we need both of these
    assert(tracks.size() == 1 || tracks.size() == 2);  // for the moment I limit to *two* dimensions (tracks), since there doesn't seem to exist a good library for n-dim matrices with n unknown until runtime
    assert(logProb->size() == 1);  // requiring order equal zero a.t.m.
    log_prob_matrix.resize(tracks[0]->getAlphaSize(), (tracks.size()==1) ? 1 : tracks[1]->getAlphaSize());
    for (size_t irow=0; irow<log_prob_matrix.rows(); ++irow) {
      for (size_t icol=0; icol<log_prob_matrix.cols(); ++icol) {
	size_t logProb_col = icol + irow*log_prob_matrix.cols();  // column index in logProb
	log_prob_matrix(irow, icol) = (*logProb)[0][logProb_col];
	// NOTE perhaps I should just use logProb? it would effectively just move this index manipulation from here up to lexicalTable::getValue
	// hmm. I like it like this. Probably just move this code to the place where the text file gets parsed
      }
    }
    // print();
  }

  // ----------------------------------------------------------------------------------------
  void lexicalTable::print() {
    std::cout << stringify() << std::endl;

    // print new matrix
    std::cout << std::setw(12) << "";
    if (tracks.size() == 2)
      for (size_t iter=0; iter<tracks[1]->getAlphaSize(); iter++)
	std::cout << std::setw(12) << tracks[1]->getAlpha(iter);
    std::cout << std::endl;
    for (size_t irow=0; irow<log_prob_matrix.rows(); ++irow) {
      std::cout << std::setw(12) << tracks[0]->getAlpha(irow);
      for (size_t icol=0; icol<log_prob_matrix.cols(); ++icol) {
	std::cout << std::setw(12) << log_prob_matrix(irow,icol);
      }
      std::cout << "" << std::endl;
    }
  }
    
  // ----------------------------------------------------------------------------------------
  std::vector<std::vector<double> >* lexicalTable::getCountsTable(){
    if (counts==NULL){
      counts = new(std::nothrow) std::vector<std::vector<double> >;
            
      if (counts==NULL){
	std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
	exit(1);
      }
                        
    }
    return counts;
  }
    
  // ----------------------------------------------------------------------------------------
  std::vector<std::vector<double> >* lexicalTable::getProbabilityTable(){
    if (prob==NULL){
      prob = new(std::nothrow) std::vector<std::vector<double> >;
            
      if (prob==NULL){
	std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
	exit(1);
      }
                        
    }
    return prob;
  }
    
  // ----------------------------------------------------------------------------------------
  std::vector<std::vector<double> >* lexicalTable::getLogProbabilityTable(){
    if (logProb==NULL){
      logProb = new(std::nothrow) std::vector<std::vector<double> >;
            
      if (logProb==NULL){
	std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
	exit(1);
      }
    }
    return logProb;
  }

  // ----------------------------------------------------------------------------------------
  double lexicalTable::getValue(sequences& seqs, size_t pos){
    assert(seqs.size() == tracks.size());
    assert(pos <= seqs[0].size());  // could also make sure seqs are all the same length. And that seqs are in the order that we expect based on <tracks>
    for (size_t iseq=0; iseq<seqs.size(); ++iseq) {
      assert(seqs[iseq].getTrack() == tracks[iseq]);
      assert(orders[iseq] == 0);
    }
    return log_prob_matrix(seqs[0][pos], seqs.size()==1 ? 0 : seqs[1][pos]);
  }
        
  // ----------------------------------------------------------------------------------------
  double lexicalTable::getValue(sequence& seq, size_t pos){
    assert(0); // not implemented with new matrix scheme
  }
    
  //!Add a track to an emission
  //!\param trk Pointer to track
  //!\param orderValue order of emission from track
  void lexicalTable::addTrack(track* trk,int orderValue){
    tracks.push_back(trk);
    orders.push_back(orderValue);
    if (orderValue>max_order){
      max_order=orderValue;
    }
        
  }
        
  // ----------------------------------------------------------------------------------------
  std::string lexicalTable::stringify(){
    std::string tbl("");
    size_t tracks_size = tracks.size();
        
    if(tracks_size==0){
      std::cerr << "Can't print out table without track and order being set for lexicalTable\n";
      exit(1);
    }
        
    //Output Column Headers
    size_t columns(1);
    std::vector<size_t> alphaSizes;
    //alphaSizes.push_back(0);
        
    for(size_t i = 0;i<tracks.size();i++){
      size_t alphaSz = tracks[i]->getAlphaSize();
      columns*=alphaSz;
      alphaSizes.push_back(alphaSz);
    }
        
        
    reverse(alphaSizes.begin(),alphaSizes.end());
        
    std::string colHeader("@");
        
    for(size_t i = 0;i<columns;i++){
      size_t indexValue = i;
      for(size_t tr=0;tr<tracks.size();tr++){
                
	if (tr>0){
	  colHeader+= "|";
	}
                
	size_t val(0);
	if (tr<tracks.size()-1){
	  val= floor(indexValue/alphaSizes[tr]);
	  indexValue-=val*alphaSizes[tr];
	}
	else{
	  val = indexValue;
	}
                
	colHeader+=tracks[tr]->convertIndexToWord(val, 1);
      }
      colHeader+="\t";
    }
        
    tbl+=colHeader + "\n";
        
    std::vector<std::vector<double> >* temp;
        
    if (logProb!=NULL){
      temp=logProb;
    }
    else if (prob!=NULL){
      temp=prob;
    }
    else if (counts!=NULL){
      temp=counts;
    }
    else{
      std::cerr << "No table is defined\n";
      return "";
    }
        
    //TODO: Fix row header for other orders
    bool rowHeader = (temp->size()>1) ? true : false;
    for(size_t i=0;i<temp->size();i++){
      std::string header("");
            
      if (rowHeader){
	size_t indexValue = i;
                
	for(size_t tr=0;tr<tracks.size();tr++){
                    
	  if (tr>0 && orders[tr]>0){
	    header+= "|";
	  }
                    
                    
	  size_t val(0);
	  if (tr<tracks.size()-1){
	    double pwr = POWER[orders[tr+1]][tracks[tr+1]->getAlphaSize()-1];
	    val= floor(indexValue/pwr);
	    indexValue-=val*pwr;
	  }
	  else{
	    val = indexValue;
	  }
                    
	  header+=tracks[tr]->convertIndexToWord(val, orders[tr]);
	}
	tbl+="@" + header + "\t";
                
      }
            
      for(size_t j=0;j<(*temp)[i].size();j++){
	if (j>0){
	  tbl+="\t";
	}
	tbl+=double_to_string((*temp)[i][j]);
      }
      tbl+="\n";
    }
    return tbl;

  }

}
