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
    log_emission = NULL;
    x_subarray=NULL;
    y_subarray=NULL;
        
        
    unknownScoreType=NO_SCORE;
    unknownDefinedScore=-INFINITY;
        
    return;
  }
    
  lexicalTable::~lexicalTable(){
    delete logProb;
    delete prob;
    delete counts;
    delete log_emission;
    delete x_subarray;
    delete y_subarray;

    assert(0);  // don't forget to delete all that new stuff you added, dude
    // also TODO:
    //    - two seqs on *different* emissions
    //    - two seqs on same track (?)
    //    - write print function for new matrices
    //    - implement reduced order
        
    logProb=NULL;
    prob=NULL;
    counts=NULL;
    log_emission = NULL;
    x_subarray=NULL;
    y_subarray=NULL;
  }
    
  void lexicalTable::createTable(int rows, int columns, int pseudocount, valueType typ){
    if (typ==COUNTS){
      if (counts!=NULL){
	delete counts;
      }
      std::vector<double> temp_columns(columns,pseudocount);
      counts=new(std::nothrow) std::vector<std::vector<double> > (rows,temp_columns);
            
      if (counts==NULL){
	std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
	exit(1);
      }
    }
    return;
  }
    
  void lexicalTable::print() {
    std::cout << stringify() << std::endl;
  }
    
    
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

  //! \return indices we need to average over (when pos < order)
  vector<size_t> lexicalTable::getIndicesToAverage(sequence seq, size_t order, size_t pos, std::vector<uint8_t> pre_word) {
    // build all possible words which end in <pre_word>
    // eg if <pre_word> is CT and this is a third order track, we want [ACT, CCT, GCT, TCT]
    // To do this construct a vector<vector> with all the posibilities for *individual* positions:
    //    0 1 2 3
    // 0  A C G T
    // 1  C
    // 2  T
    // Then loop over this vector<vector> to push back all possibilities
    vector<vector<size_t> > tmp_no_name_vector;  // no, I don't know what to fuckin call it
    for (size_t iter=0; iter<order; ++it) {
      tmp_no_name_vector.push_back(vector<size_t>());
      if (iter < order-pos) {  // push back all the characters in the alphabet -- we need to average over all possible ones
	for (size_t ichar=0; ichar<seq->track()->getAlphaSize(); ++ichar)
	  tmp_no_name_vector[iter].push_back(ichar);
      } else {  // just push back the actual sequence character
	tmp_no_name_vector[iter].push_back(pre_word[iter-order+pos]);
      }
    }

    for (size_t ipos=0; ipos<tmp_no_name_vector.size(); ++ipos) {
      for (size_t ichar=0; ichar<tmp_no_name_vector[ipos].size(); ++ichar) {
	std::cout << " " << tmp_no_name_vector[ipos][ichar];
      }
      std::cout << "" << std::endl;
    }
    return 
  }
    
  //! \return log emission probability
  double lexicalTable::getValue(sequences& seqs, size_t pos){
    assert(seqs.size() == trcks.size());
    assert(pos <= seqs[0].size());  // could also make sure seqs are all the same length. And that seqs are in the order that we expect based on <trcks>
    std::vector<std::vector<size_t> > iwords;  // indices of the words for each sequence (second dimension is for all the ones we have to average over when pos < order[iseq])
    std::vector<std::vector<size_t> > iemissions;  // indices of the emissions for each sequence
    for (size_t iseq=0; iseq<seqs.size(); ++iseq) {
      assert(seqs[iseq].getTrack() == trcks[iseq]);
      std::vector<size_t> tmp_words;
      std::vector<size_t> tmp_emissions;
      std::vector<uint8_t> pre_word = seqs[iseq].getDigitalSubSeq(max(pos-order[iseq], 0), pos);  // the preceding word, i.e. the part of the sequence on which this emission depends
      if (pre_word.size() < order[iseq]) {  // still near the start of the sequence
	assert(pos-order[iseq] < 0);  // don't need this check
	
	
	// wait, wait -- average the probs, or the log probs?
      size_t matrix_index = trcks[iseq]->convertDigiWordToIndex(pre_word);  // index of the matrix we want in matrix_ptrs

      
      iwords.push_back(matrix_index);
      iemissions.push_back(seqs[iseq][pos]);
    }
    return (*matrix_ptrs(iwords[0], iwords[1]))(iemissions[0], iemissions[1]);
  }
        
  //Return emission probability of sequences
  double lexicalTable::getValue(sequence& seq, size_t pos){
    assert(0); // not implemented with new matrix scheme
        
    if (max_order>pos){
      return getReducedOrder(seq, pos);
    }
                
    size_t index(seq[pos - subarray_position[0]] * subarray_value[0]);
                
    for(size_t i=1;i<dimensions;i++){
      index += seq[pos - subarray_position[i]] * subarray_value[i];
    }
                
    if (index > array_size){
      std::cerr << "Index is out of range of lookup table in lexicalTable" << std::endl;
      exit(2);
    }
                
    return (*log_emission)[index];
  }
    
    
  //!Add a track to an emission
  //!\param trk Pointer to track
  //!\param orderValue order of emission from track
  void lexicalTable::addTrack(track* trk,int orderValue){
    trcks.push_back(trk);
    alphabets.push_back(trk->getAlphaSize());
    order.push_back(orderValue);
    if (orderValue>max_order){
      max_order=orderValue;
    }
        
  }
        
  //!Set the type of counts in the emission 2D table provided by the user
  //!\param temp vector of vectors of doubles
  //!\param emmType Type of value (COUNTS, PROBABILITY, LOG_PROB)
  void lexicalTable::assignTable(std::vector<std::vector<double > >* temp, valueType emmType){
    if (emmType==COUNTS){
      counts=temp;
    }
    else if (emmType == PROBABILITY){
      prob=temp;
    }
    else if (emmType == LOG_PROB){
      logProb=temp;
      initialize_emission_table();
    }
  }
        
        
  std::string lexicalTable::stringifyAmbig(){
    std::string tbl("");
                
    size_t tracks_size = trcks.size();
        
    if(tracks_size==0){
      std::cerr << "Can't print out table without track and order being set for lexicalTable\n";
      exit(1);
    }
        
    //Calculate the Number of Column Headers and get alphabet for each track
    //This is the complete unambiguous and ambiguous
                
    std::vector<std::vector<std::string> > complete_alphabet(tracks_size, std::vector<std::string>());
                
    size_t columns(1);
    std::vector<size_t> alphaSizes;
        
    //Calculate the columns size 
    for(size_t i = 0;i<trcks.size();i++){
      size_t alphaSz = (trcks[i]->isAmbiguousSet()) ? trcks[i]->getMaxAmbiguous()+1 : trcks[i]->getAlphaSize();
      columns*=alphaSz;
      alphaSizes.push_back(alphaSz);
                        
      //Get complete alphabet for each track
      for(size_t j=0; j < alphaSz; ++j){
	complete_alphabet[i].push_back(trcks[i]->getAlpha(j));
      }
    }
        
    reverse(alphaSizes.begin(),alphaSizes.end());
        
    std::string colHeader("@");
        
    //Compose column heading
    for(size_t i = 0; i < columns; ++i){
      size_t indexValue = i;
      for(size_t tr=0;tr<trcks.size();tr++){
                
	if (tr>0){
	  colHeader+= "|";
	}
                
	size_t val(0);
	if (tr<trcks.size()-1){
	  val= floor(indexValue/alphaSizes[tr]);
	  indexValue-=val*alphaSizes[tr];
	}
	else{
	  val = indexValue;
	}
                
	colHeader+=complete_alphabet[tr][val];
      }
      colHeader+="\t";
    }
        
    tbl+=colHeader + "\n";
                
                
    //              for (size_t i=0; i< log_emission->size();i++){
    //                      std::cout << (*log_emission)[i] << std::endl;
    //              }
                
                
    //       bool rowHeader = (temp->size()>1) ? true : false;
    //       
    //              for(size_t i=0;i<temp->size();i++){
    //            std::string header("");
    //            
    //            if (rowHeader){
    //                size_t indexValue = i;
    //                
    //                for(size_t tr=0;tr<trcks.size();tr++){
    //                    
    //                    if (tr>0 && order[tr]>0){
    //                        header+= "|";
    //                    }
    //                    
    //                    
    //                    size_t val(0);
    //                    if (tr<trcks.size()-1){
    //                        double pwr = POWER[order[tr+1]][trcks[tr+1]->getAlphaSize()-1];
    //                        val= floor(indexValue/pwr);
    //                        indexValue-=val*pwr;
    //                    }
    //                    else{
    //                        val = indexValue;
    //                    }
    //                    
    //                    header+=trcks[tr]->convertIndexToWord(val, order[tr]);
    //                }
    //                tbl+="@" + header + "\t";
    //                
    //            }
    //            
    //            for(size_t j=0;j<(*temp)[i].size();j++){
    //                if (j>0){
    //                    tbl+="\t";
    //                }
    //                tbl+=double_to_string((*temp)[i][j]);
    //            }
    //            tbl+="\n";
    //        }

                
    return tbl;
  }
        
        
        
    
  std::string lexicalTable::stringify(){
    std::string tbl("");
    size_t tracks_size = trcks.size();
        
    if(tracks_size==0){
      std::cerr << "Can't print out table without track and order being set for lexicalTable\n";
      exit(1);
    }
        
    //Output Column Headers
    size_t columns(1);
    std::vector<size_t> alphaSizes;
    //alphaSizes.push_back(0);
        
    for(size_t i = 0;i<trcks.size();i++){
      size_t alphaSz = trcks[i]->getAlphaSize();
      columns*=alphaSz;
      alphaSizes.push_back(alphaSz);
    }
        
        
    reverse(alphaSizes.begin(),alphaSizes.end());
        
    std::string colHeader("@");
        
    for(size_t i = 0;i<columns;i++){
      size_t indexValue = i;
      for(size_t tr=0;tr<trcks.size();tr++){
                
	if (tr>0){
	  colHeader+= "|";
	}
                
	size_t val(0);
	if (tr<trcks.size()-1){
	  val= floor(indexValue/alphaSizes[tr]);
	  indexValue-=val*alphaSizes[tr];
	}
	else{
	  val = indexValue;
	}
                
	colHeader+=trcks[tr]->convertIndexToWord(val, 1);
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
                
	for(size_t tr=0;tr<trcks.size();tr++){
                    
	  if (tr>0 && order[tr]>0){
	    header+= "|";
	  }
                    
                    
	  size_t val(0);
	  if (tr<trcks.size()-1){
	    double pwr = POWER[order[tr+1]][trcks[tr+1]->getAlphaSize()-1];
	    val= floor(indexValue/pwr);
	    indexValue-=val*pwr;
	  }
	  else{
	    val = indexValue;
	  }
                    
	  header+=trcks[tr]->convertIndexToWord(val, order[tr]);
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
        
        
        
        
        
  void lexicalTable::init_table_dimension_values(){
    //
    // initialize the arrays of size_t x_subarray and y_subarray
    // for, e.g. this hmm:
    //      track     order     alphabet (size)
    //        0         2         AC (2)
    //        1         1         GT (2)
    //        2         1         MF (2)
    // we have
    //   x_subarray = [4, 2, 1]
    //   y_subarray = [8, 4, 2, 1]

    number_of_tracks = trcks.size();
    y_dim = sumVector(order);  // <order> is a vector containing the order of each track
                
                
    //Calculate subarray dimensions for logProb and new log_emission table (includes ambiguous character)
    x_subarray = new size_t[number_of_tracks];
    y_subarray = new size_t[y_dim];
                
    //Calculate Old subarray x_dimension values
    for(size_t i=0;i<number_of_tracks;++i){
      size_t value(1);
      for(size_t j=i+1;j<number_of_tracks;++j){
	value*=alphabets[j];  // multiply by the size of this track's alphabet 
      }
      x_subarray[i]=value;
    }
                
    //Calcuate Old subarray y_dimension values
    std::vector<size_t> index(order[0],alphabets[0]);  // initialize <index> with a few entries
    for(size_t i=1;i<order.size();i++){  // then push a few more onto the tail end of it
      for(size_t j=0;j<order[i];j++){
	index.push_back(alphabets[i]);  // so final length of <index> is y_dim = (sum of all orders for all emissions)
      }
    }

    for (size_t i=0;i<y_dim;i++){
      size_t val = 1;
      for(size_t j=i+1;j<y_dim;j++){
	val*=index[j];
      }
      y_subarray[i]=val;
    }
                
    return;
  }
        
  //TODO:  NEED TO COMPLETE THIS FUNCTION
  //      size_t lexicalTable::convertIndex(size_t old_row, size_t old_column){
  //              //Decompose previous value from indices to digital letter value
  //              
  //      }
        
        
  void lexicalTable::init_array_dimension_values(){
    // 
    //
    //
    //

    // NOTE that number_of_tracks > 0 only for joint emissions (as far as I can tell)

    // RECALL y_dim = (sum of the orders of each track)
    dimensions = y_dim + number_of_tracks;
                                            
    subarray_sequence.resize(dimensions);  // EG for the model above <subarray_sequence> is [0 0 0 1 1 2 2]
    subarray_value.resize(dimensions);
    subarray_value.resize(dimensions);  // NOTE I am pretty sure this second call does nothing

    // NOTE these are only used when we have ambiguous characters
    decompose_values.resize(dimensions);
    decompose_sequence.resize(dimensions);
                
    //Calculate total size of emission table needed with ambiguous characters
    array_size = 1;
    std::vector<size_t> complete_alphabet_size;
    size_t current_dim(0);  // index to keep track of how much we've filled up subarray_sequence
    for(size_t i=0;i<number_of_tracks;i++){
                        
      //Get alphabet size and store it
      size_t alpha_size = trcks[i]->getTotalAlphabetSize();  // NOTE <alpha_size> is always equal to alphabets[i] so I think we don't need this
      complete_alphabet_size.push_back(alpha_size);
      assert(complete_alphabet_size[i] == alphabets[i]); // double check to make sure that complete_alphabet_size is redundant
      array_size*=integerPower(alpha_size, (size_t) order[i]+1);  // multiply <array_size> by (alphabets[i])^(order[i]+1). NOTE integerPower(base, exponent) = base^exponent
                        
      for(size_t j=0;j<=order[i];++j){  // NOTE the equals sign
	subarray_sequence[current_dim]=i;
	current_dim++;
      }
    }
                
    //Calculate the Sequence positions
    // EG for model above subarray_position = [2 1 0 1 0 1 0]
    for (size_t i=0;i<number_of_tracks;i++){
      for (size_t j=0;j<=order[i];++j){
	subarray_position.push_back(order[i]-j);
      }
    }
                

    //Compute the decomposing values
    //Used to convert from index to word
    std::vector<size_t> decompose_index;
    for(size_t i=0;i<number_of_tracks;++i){
      for(size_t j=0;j<order[i];++j){
	decompose_index.push_back(complete_alphabet_size[i]);
      }
    }
    for(size_t i=0;i<number_of_tracks;++i){
      decompose_index.push_back(complete_alphabet_size[i]);
    }
                
    //Calculate subarray values
    for(size_t i=0;i<dimensions;++i){
      decompose_values[i]=1;
      for(size_t j=i+1;j<dimensions;++j){
	decompose_values[i]*=decompose_index[j];
      }
    }
                
                
    //Rearrange decompose values for subarray values;
    //Final values  in subarray_value will correspond to sequence AAA(A)B(B).
    //Where A is 3rd order and B is 1st order;
    size_t array_index(0);
    size_t index(0);
    for(size_t i=0;i<number_of_tracks;++i){
      for(size_t j=0;j<order[i];++j){
	subarray_value[array_index] = decompose_values[index];
	array_index++;
	index++;
      }
      subarray_value[array_index] = decompose_values[dimensions-number_of_tracks-i];
      array_index++;
    }
    std::cout << "---" << std::endl;
    for (size_t i=0; i<subarray_value.size(); ++i) std::cout << "  " << i << " " << subarray_value[i] << std::endl;
                
                
    //Finalize decompose_sequence
    std::vector<size_t> temp;
    for(size_t i=0;i<number_of_tracks;++i){
      for(size_t j=0;j<order[i];++j){
	temp.push_back(i);
      }
    }
    for(size_t i=0;i<number_of_tracks;++i){
      temp.push_back(i);
    }
                
    for(size_t i=0;i<dimensions;++i){
      decompose_sequence[i]=temp[i];
    }

    return;
  }
        
  //Transfer values from 2d table to array and also compute the ambiguous score
  void lexicalTable::transferValues(std::vector<bool>& transferred){
                
    //Transfer unambiguous scores
    std::cout << "====" << std::endl;
    for(size_t row=0;row<logProb->size();++row){
      for(size_t column=0;column<(*logProb)[row].size();++column){
	std::vector<uint8_t> alphabet;  // vector which is used to calculate the index in the (linear) array <log_emission> which corresponds to the entry [row][column] in the table <logProb>
	decompose(row, column, alphabet);  // push values into <alphabet>. *why* these particular values, you say? I don't know!
	size_t index = calculateArrayIndex(alphabet);
	(*log_emission)[index] = (*logProb)[row][column];  // NOTE as far as I can tell, this sometimes maps *different* row:column pairs to the *same* index
	std::cout << "  " << row << " " << column << " --> " << index << std::endl;
	transferred[index] = true;
      }
    }
    std::cout << "====" << std::endl;
                
    //Compute all ambiguous scores
    for(size_t i=0;i<transferred.size();i++){
      if (transferred[i]){
	continue;
      }
                        
      if (unknownScoreType == DEFINED_SCORE){
	(*log_emission)[i] = unknownDefinedScore;
	continue;
      }
      else if (unknownScoreType == NO_SCORE){
	continue;
      }
                        
      //Get the letters for index
      std::vector<uint8_t> letters;
      decompose(i,letters);
                        
      //Expand the unambiguous words and get all the possible values
      //std::vector<double> expanded;
      //expand_ambiguous(letters,expanded);
                        
      (*log_emission)[i] = getAmbiguousScore(letters);
                        
      //                      //Assign the values accordint the Score Type
      //                      if (unknownScoreType == AVERAGE_SCORE){
      //                              (*log_emission)[i] = avgLogVector(expanded);
      //                      }
      //                      else if (unknownScoreType == LOWEST_SCORE){
      //                              (*log_emission)[i] = minVector(expanded);
      //                      }
      //                      else if (unknownScoreType == HIGHEST_SCORE){
      //                              (*log_emission)[i] = maxVector(expanded);
      //                      }
    }
                
    //              for (size_t i=0;i<log_emission->size();++i){
    //                      std::vector<uint8_t> letters;
    //                      decompose(i,letters);
    //                      for (size_t j = 0; j< letters.size();j++){
    //                              std::cout << (int) letters[j];
    //                      }
    //                      std::cout << "\t" ;
    //                      std::cout << (*log_emission)[i] << std::endl;
    //              }
    return;
  }
        
  //Calculate lower order emission from the current table values
  //Given order and position/sequence
  //Calculate the values using Index and [all alphabets] for higher orders
  double lexicalTable::getReducedOrder(sequences& seqs, size_t position){
    Index indices;
    for(size_t i=0;i<dimensions;i++){
      Index subtotal;
      size_t seq = subarray_sequence[i];
      size_t pos = subarray_position[i];
                        
      if (subarray_position[i] > position){
	subtotal.setAmbiguous(trcks[seq]->getUnambiguousSet());
      }
      else if (seqs[seq][position - pos] > max_unambiguous[seq]){
	subtotal.setAmbiguous(trcks[seq]->getAmbiguousSet(seqs[seq][position-pos]));
      }
      else{
	subtotal+=seqs[seq][position-pos];
      }
                        
      subtotal *= subarray_value[i];
      indices  += subtotal;
    }
                
    if (unknownScoreType == AVERAGE_SCORE || unknownScoreType == NO_SCORE){
      double temp(0);
      for(size_t i=0;i<indices.size();i++){
	temp+=exp((*log_emission)[indices[i]]);
      }
      temp /= indices.size();
      temp = log(temp);
      return temp;
    }
    else if (unknownScoreType == LOWEST_SCORE){
      double temp(INFINITY);
      for(size_t i=0;i<indices.size();i++){
	double val = (*log_emission)[indices[i]];
	if (val < temp){
	  temp = val;
	}
      }
      return temp;
    }
    else if (unknownScoreType == HIGHEST_SCORE){
      double temp(-INFINITY);
      for(size_t i=0;i<indices.size();i++){
	double val = (*log_emission)[indices[i]];
	if (val > temp){
	  temp = val;
	}
      }
      return temp;
    }
    return -INFINITY;
  }
        
        
  //Calculate lower order emission from the current table values
  //Given order and position/sequence
  //Calculate the values using Index and [all alphabets] for higher orders
  double lexicalTable::getReducedOrder(sequence& seq, size_t position){
    Index indices;
    for(size_t i=0;i<dimensions;i++){
      Index subtotal;
      size_t sq = subarray_sequence[i];
      size_t pos = subarray_position[i];
                        
      if (subarray_position[i] > position){
	subtotal.setAmbiguous(trcks[sq]->getUnambiguousSet());
      }
      else if (seq[position - pos] > max_unambiguous[sq]){
	subtotal.setAmbiguous(trcks[sq]->getAmbiguousSet(seq[position-pos]));
      }
      else{
	subtotal+=seq[position-pos];
      }
                        
      subtotal *= subarray_value[i];
      indices  += subtotal;
    }
                
    if (unknownScoreType == AVERAGE_SCORE || unknownScoreType == NO_SCORE){
      double temp(0);
      for(size_t i=0;i<indices.size();i++){
	temp+=exp((*log_emission)[indices[i]]);
      }
      temp /= indices.size();
      temp = log(temp);
      return temp;
    }
    else if (unknownScoreType == LOWEST_SCORE){
      double temp(INFINITY);
      for(size_t i=0;i<indices.size();i++){
	double val = (*log_emission)[indices[i]];
	if (val < temp){
	  temp = val;
	}
      }
      return temp;
    }
    else if (unknownScoreType == HIGHEST_SCORE){
      double temp(-INFINITY);
      for(size_t i=0;i<indices.size();i++){
	double val = (*log_emission)[indices[i]];
	if (val > temp){
	  temp = val;
	}
      }
      return temp;
    }
    return -INFINITY;
  }

        
  double lexicalTable::getAmbiguousScore(std::vector<uint8_t>& letters){
    Index indices;
    for(size_t i=0;i<dimensions;++i){
      Index subtotal;
      if (letters[i]>max_unambiguous[decompose_sequence[i]]){
	subtotal.setAmbiguous(trcks[decompose_sequence[i]]->getAmbiguousSet(letters[i]));
      }
      else{
	subtotal+= letters[i];
      }
                        
      subtotal *= decompose_values[i];
      indices += subtotal;
    }
                
    if (unknownScoreType == AVERAGE_SCORE){
      double temp(0);
      for(size_t i=0;i<indices.size();i++){
	temp+=exp((*log_emission)[indices[i]]);
      }
      temp /= indices.size();
      temp = log(temp);
                                                
      return temp;
    }
    else if (unknownScoreType == LOWEST_SCORE){
      double temp(INFINITY);
      for(size_t i=0;i<indices.size();i++){
	double val = (*log_emission)[indices[i]];
	if (val < temp){
	  temp = val;
	}
      }
      return temp;
    }
    else if (unknownScoreType == HIGHEST_SCORE){
      double temp(-INFINITY);
      for(size_t i=0;i<indices.size();i++){
	double val = (*log_emission)[indices[i]];
	if (val > temp){
	  temp = val;
	}
      }
      return temp;
    }
                
    return -INFINITY;
  }
        
        
  //Use index instead much faster than expanding the letters and then reverse computing 
  //      void lexicalTable::expand_ambiguous(std::vector<uint8_t>& letters, std::vector<double>& expanded){
  //              
  //              std::vector<std::vector<uint8_t> >* temp_words = new std::vector<std::vector<uint8_t> >;
  //              temp_words->push_back(letters);
  //              for(size_t i=0;i<dimensions;i++){
  //                      temp_words = expand_ambiguous(temp_words, i);
  //              }
  //              
  //              for(size_t i=0;i<temp_words->size();i++){
  //                      size_t index = calculateIndexFromDecomposed((*temp_words)[i]);
  //                      expanded.push_back((*log_emission)[index]);
  //              }
  //              
  //              return;
  //      }
  //      
  //      std::vector<std::vector<uint8_t> >* lexicalTable::expand_ambiguous(std::vector<std::vector<uint8_t> >* words, size_t letter){
  //              
  //              std::vector<std::vector<uint8_t> >* temp_words= new std::vector<std::vector<uint8_t> >;
  //              for(size_t i=0; i < words->size(); ++i){
  //                      if ((*words)[i][letter] > max_unambiguous[decompose_sequence[letter]]){
  //                              std::vector<uint8_t>& set = trcks[decompose_sequence[letter]]->getAmbiguousSet((*words)[i][letter]);
  //                              for(size_t j=0;j<set.size();++j){
  //                                      (*words)[i][letter] = set[j];
  //                                      temp_words->push_back((*words)[i]);
  //                              }
  //                      }
  //                      else{
  //                              temp_words->push_back((*words)[i]);
  //                      }
  //              }
  //              
  //              delete words;
  //              words = NULL;
  //              return temp_words;
  //      }
        
  size_t lexicalTable::calculateIndexFromDecomposed(std::vector<uint8_t>& word){
    size_t index(0);
    for(size_t i=0;i<dimensions;++i){
      index += decompose_values[i] * word[i];
    }
    return index;
  }
        
  size_t lexicalTable::calculateArrayIndex(std::vector<uint8_t>& kmer){
    size_t index(0);
    for(size_t i=0;i<dimensions;++i){  // RECALL <dimensions> = (sum of orders) + <number_of_tracks>
      index += subarray_value[i] * kmer[i];
    }
    return index;
  }
  // using namespace std;
  void lexicalTable::decompose(size_t row, size_t column, std::vector<uint8_t>& letters){
    // std::cout << "decomposing"
    // 	      << " row " << row
    // 	      << " column " << column
    // 	      << std::endl;
    // std::cout
    //   << setw(12) << "i"
    //   << setw(12) << "row"
    //   << setw(18) << "y_subarray[i]"
    //   << setw(12) << "val"
    //   << std::endl;
    //Decompose the row into the preceding letters
    for(size_t i=0;i<y_dim;++i){  // RECALL y_dim is the sum of the orders of each emission
      size_t val = floor(row/y_subarray[i]);
      row -= val*y_subarray[i];
      letters.push_back(val);
      // std::cout
      // 	<< setw(12) << i
      // 	<< setw(12) << row
      // 	<< setw(18) << y_subarray[i]
      // 	<< setw(12) << val
      // 	<< std::endl;
    }
                
    //Decompose the column into the emitted letters
    for(size_t i=0;i<number_of_tracks;++i){
      size_t val = floor(column/x_subarray[i]);
      column-=val*x_subarray[i];
      letters.push_back(val);
    }

    // std::cout << "returning letters ";
    // for (size_t il=0; il<letters.size(); ++il) std::cout << " " << (int)(letters[il]);
    // std::cout << "" << std::endl;
    return;
  }
        
        
  void lexicalTable::decompose(size_t index, std::vector<uint8_t>& letters){
                
    //Decompose the index into the emitted letters
    for(size_t i=0;i<dimensions;++i){
      size_t val = floor(index/subarray_value[i]);
      index-=val*subarray_value[i];
      letters.push_back(val);
    }
                
    return;
  }

  //Todo
  //Convert table to simpleNtable compatible format
  //with ambiguous characters
  void lexicalTable::initialize_emission_table(){
    if (logProb == NULL){
      std::cerr << "Cannot initialize emission table until after the tables have been assigned";
      exit(2);
    }
    init_table_dimension_values();
    init_array_dimension_values();
                
    for(size_t i = 0; i < number_of_tracks ; ++i){
      max_unambiguous.push_back(trcks[i]->getMaxUnambiguous());
    }
                
    if(unknownDefinedScore == DEFINED_SCORE){
      log_emission = new std::vector<double> (array_size,unknownDefinedScore);
    }
    else{
      log_emission = new std::vector<double> (array_size,-INFINITY);
    }
                
    //Transfer values to emission_table
    std::vector<bool> transferred (array_size,false);
    transferValues(transferred);

    // ----------------------------------------------------------------------------------------
    assert(number_of_tracks == trcks.size());  // not sure why we need both of these
    assert(number_of_tracks == 1 || number_of_tracks == 2);  // for the moment I limit to *two* dimensions (tracks), since there doesn't seem to exist a good library for n-dim matrices with n unknown until runtime
    size_t n_rows = pow(alphabets[0], order[0]);
    size_t n_columns = 1;
    if (number_of_tracks == 2)
      n_columns = POWER[order[1]][alphabets[1]-1];
    assert(n_rows*n_columns == (*logProb).size());  // number or entries in matrix_ptrs should equal the number of rows in logProb
    matrix_ptrs.resize(n_rows, n_columns);
    std::cout << "init matrix with " << matrix_ptrs.rows() << " rows and " << matrix_ptrs.cols() << " columns" << std::endl;
    for (size_t ir=0; ir<matrix_ptrs.rows(); ++ir) {  // loop over words (of length <order>) on track 1
      std::string track_1_word = trcks[0]->convertIndexToWord(ir, order[0]);
      assert(trcks[0]->convertAlphaWordToIndex(track_1_word) == ir);  // double check for closure
      for (size_t ic=0; ic<matrix_ptrs.cols(); ++ic) {  // loop over words on track 2
	std::string track_2_word = trcks[1]->convertIndexToWord(ic, order[1]);
	assert(trcks[1]->convertAlphaWordToIndex(track_2_word) == ic);  // double check for closure
	
	matrix_ptrs(ir,ic) = new Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
	matrix_ptrs(ir,ic)->resize(alphabets[0], (number_of_tracks==1) ? 1 : alphabets[1]);
	for (size_t isubrow=0; isubrow<matrix_ptrs(ir,ic)->rows(); ++isubrow) {  // loop over track 1 emissions for this word
	  for (size_t isubcol=0; isubcol<matrix_ptrs(ir,ic)->cols(); ++isubcol) {  // loop over track 2 emissions for this word
	    // ir: track 1 word, ic: track2 word, isubrow: track 1 emission, isubcol: track 2 emission
	    // now convert from our four-index set to the indces in the 2d matrix logProb, eg for two tracks [HT] [01], both second order
	    //     00 01 10 11
	    //  HH                                                              0 1
	    //  HT                                                           H
	    //  TH    x          --> x is a pointer to the log probs here:   T
	    //  TT
	    //
	    // whereas logProb represents this as
	    //   emission: H0 H1 T0 T1
	    //  word:
	    //    HH,00
	    //    HH,01
	    //    HH,10
	    //    HH,11
	    //    HT,00
	    //    ...
	    // index of this word in track 2 (ic), plus the number of times we had to cycle through all words in track 1 (ir) to get here, times the number of possible track 1 words (matrix_ptrs.cols())
	    size_t logProb_row = ic + ir*matrix_ptrs.cols();
	    // index of this character (emission) in track 2 (isubcol), plus the number of times we had to cycle through all characters in track 1 (isubrow) to get here, times the size of the track 1 alphabet (matrix_ptrs(ir,ic)->cols())
	    size_t logProb_col = isubcol + isubrow*matrix_ptrs(ir,ic)->cols();
	    (*matrix_ptrs(ir,ic))(isubrow, isubcol) = (*logProb)[logProb_row][logProb_col];
	    // NOTE perhaps I should just use logProb? it would effectively just move this index manipulation from here up to lexicalTable::getValue
	    // hmm. I like it like this. Probably just move this code to the place where the text file gets parsed
	  }
	}
      }
    }
    print();
  }
}
