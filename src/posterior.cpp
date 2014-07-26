//
//  posterior.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "trellis.h"

namespace StochHMM {
// ----------------------------------------------------------------------------------------
void trellis::posterior(model *h, sequences* sqs){
  hmm = h;
  seqs = sqs;
  seq_size                = seqs->GetSequenceLength();
  state_size              = hmm->state_size();
  exDef_defined   = seqs->exDefDefined();
              
  if (posterior_score!=NULL){
    delete posterior_score;
    posterior_score = NULL;
  }
  ending_backward_prob = -INFINITY;
  ending_forward_prob  = -INFINITY;
              
  posterior();
}
      
// ----------------------------------------------------------------------------------------
void trellis::posterior() {
  posterior_score = new (std::nothrow) double_2D(seq_size, std::vector<double>(state_size, -INFINITY));
  scoring_current = new (std::nothrow) std::vector<double> (state_size, -INFINITY);
  scoring_previous= new (std::nothrow) std::vector<double> (state_size, -INFINITY);
  if (scoring_current == NULL || scoring_previous == NULL || posterior_score == NULL) {
    std::cerr << "Can't allocate Posterior score table. OUT OF MEMORY" << std::endl;
    exit(2);
  }
              
  std::bitset<STATE_MAX> next_states;
  std::bitset<STATE_MAX> current_states;
              
  double  forward_temp(-INFINITY);
  double  emission(-INFINITY);
  state* init = hmm->getInitial();
  std::bitset<STATE_MAX>* initial_to = hmm->getInitialTo();
  std::bitset<STATE_MAX>* from_trans(NULL);
              
  // calculate forward scores from INIT state, and initialize next_states
  for(size_t i=0; i<state_size; ++i) {
    if ((*initial_to)[i]) {  // if the bitset is set (meaning there is a transition to this state), calculate the viterbi
      forward_temp = (*hmm)[i]->get_emission_prob(*seqs,0) +  getTransition(init, i, 0);
      if (forward_temp > -INFINITY) {
        (*scoring_current)[i] = forward_temp;
        next_states |= (*(*hmm)[i]->getTo());
      }
    }
  }

  // calculate the rest of the forward scores
  for(size_t position=1; position<seq_size; ++position) {
    (*posterior_score)[position-1].assign(scoring_current->begin(), scoring_current->end());
          
    // swap current and previous viterbi scores
    scoring_previous->assign(state_size, -INFINITY);
    swap_ptr = scoring_previous;
    scoring_previous = scoring_current;
    scoring_current = swap_ptr;
                      
    //Swap current_states and next states sets
    current_states.reset();
    current_states |= next_states;
    next_states.reset();
          
    for (size_t current=0; current<state_size; ++current) { //i is current state that emits value
      if (!current_states[current])
        continue;
              
      emission = (*hmm)[current]->get_emission_prob(*seqs, position);
      if (emission == -INFINITY)
        continue;
              
      from_trans = (*hmm)[current]->getFrom();
                              
      for (size_t previous=0; previous<state_size; ++previous) {  //j is previous state
        if (!(*from_trans)[previous])
          continue;
                                      
        if ((*scoring_previous)[previous] != -INFINITY) {
          forward_temp = (*scoring_previous)[previous] + emission + getTransition((*hmm)[previous], current , position);
                                              
          if ((*scoring_current)[current] == -INFINITY) {
            (*scoring_current)[current] = forward_temp;
          } else {
            (*scoring_current)[current] = AddInLogSpace(forward_temp, (*scoring_current)[current]);  // add the unloged vals and relog, i.e. use OR
          }
          next_states |= (*(*hmm)[current]->getTo());
        }
      }
    }
  }
              
  (*posterior_score)[seq_size-1].assign(scoring_current->begin(), scoring_current->end());
              
  // Swap current and previous scores
  scoring_previous->assign(state_size, -INFINITY);
  swap_ptr = scoring_previous;
  scoring_previous = scoring_current;
  scoring_current = swap_ptr;

  // set ending_forward_prob, the total prob of all paths from INIT to END
  ending_forward_prob = -INFINITY;
  for(size_t i=0; i<state_size ;++i) {
    if ((*scoring_previous)[i] != -INFINITY) {
      forward_temp = (*scoring_previous)[i] + (*hmm)[i]->getEndTrans();
      if (forward_temp > -INFINITY) {
        if (ending_forward_prob == -INFINITY) {
          ending_forward_prob = forward_temp;
        } else {
          ending_forward_prob = AddInLogSpace(forward_temp, ending_forward_prob);
        }
      }
    }
  }

  // Perform Backward Algorithm
  double  backward_temp(-INFINITY);               
  scoring_previous->assign(state_size, -INFINITY);
  scoring_current->assign(state_size, -INFINITY);
  std::vector<double> posterior_sum(seq_size, -INFINITY);
  std::bitset<STATE_MAX>* ending_from = hmm->getEndingFrom();

  // calculate backward scores from INIT state, and initialize next_states
  for(size_t st_current=0; st_current<state_size; ++st_current) {
    if ((*ending_from)[st_current]) {  //if the bitset is set (meaning there is a transition to this state)
      backward_temp = (*hmm)[st_current]->getEndTrans();
      if (backward_temp > -INFINITY) {
        (*scoring_current)[st_current] = backward_temp;
        next_states[st_current] = 1;
      }
    }
  }

  for(size_t position=seq_size-2; position!=SIZE_MAX ; --position ) {
    //Swap current_states and next states sets
    current_states.reset();
    current_states |= next_states;
    next_states.reset();

    for (size_t i=0; i<state_size; ++i) {
      if ((*posterior_score)[position+1][i] != -INFINITY && (*scoring_current)[i]!= -INFINITY) {
        (*posterior_score)[position+1][i] = ((double)(*posterior_score)[position+1][i] + (double)(*scoring_current)[i]) - ending_forward_prob;
        if ((*posterior_score)[position+1][i] > -7.6009) {  //Above significant value;
          posterior_sum[position+1] = AddInLogSpace(posterior_sum[position+1], (*posterior_score)[position+1][i]);
        }
      }
    }

    //Swap current and previous viterbi scores
    scoring_previous->assign(state_size,-INFINITY);
    swap_ptr = scoring_previous;
    scoring_previous = scoring_current;
    scoring_current = swap_ptr;

    for (size_t st_previous=0; st_previous<state_size; ++st_previous) { //i is current state that emits value
      if (!current_states[st_previous])
        continue;

      emission = (*hmm)[st_previous]->get_emission_prob(*seqs, position+1);
      if (emission == -INFINITY)
        continue;

      from_trans = (*hmm)[st_previous]->getFrom();
      for (size_t st_current=0; st_current<state_size; ++st_current) {  //j is previous state
        if (!(*from_trans)[st_current])
          continue;

        if ((*scoring_previous)[st_previous] != -INFINITY) {
          backward_temp = (*scoring_previous)[st_previous] + emission + getTransition((*hmm)[st_current], st_previous , position+1);
          if ((*scoring_current)[st_current] == -INFINITY) {
            (*scoring_current)[st_current] = backward_temp;
          } else {
            (*scoring_current)[st_current] = AddInLogSpace(backward_temp, (*scoring_current)[st_current]);
          }
          next_states[st_current] = 1;
        }
      }
    }
  }

  for (size_t i=0; i<state_size; ++i) {
    (*posterior_score)[0][i] = ((double)(*posterior_score)[0][i] + (double)(*scoring_current)[i]) - ending_forward_prob;
    posterior_sum[0] = AddInLogSpace(posterior_sum[0], (*posterior_score)[0][i]);
  }

  ending_backward_prob = -INFINITY;
  init = hmm->getInitial();
  for(size_t i=0; i<state_size; ++i) {
    if ((*scoring_current)[i] != -INFINITY) {
      backward_temp = (*scoring_current)[i] + (*hmm)[i]->get_emission_prob(*seqs,0) + getTransition(init, i, 0);
      if (backward_temp > -INFINITY) {
        if (ending_backward_prob == -INFINITY) {
          ending_backward_prob = backward_temp;
        } else {
          ending_backward_prob = AddInLogSpace(backward_temp, ending_backward_prob);
        }
      }
    }
  }
              
  if (abs(ending_backward_prob - ending_forward_prob) > 0.0000001) {
    std::cerr << "Ending sequence probabilities calculated by Forward and Backward algorithm are different.  They should be the same.\t" << __FUNCTION__ << std::endl;
  }
              
  for (size_t i=0; i<seq_size; i++) {
    for(size_t j=0; j<state_size; j++) {
      if ((*posterior_score)[i][j] == -INFINITY)
        continue;
      if ((*posterior_score)[i][j] > -7.6009) {  //Above significant value;
        (*posterior_score)[i][j] -= posterior_sum[i];
      } else {
        (*posterior_score)[i][j] = -INFINITY;
      }
    }
  }
              
  delete scoring_previous;
  delete scoring_current;
  scoring_previous = NULL;
  scoring_current = NULL;
}
      
// ----------------------------------------------------------------------------------------
void trellis::traceback_posterior(traceback_path& path) {
  if (posterior_score == NULL) {
    std::cerr << __FUNCTION__ << " called before trellis::posterior was completed\n";
    exit(2);
  }
              
  double max(-INFINITY);
  int16_t max_ptr(-1);
  for(size_t position=seq_size-1; position != SIZE_MAX; position--) {
    max = -INFINITY;
    max_ptr = -1;
    for (size_t st=0; st < state_size; st++) {
      if ((*posterior_score)[position][st] > max) {
        max = (*posterior_score)[position][st];
        max_ptr = st;
      }
    }
    path.push_back(max_ptr);
  }
  return;
}
      
// ----------------------------------------------------------------------------------------
void trellis::traceback_stoch_posterior(traceback_path& path){
  for (size_t position =seq_size-1; position != SIZE_MAX; --position){
    double random=((double)rand()/((double)(RAND_MAX)+(double)(1)));
    double cumulative_prob(0.0);
          
    for (size_t st = 0; st < state_size ; ++st){
      cumulative_prob+=exp((*posterior_score)[position][st]);
      if (random < cumulative_prob){
        path.push_back( (int16_t) st);
      }
    }
  }
  return;
}
      
void trellis::traceback_stoch_posterior(multiTraceback& paths, size_t reps){
  for (size_t i = 0; i < reps; i++){
    traceback_path path(hmm);
    traceback_stoch_posterior(path);
    paths.assign(path);
  }
  return;
}
}
