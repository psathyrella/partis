//   SMCTC: moveset.hh  
//
//   Copyright Adam Johansen, 2008.
// 
//   This file is part of SMCTC.
//
//   SMCTC is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   SMCTC is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with SMCTC.  If not, see <http://www.gnu.org/licenses/>.
//

//! \file
//! \brief Classes and functions which deal with collections of sampler proposal "moves".
//!
//! This file contains definitions of smc::moveset.
//! It deals with the collections of proposal moves (including initialisation and MCMC moves) which must be dealt with by the sampler.

#ifndef __SMC_MOVESET_HH
#define __SMC_MOVESET_HH 1.0

#include "particle.hh"
//#include "rng.hh"

namespace smc {
  /// A template class for a set of moves for use in an SMC samplers framework.

  template <class Space> class moveset
    {
    private:
      ///The number of moves which are present in the set
      long number;
      ///The function which initialises a particle.
      particle<Space> (*pfInitialise)(rng*);
      ///The function which selects a move for a given particle at a given time.
      long (*pfMoveSelect)(long , const particle<Space> &, rng*);
      ///The functions which perform actual moves.
      void (**pfMoves)(long, particle<Space> &, rng*);
      ///A Markov Chain Monte Carlo move.
      int (*pfMCMC)(long,particle<Space> &, rng*);

    public:
      ///Create a completely unspecified moveset
      moveset();
      ///Create a reduced moveset with a single move
      moveset(particle<Space> (*pfInit)(rng*),
	      void (*pfNewMoves)(long, particle<Space> &,rng*),
	      int (*pfNewMCMC)(long,particle<Space> &,rng*));
      ///Create a fully specified moveset
      moveset(particle<Space> (*pfInit)(rng*),long (*pfMoveSelector)(long , const particle<Space> &,rng*), 
	      long nMoves, void (**pfNewMoves)(long, particle<Space> &,rng*),
	      int (*pfNewMCMC)(long,particle<Space> &,rng*));
      
      ///Initialise a particle
      particle<Space> DoInit(rng * pRng) { return (*pfInitialise)(pRng);};
      ///Perform an MCMC move on a particle
      int DoMCMC(long lTime, particle<Space> & pFrom, rng* pRng);
      ///Select an appropriate move at time lTime and apply it to pFrom
      void DoMove(long lTime, particle<Space> & pFrom,rng * pRng);
      
      ///Free the memory used for the array of move pointers when deleting
      ~moveset();

      /// \brief Set the initialisation function.
      /// \param pfInit is a function which returns a particle generated according to the initial distribution 
      void SetInitialisor( particle<Space> (*pfInit)(rng*) )
      {pfInitialise = pfInit;}

      /// \brief Set the MCMC function
      /// \param pfNewMCMC  The function which performs an MCMC move
      void SetMCMCFunction(int (*pfNewMCMC)(long,particle<Space> &,rng*)) { pfMCMC = pfNewMCMC;}

      /// \brief Set the move selection function
      /// \param pfMoveSelectNew returns the index of move to perform at the specified time given a specified particle
      void SetMoveSelectionFunction(long (*pfMoveSelectNew)(long , const particle<Space> &, rng*))
	{pfMoveSelect = pfMoveSelectNew; };

      ///Set the individual move functions to the supplied array of such functions
      void SetMoveFunctions(long nMoves, void (**pfNewMoves)(long, particle<Space> &, rng*));
      
      ///Moveset assignment should allocate buffers and deep copy all members.
      moveset<Space> & operator= (moveset<Space> & pFrom);
    };


  /// The argument free smc::moveset constructor simply sets the number of available moves to zero and sets
  /// all of the associated function pointers to NULL.
  template <class Space>
  moveset<Space>::moveset()
  {
    number = 0;

    pfInitialise = NULL;
    pfMoveSelect = NULL;
    pfMoves = NULL;
    pfMCMC = NULL;
  }

  /// The three argument moveset constructor creates a new set of moves and sets all of the relevant function
  /// pointers to the supplied values. Only a single move should exist if this constructor is used.
  /// \param pfInit The function which should be used to initialise particles when the system is initialised
  /// \param pfNewMoves An functions which moves a particle at a specified time to a new location
  /// \param pfNewMCMC The function which should be called to apply an MCMC move (if any)
  template <class Space>
  moveset<Space>::moveset(particle<Space> (*pfInit)(rng*),
			  void (*pfNewMoves)(long, particle<Space> &,rng*),
			  int (*pfNewMCMC)(long,particle<Space> &,rng*))
  {
    SetInitialisor(pfInit);
    SetMoveSelectionFunction(NULL);
    pfMoves = NULL; //This presents an erroneous deletion by the following call
    SetMoveFunctions(1, &pfNewMoves);
    SetMCMCFunction(pfNewMCMC);
  }

  /// The five argument moveset constructor creates a new set of moves and sets all of the relevant function
  /// pointers to the supplied values.
  /// \param pfInit The function which should be used to initialise particles when the system is initialised
  /// \param pfMoveSelector The function which selects a move to apply, at a specified time, to a specified particle
  /// \param nMoves The number of moves which are defined in general
  /// \param pfNewMoves An array of functions which move a particle at a specified time to a new location
  /// \param pfNewMCMC The function which should be called to apply an MCMC move (if any)
  template <class Space>
  moveset<Space>::moveset(particle<Space> (*pfInit)(rng*),long (*pfMoveSelector)(long ,const particle<Space> &,rng*), 
			  long nMoves, void (**pfNewMoves)(long, particle<Space> &,rng*),
			  int (*pfNewMCMC)(long,particle<Space> &,rng*))
  {
    SetInitialisor(pfInit);
    SetMoveSelectionFunction(pfMoveSelector);
    pfMoves = NULL; //This presents an erroneous deletion by the following call
    SetMoveFunctions(nMoves, pfNewMoves);
    SetMCMCFunction(pfNewMCMC);
  }

  template <class Space>
  moveset<Space>::~moveset()
  {
    if(pfMoves)
      delete [] pfMoves;
  }

  template <class Space>
  int moveset<Space>::DoMCMC(long lTime, particle<Space> & pFrom, rng *pRng)
  {
    if(pfMCMC) {
      return pfMCMC(lTime,pFrom,pRng);
    }
    else {
      return false;
    }
  }
  
  template <class Space>
  void moveset<Space>::DoMove(long lTime, particle<Space> & pFrom, rng *pRng)
    {
      if(number > 1)
	pfMoves[pfMoveSelect(lTime,pFrom,pRng)](lTime,pFrom,pRng);
      else
	pfMoves[0](lTime,pFrom,pRng);
    }

  /// \param nMoves The number of moves which are defined in general.
  /// \param pfNewMoves An array of functions which move a particle at a specified time to a new location
  ///
  /// The move functions accept two arguments, the first of which corresponds to the system evolution time and the
  /// second to an initial particle position and the second to a weighted starting position. It returns a new 
  /// weighted position corresponding to the moved particle.
  template <class Space>
   void moveset<Space>::SetMoveFunctions(long nMoves,void (**pfNewMoves)(long,particle<Space> &,rng*))
    {
      number = nMoves;
      if(pfMoves)
	delete [] pfMoves;
      pfMoves = (void (**)(long int, smc::particle<Space> &,rng*)) new void* [nMoves];
      for(int i = 0; i < nMoves; i++)
	pfMoves[i] = pfNewMoves[i];
      return;
    }

  template <class Space>
  moveset<Space> & moveset<Space>::operator= (moveset<Space> & pFrom)
  {
    SetInitialisor(pFrom.pfInitialise);
    SetMCMCFunction(pFrom.pfMCMC);
    SetMoveSelectionFunction(pFrom.pfMoveSelect);
    SetMoveFunctions(pFrom.number, pFrom.pfMoves);        
    
    return *this;
  }
}
#endif
