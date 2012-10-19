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
#include <functional>
#include <assert.h>
//#include "rng.hh"

namespace smc
{
/// A template class for a set of moves for use in an SMC samplers framework.

template <class Space> class moveset
{
public:
    /// Callbacks used
    typedef std::function<particle<Space>(rng*)> init_fn;
    typedef std::function<long(long, const particle<Space>&, rng*)> move_select_fn;
    typedef std::function<int(long, particle<Space> &, rng*)> mcmc_fn;
    typedef std::function<void(long, particle<Space>&, rng*)> move_fn;

private:
    ///The function which initialises a particle.
    init_fn pfInitialise;
    ///The function which selects a move for a given particle at a given time.
    move_select_fn pfMoveSelect;
    ///The functions which perform actual moves.
    std::vector<move_fn> pfMoves;
    ///The Markov Chain Monte Carlo moves to use.
    std::vector<mcmc_fn> pfMCMC;

public:
    ///Create a completely unspecified moveset
    moveset();
    ///Create a reduced moveset with a single move moveset(init_fn pfInit,
    moveset(init_fn pfInit,
            move_fn newMoves,
            mcmc_fn pfNewMCMC);
    ///Create a fully specified moveset
    moveset(init_fn pfInit, move_select_fn pfMoveSelect, std::vector<move_fn> pfNewMoves, mcmc_fn pfNewMCMC);

    ///Initialise a particle
    particle<Space> DoInit(rng * pRng) { return pfInitialise(pRng);};
    ///Perform an MCMC move on a particle
    int DoMCMC(long lTime, particle<Space> & pFrom, rng* pRng);
    ///Select an appropriate move at time lTime and apply it to pFrom
    void DoMove(long lTime, particle<Space> & pFrom, rng * pRng);

    ///Free the memory used for the array of move pointers when deleting
    ~moveset();

    /// \brief Set the initialisation function.
    /// \param pfInit is a function which returns a particle generated according to the initial distribution
    void SetInitialisor(init_fn pfInit)
    {pfInitialise = pfInit;}

    /// \brief Set the MCMC function
    /// \param pfNewMCMC  The function which performs an MCMC move
    void SetMCMCFunction(mcmc_fn pfNewMCMC)
    { pfMCMC = std::vector<mcmc_fn>(1, pfNewMCMC);}

    /// \brief Set the MCMC function
    /// \param pfNewMCMC  Vector of functions which supply the new move
    void SetMCMCFunction(std::vector<mcmc_fn> pfNewMCMC)
    { pfMCMC = pfNewMCMC;}

    /// \brief Set the move selection function
    /// \param pfMoveSelectNew returns the index of move to perform at the specified time given a specified particle
    void SetMoveSelectionFunction(move_select_fn pfMoveSelectNew)
    {pfMoveSelect = pfMoveSelectNew; };

    ///Set the individual move functions to the supplied array of such functions
    void SetMoveFunctions(std::vector<move_fn> moves);

    ///Moveset assignment should allocate buffers and deep copy all members.
    moveset<Space> & operator= (moveset<Space> & pFrom);
};


/// The argument free smc::moveset constructor simply sets the number of available moves to zero and sets
/// all of the associated function pointers to NULL.
template <class Space>
moveset<Space>::moveset()
{
    pfInitialise = NULL;
    pfMoveSelect = NULL;
}

/// The three argument moveset constructor creates a new set of moves and sets all of the relevant function
/// pointers to the supplied values. Only a single move should exist if this constructor is used.
/// \param pfInit The function which should be used to initialise particles when the system is initialised
/// \param pfNewMoves An functions which moves a particle at a specified time to a new location
/// \param pfNewMCMC The function which should be called to apply an MCMC move (if any)
template <class Space>
moveset<Space>::moveset(init_fn pfInit,
                        move_fn newMoves,
                        mcmc_fn pfNewMCMC)
{
    SetInitialisor(pfInit);
    SetMoveSelectionFunction(NULL);
    SetMoveFunctions(std::vector<move_fn>(1, newMoves));
    assert(pfMoves.size() == 1);
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
moveset<Space>::moveset(init_fn pfInit, move_select_fn pfMoveSelector, std::vector<move_fn> pfNewMoves, mcmc_fn pfNewMCMC)
{
    SetInitialisor(pfInit);
    SetMoveSelectionFunction(pfMoveSelector);
    SetMoveFunctions(pfNewMoves);
    SetMCMCFunction(pfNewMCMC);
}

template <class Space>
moveset<Space>::~moveset()
{
}

template <class Space>
int moveset<Space>::DoMCMC(long lTime, particle<Space> & pFrom, rng *pRng)
{
    bool any_accepted = false;
    typename std::vector<mcmc_fn>::const_iterator fn = pfMCMC.begin(), end = pfMCMC.end();
    for(; fn != end; ++fn) {
        if((*fn)(lTime, pFrom, pRng)) {// move accepted
            any_accepted = true;
        }
    }

    return any_accepted;
}

template <class Space>
void moveset<Space>::DoMove(long lTime, particle<Space> & pFrom, rng *pRng)
{
    if(pfMoves.size() > 1)
        pfMoves[pfMoveSelect(lTime, pFrom, pRng)](lTime, pFrom, pRng);
    else
        pfMoves[0](lTime, pFrom, pRng);
}

/// \param nMoves The number of moves which are defined in general.
/// \param pfNewMoves An array of functions which move a particle at a specified time to a new location
///
/// The move functions accept two arguments, the first of which corresponds to the system evolution time and the
/// second to an initial particle position and the second to a weighted starting position. It returns a new
/// weighted position corresponding to the moved particle.
template <class Space>
void moveset<Space>::SetMoveFunctions(std::vector<move_fn> newMoves)
{
    pfMoves = std::vector<move_fn>(newMoves.begin(), newMoves.end());
    return;
}

template <class Space>
moveset<Space> & moveset<Space>::operator= (moveset<Space> & pFrom)
{
    SetInitialisor(pFrom.pfInitialise);
    SetMCMCFunction(pFrom.pfMCMC);
    SetMoveSelectionFunction(pFrom.pfMoveSelect);
    SetMoveFunctions(pFrom.pfMoves);

    return *this;
}
}
#endif
