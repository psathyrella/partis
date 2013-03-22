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

#include "mcmc.hh"
#include "particle.hh"
#include <cassert>
#include <functional>
#include <vector>
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
    typedef std::function<void(long, particle<Space>&, rng*)> move_fn;

private:
    ///The function which initialises a particle.
    init_fn pfInitialise;
    ///The function which selects a move for a given particle at a given time.
    move_select_fn pfMoveSelect;
    ///The functions which perform actual moves.
    std::vector<move_fn> pfMoves;
    ///The Markov Chain Monte Carlo moves to use.
    mcmc_moves<Space> pfMCMC;
    ///Number of MCMC moves to make
    std::size_t nMCMC;

public:
    ///Create a completely unspecified moveset
    moveset();
    ///Create a reduced moveset with a single move moveset(init_fn pfInit,
    moveset(init_fn pfInit,
            move_fn newMoves);
    ///Create a fully specified moveset
    moveset(init_fn pfInit, move_select_fn pfMoveSelect,
            std::vector<move_fn> pfNewMoves,
            mcmc_moves<Space> selector);

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

    /// \brief Set the MCMC function, sets the number of moves to pfNewMCMC.Count()
    /// \param pfNewMCMC  The function which performs an MCMC move
    void SetMCMCSelector(mcmc_moves<Space> pfNewMCMC)
    { pfMCMC = pfNewMCMC; SetNumberOfMCMCMoves(pfNewMCMC.Count()); }

    /// \brief Set the number of MCMC moves to make
    /// \param nMCMC Number of moves to make
    void SetNumberOfMCMCMoves(const size_t n)
    { this->nMCMC = n; };

    /// \brief Set the move selection function
    /// \param pfMoveSelectNew returns the index of move to perform at the specified time given a specified particle
    void SetMoveSelectionFunction(move_select_fn pfMoveSelectNew)
    {pfMoveSelect = pfMoveSelectNew; };

    ///Set the individual move functions to the supplied array of such functions
    void SetMoveFunctions(const std::vector<move_fn>& moves);

    ///Moveset assignment should allocate buffers and deep copy all members.
    moveset<Space> & operator= (moveset<Space> & pFrom);
};


/// The argument free smc::moveset constructor simply sets the number of available moves to zero and sets
/// all of the associated function pointers to NULL.
template <class Space>
moveset<Space>::moveset() :
    pfInitialise(nullptr),
    pfMoveSelect(nullptr),
    nMCMC(0)
{
}

/// The three argument moveset constructor creates a new set of moves and sets all of the relevant function
/// pointers to the supplied values.
/// \param pfInit The function which should be used to initialise particles when the system is initialised
/// \param pfNewMoves An functions which moves a particle at a specified time to a new location
template <class Space>
moveset<Space>::moveset(init_fn pfInit,
                        move_fn newMoves) :
    nMCMC(0)
{
    SetInitialisor(pfInit);
    SetMoveSelectionFunction(nullptr);
    SetMoveFunctions(std::vector<move_fn>(1, newMoves));
    assert(pfMoves.size() == 1);
}

/// The five argument moveset constructor creates a new set of moves and sets all of the relevant function
/// pointers to the supplied values.
/// \param pfInit The function which should be used to initialise particles when the system is initialised
/// \param pfMoveSelector The function which selects a move to apply, at a specified time, to a specified particle
/// \param nMoves The number of moves which are defined in general
/// \param pfNewMoves An array of functions which move a particle at a specified time to a new location
/// \param selector The function which should be called to apply an MCMC move (if any)
template <class Space>
moveset<Space>::moveset(init_fn pfInit, move_select_fn pfMoveSelector, std::vector<move_fn> pfNewMoves,
                        mcmc_moves<Space> selector)
    : nMCMC(selector.Count())
{
    SetInitialisor(pfInit);
    SetMoveSelectionFunction(pfMoveSelector);
    SetMoveFunctions(pfNewMoves);
    SetMCMCSelector(selector);
}

template <class Space>
moveset<Space>::~moveset()
{
}

template <class Space>
int moveset<Space>::DoMCMC(long lTime, particle<Space> & pFrom, rng *pRng)
{
    assert(pfMCMC.Count() > 0 || nMCMC == 0);
    bool any_accepted = false;
    auto moves = pfMCMC.SelectMoves(pRng, nMCMC);
    for(auto move_fn : moves) {
        if((*move_fn)(lTime, pFrom, pRng))
            any_accepted = true;
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
void moveset<Space>::SetMoveFunctions(const std::vector<move_fn>& newMoves)
{
    pfMoves = std::vector<move_fn>(newMoves.begin(), newMoves.end());
    return;
}

template <class Space>
moveset<Space> & moveset<Space>::operator= (moveset<Space> & pFrom)
{
    SetInitialisor(pFrom.pfInitialise);
    SetMCMCSelector(pFrom.pfMCMC);
    SetNumberOfMCMCMoves(pFrom.nMCMC);
    SetMoveSelectionFunction(pFrom.pfMoveSelect);
    SetMoveFunctions(pFrom.pfMoves);

    return *this;
}
}
#endif
