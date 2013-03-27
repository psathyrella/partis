//   SMCTC: mcmc.hh
//
//   Copyright Connor McCoy, 2012
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
//! \brief Classes and functions which deal with collections of sampler proposal MCMC "moves".

#ifndef __SMC_MCMC_HH
#define __SMC_MCMC_HH 1.0

#include <cassert>
#include <functional>
#include <vector>

#include "particle.hh"
#include "rng.hh"

namespace smc
{

/// An MCMC move

template <typename Space>
class mcmc_moves
{
public:
    typedef std::function<int(long, particle<Space> &, rng*)> mcmc_fn;

    mcmc_moves() :
        uniform_weights(true) {};

    /// \brief Initialize with moves and weights
    mcmc_moves(const std::vector<mcmc_fn>& moves, const std::vector<double>& weights) :
        moves(moves),
        weights(weights),
        uniform_weights(AreWeightsUniform()) { };

    /// \brief Initialize with some moves
    ///
    /// All moves are weighted equally.
    mcmc_moves(const std::vector<mcmc_fn>& moves) :
        moves(moves),
        weights(std::vector<double>(moves.size(), 1.0)),
        uniform_weights(AreWeightsUniform()) { };

    /// Add an MCMC move with specied weight
    void AddMove(mcmc_fn move, double weight = 1.0);
    /// \brief Select a move randomly from a multinomial distribution using MCMC move weights as probabilities

    /// \returns A pointer to an MCMC move function
    mcmc_fn* SelectMove(smc::rng* rng);
    /// \brief Select moves randomly using weights as probabilities

    /// \returns A vector of functions to apply (may contain duplicates)
    std::vector<mcmc_fn*> SelectMoves(smc::rng* r, unsigned n);

    inline size_t Count() const { return moves.size(); };
private:
    bool AreWeightsUniform() const;

    std::vector<mcmc_fn> moves;
    std::vector<double> weights;
    bool uniform_weights;
};

// Implementation
template <typename Space>
void mcmc_moves<Space>::AddMove(mcmc_fn move, double weight)
{
    moves.push_back(move);
    weights.push_back(weight);
    uniform_weights = AreWeightsUniform();
}

template <typename Space>
std::vector<std::function<int(long, particle<Space> &, rng*)>*> mcmc_moves<Space>::SelectMoves(smc::rng* r, unsigned n)
{
    std::vector<mcmc_moves<Space>::mcmc_fn*> result;
    result.reserve(n);

    assert(moves.size() >= 0);
    assert(moves.size() == weights.size());
    if(moves.size() == 1) return std::vector<mcmc_moves<Space>::mcmc_fn*>(n, &(moves[0]));

    // Unequal weights: sample from multinomial
    std::vector<unsigned> counts(moves.size(), 0);
    r->Multinomial(n, moves.size(), weights.data(), counts.data());
    for(size_t i = 0; i < counts.size(); i++) {
        for(size_t j = 0; j < counts[i]; j++)
            result.push_back(&moves[i]);
    }
    assert(result.size() == n);
    return result;
}

template <typename Space>
std::function<int(long, particle<Space> &, rng*)>* mcmc_moves<Space>::SelectMove(smc::rng* rng)
{
    return SelectMoves(rng, 1)[0];
}

template <typename Space>
bool mcmc_moves<Space>::AreWeightsUniform() const
{
    if(weights.size() < 2) {
        return true;
    }
    double initial = weights[0];
    for(std::vector<double>::const_iterator it = weights.begin() + 1; it != weights.end(); ++it) {
        if(*it != initial)
            return false;
    }
    return true;
}

} // namespace smc

#endif
