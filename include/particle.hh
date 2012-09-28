//   SMCTC: particle.hh
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
//! \brief Class used to store and manipulate a single particle.
//!
//! This file contains the smc::particle class which is used internally and passed to move functions.

#ifndef __SMC_PARTICLE_HH
#define __SMC_PARTICLE_HH 1.0

#include <float.h>
#include <limits>
#include <cmath>

namespace smc {
  /// A template class for the particles of an SMC algorithm
  template <class Space> class particle
    {
    private:
      /// Value of this particle
      Space    value;
      /// Natural logarithm of this particle's weight.
      double   logweight;

    public:
      particle();
      /// Constructor which initialises the particles value and weight.
      particle(Space sInit,double dLogWeight);
      /// The copy constructor performs a shallow copy.
      particle(const particle<Space> & pFrom);
      /// The assignment operator performs a shallow copy.
      particle<Space> & operator= (const particle<Space> & pFrom);

      ~particle();

      /// Returns the particle's value 
      Space const & GetValue(void) const {return value;}
      /// Returns a pointer to the value to allow for more efficient changes
      Space* GetValuePointer(void) {return &value;}
      /// Returns the particle's log weight.
      double GetLogWeight(void) const {return logweight;}
      /// Returns the particle's unnormalised weight.
      double GetWeight(void) const {return exp(logweight);}
      
      /// \brief Sets the particle's value and weight explicitly
      ///
      /// \param sValue The particle value to use 
      /// \param dLogWeight The natural logarithm of the new particle weight
      void Set(Space sValue,double dLogWeight){value = sValue; logweight = dLogWeight;}
      /// \brief Sets the particle's value explicitly
      ///
      /// \param sValue The particle value to use
      void SetValue(const Space & sValue){value = sValue;}
      /// \brief Sets the particle's log weight explicitly
      ///
      /// \param dLogWeight The natural logarithm of the new particle weight
      void SetLogWeight(const double & dLogWeight) {logweight = dLogWeight;}
      /// \brief Sets the particles weight explicitly
      ///
      /// \param dWeight The new particle weight
      void SetWeight(const double & dWeight) {logweight = log(dWeight);}

      /// \brief Increase the log weight by a specified amount
      ///
      /// \param dIncrement The amount to add to the log weight.
      void AddToLogWeight(double dIncrement) { logweight += dIncrement;}
      /// \brief Multiply the weight by a specified factor
      ///
      /// \param dMultiplier The factor to multiply the weight by.
      void MultiplyWeightBy(double dMultiplier) { logweight += log(dMultiplier);}
  };


  /// Create a particle with undefined value and weight NAN
  template <class Space>
    particle<Space>::particle()
    {
        logweight = std::numeric_limits<double>::quiet_NaN();
    }

  ///Copy constructor
  template <class Space>
  particle<Space>::particle(const particle<Space> & pFrom)
  {
    value = pFrom.value;
    logweight = pFrom.logweight;
  }
  
  /// Create a particle with value sInit and log weight dLogWeight 
  /// \param sInit The initial value of the particle
  /// \param dLogWeight The initial value of the natural logarithm of the particle weight
  template <class Space>
    particle<Space>::particle(Space sInit, double dLogWeight)
    {
      value = sInit;
      logweight =dLogWeight;
    }

  /// Dispose of a particle which is no longer required
  template <class Space>
  particle<Space>::~particle()
  {
  }

  /// Copy the values of pFrom to the values of this to set this particle identical to pFrom in a deep
  /// copy sense.
  template <class Space>
  particle<Space> & particle<Space>::operator= (const particle<Space> & pFrom)
  {
    this->value = pFrom.value;
    this->logweight = pFrom.logweight;

    return *this;
  }
}

namespace std {
  /// Produce a human readable display of an smc::particle class using the standard stream operators

  /// \param os The output stream to which the display should be made.
  /// \param p  The particle which is to be displayed.
  template <class Space>
  std::ostream & operator << (std::ostream & os, smc::particle<Space> & p)
  {
    Space val = p.GetValue();
     os << val << "," << exp(p.GetLogWeight());
    return os;
  }
}
#endif
