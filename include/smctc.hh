//
//     smctc.h
//
//     The principle header file for the SMC template class.
//
//     Copyright Adam Johansen, 2008
//
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

/// \mainpage smctc -- A Template Class Library for SMC Simulation
///  
/// \version 1.0
/// \author Adam M. Johansen
///
/// \section intro_sec Summary
///
/// The SMC template class library (SMCTC) is intended to be used to implement SMC 
/// algorithms ranging from simple particle filter to complicated SMC algorithms
/// from simple particle filters to the SMC samplers of Del Moral, Doucet and Jasra (2006) 
/// within a generic framework.
///
/// \section license License
///
/// SMCTC is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// SMCTC is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with SMCTC.  If not, see <http://www.gnu.org/licenses/>.
/// 
/// The software is still in development but is thought to be sufficiently fully-featured and
/// stable for use.
///
/// \section using_sec Using the Template Library
///
/// In order to use the template libary it is necessary to place the set of header files in 
/// a suitable include directory and then to instantiate the various classes with the 
/// appropriate type for your sampler.
///
/// For implementation details, please see the accompanying user guide.
/// 
///
/// \section example_sec Examples
///
/// \subsection A Simple Particle Filter
///
/// This example provides a simple particle filter for the approximate filtering of
/// nonlinear, nongaussian state space models.
///
/// \subsection SMC Samplers for Rare Event Simulation
/// 
/// The main.cc file provides an example of a use of the template library to estimate 
/// Gaussian tail probabilities within an SMC samplers framework by making use of a sequence 
/// of intermediate distributions defined by introducing a potential function which
/// varies from flat to being very close to an indicator function on the tail event.
///
/// This example is taken from from Johansen, Del Moral and Doucet (2006).
///
/// \subsection Acknowledgements
///
/// Thanks are due to Edmund Jackson and Mark Briers for testing SMCTC on a variety of platforms.
/// The Visual C++ project and solution files are kindly provided by Mark Briers.
///

//! \file
//! \brief The main header file for SMCTC.
//!
//! This file serves as an interface between user-space programs and the SMCTC library.
//! This is the only header file which applications need to include to make use of the library (it includes
//! such additional files as are necessary).

#ifndef __SMC_TDSMC_HH

#define __SMC_TDSMC_HH 1.0

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <gsl/gsl_rng.h>

#include "smc-exception.hh"
#include "sampler.hh"

/// The Sequential Monte Carlo namespace

///
///    The classes and functions within this namespace are intended to be used for producing
///    implemenetations of SMC samplers and related simulation techniques.

namespace smc {}

/// The standard namespace 

///     The classes provided within the standard libraries reside within this namespace and 
///     the TDSMC class library adds a number of additional operator overloads to some of
///     the standard classes to allow them to deal with our structures.

namespace std {}
#endif
