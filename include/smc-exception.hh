//   SMCTC: smc-exception.hh
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

//! \file
//! \brief This file defines exception-handling facilities.
//!
//! The smc::exception class, which is used for exception handling by SMCTC, is defined.

#ifndef __SMC_EXCEPT_HH
#define __SMC_EXCEPT_HH 1.0

#include <iostream>

///A macro which autocompletes the housekeeping components of an smc::exception
#define SMC_EXCEPTION(code,error) smc::exception(__FILE__, __LINE__, code, error)

///Exception thrown if a file cannot be accessed.
#define SMCX_FILE_NOT_FOUND 0x0020
///Exception thrown if the sampler attempts to access history data which wasn't stored.
#define SMCX_MISSING_HISTORY 0x0010
///Exception thrown if an attempt is made to instantiate a class of which a single instance is permitted more than once.
#define SMCX_MULTIPLE_INSTANTIATION 0x1000

namespace smc {
  ///SMC Exception class

  /// This class holds details of unrecoverable errors which occur within the SMCTC library.
  /// An instance of it is thrown whenever such an error occurs.
  class exception {
  public:
    char const * szFile; //!< The source file from which the code generating the exception was generated.
    long lLine;   //!< The line of that source file which generates the exception.
    long lCode;   //!< A numerical code indicating the nature of the exception generated.
    char const * szMessage; //!< A human-readable explanation of the cause of the exception.   

    exception(char const *, long, long, char const *);
  };
}

namespace std {
  /// Produce a human-readable display of the state of an smc::exception class using the stream operator.

  /// Produce a human-readable display of the state of an smc::exception class using the stream operator.
  /// \param os The output stream to which the display should be made.
  /// \param e  The exception which is to be displayed.
  std::ostream & operator<< (std::ostream &, smc::exception &);
}

#endif
