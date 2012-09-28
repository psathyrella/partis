#include "smc-exception.hh"

//! \file
//! \brief The untemplated smc::exception class is implemented here.

namespace smc
{
  //! Generate an SMCTC Exception class with the specified initialisation.

  //! This constructor fills the four elements of the class with their specified values.
  //! It is used to allow a single-line command to create and throw an exception.
  //!
  //! \param szN The name of the source file generating the exception.
  //! \param lL The line in that file responsible for the exception.
  //! \param lC The numerical code identifying the exception.
  //! \param szM An textual explanation of the problem.
  exception::exception(char const * szN, long lL, long lC, char const * szM)
  {
    szFile = szN; 
    lLine = lL; 
    lCode = lC; 
    szMessage = szM;
  }
}

namespace std {
  ///Display a human-readable version of an SMC exception.

  /// \param os The stream to write to.
  /// \param e The exception class to display.
  /// \return os
  std::ostream & operator<< (std::ostream & os, smc::exception & e)
  {
    os << "SMC Exception: #" << e.lCode << endl;
    os << "Error occured in file " << e.szFile << " at line " << e.lLine << "." << endl;
    os << "Details: " << endl << '\t' << e.szMessage << endl;
    return os;
  }
}

