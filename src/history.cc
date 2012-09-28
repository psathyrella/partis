#include "smctc.hh"

//! \file
//! \brief This file contains the untemplated functions used for storing the history of the system.

namespace smc {
  /// This constructor produces an initialised historyflags instance.
  ///
  /// \param wasResampled An indicator which should be nonzero if the particle 
  /// system was resampled during the iteration being described

  historyflags::historyflags(int wasResampled)
  {
    if(wasResampled)
      Resampled = 1;
    else
      Resampled = 0;
  }
}
