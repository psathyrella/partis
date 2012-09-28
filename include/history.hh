//   SMCTC: history.hh  
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
//! \brief Classes and function related to the history of the sampler.
//!
//! This file contains template definitions for the classes used to store the history of an SMCTC sampler.
//! It defines smc::history, smc::historyelement and smc::history.

#ifndef __SMC_HISTORY_HH
#define __SMC_HISTORY_HH 1.0

namespace smc {
  /// The historyflags class holds a set of flags which describe various properties of the particle system at a given time.
  class historyflags
  {
  private:
    /// true if the particle system was resampled during the described iteration.
    unsigned int Resampled : 1;
  public:
    ///Create a new set of history flags corresponding to the specified properties
    historyflags(int wasResampled);

    ///This function returns true if the flag set indicates that the ensemble was resampled during the described iteration.
    int WasResampled(void) {return Resampled;}
  };

  /// A template class for the elements of a linked list to be used for the history of the sampler.
  template <class Particle>class historyelement
  {
  private:
    long number; //!< The number of particles (presently redundant as this is not a function of iteration)
    int nAccepted; //!< Number of MCMC moves accepted during this iteration.
    Particle* value; //!< The particles themselves (values and weights)
    historyflags flags; //!< Flags associated with this iteration.
    historyelement<Particle> *pLast; //!< The parent of this node.
    historyelement<Particle> *pNext; //!< The child of this node.
    
  public:
    /// The null constructor creates an empty history element.
    historyelement();
    /// A constructor with four arguments initialises the particle set.
    historyelement(long lNumber, Particle* pNew, int nAccepts, historyflags hf);
    /// A constructor with six arguments initialises the whole structure.
    historyelement(long lNumber, Particle* pNew, historyelement<Particle>* pParent, historyelement<Particle>* pChild, int nAccepts, historyflags  hf);

    /// The destructor tidies up.
    ~historyelement();

    /// Returns the effective sample size of this particle generation.
    double GetESS(void) const;
    /// Returns the flags
    historyflags GetFlags(void) const{return flags;}
    /// Returns the parent of this element.
    historyelement<Particle> * GetLast(void) const { return pLast; }
    /// Returns the child of this element.
    historyelement<Particle> * GetNext(void) const { return pNext; }
    /// Returns the number of particles present.
    long GetNumber(void) const {return number;} 
    /// Returns a pointer to the current particle set.
    Particle * GetValues(void) const { return value; }
    /// Add a new history element with the specified values after the current one and maintain the list structure.
    void InsertAfter(long lNumber, Particle * pNew);
    /// Integrate the supplied function according to the empirical measure of the particle ensemble.
    long double Integrate(long lTime, double (*pIntegrand)(long,const Particle&,void*), void* pAuxiliary);
    /// Sets the elements parent.
    void SetLast(historyelement<Particle>* pParent) {pLast = pParent; }
    /// Sets the elements child.
    void SetNext(historyelement<Particle>* pChild) {pNext = pChild; }
    /// Sets the particle set to the specified values.    
    void SetValue(long lNumber, Particle * pNew);

    /// Returns the number of MCMC moves accepted during this iteration.
    int AcceptCount(void) {return nAccepted; }
    /// Returns true if the particle set 
    int WasResampled(void) {return flags.WasResampled(); }
    /// \brief The left shift operator returns the element a number of positions prior to this one.
    ///
    /// \param ulCount The number of positions by which to shift.
    /// \return The element a number of positions before this one.

    historyelement<Particle> operator<<(unsigned long ulCount)
    {
      if(ulCount)
	return this->pLast << (--ulCount);
	else
	  return *this;
    }

    ///\brief The right shift operator returns the element a number of positions after this one.
    ///
    /// \param ulCount The number of positions by which to shift.
    /// \return The right shift operator returns the element a number of positions after this one.
    historyelement<Particle> operator>>(unsigned long ulCount)
    {
      if(ulCount)
	return this->pNext << (--ulCount);
      else
	return *this;
    }

  };

  template <class Particle>
  historyelement<Particle>::historyelement()
  {
    number = 0;
    value = NULL;
    nAccepted = 0;
    pLast = NULL;
    pNext = NULL;
  }


  /// \param lNumber The number of particles present in the particle generation
  /// \param pNew    The array of particles which are present in the particle generation
  /// \param nAccepts The number of MCMC moves that were accepted during this particle generation
  /// \param hf      The historyflags associated with the particle generation

  template <class Particle>
  historyelement<Particle>::historyelement(long lNumber, Particle* pNew, int nAccepts, historyflags hf) :
    flags(hf)
  {
    number = lNumber;
    value = new Particle[number];
    for(int i = 0; i < number; i++)
      value[i] = pNew[i];

    nAccepted = nAccepts;
    pLast = NULL;
    pNext = NULL;
  }

  /// \param lNumber The number of particles present in the particle generation
  /// \param pNew    The array of particles which are present in the particle generation
  /// \param pParent A pointer to the previous particle generation
  /// \param pChild  A pointer to the next particle generation
  /// \param nAccepts The number of MCMC moves that were accepted during this particle generation
  /// \param hf      The historyflags associated with the particle generation
  template <class Particle>
  historyelement<Particle>::historyelement(long lNumber, Particle* pNew,
					   historyelement<Particle>* pParent, historyelement<Particle>* pChild,
					   int nAccepts, historyflags hf) :
    flags(hf)
  {
    number = lNumber;
    value = new Particle[number];
    for(int i = 0; i < number; i++)
      value[i] = pNew[i];

    nAccepted = nAccepts;
    pLast = pParent;
    pNext = pChild;
  }

  template <class Particle>
  historyelement<Particle>::~historyelement(void)
  {
    if(value)
      delete [] value;
  }

  template <class Particle>
  double historyelement<Particle>::GetESS(void) const
  {
    double sum = 0;
    double sumsq = 0;

    for(int i = 0; i < number; i++) {
      sum += exp(value[i].GetLogWeight());
      sumsq += exp(2.0*value[i].GetLogWeight());
    }
    return (sum*sum)/sumsq;
  }

  /// \param lTime The timestep at which the integration is to be carried out
  /// \param pIntegrand The function which is to be integrated
  /// \param pAuxiliary A pointer to additional information which is passed to the integrand function

  template <class Particle>
  long double historyelement<Particle>::Integrate(long lTime, double (*pIntegrand)(long,const Particle&,void*), void* pAuxiliary)
  {
    long double rValue = 0;
    long double wSum = 0;

    for(int i =0; i < number; i++)
      {
	rValue += expl(value[i].GetLogWeight()) * (long double)pIntegrand(lTime, value[i], pAuxiliary);
	wSum  += expl(value[i].GetLogWeight());
      }


    rValue /= wSum;

    return rValue;
  }

  /// \param lNumber The number of particles in the generation to be inserted
  /// \param pNew The value of the particle generation to be inserted

  template <class Particle>
  void historyelement<Particle>::InsertAfter(long lNumber, Particle * pNew)
  {
    pNext = new historyelement<Particle>(lNumber, pNew, this, pNext);    
  }

  /// A template class for the history associated with a particle system evolving in SMC.

  ///  The history is a template class which should have an associated class type corresponding to
  ///    a _particle_ of the desired type, not the type itself.
  ///
  ///    Essentially, this is implemented as a doubly linked list. 


  template <class Particle> class history
    {
    private:
      ///The length of the history in time steps
      long  lLength;
      ///The first time step
      historyelement<Particle> *pRoot;
      ///The most recent time step
      historyelement<Particle> *pLeaf;

    public:
      ///The argument free constructor creates an empty list.
      history();

      ///This function returns a pointer to the first element of the history.
      const historyelement<Particle > * GetElement(void) const {return pRoot; }

      /// Returns the effective sample size of the specified particle generation.
      double GetESS(long lGeneration) const;
      ///Returns the number of generations recorded within the history.
      long GetLength (void) const { return lLength; }
      ///Integrate the supplied function over the path of the particle ensemble.
      double IntegratePathSampling(double (*pIntegrand)(long,const Particle&,void*), double (*pWidth)(long,void*), void* pAuxiliary);
      double IntegratePathSamplingFinalStep(double (*pIntegrand)(long,const Particle&,void*), void* pAuxiliary) const;

      ///Output a vector indicating the number of accepted MCMC moves at each time instance
      void OstreamMCMCRecordToStream(std::ostream &os) const;
      ///Output a 0-1 value vector indicating the times at which resampling occured to an output stream
      void OstreamResamplingRecordToStream(std::ostream &os) const;

      ///Remove the terminal particle generation from the list and return that particle.
      Particle * Pop(void);
      ///Remove the terminal particle generation from the list and stick it in the supplied data structures
      void Pop(long* plNumber, Particle** ppNew, int* pnAccept, historyflags * phf);
      ///Append the supplied particle generation to the end of the list.
      void Push(long lNumber, Particle * pNew, int nAccept, historyflags hf);


      ///Display the list of particles in a human readable form.
      //  void StreamParticles(std::ostream & os);
    };

  /// This constructor simply sets the root and leaf pointers to NULL and the list length to zero.
  template <class Particle>
  history<Particle>::history()
  {
    pRoot = NULL;
    pLeaf = NULL;
    lLength = 0;
  }
   /// Returns the effective sample size of the specified particle generation.
  template <class Particle>
  double  history<Particle>::GetESS(long lGeneration) const
  {
    historyelement<Particle> * pCurrent = pRoot;
    for(long l = 0; l < lGeneration; l++)
      pCurrent = pCurrent->GetNext();
    return pRoot->GetESS();
  }



  /// This function records the MCMC acceptance history to the specified output stream as a list of
  /// the number of moves accepted at each time instant.
  ///
  /// \param os The output stream to send the data to.
  template <class Particle>
  void history<Particle>:: OstreamMCMCRecordToStream(std::ostream &os) const
  {
    historyelement<Particle> * pCurrent = pRoot;

    while(pCurrent) {
      os << pCurrent->AcceptCount() << std::endl;
      pCurrent = pCurrent->GetNext();
    }
  }
  /// This function records the resampling history to the specified output stream as a 0-1 valued list which takes
  /// the value 1 for those time instances when resampling occured and 0 otherwise.
  ///
  /// \param os The output stream to send the data to.
  template <class Particle>
  void history<Particle>:: OstreamResamplingRecordToStream(std::ostream &os) const
  {
    historyelement<Particle> * pCurrent = pRoot;

    while(pCurrent) {
      if(pCurrent->WasResampled())
	os << "1\t";
      else
	os << "0\t";

      os << pCurrent->GetESS() << std::endl;

      pCurrent = pCurrent->GetNext();
    }
  }

  /// This function performs a trapezoidal integration of the type which is useful when using path sampling to estimate the
  /// normalising constant of a potential function in those cases where a sequence of distributions is produced by deforming
  /// the initial distribution by a sequence of progressively more influential potential functions which are proportional
  /// to the density of some other distribution with respect to the starting distribution.
  ///
  /// The integrand is integrated at every time point in the particle system history. The results of this integration are
  /// taken to be point-evaluations of the path sampling integrand which are spaced on a grid of intervals given by the
  /// width function. The path sampling integral is then calculated by performing a suitable trapezoidal integration and
  /// the results of this integration is returned.
  ///
  /// pAuxiliary is passed to both of the user specified functions to allow the user to pass additional data to either or
  /// both of these functions in a convenient manner. It is safe to use NULL if no such data is used by either function.
  ///
  /// \param pIntegrand  The function to integrated. The first argument is evolution time, the second a particle at which the function is to be evaluated and the final argument is always pAuxiliary.
  /// \param pWidth      The function which returns the width of the path sampling grid at the specified evolution time. The final argument is always pAuxiliary
  /// \param pAuxiliary  A pointer to auxiliary data to pass to both of the above functions

  template <class Particle>
  double history<Particle>::IntegratePathSampling(double (*pIntegrand)(long,const Particle&,void*), double (*pWidth)(long,void*), void* pAuxiliary)
  {
    long lTime = 0;
    long double rValue = 0.0;
    
    historyelement<Particle> * pCurrent = pRoot;

    lTime = 1;
    pCurrent=pCurrent->GetNext();
    while(pCurrent) 
      {
	rValue += pCurrent->Integrate(lTime, pIntegrand, pAuxiliary) * (long double)pWidth(lTime,pAuxiliary);
	pCurrent = pCurrent->GetNext();
	lTime++; 
      }
    return ((double)rValue);
  }

  template <class Particle>
  double history<Particle>::IntegratePathSamplingFinalStep(double (*pIntegrand)(long,const Particle&,void*), void* pAuxiliary) const
  {
    return pLeaf->Integrate(lLength-1,pIntegrand,pAuxiliary);
  }


  /// Pop() operates just as the standard stack operation does. It removes the particle generation currently occupying
  /// the terminal position in the chain, decrements the length counter and returns the particle set as an array.
  template <class Particle>
  Particle * history<Particle>::Pop(void)
  {
    if(lLength == 0)
      return NULL;

    Particle * rValue = pLeaf->GetValues();

    lLength--;

    if(lLength == 0)
      pRoot = pLeaf = 0;
    else {
      pLeaf = pLeaf->GetLast();
      delete pLeaf->GetNext();
      pLeaf->SetNext(NULL);

    }
    return rValue;
  }

  /// Pop operates as the usual stack operation
  ///
  /// If called with four pointers of this sort, it removes the last particle generation from the history stack and
  /// places them in the reference objects.
  template <class Particle>
  void history<Particle>::Pop(long* plNumber, Particle** ppNew, int* pnAccept, historyflags * phf)
  {
    if(plNumber)
      (*plNumber) = pLeaf->GetNumber();
    if(ppNew) {
      for(long l = 0; l < *plNumber; l++)
      (*ppNew)[l]    = pLeaf->GetValues()[l];
    }
    if(pnAccept)
      (*pnAccept) = pLeaf->AcceptCount();
    if(phf)
      (*phf)      = pLeaf->GetFlags();

    if(lLength <= 1) {
      pRoot = NULL;
      pLeaf = NULL;
    }
    else {
      pLeaf = pLeaf->GetLast();
      //      delete pLeaf->GetNext();
      pLeaf->SetNext(NULL);

    }

    lLength--;

    return;

  }

  /// Push operates just like the standard stack operation: it adds the specified particle set generation to the history
  /// of the particle set and increments the stack counter.
  ///
  /// \param lNumber The number of particles present in this generation of the system.
  /// \param pNew    An array containing the particles present in this generation of the system.
  /// \param nAccepts The number of accepted MCMC moves during this iteration of the system
  /// \param hf      The historyflags associated with this generation of the system.

  template <class Particle>
  void history<Particle>::Push(long lNumber, Particle * pNew, int nAccepts, historyflags hf)
  {
    if(lLength == 0) {
      pRoot = new historyelement<Particle>(lNumber, pNew, nAccepts, hf);
      pLeaf = pRoot;
    }
    else {
      pLeaf = new historyelement<Particle>(lNumber, pNew, pLeaf, NULL, nAccepts, hf);
      pLeaf->GetLast()->SetNext(pLeaf);
    } 
    lLength++;
  }
}



namespace std {
  /// This function will ultimately allow the standard stream operators to be used to display a particle history in a human readable form.

  /// It isn't yet fully implemented.
  template <class Particle>
  ostream & operator<<(ostream & os, smc::history<Particle> h)
  {
    h.StreamParticles(os);
    return os;
  }
}
#endif
