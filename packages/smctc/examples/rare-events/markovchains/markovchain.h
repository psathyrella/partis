#ifndef __MARKOVCHAIN_H
#define __MARKOVCHAIN_H

#include <cstdlib>
#include <iostream>

using namespace std;

//!A template class for the elements of a Markov chain.
template <class Space> class mElement
{
 public:
  //! Pointer to the previous element in the chain.
  mElement *pPrev;
  //! Pointer to the next element in the chain.
  mElement *pNext;
  //! Value of this element of the chain.
  Space    value;

  mElement();
  mElement(Space);
  mElement(mElement<Space>*,Space,mElement<Space>*);

  mElement<Space> operator=(mElement<Space> & mE) { this->value = mE.value; }
  mElement<Space> operator=(mElement<Space> mE) { this->value = mE.value; }
};


template <class Space>
mElement<Space>::mElement()
{
  pPrev = NULL;
  pNext = NULL;
}

template <class Space>
mElement<Space>::mElement(Space sInit)
{
  value = sInit;

  pPrev = NULL;
  pNext = NULL;
}

template <class Space>
mElement<Space>::mElement(mElement<Space>* pPInit, Space sInit, mElement<Space>* pNInit)
{
  value = sInit;
 
  pPrev = pPInit;
  pNext = pNInit;
}

//!A template class for a Markov Chain itself.
template <class Space> class mChain
{
 private:
  long nLength;
  mElement<Space> *pRoot;
  mElement<Space> *pLeaf;

  Space (*pfKernel) (Space);
 public:
  mChain();
  mChain(long nInit, Space* pElements);
  mChain(mChain<Space> *pMci);
  mChain(const mChain<Space> &mci);

  ~mChain();

  mChain<Space> operator= (const mChain<Space>&);
  mChain<Space> operator+ (Space const &) const;
  mChain<Space>& operator+=(Space const &);
  mChain<Space> operator- (Space const &) const;
  mChain<Space>& operator-=(Space const &);
  const Space & AddToElement(long nIndex, const Space &);
  int AppendElement(const Space & );
  int CheckSanity(void) const;
  int DeleteElement(long);
  int DeleteRange(long,long);
  void Empty(void);
  int KernelElement(void);
  int KernelRange(long);
  int InsertElement(long, Space);
  int PrependElement(Space);

  mElement<Space> *GetElement(long ind) const {if(ind > nLength) return NULL; mElement<Space> *pCount = pRoot; for(int i=0; i < ind; i++) pCount = pCount->pNext; return pCount; }
  mElement<Space> const *GetTerminal(void) const {return pLeaf;} 
  long GetLength(void) const { return nLength; }

  Space (*GetKernel(void))(Space) { return pfKernel; }
  void SetKernel(Space (*pfKernelNew) (Space)) { pfKernel = pfKernelNew; }
};

//! Constructor for an empty Markov Chain.
template <class Space> 
mChain<Space>::mChain()
{
    nLength = 0; 
    pRoot = NULL;
    pLeaf = NULL;

    pfKernel = NULL;

}

//! Constructor for an initialised Markov Chain.
template <class Space>
mChain<Space>::mChain(long nInit, Space* pElements)
{
  nLength = nInit;

  if(nInit) {
    pRoot = new mElement<Space>(pElements[0]);   
    
    mElement<Space> *pCurrent = pRoot;
    
    for(int i = 1; i < nInit; i++) {
      pCurrent->pNext = new mElement<Space>(pCurrent, pElements[i], NULL);
      pCurrent = pCurrent->pNext;
    }

    pLeaf = pCurrent;
  }
  else {
    pRoot = NULL;
    pLeaf = NULL;
  }

  pfKernel = NULL;
}

//! Constructor for an initialisd Markov Chain
template <class Space>
mChain<Space>::mChain(mChain<Space> *pMci)
{
  nLength = pMci->mLength;
  if(pMci->nLength == 0) {
    //    nLength = 0;
    pRoot = pLeaf = NULL;
  }
  else if(pMci->nLength == 1) {
    pRoot = new mElement<Space>(NULL, pMci->GetElement(0)->value,NULL);
    pLeaf = pRoot;
    // nLength = 1;
  }
  else {
    mElement<Space> * mEN = pMci->GetElement(0);
    pRoot = new mElement<Space>(NULL,mEN->value,NULL);
    mElement<Space> * mCurrent = pRoot;
    for(int i = 1; i < nLength; i++) {
      mCurrent->pNext = new mElement<Space>(mCurrent, mEN->value, NULL);
      mEN = mEN->pNext;
      mCurrent = mCurrent->pNext;
    }    
    pLeaf = mCurrent;
  }

}

template <class Space>
mChain<Space>::mChain(const mChain<Space> &Mci)
{
  /*
   nLength = 0;
   pRoot = pLeaf = NULL;
  
   for(int i = 0; i < Mci.nLength; i++) 
     AppendElement(Mci.GetElement(i)->value);
  */

  nLength = Mci.nLength;
  if(Mci.nLength == 0) {
    pRoot = pLeaf = NULL;
  }
  else  {
    if(Mci.nLength == 1) {
      pRoot = new mElement<Space>(NULL, Mci.GetElement(0)->value,NULL);
      pLeaf = pRoot;
    }
    else {
      mElement<Space> * mEN = Mci.GetElement(0);
      pRoot = new mElement<Space>(NULL,mEN->value,NULL);
      mElement<Space> * mCurrent = pRoot;
      for(int i = 1; i < nLength; i++) {
	mCurrent->pNext = new mElement<Space>(mCurrent, mEN->pNext->value, NULL);
	mEN = mEN->pNext;
	mCurrent = mCurrent->pNext;
      }    
      pLeaf = mCurrent;
    }
  }

  pfKernel = Mci.pfKernel;
}

template <class Space>
mChain<Space>::~mChain()
{
  //It's wise to have something here or things leak somewhat....

 if(nLength > 0) {
    mElement<Space>* pCurrent = pLeaf;
    while(pCurrent->pPrev) {
      pCurrent = pCurrent->pPrev;
      pCurrent->pNext->pPrev = NULL;
      delete pCurrent->pNext;
      pCurrent->pNext=NULL;
    }
    delete pCurrent;
  }
  pRoot=pLeaf=NULL;
  nLength = 0;
}

/// Add the specified value to the specified element of the Markov Chain

template <class Space>
const Space & mChain<Space>::AddToElement(long nElement, const Space & sAdd)
{
  mElement<Space>* pCurrent = pRoot;

  for(int n = 0; n < nElement; n++)
    {
      pCurrent=pCurrent->pNext;
      if(!pCurrent)
	return sAdd;
    }

  return pCurrent->value = pCurrent->value + sAdd;
}

//! Add a specified element to the end of a Markov Chain.
template <class Space>
int inline mChain<Space>::AppendElement(const Space &snew)
{
  if(nLength) {
    pLeaf->pNext = new mElement<Space>(pLeaf,snew,NULL);
    pLeaf = pLeaf->pNext;

    nLength++;

    return 0;
  }
  else {
    pRoot = new mElement<Space>(NULL,snew,NULL);
    pLeaf = pRoot;
    nLength  = 1;
    return 0;
  }
}

//! Check that the doubly-linked list is sane
template <class Space>
int mChain<Space>::CheckSanity(void) const
{
  int nLengthCheck = 0, nLengthCheck2 = 0;
  int nLinkFails   = 0;

  mElement<Space> *pHead = pRoot;

  while(pHead) {
    nLengthCheck++;
    
    if(pHead->pNext) {
      if(pHead->pNext->pPrev != pHead)
	nLinkFails++;
    }

    pHead=pHead->pNext;
  }

  pHead = pLeaf;
  while(pHead) {
    nLengthCheck2++;
    
    if(pHead->pPrev) {
      if(pHead->pPrev->pNext != pHead)
	nLinkFails++;
    }

    pHead=pHead->pPrev;
  }

  if(nLengthCheck != nLength ||  nLengthCheck2 != nLength)
    cerr << "LengthCheck " << nLengthCheck << "," << nLengthCheck2 << "(" << nLength << ")" << endl;
  if(nLinkFails)
    cerr << "Link Failure count..." << nLinkFails << endl;

  if(nLengthCheck != nLength)
    return -1;
  return nLinkFails;
}

//! Remove an element from a specified position in a Markov Chain.
template <class Space>
int mChain<Space>::DeleteElement(long nPos)
{
  if(nPos < 0 || nPos > nLength)
    return 1;
  
  if(nPos == 0) {
    mElement<Space> *pRootN = pRoot->pNext;
    delete pRoot;
    pRoot=pRootN;
    if(pRoot)
      pRoot->pPrev = NULL;

    nLength--;

    if(nLength < 2)
      pLeaf = pRoot;

    return 0;
  }
  if(nPos == nLength) {
    mElement<Space> *pLeafN = pLeaf->pPrev;
    delete pLeaf;
    pLeaf=pLeafN;
    if(pLeaf)
      pLeaf->pNext = NULL;

    nLength--;

    if(nLength < 2)
      pRoot = pLeaf;

    return 0;
  }

   mElement<Space> *pPoint = pRoot, *pPointN, *pPointP;
   for(int i = 1; i < nPos; i++)
     pPoint = pPoint->pNext;
   pPointP = pPoint->pPrev;
   pPointN = pPoint->pNext;

   delete pPoint;
   pPointP->pNext = pPointN;
   pPointN->pPrev = pPointP;
   
   nLength--;

   return 0;
}

//! Delete a contiguous region from within a Markov Chain.
template <class Space>
int mChain<Space>::DeleteRange(long nPos1, long nPos2)
{
  for(int i = nPos1; i < nPos2; i++)
    if(DeleteElement(nPos1))
      return 1;

  return 0;


  //This should be replaced wite an efficient implementation at some point!
}

//! Remove all of the elements from a Markov chain, leaving only an empty container
template <class Space>
void inline mChain<Space>::Empty(void)
{
 if(nLength > 0) {
    mElement<Space>* pCurrent = pLeaf;
    while(pCurrent->pPrev) {
      pCurrent = pCurrent->pPrev;
      pCurrent->pNext->pPrev = NULL;
      delete pCurrent->pNext;
      pCurrent->pNext=NULL;
    }

    delete pCurrent;
  }
  pRoot=pLeaf=NULL;
  nLength = 0;
}

//! Insert a specified element at a specified position in a Markov Chain.
template <class Space>
int  mChain<Space>::InsertElement(long nPos, Space snew)
{
  //If the insertion position is zero then insert the new element at the root of the chain
  if(nPos == 0) {
    return PrependElement(snew);
  }
  if(nPos == nLength)
    return AppendElement(snew);

  if(nPos > 0 && nPos < nLength) {
    mElement<Space> *pPoint = pRoot, *pPointN;
    for(int i = 1; i < nPos; i++)
      pPoint = pPoint->pNext;
    pPointN = pPoint->pNext;
    pPoint->pNext = new mElement<Space>(pPoint,snew,pPointN);
    pPointN->pPrev= pPoint->pNext;

    nLength++;
    return 0;
  }
  
  //If the insertion position does not lie between 0 and the length of the chain then we have a problem
  return 1;
}

//! Extend the chain by a single element obtained from the associated Markov Kernel
template <class Space>
int mChain<Space>::KernelElement(void)
{
  if(nLength) {
    pLeaf->pNext = new mElement<Space>(pLeaf, (*pfKernel)(pLeaf->value),NULL);
    pLeaf = pLeaf->pNext;

    nLength++;
    return 0;
  }
  else {
    return 1; //We can't very well extend a non-existant chain according to a kernel
  }
}

//! Extend the chain by a specified number of elements obtained by iterative application of the associated Markov Kernel.
template <class Space>
int mChain<Space>::KernelRange(long n)
{
  for(int i = 0; i < n; i++) {
    if(KernelElement())
      return 1+i;
  }

  return 0;    
}

//! Insert a specified element at the beginning of a Markov Chain.
template <class Space>
int mChain<Space>::PrependElement(Space snew)
{
  if(nLength) {
    pRoot->pPrev = new mElement<Space>(NULL,snew,pRoot);
    pRoot=pRoot->pPrev;

    nLength++;
    return 0;
  }
  else {
    pRoot = new mElement<Space>(NULL,snew,NULL);
    pLeaf = pRoot;
    nLength = 1;
    return 0;
  }
}


//Overloaded operators for use with these template classes

//! Assigning a Markov Chain performs a deep copy
/*
template <class Space> 
mChain<Space> & mChain<Space>::operator= (mChain<Space> & mci)
{
 nLength = mci.nLength;
  
 //  mChain<Space> *that = new mChain<Space>(&mci);
  //
  mElement<Space> *pi, *po;

  pi = mci.pRoot;
  
  Empty();

  pRoot = new mElement<Space>(NULL,pi->value,NULL);
  po = pRoot;
  //  po->value = pi->value;

  while(pi->pNext) {
    po->pNext = new mElement<Space>(po,pi->pNext->value,NULL);
    //    po->pNext->value = pi->pNext->value;
    //    po->pNext->pPrev = po;
    po=po->pNext;
    pi=pi->pNext;
  }
  pLeaf = pi;

  return (*this);

  //return *that;
}
*/
  //! Assigning a Markov Chain performs a deep copy
template <class Space> 
mChain<Space> mChain<Space>::operator= (const mChain<Space> &mci)
{
  // First, get rid of the old stuff
  Empty();

  // Now assign some new stuff

 nLength = mci.nLength;

 if(!nLength)
   return *this;
  
 //  mChain<Space> *that = new mChain<Space>(&mci);
  //
  mElement<Space> *pi, *po;

  pi = mci.pRoot;

  pRoot = new mElement<Space>(NULL,pi->value,NULL);
  po = pRoot;

  while(pi->pNext) {
    pi=pi->pNext;
    po->pNext = new mElement<Space>(po,pi->value,NULL);
    po=po->pNext;
  }
  pLeaf = po;

  return (*this);

  //return *that;
}

// Provide shortcuts for shifting the entire chain
template <class Space>
mChain<Space> mChain<Space>::operator+ (Space const & sInc) const
{
  static mChain<Space> that; // = new mChain<Space>;

  that = *this;

  mElement<Space> *pCurrent;
  pCurrent = that.pRoot;
  do {
    pCurrent->value = pCurrent->value + sInc;
    pCurrent = pCurrent->pNext;
  } while (pCurrent);
  
  return that;
}

template <class Space>
mChain<Space>& mChain<Space>::operator+=(Space const & sInc)
{
  *this = *this + sInc;

  return *this;
}

template <class Space>
mChain<Space> mChain<Space>::operator- (Space const & sInc) const
{
  static  mChain<Space> that; // = new mChain<Space>;

  that = *this;

  mElement<Space> *pCurrent;
  pCurrent = that.pRoot;
  do {
    pCurrent->value = pCurrent->value - sInc;
    pCurrent = pCurrent->pNext;
  } while (pCurrent);

  return that;
}

template <class Space>
mChain<Space>& mChain<Space>::operator-=(Space const & sInc)
{
  *this = *this - sInc;
  
  return *this;
}
//Overloaded stream operators for use with these template classes
//! Allow the use of output streams with Markov Chains.
template <class Space>
ostream & operator<< (ostream& os, mChain<Space>& mc)
{
  if(!mc.GetLength())
    return os;

  mElement<Space> *me = mc.GetElement(0);

  os << "[" << me->value;
  while(me->pNext) {
    me = me->pNext;
    os << "," << me->value;
  }
  os << "]";
  return os;
}

#endif
