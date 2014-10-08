//
// File: MixedSubstitutionModelSet.h
// Created by: Laurent Guéguen
// Created on: mercredi 25 mai 2011, à 22h 04
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _MIXEDSUBSTITUTIONMODELSET_H_
#define _MIXEDSUBSTITUTIONMODELSET_H_


#include "SubstitutionModelSet.h"
#include "AbstractSubstitutionModel.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/VectorTools.h>

namespace bpp
{
  /**
   * @brief Substitution models manager for Mixed Substitution Models.
   * This class inherits from SubstitutionModelSet.
   *
   * This class is done to handle specific cases of choices among the
   * submodels of mixed substitution models. Each branch of the tree
   * is labelled by a mixed model, and a site may be restricted to a
   * set of submodels it is allowed to follow. These sets are defined
   * through an hypergrap, ie a list of hypernodes.
   *
   * For example, suppose there are 3 mixed models (M1,M2 and M3),
   * with 2, 3, 4 submodels (S1, S2, ...) each.
   *
   * If the sites are allowed to follow any combination of submodels
   * (12 combinations) the corresponding hypergraph has only one
   * hypernode: (<1,2>,<1,2,3>,<1,2,3,4>).
   *
   * The hypergraph with hypernodes
   * ((<1>,<1,2>,<1,2>),(<2>,<3>,<3,4>)) means that a site either
   * follows 6 combinations:
   *
   * M1:S1, M2:S1 or M2:S2, and M3:S1 or M3:S2.
   *
   * or
   *
   * M1:S2, M2:S3, and M3:S3 or M3:S4.
   *
   *
   * Actually, additional coordinates are set when there are non mixed
   * models, with no value in them, and not used in practice.
   *
   * An hypernode is valid only if each mixed model is represented at
   * least by one submodel.
   *
   * Dependency of the submodels entails constraints in the
   * probabilities of the submodels, and definition of the hypernodes
   * must be taken with care for the whole modelling to be possible.
   *
   *
   *
   * In this implementation, for sake of simplification (and for
   * reason of time), all the submodels must belong to exactly one
   * given hypernode, but in theory more complex dependencies are 
   * possible.
   *
   * Concerning the probabilities of the submodels in each hypernode,
   * the first coordinate (ie set of submodels inside a mixed model)
   * in the list defines the probability of each hypernode. For each
   * coordinate (the first included), when there are several
   * submodels, their probabilities are conditional probabilities,
   * which means that they sum 1 and their ratio are unchanged.
   *
   * For instance, for hypergraph ((<1>,<1,2>,<1,2>),(<2>,<3>,<3,4>)),
   * the probabilities of hypernodes are the probabilities of M1:S1
   * and M1:S2. In the first hypernode, the probabilities of M2:S1 and
   * M2:S2 are P(M2:S1)/(P(M2:S1)+P(M2:S2)) and
   * P(M2:S2)/(P(M2:S1)+P(M2:S2)).
   *
   * We do not certify that the probability parameters of the mixed
   * models are all useful, and then identifiability problems may be
   * encountered. 
   *
   * There is a method ("complete") that creates an additional
   * hypernode to ensure that all submodels belong to at least an
   * hypernode.
   * 
   *
   */
  
  class MixedSubstitutionModelSet :
    public SubstitutionModelSet
  {
  public:
     class HyperNode
     {
     public: 
       class Node {

         /**
          * A vector<int> where all elements are different and in
          * increasing order.
          *
          */

         Vint vNumb_;

       public:
         Node() : vNumb_() {};
         Node(const Node& n): vNumb_(n.vNumb_){};
         Node& operator=(const Node& n){
           vNumb_=n.vNumb_;
           return *this;
         }
         
         ~Node(){};

         Node& operator=(const Vint& n) {
           vNumb_=n;
           return *this;
         }

         void insertN(const Vint& vn);

         size_t size() const{
           return vNumb_.size();
         }

         /*
          *@brief Cumulates the elements of the given Node into this one.
          *
          */

         Node& operator+=(const Node&);
         
         /*
          *@brief checks if this Node is included in another one.
          *
          */

         bool operator<=(const Node&) const;

         /*
          *@brief checks if this HyperNode includes another one.
          *
          */

         bool operator>=(const Node&) const;

         /*
          *@brief checks if this Node intersects another one.
          *
          */
       
         bool intersects(const Node&) const;

         int operator[](size_t i) const { return vNumb_[i]; }

       };

     private:
       
       std::vector<Node> vNumbers_;

       /*
        *@brief the coordinates of the Nodes that are not used.
        *
        */
       
       Vint vUnused_;

       /*
        *@brief probability of this HyperNode.
        *
        */

       double proba_;
       
     public:
       HyperNode(const MixedSubstitutionModelSet*);
       HyperNode(const HyperNode&);
       HyperNode& operator=(const HyperNode&);
       ~HyperNode(){};

       /*
        *@brief sets submodel numbers in the nMth mixed model. Checks
        *  if all the numbers are valid.
        *
        *@param nM number of the mixed model
        *@param vnS vector of numbers of the submodel
        */
    
       void setModel(size_t nM, const Vint& vnS);

       /*
        *@brief adds submodel numbers to the nMth mixed model. Checks
        *  if all the numbers are valid.
        *
        *@param nM number of the mixed model
        *@param vnS vector of numbers of the submodel
        */
    
       void addToModel(size_t nM, const Vint& vnS);
       /*
        *@brief Cumulates the Nodes of the given HyperNode into this one.
        *
        */

       HyperNode& operator+=(const HyperNode&);
         
       /*
        *@brief checks if this HyperNode is included in another one.
        *
        */
       
       bool operator<=(const HyperNode&) const;

       /*
        *@brief checks if this HyperNode includes at least a submodel of each mixed model
        *
        */
       bool isComplete() const;
       /*
        *@brief checks if this HyperNode includes another one.
        *
        */
       bool operator>=(const HyperNode&) const;

       /*
        *@brief checks if this HyperNode intersects another one.
        *
        */       
       bool intersects(const HyperNode&) const;

       /*
        *@brief returns the probability
        *
        */

       double getProbability() const {return proba_;}

       /*
        *@brief sets the probability
        *
        */

       void setProbability(double x) { proba_ = x; }
       
       const Node& getNode(size_t i) const { return vNumbers_[i]; }

     };

  private:
    
    std::vector<HyperNode*> vpHyperNodes_;

  public:
    /**
     * @brief Create a model set according to the specified alphabet.
     *
     * @param alpha The alphabet to use for this set.
     */
    MixedSubstitutionModelSet(const Alphabet* alpha): 
      SubstitutionModelSet(alpha), vpHyperNodes_() {}
    
    ~MixedSubstitutionModelSet();
    
    MixedSubstitutionModelSet(const MixedSubstitutionModelSet& set);
    
    MixedSubstitutionModelSet(const SubstitutionModelSet& set);

    MixedSubstitutionModelSet& operator=(const MixedSubstitutionModelSet& set);
    
#ifndef NO_VIRTUAL_COV
    MixedSubstitutionModelSet*
#else
    Clonable*
#endif
    clone() const { return new MixedSubstitutionModelSet(*this); }

    /**
     * @brief Resets the list of the HyperNodes
     */
    
    void clear();

    /*
     *@brief adds a new empty HyperNode to the end of the HyperNodes
     * list.
     */
    
    void addEmptyHyperNode();

    /*
     *@brief adds the copy of an HyperNode to the end of the
     * HyperNodes list.
     */
    
    void addHyperNode(const HyperNode& hn);

    /*
     *@brief If necessary, adds a new HyperNode such that all
     *       submodels of the mixture models are at least in an
     *       HyperNode.
     *
     * Returns true iff a new path has been built.
     * 
     */

    bool complete();
    
    /*
     *@brief adds a submodel number to the nMth mixed model of the
     *  nHth HyperNode of the list (default nH=0). Checks if all the
     *  numbers are valid.
     *
     *@param nM number of the mixed model
     *@param vnS number of the submodel
     *@param nH number of the concerned HyperNode (default the last element of
     *     the list)
     */
    
    void addToHyperNode(size_t nM, const Vint& vnS, int nH = -1);

    size_t getNumberOfHyperNodes() const{ return vpHyperNodes_.size();}

    HyperNode& getHyperNode(size_t i) {return *vpHyperNodes_[i];} 

    const HyperNode& getHyperNode(size_t i) const {return *vpHyperNodes_[i];} 

    /*
     *@brief Checks if all the path (ie hypernodes) are exclusive.
     *
     */
    
    bool hasExclusivePaths() const;

    void fireParameterChanged(const ParameterList& parameters);

    /*
     *@brief compute the probabilities in all the HyperNodes
     *
     */

    void computeHyperNodesProbabilities();

    /*
     *@brief computes the probability of an HyperNode, given
     *     the conditional probabilities of the submodels computed
     *     from the hypernodes of this MixedSubstitutionModelSet
     *     object. If the HyperNode does not match the structure of
     *     allowed by this MixedSubstitutionModelSet, an Exception
     *     is thrown.
     *
     *     The probability of an HyperNode is the product -- on the
     *     set of the mixed models -- of the sums of the
     *     conditional probabilities of the submodels that belon to
     *     this hypernode for each mixed model.
     *
     *@param hn the HyperNode which conditional probability is computed.
     */

    double getHyperNodeProbability(const HyperNode& hn) const;

  };
  
} // end of namespace bpp.

#endif // _MIXEDSUBSTITUTIONMODELSET_H_

