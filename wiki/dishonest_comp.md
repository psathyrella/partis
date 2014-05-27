#Comparison of Dishonest Casino Models in StochHMM, R-HMM, Mamot, and HMMoc

##StochHMM Model
```
#STOCHHMM MODEL FILE
MODEL INFORMATION
======================================================
MODEL_NAME:	CASINO DICE MODEL
MODEL_DESCRIPTION:	Taken from CH3 Durbin/Eddy
MODEL_CREATION_DATE:	August 28,2009

TRACK SYMBOL DEFINITIONS
======================================================
DICE:	1,2,3,4,5,6

STATE DEFINITIONS
#############################################
STATE:	
	NAME:	INIT
TRANSITION:	STANDARD: P(X)
	FAIR:	0.5
	LOADED:	0.5
#############################################
STATE:	
	NAME:	FAIR
	PATH_LABEL:	F
	GFF_DESC:	FAIR
TRANSITION:	STANDARD: P(X)
	FAIR:	0.95
	LOADED:	0.05
	END:	1
EMISSION:	DICE: P(X)
	ORDER:	0
@1	2	3	4	5	6
0.167	0.167	0.167	0.167	0.167	0.167
#############################################
STATE:
	NAME:	LOADED
	PATH_LABEL:	L
	GFF_DESC:	LOADED
TRANSITION:	STANDARD: P(X)
	FAIR:	0.1
	LOADED:	0.9
	END:	1
EMISSION:	DICE: P(X)
	ORDER:	0
@1	2	3	4	5	6	
0.1	0.1	0.1	0.1	0.1	0.5
#############################################
//END
```

##R HMM Model  Adapted from [R HMM Library](http://cran.r-project.org/web/packages/HMM/HMM.pdf)
```
library(HMM)
args = commandArgs(trailingOnly = TRUE);

#Name of dice roll file (which is comma separated)
seqFile = args[1];

#Flag for Viterbi (if 1 then perform viterbi, otherwise perform posterior)
viterbi = as.numeric(args[2])

alphabet = c(1, 2, 3, 4, 5, 6)
states = c("Fair", "Loaded")

#State Emissions
fair = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
loaded = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.5)

#Transition matrix
transProbs = matrix(c(0.99, 0.01, 0.02, 0.98), c(2,2), byrow = TRUE)

#Emission Matrix
emissionProbs = matrix(c(fair, loaded), c(2, 2), byrow = TRUE)

hmm = initHMM(states, alphabet, transProbs = transProbs, emissionProbs = emissionProbs)

observations = scan(file = seqFile,sep=",")

if (viterbi == 1){
    vit = viterbi(hmm, observations)
}

if (viterbi == 0){
    f = forward(hmm, observations)
    b = backward(hmm, observations)
    i <- f[1, length(observations)]
    j <- f[2, length(observations)]
    probObservations = (i + log(1 + exp(j - i)))
    posterior = exp((f + b) - probObservations)
}

gc()
```


##[Mamot Model](http://bcf.isb-sib.ch/mamot/)
```
# HMM for the "dishonest casino" presented by Durbin, Eddy, Krogh, Mitchison 
# in "Biological sequence analysis" (Cambridge University Press, 1998).
# Exercise 1.1, p. 6 and p. 54, 56, 59, 61, 65

Alphabet: abcdef
################################################
State BEGIN
	E: 0  0  0  0  0  0
	T: F 1
# Fair
State F
	E: 0.1666  0.1666  0.1666  0.1666  0.1666  0.1666
	T: F 0.95  L 0.049  END 0.001
# not fair
State L 
	E: 0.1 0.1 0.1 0.1 0.1 0.5
	T: F 0.1  L 0.9
State END  
	E: 0  0  0  0  0  0
	T: END 0
################################################
NULLMODEL:  0.1666  0.1666  0.1666  0.1666  0.1666  0.1666
```


#[HMMoc](http://genserv.anat.ox.ac.uk/downloads/software/hmmoc/)
```
<?xml version="1.0"?>
<!--
      This file is part of HMMoC 1.3, a hidden Markov model compiler.
      Copyright (C) 2007 by Gerton Lunter, Oxford University.
  
      HMMoC is free software; you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation; either version 2 of the License, or
      (at your option) any later version.
  
      HMMOC is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.
  
      You should have received a copy of the GNU General Public License
      along with HMMoC; if not, write to the Free Software
      Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
-->



<hml debug="true">



<author>Gerton Lunter</author>



<alphabet id="dice">
 123456
</alphabet>



<output id="sequence">
  <alphabet idref="dice"/>
  <identifier type="length" value="iLen"/>
  <identifier type="sequence" value="aSeq"/>
  <code type="parameter" value="char *aSeq"/>
  <code type="parameter" value="int iLen"/>
</output>


<hmm id="Casino">

 <description>  The occasionally dishonest casino  </description>

 <outputs id="casinooutputs">
  <output idref="sequence"/>
 </outputs>


 <clique id="block1">
  <state id="start"/>
 </clique>

 <clique id="block2">
  <state id="honest"/>
  <state id="dishonest"/>
 </clique>

 <clique id="block3">
  <state id="end"/>
 </clique>


 <graph>
  <clique idref="block1"/>
  <clique idref="block2"/>
  <clique idref="block3"/>
 </graph>


 <transitions>
  <transition from="start" to="honest" probability="one" emission="emitHonest"/>
  <transition from="honest" to="honest" probability="stayHonest" emission="emitHonest"/>
  <transition from="honest" to="dishonest" probability="goDishonest" emission="emitDishonest"/>
  <transition from="dishonest" to="dishonest" probability="stayDishonest" emission="emitDishonest"/>
  <transition from="dishonest" to="honest" probability="goHonest" emission="emitHonest"/>
  <transition from="honest" to="end" probability="goStop" emission="empty"/>
  <transition from="dishonest" to="end" probability="goStop" emission="empty"/>
 </transitions>


 <code id="paramsClassDef" where="classdefinitions">
   <![CDATA[
     struct Params {
       double iGoHonest;
       double iGoDishonest;
       double iGoStop;
       double aEmitDishonest[6];
     };
   ]]>
  </code>


  <emission id="empty">
   <probability>
    <code type="expression"> 1.0 </code>
   </probability>
  </emission>


  <emission id="emitHonest">
   <output idref="sequence"/>
   <probability>
    <code type="statement">
     <identifier output="sequence" value="iEmission"/>
     <identifier type="result" value="iProb"/>
     <![CDATA[
  
       iProb = 1/6.0;
       /* This probability does not depend on the symbol.  HMMoC warns if it does not see the label 'iEmission'
          somewhere in the code -- its appearance in this comment stops it from warning */

     ]]>
    </code>
   </probability>
  </emission>


  <emission id="emitDishonest">
   <output idref="sequence"/>
   <probability>
    <code type="statement">
     <identifier output="sequence" value="iEmission"/>
     <identifier type="result" value="iProb"/>
     <!--  Here goes the code computing the probability -->
     <![CDATA[
  
       iProb = iPar.aEmitDishonest[ iEmission - '1' ];

     ]]>
    </code>
   </probability>
  </emission>


  <probability id="one"><code> 1.0 </code></probability>


  <probability id="goDishonest">
    <code>
      <!--  Tell HMMoC that this code requires an input parameter, which itself need a definition to make sense -->
      <code type="parameter" init="paramsClassDef" value="Params iPar"/>
      <!-- The actual code for this probability follows (no need to quote this) -->

        iPar.iGoDishonest 

    </code>
  </probability>

  <probability id="goHonest"><code> iPar.iGoHonest </code></probability>
  <probability id="goStop"><code> iPar.iGoStop </code></probability>
  <probability id="stayHonest"><code> 1.0 - iPar.iGoDishonest - iPar.iGoStop </code></probability>
  <probability id="stayDishonest"><code> 1.0 - iPar.iGoHonest - iPar.iGoStop </code></probability>

</hmm>

<!-- Code generation -->


<forward  outputTable="yes" name="Forward" id="fw">
  <!-- Specify HMM to make code for -->
  <hmm idref="Casino"/>
</forward>

<backward  outputTable="yes" baumWelch="yes" name="Backward" id="bw">
  <!-- Specify HMM to make code for -->
  <hmm idref="Casino"/>
</backward>

<viterbi  name="Viterbi" id="vit">
  <hmm idref="Casino"/>
</viterbi>

<codeGeneration realtype="bfloat" file="casino.cc" header="casino.h" language="C++">

  <forward idref="fw"/>
  <backward idref="bw"/>
  <viterbi idref="vit"/>

</codeGeneration>

 
</hml>
```

##HMMoc Viterbi Code
```
/*
 *    This file is part of HMMoC 1.3, a hidden Markov model compiler.
 *    Copyright (C) 2007 by Gerton Lunter, Oxford University.
 *
 *    HMMoC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    HMMOC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HMMoC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 \*/
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "casino.h"


void import_fasta(char* file,std::string& seq){
	std::ifstream fa;
	fa.open(file);
	
	if (!fa.good()){
		return;
	}
	
	std::string header;
	getline(fa,header,'\n');
	
	std::string line;
	while(getline(fa,line,'\n')){
		if (line[0] == ' ' || line[0]=='>'){
			continue;
		}
		
		seq += line;
	}

	return;
}


int main(int argc, char** argv) {
	
	// The parameters of the model
	Params iPar;
	
	// Fill it with some values
	iPar.iGoDishonest = 0.05;     // probability of going from Honest to the Dishonest state
	iPar.iGoHonest = 0.02;        // probability of going from Dishonest to the Honest state
	iPar.iGoStop = 0.00001;       // probability of going from either to the End state
	
	iPar.aEmitDishonest[0] = 0.1;
	iPar.aEmitDishonest[1] = 0.1;
	iPar.aEmitDishonest[2] = 0.1;
	iPar.aEmitDishonest[3] = 0.1;
	iPar.aEmitDishonest[4] = 0.1;
	iPar.aEmitDishonest[5] = 0.5; // Probability of throwing a 6 is heavily favoured in the Dishonest state
	
//Import Sequence
	
	std::string aSequence;
	import_fasta(argv[1],aSequence);
	int iPathLength=aSequence.length();
	char* seq = new char[iPathLength+1];
	strcpy(seq, aSequence.c_str());
	aSequence="";

	// Decode the emission sequence using Viterbi, and compute posteriors and Baum Welch counts using Forward and Backward
	CasinoDPTable *pViterbiDP;
	
	cout << "Calculating Viterbi probability..." << endl;
	bfloat iVitProb = Viterbi_recurse(&pViterbiDP, iPar, seq, iPathLength );
	cout << "Viterbi: "<<iVitProb<<endl;
	
	cout << "Calculating Viterbi path..." << endl;
	Path& iViterbiPath = Viterbi_trace(pViterbiDP, iPar, seq, iPathLength );
	
	// Compare the true and Viterbi paths, and print the posterior probability of being in the honest state
	int iVHonest = pViterbiDP->getId("honest");
	
	for (int i=0; i<iPathLength; i++) {
		
		if (iViterbiPath.toState(i) == iVHonest) {
			cout << "H";
		} else {
			cout << "D";
		}
	}
	std::cout << std::endl;
	
}
```


##HMMoc Posterior Code
```
/*
 *    This file is part of HMMoC 1.3, a hidden Markov model compiler.
 *    Copyright (C) 2007 by Gerton Lunter, Oxford University.
 *
 *    HMMoC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    HMMOC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HMMoC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 \*/
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "casino.h"

void import_fasta(char* file,std::string& seq){
	std::ifstream fa;
	fa.open(file);
	
	if (!fa.good()){
		return;
	}
	
	std::string header;
	getline(fa,header,'\n');
	
	std::string line;
	while(getline(fa,line,'\n')){
		if (line[0] == ' ' || line[0]=='>'){
			continue;
		}
		
		seq += line;
	}
	
	return;
}


int main(int argc, char** argv) {
	
	// The parameters of the model
	Params iPar;
	
	// Fill it with some values
	iPar.iGoDishonest = 0.05;     // probability of going from Honest to the Dishonest state
	iPar.iGoHonest = 0.02;        // probability of going from Dishonest to the Honest state
	iPar.iGoStop = 0.00001;       // probability of going from either to the End state
	
	iPar.aEmitDishonest[0] = 0.1;
	iPar.aEmitDishonest[1] = 0.1;
	iPar.aEmitDishonest[2] = 0.1;
	iPar.aEmitDishonest[3] = 0.1;
	iPar.aEmitDishonest[4] = 0.1;
	iPar.aEmitDishonest[5] = 0.5; // Probability of throwing a 6 is heavily favoured in the Dishonest state

	
	std::string aSequence;
	import_fasta(argv[1],aSequence);
	int iPathLength=aSequence.length();
	char* seq = new char[iPathLength+1];
	strcpy(seq, aSequence.c_str());
	aSequence="";
		
	// Decode the emission sequence using Viterbi, and compute posteriors and Baum Welch counts using Forward and Backward
	CasinoDPTable *pViterbiDP, *pFWDP, *pBWDP;
	CasinoBaumWelch bw;
	
	cout << "Calculating Forward probability..." << endl;
	bfloat iFWProb = Forward(&pFWDP, iPar, seq, iPathLength );
	cout << "Forward: "<<iFWProb<<endl;
	
	cout << "Calculating Backward probability..." << endl;
	bfloat iBWProb = Backward(bw, pFWDP, &pBWDP, iPar, seq, iPathLength );
	cout << "Backward:"<<iBWProb<<endl;
	
	
	// Compare the true and Viterbi paths, and print the posterior probability of being in the honest state
	int iVHonest = pViterbiDP->getId("honest");
	
	for (int i=0; i<iPathLength; i++) {

		double iPosterior = pFWDP->getProb("honest",i+1)*pBWDP->getProb("honest",i+1)/iFWProb;
		cout << " " << iPosterior << endl;
		
	}
}
```

