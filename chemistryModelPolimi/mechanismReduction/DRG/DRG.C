/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "DRG.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::DRG<CompType,ThermoType>::DRG
(
    const Foam::dictionary& dict,
    Foam::TDACChemistryModel<CompType,ThermoType>& chemistry
)
:
    mechanismReduction<CompType,ThermoType>(dict, chemistry)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::DRG<CompType,ThermoType>::~DRG()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::DRG<CompType,ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
) 
{
    scalarField& completeC(this->chemistry_.completeC());
    scalarField c1(this->nSpecie_+2, 0.0);
    
    for(label i=0; i<this->nSpecie_; i++)
    {
        c1[i] = c[i];
	completeC[i] = c[i];
    }

    c1[this->nSpecie_] = T;
    c1[this->nSpecie_+1] = p;
    
    //Compute the rAB matrix
    RectangularMatrix<scalar> rABNum(this->nSpecie_,this->nSpecie_,0.0);
    scalarField rABDen(this->nSpecie_,0.0);
    
    //Number of initialized rAB for each lines
    Field<label> NbrABInit(this->nSpecie_,0);
    //Position of the initialized rAB, -1 when not initialized
    RectangularMatrix<label> rABPos(this->nSpecie_, this->nSpecie_, -1);
    //Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec(this->nSpecie_, this->nSpecie_, -1);
	
    scalar pf,cf,pr,cr;
    label lRef, rRef;
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
	//for each reaction compute omegai        
	scalar omegai = this->chemistry_.omega
        (
	    R, c1, T, p, pf, cf, lRef, pr, cr, rRef
	);


	//then for each pair of species composing this reaction,
	//compute the rAB matrix (separate the numerator and 
	//denominator)
        /*
            While computing the rAB for all the species involved in the reaction
            it should be considered that one can write a reaction A+B=2C as A+B=C+C
            In this case, the following algorithm only take once the effect 
            of the species. It stores the species encountered in the reaction but
            use another list to see if this species has already been used
        */
        
        DynamicList<scalar> wA(R.lhs().size()+R.rhs().size()); 
        DynamicList<label> wAID(R.lhs().size()+R.rhs().size());
        
        forAll(R.lhs(), s)
        {
	    label ss = R.lhs()[s].index;
	    scalar sl = -R.lhs()[s].stoichCoeff;
            bool found(false);
            forAll(wAID,id)
            {
                if(ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
	}
	forAll(R.rhs(), s)
        {
            label ss = R.rhs()[s].index;
            scalar sl = R.rhs()[s].stoichCoeff;
            bool found(false);
            forAll(wAID,id)
            {
                if(ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }
        
        //Now that all nuAi*wi are computed, without counting twice species
        //present in both rhs and lhs, we can update the denominator and numerator
        // for the rAB
        wAID.shrink();
        forAll(wAID,id)
        {
            label curID = wAID[id];
            scalar curwA = ((wA[id]>=0) ? wA[id] : -wA[id]); //absolute value of aggregated value
            List<bool> deltaBi(this->nSpecie_, false);
            FIFOStack<label> usedIndex;
            forAll(R.lhs(), j)
            {
                label sj = R.lhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj]=true;
            }
            forAll(R.rhs(), j)
            {
                label sj = R.rhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj]=true;
            }

            deltaBi[curID]=false; //disable for self reference (by definition rAA=0)
            
            while(!usedIndex.empty())
            {
                label curIndex = usedIndex.pop();
                
                if(deltaBi[curIndex])
                {
                    deltaBi[curIndex]=false; //disable to avoid counting it more than once
                    if (rABPos[curID][curIndex]==-1)//test if this rAB is not initialized
                    {
                        rABPos[curID][curIndex] = NbrABInit[curID];//it starts at rABPos[ss][sj]=0
                        NbrABInit[curID]++;
                        rABNum[curID][rABPos[curID][curIndex]] = curwA;//to avoid overflow
                        rABOtherSpec[curID][rABPos[curID][curIndex]] = curIndex; //store the other specie involved					
                    }
                    else
                        rABNum[curID][rABPos[curID][curIndex]] += curwA;
                }
            }
            if (rABDen[curID] == 0.0)
            {
                rABDen[curID] = curwA;
            }
            else 
            {
                rABDen[curID] +=curwA;
            }
        }
    }//end forAll reactions
     //rii = 0.0 by definition
		
    label speciesNumber = 0;

    //set all species to inactive and activate them according
    //to rAB and initial set
    for (label i=0; i<this->nSpecie_; i++)
	this->activeSpecies_[i] = false;

    const labelList& SIS(this->searchInitSet());
    FIFOStack<label> Q;
    //Initialize the list of active species with the search initiating set (SIS)
    for (label i=0; i<SIS.size(); i++)
    {
	label q = SIS[i];
        this->activeSpecies_[q] = true;
        speciesNumber++;
        Q.push(q);
    }

    //Depth first search with rAB
    while (!Q.empty())
    {
        label u = Q.pop();
        scalar Den = rABDen[u];

        if (Den > VSMALL)
        {
            for (label v=0; v<NbrABInit[u]; v++)
            {
                label otherSpec = rABOtherSpec[u][v];                
                scalar rAB = rABNum[u][v]/Den;

                if(rAB>1)
                {
                    Info << "Badly Conditioned rAB : " << rAB << "species involved : "<<u << "," << otherSpec << endl;
                    rAB=1.0;
                }
                //do a DFS on B only if rAB is above the tolerance and if the species was not searched before
                if (rAB >= this->epsDAC() && !this->activeSpecies_[otherSpec])
                {
                    Q.push(otherSpec);
                    this->activeSpecies_[otherSpec] = true;
                    speciesNumber++;
                }
            }
        }
    }//end of Q.empty()
    
    //Put a flag on the reactions containing at least one removed species
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
		this->chemistry_.reactionsDisabled()[i]=false;
        
        forAll(R.lhs(), s)
        {
            label ss = R.lhs()[s].index;
            if (!this->activeSpecies_[ss]) //the species is inactive then the reaction is removed
            {		
                this->chemistry_.reactionsDisabled()[i]=true;//flag the reaction to disable it
                break; //further search is not needed
            }	
        }
        if (!this->chemistry_.reactionsDisabled()[i])//if the reaction has not been disabled yet
        {
            forAll(R.rhs(), s)
            {
                label ss = R.rhs()[s].index;
                if (!this->activeSpecies_[ss]) //the species is inactive then the reaction is removed
                {
                    this->chemistry_.reactionsDisabled()[i]=true;//flag the reaction to disable it
                    break; //further search is not needed
                }
            }
        }
    }//end of loop over reactions
    this->NsSimp_ = speciesNumber;
    scalarField& simplifiedC(this->chemistry_.simplifiedC());
    simplifiedC.setSize(this->NsSimp_+2);
    DynamicList<label>& s2c(this->chemistry_.simplifiedToCompleteIndex());
    s2c.setSize(this->NsSimp_);
    Field<label>& c2s(this->chemistry_.completeToSimplifiedIndex());

    label j = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        if (this->activeSpecies_[i])
        {
            s2c[j] = i;
            simplifiedC[j] = c[i];
            c2s[i] = j++;
	    if (!this->chemistry_.isActive(i))
            	this->chemistry_.setActive(i);            
        }
        else 
        {
            c2s[i] = -1;
        }
    }
            
    simplifiedC[this->NsSimp_] = T;
    simplifiedC[this->NsSimp_+1] = p;
    this->chemistry_.NsDAC(this->NsSimp_);
    //change temporary Ns in chemistryModel
    //to make the function nEqns working
    this->chemistry_.nSpecie() = this->NsSimp_;
}


// ************************************************************************* //
