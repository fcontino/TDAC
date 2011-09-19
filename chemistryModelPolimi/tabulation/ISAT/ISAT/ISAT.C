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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "ISAT.H"
#include "error.H"
#include "TDACChemistryModel.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "SLList.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
template<class CompType, class ThermoType>
Foam::ISAT<CompType, ThermoType>::ISAT
(
    const dictionary& chemistryProperties,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    tabulation<CompType,ThermoType>(chemistryProperties, chemistry),
    chemistry_(chemistry),
    chemisTree_(chemistry,this->coeffsDict_),
    tolerance_(readScalar(this->coeffsDict_.lookup("tolerance"))),
    scaleFactor_(chemistry_.Y().size()+2,1.0),
    tauStar_(false),
    clean_(this->coeffsDict_.lookupOrDefault("cleanAll", false)),
    checkUsed_(this->coeffsDict_.lookupOrDefault("checkUsed", (runTime_->endTime().value()-runTime_->startTime().value())/runTime_->deltaT().value())),
    checkGrown_(this->coeffsDict_.lookupOrDefault("checkGrown", INT_MAX)),
    MRUSize_(this->coeffsDict_.lookupOrDefault("MRUSize", 0)),
    cleaningRequired_(false),
    toRemoveList_(),
    nFailedFirst_(0),
    totRetrieve_(0),
    MRURetrieve_(this->coeffsDict_.lookupOrDefault("MRURetrieve", false)),
    max2ndRetBalance_(this->coeffsDict_.lookupOrDefault("max2ndRetBalance",1.0)),
    maxDepthFactor_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "maxDepthFactor",
            (chemisTree_.maxElements()-1)/(std::log(chemisTree_.maxElements())/std::log(2.0))
        )
    ),
    runTime_(&chemistry.time()),
    previousTime_(runTime_->timeToUserTime(runTime_->startTime().value())),
    checkEntireTreeInterval_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "checkEntireTreeInterval",
            (runTime_->endTime().value()-runTime_->startTime().value())/runTime_->deltaT().value()
        )
    ),
    chPMaxLifeTime_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "chPMaxLifeTime",
            (runTime_->endTime().value()-runTime_->startTime().value())/runTime_->deltaT().value()
        )
    ),
    chPMaxUseInterval_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "chPMaxUseInterval",
            (runTime_->endTime().value()-runTime_->startTime().value())/runTime_->deltaT().value()
        )
    )
{

    if(this->online_)
    {
        
        dictionary scaleDict(this->coeffsDict_.subDict("scaleFactor"));
        label Ysize = chemistry_.Y().size();
        for(label i = 0; i<Ysize; i++)
        {
            if(!scaleDict.found(chemistry_.Y()[i].name()))
            {
                scaleFactor_[i] = readScalar(scaleDict.lookup("otherSpecies"));
            }
            else
            {
                scaleFactor_[i] = readScalar(scaleDict.lookup(chemistry_.Y()[i].name()));
            }
        } 
        scaleFactor_[Ysize] = readScalar(scaleDict.lookup("Temperature"));    
        scaleFactor_[Ysize+1] = readScalar(scaleDict.lookup("Pressure"));
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ISAT<CompType, ThermoType>::~ISAT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::retrieve
(
    const Foam::scalarField& phiq,
    chemPointBase*& closest
)
{
    chemisTree_.binaryTreeSearch(phiq, chemisTree_.root(),closest);
    if (!closest)
    {
	return false;
    }
    else
    {
	chemPointISAT<CompType, ThermoType>* phi0 = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(closest);
        if(phi0->inEOA(phiq))
        {	
            if(phi0->nUsed() > checkUsed()*chemistry_.Y()[0].size())
            {
                cleaningRequired_ = true;
                bool inList(false);
                forAll(toRemoveList_,tRi)
                {
                    if(toRemoveList_[tRi]==phi0)
                    {
                        inList=true;
                        break;
                    }
                }
                if(!inList)
                {
                    toRemoveList_.append(phi0);
                }    
            }
            phi0->lastTimeUsed()=runTime_->timeOutputValue();
            addToMRU(phi0);
            totRetrieve_++;
            return true;                
        }
	else if(chemistry_.exhaustiveSearch())
	{   
            //exhaustiveSearch if BT search failed
            phi0=chemisTree_.treeMin();
            while(phi0!=NULL)
            {
                if(phi0->inEOA(phiq))
                {
                    closest = phi0;
                    chemistry_.nFailBTGoodEOA()++;
                    if(phi0->nUsed() > checkUsed()*chemistry_.Y()[0].size())
                    {
                        cleaningRequired_ = true;
                        bool inList(false);
                        forAll(toRemoveList_,tRi)
                        {
                            if(toRemoveList_[tRi]==phi0)
                            {
                                inList=true;
                                break;
                            }
                        }
                        if(!inList)
                        {
                            toRemoveList_.append(phi0);
                        }    
                    }
                    phi0->lastTimeUsed()=runTime_->timeOutputValue();
                    addToMRU(phi0);
                    return true;
                }
                phi0=chemisTree_.treeSuccessor(phi0);
            }
            //if exhaustive search failed, return false
            return false;
	}
        else if(chemisTree_.secondaryBTSearch(phiq, phi0))
        {
            if(phi0->nUsed() > checkUsed()*chemistry_.Y()[0].size())
            {
                cleaningRequired_ = true;
                bool inList(false);
                forAll(toRemoveList_,tRi)
                {
                    if(toRemoveList_[tRi]==phi0)
                    {
                        inList=true;
                        break;
                    }
                }
                if(!inList)
                {
                    toRemoveList_.append(phi0);
                }    
            }
            closest = phi0;
            chemistry_.nFailBTGoodEOA()++;
            phi0->lastTimeUsed()=runTime_->timeOutputValue();
            addToMRU(phi0);
            nFailedFirst_++;
            totRetrieve_++;
            return true;
        }
        else if(MRURetrieve_)
        {
            typename SLList<chemPointISAT<CompType, ThermoType>*>::iterator iter = MRUList_.begin();
            for ( ; iter != MRUList_.end(); ++iter)
            {
                phi0=iter();
                if(phi0->inEOA(phiq))
                {
                    chemistry_.nFailBTGoodEOA()++;
                    if(phi0->nUsed() > checkUsed()*chemistry_.Y()[0].size())
                    {
                        cleaningRequired_ = true;
                        bool inList(false);
                        forAll(toRemoveList_,tRi)
                        {
                            if(toRemoveList_[tRi]==phi0)
                            {
                                inList=true;
                                break;
                            }
                        }
                        if(!inList)
                        {
                            toRemoveList_.append(phi0);
                        }    
                    }
                    phi0->lastTimeUsed()=runTime_->timeOutputValue();
                    addToMRU(phi0);
                    nFailedFirst_++;
                    totRetrieve_++;
                    return true;
                }
            }
        }
        return false;
    }
}


/*---------------------------------------------------------------------------*\
	Check if the composition of the query point phiq lies in the ellipsoid of 
	accuracy approximating the region of accuracy of the stored chemPoint phi0
	Input : phi0 the nearest chemPoint used in the linear interpolation
			phiq the composition of the query point for which we want to 
				compute the mapping
	Output: true if phiq is in the EOA, false if not
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::grow
(
	chemPointBase*& phi0Base,
	const scalarField& phiq,
	const scalarField& Rphiq
)
{
    if(!phi0Base) 
    {
	return false;
    }
	
    chemPointISAT<CompType, ThermoType>* phi0 = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(phi0Base);

    
    if (phi0->checkSolution(phiq,Rphiq))
    {
        //phi0 is only grown when checkSolution returns true
        if (phi0->nGrown() > checkGrown())
        {
            cleaningRequired_ = true;
            bool inList(false);
            forAll(toRemoveList_,tRi)
            {
                if(toRemoveList_[tRi]==phi0)
                {
                    inList=true;
                    break;
                }
            }
            if(!inList)
            {
                toRemoveList_.append(phi0);
            }    
        }
	return true;
    }
    else
    {
	return false;
    }
}

/*---------------------------------------------------------------------------*\
    Compute and return the mapping of the composition phiq from stored data
	Input : phi0 the nearest chemPoint used in the linear interpolation
			phiq the composition of the query point for which we want to 
				compute the mapping
			Rphiq the mapping of the new composition point (given as empty)
	Output: void (the mapping is stored in the Rphiq array)
	Rphiq = Rphi0 + A * (phiq-phi0)
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::calcNewC
(
	chemPointBase*& phi0Base,
	const scalarField& phiq,
              scalarField& Rphiq
)
{
    label nEqns = chemistry_.nEqns();
    chemPointISAT<CompType, ThermoType>* phi0 = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(phi0Base);
    bool isDACActive = phi0->DAC();
    List<label>& completeToSimplified(phi0->completeToSimplifiedIndex());		
    Rphiq = phi0->Rphi(); //Rphiq=Rphi0
    scalarField dphi=phiq-phi0->phi();
    
    const List<List<scalar> >& Avar = phi0->A();


    //Rphiq[i]=Rphi0[i]+A[i][j]dphi[j]
    for (label i=0; i<nEqns-2; i++)
    {
        if (isDACActive)
        {
            label si=completeToSimplified[i];
            //the species is active
            if (si!=-1)
            {
                for (label j=0; j<nEqns-2; j++) 
                {
                    label sj=completeToSimplified[j];
                    if (sj!=-1)
                        Rphiq[i] += Avar[si][sj]*dphi[j];
                }
                Rphiq[i] += Avar[si][phi0->NsDAC()]*dphi[nEqns-2];
                Rphiq[i] += Avar[si][phi0->NsDAC()+1]*dphi[nEqns-1];
                //As we use an approximation of A, Rphiq should be ckeck for 
                //negative value
                Rphiq[i] = max(0.0,Rphiq[i]);
            }
            //the species is not active A[i][j] = I[i][j]
            else
            {
                Rphiq[i] += dphi[i];
                Rphiq[i] = max(0.0,Rphiq[i]);
            }
        }
        else //DAC is not active
        {
            for (label j=0; j<nEqns; j++) Rphiq[i] += Avar[i][j]*dphi[j];
            //As we use an approximation of A, Rphiq should be ckeck for 
            //negative value
            Rphiq[i] = max(0.0,Rphiq[i]);
            
        }
    }
}//end calcNewC


/*---------------------------------------------------------------------------*\
	Add a new leaf to the binary tree
	(with reference to an existing chemPoint))
	Input : phiq the new composition to store
			Rphiq the mapping of the new composition point
			A the mapping gradient matrix
			phi0 the chemPoint which is the nearest from phiq and which will
				be replaced by a node splitting the composition space between
				phi0 and phiq
			nCols the size of the matrix
	Output: void
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::add
(
	const scalarField& phiq,    
	const scalarField& Rphiq, 
              List<List<scalar> >& A,
              chemPointBase*& phi0,
	const label nCols
)
{
    if (chemisTree().isFull())
    {
        if (MRUSize_>0)
        {
            List<chemPointISAT<CompType, ThermoType>*> tempList(MRUList_);
            //create a copy of each chemPointISAT of the MRUList_
            
            chemisTree().clear();
	    toRemoveList_.clear();

            chemPointISAT<CompType, ThermoType>* nulPhi=0;
            //insert the point to add first
            chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols,nulPhi);
            
            addToMRU(chemisTree().treeMin());
            
            forAll(tempList,i)
            {
                chemisTree().insertNewLeaf
                (
                    tempList[i]->phi(),
                    tempList[i]->Rphi(),
                    tempList[i]->A(),
                    scaleFactor(),
                    tolerance(),
                    nCols,
                    nulPhi
                );
            }
        }
        else
        {
            chemisTree().clear();
	    toRemoveList_.clear();
            chemPointISAT<CompType, ThermoType>* nulPhi=0;
            chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols,nulPhi);
        }
        return true;
    }
    else
    {
        chemPointISAT<CompType, ThermoType>* phi0ISAT = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(phi0);
        chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols, phi0ISAT);
        phi0 = phi0ISAT;
        return false;
    }
}//end add


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::clear()
{

    Info<< "Clearing chemistry library" << endl;
    chemisTree_.clear();
    toRemoveList_.clear();
}


/*---------------------------------------------------------------------------*\
	Add a chemPoint to the MRU list 
	Input : cp the chemPoint to add
	Output: void
	Description: If cp is not in the list, insert it in front of it.
				 If cp is in the list, move it to the front.
	Note : tested with entry of type label			 
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::addToMRU(chemPointISAT<CompType, ThermoType>* phi0)
{
    if (MRUSize_ > 0)
    {
        //first search if the chemPoint is already in the list
        bool isInList = false;
        typename SLList<chemPointISAT<CompType, ThermoType>*>::iterator iter = MRUList_.begin();
        for ( ; iter != MRUList_.end(); ++iter)
        {
            if(iter()==phi0)
            {
                isInList = true;
                break;
            }
        }
        //if it is in the list, then move it to front
        if (isInList)
        {
            if (iter()!=MRUList_.last())
            {
                //iter hold the position of the element to move
                MRUList_.remove(iter);
                
                //insert the element in front of the list
                MRUList_.append(phi0);
            }
        }
        else //chemPoint not yet in the list
        {
            if (MRUList_.size()==MRUSize_)
            {	
                
                MRUList_.removeHead();
                MRUList_.append(phi0);
            }
            else
            {
                MRUList_.append(phi0);
            }
        }
    }
}


template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::cleanAndBalance()
{
    
    bool treeModified(false);
    //1- check if the tree should be cleaned (flag from nUsed or nGrown)
    if(cleaningRequired_)
    {
        cleaningRequired_=false;
        //2- remove the points that have raised a flag because of number of growth or used 
        //(they are stored in the toRemoveList)
        forAll(toRemoveList_,trli)
        {
            chemisTree_.deleteLeaf(toRemoveList_[trli]);
        }
        toRemoveList_.clear(); //set size to 0, the pointers have been deleted in deleteLeaf function
        treeModified=true;
    }
    
    //3- check if the entire tree should be scanned (after a given number of time-steps or for other criterion)
    if
    (
        (runTime_->timeOutputValue()-previousTime_)
        >
        (checkEntireTreeInterval_*runTime_->timeToUserTime(runTime_->deltaTValue()))
    )
    {
        previousTime_ = runTime_->timeOutputValue();
        
        //3a- remove the points that are too old or not used recently

        //scan the entire tree
        chemPointISAT<CompType, ThermoType>* x = chemisTree_.treeMin();
        while(x!=NULL)
        {
            chemPointISAT<CompType, ThermoType>* xtmp = chemisTree_.treeSuccessor(x);
            if
            (
                ((runTime_->timeOutputValue() - x->timeTag()) > (chPMaxLifeTime_*runTime_->timeToUserTime(runTime_->deltaTValue()))) 
                || 
                ((runTime_->timeOutputValue() - x->lastTimeUsed()) > (chPMaxUseInterval_*runTime_->timeToUserTime(runTime_->deltaTValue())))
            )
            {
                chemisTree_.deleteLeaf(x);
                treeModified=true;
            }
            x = xtmp;
        }
        //3b- check if the tree should be balanced according to criteria:
        //      number of secondaryRetrieve above a given threshold (portion of totRetrieve = primary+secondary retrieve)
        //      depth of the tree bigger than a*log2(size), where a is a given parameter
        if(chemisTree_.size()>0)
        {
            if
                (
                 (nFailedFirst_ > max2ndRetBalance_*totRetrieve_) 
                 || 
                 (chemisTree_.depth() > maxDepthFactor_*std::log(chemisTree_.size())/std::log(2.0))
                 )
            {
                totRetrieve_=0;
                nFailedFirst_=0;
                treeModified=chemisTree_.balance();
            }
        }
    }
    
    //return a bool to specify if the tree structure has been modified
    return treeModified;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::operator=(const Foam::ISAT<CompType, ThermoType>& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::ISAT::operator=(const Foam::ISAT&)")
            << "Attempted to assignment to self"
            << abort(FatalError);
    }
}

// ************************************************************************* //
