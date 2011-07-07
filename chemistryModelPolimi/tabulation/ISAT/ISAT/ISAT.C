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
    const TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    tabulation<CompType,ThermoType>(chemistryProperties, chemistry),
    chemistry_(chemistry),
    chemisTree_(chemistry,readLabel(this->coeffsDict_.lookup("maxElements"))),
    tolerance_(0.0),
    scaleFactor_(chemistry_.Y().size()+2,1.0),
    logT_(false),
    tauStar_(false),
    clean_(false),
    checkUsed_(INT_MAX),
    checkGrown_(INT_MAX),
    MRUSize_(0)
    
{

    tauStar_.readIfPresent("tauStar",this->coeffsDict_); 
/*    if(this->coeffsDict_.found("tauStar"))
    {
        tauStar_ = this->coeffsDict_.lookup("tauStar");
    }
  */  
    if(this->online_)
    {
        tolerance_ = readScalar(this->coeffsDict_.lookup("tolerance"));
        if (this->coeffsDict_.found("checkUsed")) //if not specified equal to 1000 (most probably not used)
            checkUsed_ = readScalar(this->coeffsDict_.lookup("checkUsed"));
        else Info << "checkUsed not specified, set to 1000 by default" << endl;
        if (this->coeffsDict_.found("checkGrown")) //if not specified equal to 1000 (most probably not used)    
            checkGrown_ = readLabel(this->coeffsDict_.lookup("checkGrown"));
        else Info << "checkGrown not specified, set to 1000 by default" << endl;
        //       logT_ = onlineDict_.lookup("logT");
        clean_.readIfPresent("cleanAll",this->coeffsDict_);
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
        
        if (this->coeffsDict_.found("MRUSize")) //if not specified, considered as 0 (and not used)
            MRUSize_ =  readLabel(this->coeffsDict_.lookup("MRUSize"));   
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
                //delete leaf from tree and compact the list of node and chemPoint
                chemisTree_.deleteLeafCompact(phi0->listIndex());
		closest = NULL;
		return false;
            }    
            else 
            {
                addToMRU(phi0->listIndex());
                return true;                
            }
        }
	else
	{
	    return false;
	}
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
    if (phi0->nGrown() > checkGrown())
    {
        chemisTree_.deleteLeafCompact(phi0->listIndex());
	phi0Base = NULL;
	return false;
    }
    else if (phi0->checkSolution(phiq,Rphiq))
    {
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
void Foam::ISAT<CompType, ThermoType>::add
(
	const scalarField& phiq,    
	const scalarField& Rphiq, 
              List<List<scalar> >& A,
              chemPointBase*& phi0,
	const label nCols
)
{
    if(!phi0) add(phiq,Rphiq,A,nCols);
    else if (chemisTree().isFull())
    {
        if (MRUSize_>0)
        {
            List<chemPointISAT<CompType, ThermoType>*> tempList(MRUList_.size());
            List<chemPointISAT<CompType, ThermoType>*>& chemPointISATListTemp = chemisTree().chemPointISATList();
            //create a copy of each chemPointISAT of the MRUList_
            SLList<label>::iterator iter = MRUList_.begin();
            label j = 0;
            
            //the list of chemPoint stored in the binaryTree will be removed
            //a temporary list is then created to add them back to the tree
            for ( ; iter != MRUList_.end(); ++iter)
            {
                label chemPointISATIndex = iter();
                /*tempList[j]=new chemPoint(chemPointListTemp[chemPointIndex]->phi(),
                 chemPointListTemp[chemPointIndex]->Rphi(),
                 chemPointListTemp[chemPointIndex]->A(),
                 scaleFactor(),
                 tolerance(),
                 nCols	
                 );
                 */
                tempList[j]=new chemPointISAT<CompType, ThermoType>(*chemPointISATListTemp[chemPointISATIndex]);
                j++;
                
            }
            
            chemisTree().clear();
            chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols);
            
            for (label i=0; i<tempList.size() ; i++)
            {
                
                chemisTree().insertNewLeaf(	tempList[i]->phi(),
                                           tempList[i]->Rphi(),
                                           tempList[i]->A(),
                                           scaleFactor(),
                                           tolerance(),
                                           nCols
                                           );
            }
            
        }
        else
        {
            chemisTree().clear();
            chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols);
        }
        
        Info << "depth : " << chemisTree().depth() << endl;
    }
    else if (chemisTree().size() ==1)
    {
        chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols);
    }
    else
    {
        chemPointISAT<CompType, ThermoType>* phi0ISAT = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(phi0);
        chemisTree().insertNewLeaf(phi0ISAT, phiq, Rphiq, A, scaleFactor(), tolerance(), nCols);
        phi0 = phi0ISAT;
    }
}//end add

/*---------------------------------------------------------------------------*\
	Add a new leaf to the binary tree 
	(without reference to an existing chemPoint)
	Input : phiq the new composition to store
			Rphiq the mapping of the new composition point
			A the mapping gradient matrix
			nCols the size of the matrix
	Output: void
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::add
(
	const scalarField& phiq,    
	const scalarField& Rphiq, 
              List<List<scalar> >& A,
	const label nCols
)
{
	chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols);
}//end add

	
/*---------------------------------------------------------------------------*\
	Replace a new leaf of the binary tree 
	Input : phi0 the chemPoint to replace in the binary tree
			phiq the new composition to store
			Rphiq the mapping of the new composition point
			A the mapping gradient matrix
			nCols the size of the matrix
	Output: void
	Description: the replacement of a leaf of the binary tree is done in
			three steps. First the index of the chemPoint in the list hold
			in the binaryTree class is saved. Then the leaf is deleted and 
			the tree is reshaped in order to maintain the current splitting
			of the composition space and allow following binary tree search.
			Finally the leaf is inserted as usual but its position in the 
			chemPoint list is specified by the deleted chemPoint.
			(Note: position in the list and in the binary tree is independant)
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::replace
(
              chemPointBase*& phi0Base,
	const scalarField& phiq,    
	const scalarField& Rphiq, 
              List<List<scalar> >& A,
	const label nCols
)
{
	chemPointISAT<CompType, ThermoType>* phi0 = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(phi0Base);
	label phi0Index = phi0->listIndex();
	label nodeIndex = chemisTree().deleteLeaf(phi0Index);

	chemisTree().insertNewLeaf(phiq, Rphiq, A, scaleFactor(), tolerance(), nCols, nodeIndex, phi0Index);
}//end add


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::clear()
{

    Info<< "Clearing chemistry library" << endl;
    chemisTree_.clear();
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
void Foam::ISAT<CompType, ThermoType>::addToMRU(label chemPointISATIndex)
{
    if (MRUSize_ > 0)
    {
        //first search if the cp is already in the list
        bool isInTheTree = false;
        SLList<label>::iterator iter = MRUList_.begin();
        for ( ; iter != MRUList_.end(); ++iter)
        {
            if(iter()==chemPointISATIndex)
            {
                isInTheTree = true;
                break;
            }
        }
        //if it is in the tree, then move it to front
        if (isInTheTree)
        {
            if (iter()!=MRUList_.last())
            {
                //iter hold the position of the element to move
                label toTail = MRUList_.remove(iter);
                
                //insert the element in front of the list
                MRUList_.append(toTail);
            }
        }
        else //cp not yet in the list
        {
            if (MRUList_.size()==MRUSize_)
            {	
                
                MRUList_.removeHead();
                MRUList_.append(chemPointISATIndex);
            }
            else
            {
                MRUList_.append(chemPointISATIndex);
            }
        }
    }
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
