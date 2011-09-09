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

#include "chemPointISAT.H"
#include "scalarField.H"
#include "binaryNode.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
binaryNode<CompType, ThermoType>::binaryNode
(
)
:
    elementLeft_(NULL),
    elementRight_(NULL),
    left_(NULL), 
    right_(NULL),
    parent_(NULL)
{}


template<class CompType, class ThermoType>
binaryNode<CompType, ThermoType>::binaryNode
(
    chemPointISAT<CompType, ThermoType>* elementLeft,
    chemPointISAT<CompType, ThermoType>* elementRight,
    binaryNode<CompType, ThermoType>* parent
)
:
    elementLeft_(elementLeft),
    elementRight_(elementRight),
    left_(NULL), 
    right_(NULL),
    parent_(parent),
    v_(elementLeft->spaceSize(),0.0)
{
	calcV(elementLeft, elementRight, v_);
	a_ = calcA(elementLeft, elementRight);
}

template<class CompType, class ThermoType>
binaryNode<CompType, ThermoType>::binaryNode
(
    binaryNode<CompType, ThermoType> *bn
)
:
    elementLeft_(bn->elementLeft()),
    elementRight_(bn->elementRight()),
    left_(bn->left()), 
    right_(bn->right()),
    parent_(bn->parent()),
    v_(bn->v()),
    a_(bn->a())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*---------------------------------------------------------------------------*\			  
    Compute vector v:
	Let E be the ellipsoid which covers the region of accuracy of 
	the left leaf (previously defined). E is described by 
	E={phi| ||L^T.(phi-phi0)|| <= 1}, (see chemPoint for more details).
	let E' be the tranformation of E in a space where E' is a hypersphere
	centered at the origin, in this space y=L^T.(phi-phi0) and then
	E'={y| ||y||<=1}
	let u be the unit vector joining the center of E' and the newly added 
	composition point in the transformed space (y2=L^T.(phiq-phi0)),u = y2/||y2| 
	Then the hyperplane separating the two points is defined as the
	perpendicular bisector of the segment linking 0 to y2
		H' = {y| u^T.(y-yh) = 0},
	where yh = y2/2.
	In the orginal composition space, the hyperplane H is defined by
		H = {y| v^T(phi-phih) = 0},
	where phih = phi0 + L^-T.yh = (phi0 + phiq) / 2 and v is
			  L.L^T (phiq-phi0)
		v = -------------------- .
			||L.L^T (phiq-phi0)||
			
	Note :  v is not normalised in this implementation since it is used
                on both side of an equality to know if one should go on the
                left or the right in the binary tree
    Input : elementLeft : chemPoint of the left element
	    elementRight: chemPoint of the right element
	    v		: empty scalar field to store v
    Output : void (v is stored in the empty scalarField)
\*---------------------------------------------------------------------------*/						
template<class CompType, class ThermoType>
void binaryNode<CompType, ThermoType>::calcV(chemPointISAT<CompType, ThermoType>*& elementLeft, chemPointISAT<CompType, ThermoType>*& elementRight, scalarField& v)
{
    //LT is the transpose of the L matrix
    List<List<scalar> >& LT = elementLeft->LT();
    label dim = elementLeft->spaceSize();
    if (elementLeft->DAC()) dim = elementLeft->NsDAC()+2;
    scalarField phiDif = elementRight->phi() - elementLeft->phi();
    const scalarField& scaleFactor = elementLeft->scaleFactor();
    const scalar epsTol = elementLeft->epsTol();
    
    for (label i=0; i<elementLeft->spaceSize(); i++)
    {
        label si = i;
        bool outOfIndexI = true;
        if(elementLeft->DAC())
        {
            if(i<elementLeft->spaceSize()-2)
            {
                si = elementLeft->completeToSimplifiedIndex(i);
                outOfIndexI = (si==-1);
            }
            else
            {
                outOfIndexI = false;
                label dif = i-(elementLeft->spaceSize()-2);
                si = elementLeft->NsDAC()+dif;
            }
        }
        if(!(elementLeft->DAC()) || (elementLeft->DAC() && !(outOfIndexI)))
        {
            v[i]=0.0;
            for (label j=0; j<elementLeft->spaceSize(); j++)
            {
                label sj = j;
                bool outOfIndexJ = true;
                if(elementLeft->DAC())
                {
                    if(j<elementLeft->spaceSize()-2)
                    {
                        sj = elementLeft->completeToSimplifiedIndex(j);	
                        outOfIndexJ = (sj==-1);
                    }
                    else
                    {
                        outOfIndexJ = false;
                        label dif = j-(elementLeft->spaceSize()-2);
                        sj = elementLeft->NsDAC()+dif;
                    }
                }			     
                if(!(elementLeft->DAC()) || (elementLeft->DAC() && !(outOfIndexJ)))
                {
                    //since L is a lower triangular matrix k=0->min(i,j)
                    for (label k=0; k<=min(si,sj); k++) v[i] += LT[k][si]*LT[k][sj]*phiDif[j];
                }
            }
        }
        else
        {
            //when it is an inactive species the diagonal element of LT is simply 1/(scaleFactor*epsTol)
            scalar div = scaleFactor[i]*scaleFactor[i]*epsTol*epsTol;
            v[i] = phiDif[i]/div;
        }
    }
}
    
/*---------------------------------------------------------------------------*\
	Let 'a' be the product v^T.phih, with phih = (phi0 + phiq)/2. 
	When travelling in the binary tree, 
	to know in which part of the composition space the query point 'phi' 
	belongs to, v^T.phi is performed. If the result is "> a" then it belongs
	to the right part (where phiq is), otherwise it belongs to the left
	part (where phi0 is).
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
scalar binaryNode<CompType, ThermoType>::calcA(chemPointISAT<CompType, ThermoType>* elementLeft, chemPointISAT<CompType, ThermoType>* elementRight)
{
	scalar a = 0.0;
	scalarField phih = (elementLeft->phi()+elementRight->phi())/2;
	label spaceSize = elementLeft->spaceSize();
	const scalarField& V = v();
	for (label i=0; i<spaceSize; i++)
	{
		a += V[i]*phih[i]; 
	}
//	Info << "a = " << a << endl;
	return a;
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
