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
    
    Node of the binary tree

\*---------------------------------------------------------------------------*/

#ifndef BINARY_NODE_H
#define BINARY_NODE_H

#include "chemPointISAT.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class binaryNode Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class binaryNode
{

public:

        //- Element on the left
        chemPointISAT<CompType, ThermoType>* elementLeft_;

        //- Element on the right
        chemPointISAT<CompType, ThermoType>* elementRight_;

        //- Node which follows on the left
        binaryNode<CompType, ThermoType>* left_;

        //- Node which follows on the right
        binaryNode<CompType, ThermoType>* right_;

        //- Parent node
        binaryNode<CompType, ThermoType>* parent_;
        

        /*---------------------------------------------------------------------------*\
	    'v' and 'a' are used for travelling the node during the binary tree
	    search algorithm. They define a hyper plane in the composition space
	    which split it between the two compositions linked to the node. They are 
	    initialized when the node is created.
		
	    Description: 
		1)Let E be the ellipsoid which covers the region of accuracy of 
		  the left leaf (previously defined). E is described by 
		  E={phi| ||L^T.(phi-phi0)|| <= 1}, (see chemPoint for more details).
		  let E' be the tranformation of E in a space where E' is a hypersphere
		  centered at the origin, in this space y=L^T.(phi-phi0) and then
		  E'={y| ||y||<=1}
		  let u be the unit vector joining the center of E' and the newly added 
		  composition point in the transformed space (y2=L^T.(phiq-phi0)),u = y2/||y2| 
		  Then the hyperplane separeting the two points is defined as the
		  perpendicular bisector of the segment linking 0 to y2
		 	H' = {y| u^T.(y-yh) = 0},
		  where yh = y2/2.
		  In the orginal composition space, the hyperplane H is defined by
		 	H = {y| v^T(phi-phih) = 0},
		  where phih = phi0 + L^-T.yh = (phi0 + phiq) / 2 and v is
			    L.L^T (phiq-phi0)
			v = -------------------- .
			   ||L.L^T (phiq-phi0)||
				
		2)Let a be the product v^T.phih. When travelling in the binary tree, 
		  to know in which part of the composition space the query point 'phi' 
		  belongs to, v^T.phi is performed. If the result is > a then it belongs
		  to the right part (where phiq is), otherwise it belongs to the left
		  part (where phi0 is).
				
		 3)As v multiply both phi and phih, in the implementation, v is not
		   normalised
	\*---------------------------------------------------------------------------*/
        scalarField v_;
        scalar a_;
        
 
	void calcV(chemPointISAT<CompType, ThermoType>*& elementLeft, chemPointISAT<CompType, ThermoType>*& elementRight, scalarField& v);
	scalar calcA(chemPointISAT<CompType, ThermoType>* elementLeft, chemPointISAT<CompType, ThermoType>* elementRight);
	    
//public:

    // Constructors
        
        //- Construct null
        binaryNode();
        
        //- Construct from components
        binaryNode
        (
            chemPointISAT<CompType, ThermoType>* elementLeft,
            chemPointISAT<CompType, ThermoType>* elementRight,
            binaryNode<CompType, ThermoType>* parent
        );
        //- Construct from another binary node
        binaryNode
        (
            binaryNode<CompType, ThermoType> *bn
        );
        

    // Member functions

        //- Access
        
        inline chemPointISAT<CompType, ThermoType>*& elementLeft()
        {
            return elementLeft_;
        }

        inline chemPointISAT<CompType, ThermoType>*& elementRight()
        {
            return elementRight_;
        }

        inline binaryNode<CompType, ThermoType>*& left()
        {
            return left_;
        }

        inline binaryNode<CompType, ThermoType>*& right()
        {
            return right_;
        }

        inline binaryNode<CompType, ThermoType>*& parent()
        {
            return parent_;
        }

        //- Topology

        inline const scalarField& v() const
        {
            return v_;
        }

        inline scalarField& v()
        {
            return v_;
        }
        
        inline const scalar& a() const
        {
            return a_;        
        }
        
        inline scalar& a()
        {
            return a_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "binaryNode.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
