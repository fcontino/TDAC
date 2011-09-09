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

    Implementation of the ISAT (In-situ adaptive tabulation), for chemistry 
    calculation. For forward information:
    
    S.B. Pope (1997) 
    ``Computationally efficient implementation of combustion chemistry using 
    in situ adaptive tabulation,'' 
    Combustion Theory and Modelling,  1, 41--63.

\*---------------------------------------------------------------------------*/

#ifndef ISAT_H
#define ISAT_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulation.H"
#include "IOstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "binaryTree.H" 
#include "Time.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ISAT Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class TDACChemistryModel;

template<class CompType, class ThermoType>
class ISAT
:
    public tabulation<CompType, ThermoType>
{

    // Private data
   
//        //- Dictionary for data lookup
//        const dictionary& onlineDict_;

private:
        //- Reference to the chemistryModel
        TDACChemistryModel<CompType, ThermoType>& chemistry_;

        //- List of the stored 'points' organized in a binary tree
        binaryTree<CompType, ThermoType> chemisTree_;
    
        //- Tolerance for the ISAT algorithm
        scalar tolerance_;

        //- List of scale factors for species, temperature and pressure
        scalarField scaleFactor_;
       
        //- how to deal with time
        Switch tauStar_;
        
        //- Clean or not to Clean?
        Switch clean_;
        
        //- Perform a check EOA every "checkUsed*meshSize" use of a leaf
        scalar checkUsed_;
        
        //- Perform a check EOA every "checkGrown*meshSize" grows of a leaf
        label checkGrown_; 
        
        //- Most Recently Used (MRU) list of chemPoint
        SLList<chemPointISAT<CompType, ThermoType>*> MRUList_;
        
        //- User defined size of the MRU list
	label MRUSize_;
        
        //- Variable that raise a flag when a cleaning operation is required
        bool cleaningRequired_;
        
        DynamicList<chemPointISAT<CompType, ThermoType>*> toRemoveList_;

        //- Number of points failed to be found by primary retrieve
        //  but found by either secondary retrieve or in the MRU list
        label nFailedFirst_;
        label totRetrieve_;
        
        //- After a failed primary retrieve, loo in the MRU list
        bool MRURetrieve_;
        
        //- Max proportion of the retrieve done by secondary retrieve before calling balance function
        scalar max2ndRetBalance_;
        scalar maxDepthFactor_;
        
        const Time* runTime_;
        scalar previousTime_;
        scalar checkEntireTreeInterval_;
        label chPMaxLifeTime_;
        label chPMaxUseInterval_;
        
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        ISAT(const ISAT&);

        //- Disallow default bitwise assignment
        void operator=(const ISAT&);
        
        //- Add to MRUList
        void addToMRU(chemPointISAT<CompType, ThermoType>* phi0);
		
	


public:

    //- Runtime type information
    TypeName("ISAT");

    // Constructors

        //- Construct from dictionary
        ISAT
        (
            const dictionary& chemistryProperties,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selectors

    // Destructor

        ~ISAT();


    // Member Functions

        // Access
        inline binaryTree<CompType, ThermoType>& chemisTree() 
        {
            return chemisTree_;
        }

        inline const scalarField& scaleFactor() const
        {
            return scaleFactor_;
        }
        
        inline const scalar& tolerance() const
        {
            return tolerance_;
        }
        
        inline const scalar& checkUsed() const
        {
            return checkUsed_;
        }
        
        inline Switch clean() const
        {
            return clean_;
        }

        // Database

        //- Return the size of the binary tree
        inline label size()
        {
            return chemisTree_.size();
        }
        

        //- Return the number of grows allowed before replacing the leaf
        inline const label& checkGrown()
        {
            return checkGrown_;
        }

	//- Return the depth of the binary tree
	inline label depth()
	{
	    return chemisTree_.depth();
	}
        
        inline bool& cleaningRequired()
        {
            return cleaningRequired_;
        }
        
        // Edit
				
        /*---------------------------------------------------------------------------*\
                Add a new leaf to the binary tree
                Input : phiq the new composition to store
                                Rphiq the mapping of the new composition point
                                A the mapping gradient matrix
                                phi0 the chemPoint which is the nearest from phiq and which will
                                        be replaced by a node splitting the composition space between
                                        phi0 and phiq
                                nCols the size of the matrix
                Output: void
        \*---------------------------------------------------------------------------*/
        bool add
        (
                const scalarField& phiq,    
                const scalarField& Rphiq, 
                      List<List<scalar> >& A,
                      chemPointBase*& phi0,
                label nCols
        );
        
        /*---------------------------------------------------------------------------*\
            Compute and return the mapping of the composition phiq from stored data
            Input : phi0 the nearest chemPoint used in the linear interpolation
                phiq the composition of the query point for which we want to 
                        compute the mapping
                nEqns the size of the composition space
                Rphiq the mapping of the new composition point (given as empty)
            Output: void (the mapping is stored in the Rphiq array)
        \*---------------------------------------------------------------------------*/


        void calcNewC
        (
            chemPointBase*& phi0,
            const scalarField& phiq,
            	  scalarField& Rphiq
        );
	
	/*---------------------------------------------------------------------------*\
	    Perform the grow operation 
	    Input : phi0 the nearest chemPoint used in the linear interpolation
	            phiq the composition of the query point for which we want to 
			compute the mapping
		    Rphiq
	    Output: true if the EOA of phi0 has been grown
	\*---------------------------------------------------------------------------*/
	bool grow
	(
	    chemPointBase*& phi0,
	    const scalarField& phiq,
	    const scalarField& Rphiq
	);
	
	
        //- Clear the library
        void clear();
        
        
        //- Find the closest stored leaf according to the specified composition
	bool retrieve 
	(
	    const Foam::scalarField& v0,
                  chemPointBase*& closest
	);
        
        
        //- Clean and balance the tree if needed
        bool cleanAndBalance();
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ISAT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
