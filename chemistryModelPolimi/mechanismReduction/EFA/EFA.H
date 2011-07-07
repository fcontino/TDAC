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

Class
    Foam::EFA

Description
    

SourceFiles
    EFA.C

\*---------------------------------------------------------------------------*/

#ifndef EFA_H
#define EFA_H

#include "mechanismReduction.H"
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class EFA
:
    public mechanismReduction<CompType, ThermoType>
{
    // Private data
    List<label> sC_,sH_,sO_,sN_;
    scalar  sortPart_;

public:

    //- Runtime type information
    TypeName("EFA");


    // Constructors

        //- Construct from components
        EFA
        (
            const dictionary& dict,
            TDACChemistryModel<CompType,ThermoType>& chemistry
        );


    // Destructor

        ~EFA();


    // Member Functions

        //- Reduce the mechanism
        void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        ) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EFA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
