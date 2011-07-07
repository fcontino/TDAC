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
    Foam::PFA

Description
    Path flux analysis

SourceFiles
    PFA.C

\*---------------------------------------------------------------------------*/

#ifndef PFA_H
#define PFA_H

#include "mechanismReduction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class PFA
:
    public mechanismReduction<CompType, ThermoType>
{
    // Private data



public:

    //- Runtime type information
    TypeName("PFA");


    // Constructors

        //- Construct from components
        PFA
        (
            const dictionary& dict,
            TDACChemistryModel<CompType,ThermoType>& chemistry
        );


    // Destructor

        ~PFA();


    // Member Functions

        //- Reduce the mechanism
        void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PFA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
