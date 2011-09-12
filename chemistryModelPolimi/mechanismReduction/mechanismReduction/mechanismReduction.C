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

#include "mechanismReduction.H"
#include "Switch.H"
#include "error.H"

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class CompType, class ThermoType>
Foam::mechanismReduction<CompType,ThermoType>::mechanismReduction
(
    const Foam::dictionary& dict,
    Foam::TDACChemistryModel<CompType,ThermoType>& chemistry
)
:
    dict_(dict),
    chemistry_(chemistry),
    activeSpecies_(chemistry.nSpecie(),false),
    NsSimp_(chemistry.nSpecie()),
    nSpecie_(chemistry.nSpecie()),
    coeffsDict_(dict.subDict("mechanismReduction")),
    epsDAC_(readScalar(coeffsDict_.lookup("epsDAC"))),
    initSet_(coeffsDict_.subDict("initialSet")),
    searchInitSet_(initSet_.size()),
    online_(coeffsDict_.lookup("online"))
{
    label j=0;
    for (label i=0; i<chemistry.nSpecie(); i++)
    {
        if(initSet_.found(chemistry.Y()[i].name())) searchInitSet_[j++]=i;
    }
    if (j<searchInitSet_.size())
    {
        FatalErrorIn("Foam::mechanismReduction::mechanismReduction(const Foam::dictionary& dict,Foam::TDACChemistryModel<CompType,ThermoType>& chemistry)")
            << "At least one species in the intial set is not in the mechanism "
            << abort(FatalError);
    }
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::mechanismReduction<CompType,ThermoType>::~mechanismReduction()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
