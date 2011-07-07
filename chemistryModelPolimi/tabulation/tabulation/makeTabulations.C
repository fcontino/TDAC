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

#include "thermoPhysicsTypes.H"
#include "tabulation.H"

#include "psiTDACChemistryModel.H"
#include "rhoTDACChemistryModel.H"

#include "ISAT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTabulation(psiTDACChemistryModel, gasThermoPhysics)
    makeTabulationType(ISAT, psiTDACChemistryModel, gasThermoPhysics)
            
    makeTabulation(psiTDACChemistryModel, icoPoly8ThermoPhysics)
    makeTabulationType(ISAT,psiTDACChemistryModel,icoPoly8ThermoPhysics)
    
    makeTabulation(rhoTDACChemistryModel, gasThermoPhysics)
    makeTabulationType(ISAT, rhoTDACChemistryModel, gasThermoPhysics)
   
    makeTabulation(rhoTDACChemistryModel, icoPoly8ThermoPhysics)
    makeTabulationType(ISAT, rhoTDACChemistryModel, icoPoly8ThermoPhysics)
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
