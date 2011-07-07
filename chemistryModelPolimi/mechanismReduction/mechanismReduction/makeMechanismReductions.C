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
#include "mechanismReduction.H"

#include "psiTDACChemistryModel.H"
#include "rhoTDACChemistryModel.H"

#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeMechanismReduction(psiTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(DAC, psiTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(DRG, psiTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(DRGEP, psiTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(EFA, psiTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(PFA, psiTDACChemistryModel, gasThermoPhysics)
                  
    makeMechanismReduction(psiTDACChemistryModel, icoPoly8ThermoPhysics)
    makeMechanismReductionType(DAC,psiTDACChemistryModel,icoPoly8ThermoPhysics)
    makeMechanismReductionType(DRG,psiTDACChemistryModel,icoPoly8ThermoPhysics)
    makeMechanismReductionType(DRGEP,psiTDACChemistryModel,icoPoly8ThermoPhysics)
    makeMechanismReductionType(EFA,psiTDACChemistryModel,icoPoly8ThermoPhysics)
    makeMechanismReductionType(PFA,psiTDACChemistryModel,icoPoly8ThermoPhysics)
    
    makeMechanismReduction(rhoTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(DAC, rhoTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(DRG, rhoTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(DRGEP, rhoTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(EFA, rhoTDACChemistryModel, gasThermoPhysics)
    makeMechanismReductionType(PFA, rhoTDACChemistryModel, gasThermoPhysics)
                   
    makeMechanismReduction(rhoTDACChemistryModel, icoPoly8ThermoPhysics)
    makeMechanismReductionType(DAC, rhoTDACChemistryModel, icoPoly8ThermoPhysics)
    makeMechanismReductionType(DRG, rhoTDACChemistryModel, icoPoly8ThermoPhysics)
    makeMechanismReductionType(DRGEP, rhoTDACChemistryModel, icoPoly8ThermoPhysics)
    makeMechanismReductionType(EFA, rhoTDACChemistryModel, icoPoly8ThermoPhysics)
    makeMechanismReductionType(PFA, rhoTDACChemistryModel, icoPoly8ThermoPhysics)                
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
