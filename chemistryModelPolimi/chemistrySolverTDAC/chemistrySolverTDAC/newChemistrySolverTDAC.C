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

#include "chemistrySolverTDAC.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::autoPtr<Foam::chemistrySolverTDAC<CompType, ThermoType> >
Foam::chemistrySolverTDAC<CompType, ThermoType>::New
(
    TDACChemistryModel<CompType, ThermoType>& model,
    const word& compTypeName,
    const word& thermoTypeName
)
{
    word modelName(model.lookup("chemistrySolverTDAC"));

    word chemistrySolverTDACType =
        modelName + '<' + compTypeName + ',' + thermoTypeName + '>';

    Info<< "Selecting chemistrySolverTDAC " << modelName << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistrySolverTDACType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        wordList models = dictionaryConstructorTablePtr_->toc();
        forAll(models, i)
        {
            models[i] = models[i].replace
            (
                '<' + compTypeName + ',' + thermoTypeName + '>',
                ""
            );
        }

        FatalErrorIn
        (
            "chemistrySolverTDAC::New"
            "("
                "const TDACChemistryModel&, "
                "const word&, "
                "const word&"
            ")"
        )   << "Unknown chemistrySolverTDAC type " << modelName
            << nl << nl << "Valid chemistrySolverTDAC types are:" << nl
            << models << nl << exit(FatalError);
    }

    return autoPtr<chemistrySolverTDAC<CompType, ThermoType> >
        (cstrIter()(model, modelName));
}


// ************************************************************************* //
