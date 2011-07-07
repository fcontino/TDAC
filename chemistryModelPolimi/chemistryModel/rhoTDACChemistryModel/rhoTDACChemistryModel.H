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
    Foam::rhoTDACChemistryModel

Description
    Chemistry model for density-based thermodynamics

SourceFiles
    rhoTDACChemistryModelI.H
    rhoTDACChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef rhoTDACChemistryModel_H
#define rhoTDACChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "hsReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class rhoTDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class rhoTDACChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoTDACChemistryModel(const rhoTDACChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const rhoTDACChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<hsReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("rhoTDACChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoTDACChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& thermoTypeName
        ),
        (mesh, compTypeName, thermoTypeName)
    );


    // Constructors

        //- Construct from mesh
        rhoTDACChemistryModel(const fvMesh& mesh, const word& thermoTypeName);


    //- Selector
    static autoPtr<rhoTDACChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~rhoTDACChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline hsReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const hsReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoTDACChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
