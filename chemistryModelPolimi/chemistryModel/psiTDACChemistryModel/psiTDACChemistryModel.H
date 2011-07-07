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
    Foam::psiTDACChemistryModel

Description
    Chemistry model for density-based thermodynamics

SourceFiles
    psiTDACChemistryModelI.H
    psiTDACChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef psiTDACChemistryModel_H
#define psiTDACChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "hsCombustionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class psiTDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class psiTDACChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        psiTDACChemistryModel(const psiTDACChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const psiTDACChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<hsCombustionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("psiTDACChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiTDACChemistryModel,
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
        psiTDACChemistryModel(const fvMesh& mesh, const word& thermoTypeName);


    //- Selector
    static autoPtr<psiTDACChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~psiTDACChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline hsCombustionThermo& thermo();

        //- Return const access to the thermo package
        inline const hsCombustionThermo& thermo() const;
        
        virtual label tabSize() = 0; 
        
        virtual label tabDepth() = 0;
        
         
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "psiTDACChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
