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
    Foam::chemistrySolverTDAC

Description
    An abstract base class for solving chemistry

SourceFiles
    chemistrySolverTDAC.C

\*---------------------------------------------------------------------------*/

#ifndef chemistrySolverTDAC_H
#define chemistrySolverTDAC_H

#include "TDACChemistryModel.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CompType, class ThermoType>
class chemistrySolverTDAC;

/*---------------------------------------------------------------------------*\
                      Class chemistrySolverTDAC Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistrySolverTDAC
{
protected:

    // Protected data

        //- Reference to the chemistry model
        TDACChemistryModel<CompType, ThermoType>& model_;

        //- Name of the chemistry solver
        const word name_;


public:

        //- Runtime type information
        TypeName("chemistrySolverTDAC");


        // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistrySolverTDAC,
            dictionary,
            (
                TDACChemistryModel<CompType, ThermoType>& model,
                const word& modelName
            ),
            (model, modelName)
        );


    // Constructors

        //- Construct from components
        chemistrySolverTDAC
        (
            TDACChemistryModel<CompType, ThermoType>& model,
            const word& modelName
        );


    //- Selector
    static autoPtr<chemistrySolverTDAC> New
    (
        TDACChemistryModel<CompType, ThermoType>& model,
        const word& compTypeName,
        const word& thermoTypeName
    );


    //- Destructor
    virtual ~chemistrySolverTDAC();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual scalar solve
        (
            scalarField &c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistrySolverTDAC(Comp, Thermo)                                 \
                                                                              \
    typedef chemistrySolverTDAC<Comp, Thermo>                                 \
        chemistrySolverTDAC##Comp##Thermo;                                    \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        chemistrySolverTDAC##Comp##Thermo,                                    \
        "chemistryModel<"#Comp","#Thermo">",                                  \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        chemistrySolverTDAC##Comp##Thermo,                                    \
        dictionary                                                            \
    );


#define makeChemistrySolverTDACType(SS, Comp, Thermo)                         \
                                                                              \
    typedef SS<Comp, Thermo> SS##Comp##Thermo;                                \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##Thermo,                                                     \
        #SS"<"#Comp","#Thermo">",                                             \
        0                                                                     \
    );                                                                        \
                                                                              \
    chemistrySolverTDAC<Comp, Thermo>::                                       \
        adddictionaryConstructorToTable<SS<Comp, Thermo> >                    \
            add##SS##Comp##Thermo##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistrySolverTDAC.C"
#   include "newChemistrySolverTDAC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
