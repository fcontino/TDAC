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
    Foam::mechanismReduction

Description
    An abstract class for reducing chemical mechanism

SourceFiles
    mechanismReduction.C

\*---------------------------------------------------------------------------*/

#ifndef mechanismReduction_H
#define mechanismReduction_H

#include "TDACChemistryModel.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class mechanismReduction
{


protected:

        const dictionary& dict_;
        TDACChemistryModel<CompType, ThermoType>& chemistry_;

        //List of active species (active = true)
        List<bool> activeSpecies_;
        
        //Number of active species
        label NsSimp_;
        
        //Number of species
        const label nSpecie_;
        
        //Dictionary that store the algorithm data
        const dictionary coeffsDict_;
        
private:

        
        //Tolerance for the mechanism reduction algorithm
        const scalar epsDAC_;
        
        //Dictionary that store the initial set
        const dictionary initSet_;
        
        //List of label for the search initiating set
        labelList searchInitSet_;
        
        //Is mechanism reduction active?
        const Switch online_;

public:

        //- Runtime type information
        TypeName("mechanismReduction");


        // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanismReduction,
            dictionary,
            (
                const dictionary& dict,
                TDACChemistryModel<CompType, ThermoType>& chemistry
            ),
            (dict, chemistry)
        );


    // Constructors

        //- Construct from components
        mechanismReduction
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selectors

        static autoPtr<mechanismReduction<CompType, ThermoType> > New
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry,
            const word& compTypeName,
            const word& thermoTypeName
        );


    // Destructor

        virtual ~mechanismReduction();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        )  = 0;
	
	//- Return the active species
	inline const List<bool>& activeSpecies() const
	{
	    return activeSpecies_;
	} 

        
        //- Return the number of active species
        inline label NsSimp()
        {
            return NsSimp_;
        }

        //- Return the initial number of species
        inline label nSpecie()
        {
            return nSpecie_;
        }        
        
        //- Return the tolerance
        inline scalar epsDAC() const
        {
            return epsDAC_;
        }
        
        inline const labelList& searchInitSet() const
        {
            return searchInitSet_;
        } 
        
        inline Switch online()
        {
            return online_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMechanismReduction(Comp, Thermo)                                  \
                                                                              \
    typedef mechanismReduction<Comp, Thermo>                                  \
        mechanismReduction##Comp##Thermo;                                     \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        mechanismReduction##Comp##Thermo,                                     \
        "chemistryModel<"#Comp","#Thermo">",                                  \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        mechanismReduction##Comp##Thermo,                                     \
        dictionary                                                            \
    );


#define makeMechanismReductionType(SS, Comp, Thermo)                          \
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
    mechanismReduction<Comp, Thermo>::                                        \
        adddictionaryConstructorToTable<SS<Comp, Thermo> >                    \
            add##SS##Comp##Thermo##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mechanismReduction.C"
#   include "newMechanismReduction.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
