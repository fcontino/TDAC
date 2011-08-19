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
    Foam::tabulation

Description
    An abstract class for tabulation

SourceFiles
    tabulation.C

\*---------------------------------------------------------------------------*/

#ifndef tabulation_H
#define tabulation_H

#include "TDACChemistryModel.H"
#include "chemPointBase.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class tabulation Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class tabulation
{


protected:

        const dictionary& dict_;
        TDACChemistryModel<CompType, ThermoType>& chemistry_;
     	//Dictionary that store the algorithm data
        const dictionary coeffsDict_;
        //Is tabulation active?
        const Switch online_;


        
private:
        
        

public:

        //- Runtime type information
        TypeName("tabulation");


        // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            tabulation,
            dictionary,
            (
                const dictionary& dict,
                TDACChemistryModel<CompType, ThermoType>& chemistry
            ),
            (dict, chemistry)
        );


    // Constructors

        //- Construct from components
        tabulation
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selectors

        static autoPtr<tabulation> New
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry,
            const word& compTypeName,
            const word& thermoTypeName
        );


    // Destructor

        virtual ~tabulation();


    // Member Functions

        inline Switch online()
        {
            return online_;
        }


    // Virtual functions 

	virtual const scalar& tolerance() const = 0;
	
	virtual const scalar& checkUsed() const = 0;

	virtual const label& checkGrown() = 0;

	virtual Switch clean() const = 0;

	virtual label size() = 0;

	virtual label depth() = 0;

	virtual void add
        (       
            const scalarField&,
            const scalarField&,
                List<List<scalar> >&,
            const label
        ) = 0;

	virtual void add
	(
	    const scalarField&,
	    const scalarField&,
		List<List<scalar> >&,
		chemPointBase*&,
	    const label
	) = 0;

        virtual void replace
        (
                chemPointBase*&,
            const scalarField&,
            const scalarField&,
                List<List<scalar> >&,
            const label
        ) = 0;

        virtual void calcNewC
        (
                chemPointBase*&,
            const scalarField&,
            	scalarField&
        ) = 0;

	virtual bool grow
        (
                chemPointBase*&,
            const scalarField&,
	    const scalarField&
        ) = 0;

	virtual void clear() = 0;

	virtual bool retrieve
        (
            const scalarField&,
                chemPointBase*&
        ) = 0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeTabulation(Comp, Thermo)                                          \
                                                                              \
    typedef tabulation<Comp, Thermo>                                          \
        tabulation##Comp##Thermo;                                             \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        tabulation##Comp##Thermo,                                             \
        "chemistryModel<"#Comp","#Thermo">",                                  \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        tabulation##Comp##Thermo,                                             \
        dictionary                                                            \
    );


#define makeTabulationType(SS, Comp, Thermo)                                  \
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
    tabulation<Comp, Thermo>::                                                \
        adddictionaryConstructorToTable<SS<Comp, Thermo> >                    \
            add##SS##Comp##Thermo##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tabulation.C"
#   include "newTabulation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
