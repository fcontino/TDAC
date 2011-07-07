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
    TDACChemistryModel

Description

SourceFiles
    TDACChemistryModelI.H
    TDACChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef TDACChemistryModel_H
#define TDACChemistryModel_H

#include "chemkinReader.H"
#include "Reaction.H"
#include "ODE.H"
#include "volFieldsFwd.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

template<class CompType, class ThermoType>
class chemistrySolverTDAC;
template<class CompType, class ThermoType>
class mechanismReduction;
template<class CompType, class ThermoType>
class tabulation;
/*---------------------------------------------------------------------------*\
                           Class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class TDACChemistryModel
:
    public CompType,
    public ODE
{

private:

    // Private data
        
        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reactions
        const PtrList<Reaction<ThermoType> >& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- Chemistry solver
        autoPtr<chemistrySolverTDAC<CompType, ThermoType> > solver_;

        //- Chemical source term [kg/m3/s]
        PtrList<scalarField> RR_;

        //- ODE coefficients
        scalarField coeffs_;

        
	const Time& runTime_;
	scalar solveChemistryCpuTime_;
	scalar reduceMechCpuTime_;
	scalar searchISATCpuTime_;
	scalar addNewLeafCpuTime_;

	//- Use the tabulation switch
	Switch isTabUsed_;
	
	//- Keep track of the number of species when the DAC algorithm is used
	label NsDAC_;
	label nNsDAC_;
	label meanNsDAC_;

	//- Number of assumed QSS species
	label Ntau_;


	//- mechanism reduction algorithm
	autoPtr<mechanismReduction<CompType, ThermoType> > mechRed_;

        //- the tabulation method
        autoPtr<tabulation<CompType, ThermoType> > tabPtr_;

        //- Points found in the online library
        label nFound_;
        
        //- Chem points grown by the library
        label nGrown_;
	
	//- List of bool to disable reactions
	Field<bool> reactionsDisabled_;
		
	//- Use DAC algorithm during solving
	Switch DAC_;
	
	//- Index in the complete mechanism of species contained in the simplified mechanism  
	DynamicList<label> simplifiedToCompleteIndex_;
	
	//- Index in the simplified mechanism of species contained in the complete mechanism  
	Field<label> completeToSimplifiedIndex_;
			
	//- Lists used to store the complete and simplified array of concentration
	//	(used when DAC is active)
	scalarField completeC_;	
	scalarField simplifiedC_;

	//- List of active species
	List<bool> activeSpecies_;
        
        List<List<chemkinReader::specieElement> > specieComp_;
        
 
    // Private Member Functions

	/*---------------------------------------------------------------------------*\
	    Function to compute the mapping gradient matrix
	    Input :	A the mapping gradient matrix (empty matrix which will contain it)
			cq the composition of the point where A is computed
			    cq[0->Ns()-1] in [mol/m3], cq[Ns()]=T, cq[Ns()+1]=p
			t0 the initial time
			dt CFD time-step
	    Output : void (the mapping gradient matrix is stored in A)
	
	    The following paper describes the ode to integrate in order
	    to calculate the mapping gradient matrix A:
	    S.B. Pope,"Computationally efficient implementation of combustion chemistry
	    using in situ adaptative tabulation", Combustion Theory Modelling, 1997
	    (see p46 for the ode)
	    The assumption that the jacobian is constant during the time-step is used
	     to integrate it explicitely:
	    !!Assumption need to be validated!!
	    A(phiq) = J(phiq)*dt + I 
	\*---------------------------------------------------------------------------*/
	void computeA
	(	
		  List<List<scalar> >& A, 
	    const scalarField& Rcq,
		  scalarField& cq, 
	    const scalar& t0, 
	    const scalar& dt,
	    const scalarField& Wi,
	    const scalar& rhoi
	);
		
	//Gauss Jordan elimination
	void gaussj
	(
	    List<List<scalar> >& A,
	    List<List<scalar> >& B,
	    label n
	);
		
	void gaussj
	(
	    List<List<scalar> >& A,
	    label n
	);
		
	//Compute the jacobian for the A matrix
	//(without computing dcdt)
	void jacobianForA
        (
            const scalar t,
            const scalarField& c,
            List<List<scalar> >& dfdc
        ) const;
		
        //- Disallow default bitwise assignment
        void operator=(const TDACChemistryModel<CompType, ThermoType>&);

protected:

    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<scalarField>& RR();

public:

    //- Runtime type information
    TypeName("TDACChemistryModel");

    // Constructors

        //- Construct from components
        TDACChemistryModel
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& thermoTypeName
        );


    // Destructor

        virtual ~TDACChemistryModel();

    // Member Functions
    
        //- Access the fractions of the species
        inline PtrList<volScalarField>& Y() const;

        //- The reactions
        inline const PtrList<Reaction<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;
        
        label& nSpecie()
        {
            return nSpecie_;
        }

        //- The number of reactions
        inline label nReaction() const;

        //- Return the chemisty solver
        inline const chemistrySolverTDAC<CompType, ThermoType>& solver() const;    
    
	//- CpuTime
	inline scalar solveChemistryCpuTime()
	{
	    return solveChemistryCpuTime_;	
	}
        
        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;
        
        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

	
        inline scalar reduceMechCpuTime()
	{
	    return reduceMechCpuTime_;
	}

	inline scalar searchISATCpuTime()
	{
	    return searchISATCpuTime_;
	}

	inline scalar addNewLeafCpuTime()
	{
	    return addNewLeafCpuTime_;
	}

	inline void  NsDAC(label newNsDAC)
        {
	    NsDAC_ = newNsDAC;
        }
	
	label NsDAC() const
	{
	    return NsDAC_;
	}
	
	inline label& Ntau()
	{
	    return Ntau_;
	}	

	inline label meanNsDAC() const
	{
	    return meanNsDAC_;
	}	

	Switch DAC() const
	{
	    return DAC_;
	}
		
	inline label& simplifiedToCompleteIndex(label i)
	{
	    return simplifiedToCompleteIndex_[i];
	}
		
	inline DynamicList<label>& simplifiedToCompleteIndex()
	{
	    return simplifiedToCompleteIndex_;
	}

        inline Field<label>& completeToSimplifiedIndex()
        {
            return completeToSimplifiedIndex_;
        }
		
	inline label& completeToSimplifiedIndex(label i)
	{
	    return completeToSimplifiedIndex_[i];
	}
	
	inline const label& simplifiedToCompleteIndex(label i) const
	{
	    return simplifiedToCompleteIndex_[i];
	}
		
	inline const label& completeToSimplifiedIndex(label i) const
	{
	    return completeToSimplifiedIndex_[i];
	}

	inline const Field<label>& completeToSimplifiedIndex() const
	{
	    return completeToSimplifiedIndex_;
	}
	
	inline Field<bool>& reactionsDisabled() 
	{
	    return reactionsDisabled_;
	}
	
        inline scalarField& completeC()
	{
	    return completeC_;
	}
	
	inline scalarField& simplifiedC()
	{
	    return simplifiedC_;
	}
	
        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms
            inline tmp<volScalarField> RR(const label i) const;

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            //- Return the heat release, i.e. enthalpy/sec [m2/s3]
            virtual tmp<volScalarField> dQ() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const;

            //- Return access to coefficients
            virtual scalarField& coeffs();

            //- Return reference to coefficients
            virtual const scalarField& coeffs() const;
                
                        	
        virtual void derivatives
        (
            const scalar t,
            const scalarField& c,
            scalarField& dcdt
        ) const;

        virtual void jacobian
        (
            const scalar t,
            const scalarField& c,
            scalarField& dcdt,
            scalarSquareMatrix& dfdc
        ) const;

        /*---------------------------------------------------------------------------*\
	    Solve function
	    Compute the Rates of Reaction (RR) of the species			  
	    Input : initial time and time step							  
	    Output: characteristic time
   	    The RR are stored in the RR_ variable returned by RR(i) funct.[kg/(m3.s)]
	\*---------------------------------------------------------------------------*/
        virtual scalar solve
        (
            const scalar t0,
            const scalar deltaT
        );

        const label& nFound() const
        {
            return nFound_;
        }
        
        const label& nGrown() const
        {
            return nGrown_;
        }

        label tabSize();
	label tabDepth();
	
	//set species Y[i] to active
	void setActive(label i);
	
	inline bool isActive(label i)
	{
	    return activeSpecies_[i];
	}

        inline List<List<chemkinReader::specieElement> >& specieComp()
        {
            return specieComp_;
        }
        
        //- Update ODE after the solution, advancing by delta
        virtual void update(const scalar delta)
        {}
    
    
	
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TDACChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TDACChemistryModel.C"
#   include "TDACChemistryModelSolve.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
