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

#include "TDACChemistryModel.H"
#include "chemistrySolverTDAC.H"
#include "mechanismReduction.H"
#include "tabulation.H"
#include <sys/time.h>
#include "Random.H"
#include "reactingMixture.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::TDACChemistryModel<CompType, ThermoType>::TDACChemistryModel
(
    const fvMesh& mesh,
    const word& compTypeName,
    const word& thermoTypeName
)
:
    CompType(mesh, thermoTypeName),

    ODE(),

    Y_(this->thermo().composition().Y()),

    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),
    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    solver_
    (
        chemistrySolverTDAC<CompType, ThermoType>::New
        (
            *this,
            compTypeName,
            thermoTypeName
        )
    ),
    RR_(nSpecie_),
    coeffs_(nSpecie_ + 2),
    runTime_(mesh.time()),
    solveChemistryCpuTime_(0.0),
    reduceMechCpuTime_(0.0),
    searchISATCpuTime_(0.0),
    addNewLeafCpuTime_(0.0),
    isTabUsed_(false),
    NsDAC_(nSpecie_), //NsDAC is initialized to nSpecie()
    nNsDAC_(0),
    meanNsDAC_(nSpecie_),
    Ntau_(0),
    mechRed_(NULL),
    tabPtr_(NULL),
    nFound_(0),
    nGrown_(0),
    nFailBTGoodEOA_(0),
    nCellsVisited_(0),
    //by default, the solve function will check the tabulation every 1000 time-steps
    //note: this is an approximation since it use meshSize to allow the use of floating point value
    checkTab_(this->subDict("tabulation").lookupOrDefault("checkTab",1000)),
    //by default the size of the maxToComputeList corresponds to a direct treatment of not in EOA points
    maxToComputeList_(this->subDict("tabulation").lookupOrDefault("maxToComputeList",1)),
    reactionsDisabled_(this->nReaction(), false),
    DAC_(false),
    simplifiedToCompleteIndex_(nSpecie_),//maximum size of the DynamicList
    completeToSimplifiedIndex_(nSpecie_,-1), //by default it doesn't point to anything
    completeC_(nSpecie_,0.0),
    activeSpecies_(nSpecie_,false),
    specieComp_(nSpecie_),
    fuelSpecies_(),
    fuelSpeciesID_(),
    speciesNotInEOA_(),
    speciesImpact_(),
    notInEOAToGrow_(),
    notInEOAToAdd_(),
    previousTime_(mesh.time().value()),
    timeBin_(1e-5),
    curTimeBinIndex_(0),
    growOrAddImpact_(),
    growOrAddNotInEOA_(),
    analyzeTab_(this->subDict("tabulation").lookupOrDefault("analyzeTab",false)),
    exhaustiveSearch_(false),
    nbCellsVisited_(0)
{

    // create the fields for the chemistry sources
    forAll(RR_, fieldI)
    {
        RR_.set
        (
            fieldI,
            new scalarField(mesh.nCells(), 0.0)
        );
    }
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        
        // Store the species composition according to the species index
        chemkinReader tchemRead(thermoDict);
        
        const HashTable<List<chemkinReader::specieElement> >& specComp(tchemRead.specieComposition());
        forAll(specieComp_,i)
        {
            specieComp_[i] = specComp[this->Y()[i].name()];
        }
    }

    if(this->found("tabulation"))
    {
	tabPtr_ = tabulation<CompType, ThermoType>::New(*this, *this, compTypeName,thermoTypeName);
	isTabUsed_ = tabPtr_->online();
        exhaustiveSearch_.readIfPresent("exhaustiveSearch",this->subDict("tabulation"));
    }    
    else
    {
    	isTabUsed_ = false;
    }
    
    if(this->found("mechanismReduction"))
    {
	mechRed_ = mechanismReduction<CompType, ThermoType>::New(*this,*this, compTypeName,thermoTypeName);
	DAC_ = mechRed_->online();
    }
    Info<< "chemistryModel::chemistryModel: Number of species = " << nSpecie()
        << " and reactions = " << nReaction() << endl;

    //find active species

    if(DAC_)
    {
        volScalarField Ydefault
        (
            IOobject
            (
                "Ydefault",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        bool noWrite(Ydefault[0] == 0.0); 
        forAll(this->Y(), i)
        {
            IOobject header
            (
             this->Y()[i].name(),
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ
             );
            
            // check if the species file is specified if not and zero default value for species
            // we can avoid writing this species since its value is null until it is activated
            // because it is already created in basicMultiComponentMixture, we only specify NO_WRITE
            if (!header.headerOk())
            {
                if (noWrite)
                {
                    Y_[i].writeOpt()=IOobject::NO_WRITE;
                }
            }
            else 
            {
                activeSpecies_[i]=true;
            }

        }
    }
    else
    {
        forAll(this->Y(),i)
        {
            activeSpecies_[i]=true;
        }
    }  

    OFstream speciesName_(mesh.time().path()+"/speciesName.out");
    forAll(this->Y(),i)
    {
        speciesName_ << i << "    " << this->Y()[i].name() << endl;
    }
 
    //used if one would like to use the fuel to order the cell visiting order
    if(this->subDict("mechanismReduction").found("fuelSpecies"))
    {
        fuelSpecies_ = this->subDict("mechanismReduction").subDict("fuelSpecies").toc();
    }
    else 
    {
        fuelSpecies_.setSize(1);
        fuelSpecies_[0]="NC7H16";
    }

    fuelSpeciesID_.setSize(fuelSpecies_.size());
    fuelSpeciesID_[0]=0;//in case of no fuel specie specified and no nheptane
    forAll(fuelSpecies_, i)
    {
        for (label j=0; j<this->nSpecie_; j++)
        {
            if(this->Y()[j].name() == fuelSpecies_[i])
            {
                fuelSpeciesID_[i]=j;
                break;
            }
        }
    }

    if(analyzeTab_)
    {
        //initialize all variables related to ISAT analysis
        speciesNotInEOA_.append(new List<label>(nSpecie_+2,0));
        speciesImpact_.append(new List<label>(nSpecie_+2,0));
        notInEOAToGrow_.append(new List<label>(nSpecie_+2,0));
        notInEOAToAdd_.append(new List<label>(nSpecie_+2,0));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::TDACChemistryModel<CompType, ThermoType>::~TDACChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalarField Foam::TDACChemistryModel<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf,cf,pr,cr;
    label lRef, rRef;
    label omegaSize;

    //when the set of species is reduced by the DAC algorithm,
    //the size of the omega field is not equal to nEqns
    if(DAC_) omegaSize = NsDAC_+2;
    else	 omegaSize = this->nEqns();
    scalarField om(omegaSize, 0.0);

    scalarField c2(completeC_.size(), 0.0);
    if(DAC_)
    {
        //when using DAC, the ODE solver submit a reduced set of species
        //but in order to model third-body reactions properly the complete
        //set of species  is used and only the species in the simplified
        //mechanism are updated
        c2 = completeC_;
        //update the concentration of the species in the simplified mechanism
        //the other species remain the same and are used only for third-body efficiencies
        for(label i=0; i<NsDAC_; i++)
        {
            c2[simplifiedToCompleteIndex(i)] = max(0.0, c[i]);
        }
    }
    else
    {
        for(label i=0; i<this->nSpecie(); i++)
        {
            c2[i] = max(0.0, c[i]);
        }
    }

    forAll(this->reactions(), i)
    {
        if (!reactionsDisabled_[i])
        {
            const Reaction<ThermoType>& R = this->reactions()[i];
            
            scalar omegai = this->omega
            (
             R, c2, T, p, pf, cf, lRef, pr, cr, rRef
             );
            
            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                if (DAC_) si = completeToSimplifiedIndex(si);
                scalar sl = R.lhs()[s].stoichCoeff;
                om[si] -= sl*omegai;
            }
            
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                if (DAC_) si = completeToSimplifiedIndex(si);
                scalar sr = R.rhs()[s].stoichCoeff;
                om[si] += sr*omegai;
            }
        }
    } 
    return om;
} // end omega

template<class CompType, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c2(completeC_.size(), 0.0);
    for (label i=0; i<completeC_.size(); i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(T, p, c2);
    scalar kr = R.kr(kf, T, p, c2);

    pf = 1.0;
    pr = 1.0;

    label Nl = R.lhs().size();
    label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        scalar exp = R.lhs()[slRef].exponent;
        if (exp<1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        scalar exp = R.rhs()[srRef].exponent;
        if (exp<1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::TDACChemistryModel<CompType, ThermoType>::tc() const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, SMALL),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc();

    label nReaction = reactions_.size();

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            scalar rhoi = rho[celli];
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];
            scalarField c(nSpecie_);
            scalar cSum = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                cSum += c[i];
            }

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];

                omega
                (
                    R, c, Ti, pi, pf, cf, lRef, pr, cr, rRef
                );

                forAll(R.rhs(), s)
                {
                    scalar sr = R.rhs()[s].stoichCoeff;
                    tc[celli] += sr*pf*cf;
                }
            }
            tc[celli] = nReaction*cSum/tc[celli];
        }
    }


    ttc().correctBoundaryConditions();

    return ttc;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::TDACChemistryModel<CompType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Sh = tSh();

        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                scalar hi = specieThermo_[i].Hc();
                Sh[cellI] -= hi*RR_[i][cellI];
            }
        }
    }

    return tSh;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::TDACChemistryModel<CompType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class ThermoType>
Foam::label Foam::TDACChemistryModel<CompType, ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::calculate()
{
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    for (label i=0; i<nSpecie_; i++)
    {
        RR_[i].setSize(rho.size());
    }

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0.0;
            }

            scalar rhoi = rho[celli];
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];

            scalarField c(nSpecie_);
            scalarField dcdt(nEqns(), 0.0);

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
            }

            dcdt = omega(c, Ti, pi);

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = dcdt[i]*specieThermo_[i].W();
            }
        }
    }
}

template<class CompType, class ThermoType>
inline Foam::scalarField&
Foam::TDACChemistryModel<CompType, ThermoType>::coeffs()
{
    coeffs_.setSize(nEqns());
    return coeffs_;
}


template<class CompType, class ThermoType>
inline const Foam::scalarField&
Foam::TDACChemistryModel<CompType, ThermoType>::coeffs() const
{
    return coeffs_;
}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{

    scalar T = c[this->nSpecie()];
    scalar p = c[this->nSpecie() + 1];
    //the size of dcdt is c.size() (i.e. speciesNumber+2)
    scalarField tdcdt(omega(c, T, p));
    forAll(tdcdt, i)
    {
        dcdt[i] = tdcdt[i];
    }
    scalarField c2(completeC_.size(), 0.0);
    if(DAC_)
    {
        //when using DAC, the ODE solver submit a reduced set of species
	//the complete set is used and only the species in the simplified 
	//mechanism are updated
	c2 = completeC_;
		
	//update the concentration of the species in the simplified mechanism
	//the other species remain the same and are used only for third-body efficiencies
	for(label i=0; i<NsDAC_; i++)
	{
	    c2[simplifiedToCompleteIndex(i)] = max(0.0, c[i]);
	}
    }
    else
    {
	for(label i=0; i<this->nSpecie(); i++)
	{
	    c2[i] = max(0.0, c[i]);
	}
    }	

    // constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    scalar cSum = 0.0;
    for(label i=0; i<c2.size(); i++)
    {
        scalar W = this->specieThermo()[i].W();
        cSum += c2[i];
        rho += W*c2[i];
    }
    scalar mw = rho/cSum;
    scalar cp = 0.0;
	//cp is computed according to the entire set of species
	//composing the mixture
    for(label i=0; i<c2.size(); i++)
    {
        scalar cpi = this->specieThermo()[i].cp(T);
        scalar Xi = c2[i]/rho;
        cp += Xi*cpi;
    }
    cp /= mw;
    scalar dT = 0.0;
    //dT is computed on speciesNumber and not Ns since dcdt is null
    //for species not involved in the simplified mechanism
    //without DAC speciesNumber=Ns
    for(label i=0; i<this->nSpecie(); i++)
    {
	label si;
	if (DAC_) si = simplifiedToCompleteIndex(i);
	else si = i;
        scalar hi = this->specieThermo()[si].h(T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    // limit the time-derivative, this is more stable for the ODE
    // solver when calculating the allowed time step
    scalar dtMag = min(500.0, mag(dT));
    dcdt[this->nSpecie()] = -dT*dtMag/(mag(dT) + 1.0e-10);

    // dp/dt = ...
    dcdt[this->nSpecie()+1] = 0.0;

}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
	
    //if the DAC algorithm is used, the computed Jacobian
    //is compact (size of the reduced set of species)
    //but according to the informations of the complete set
    //(i.e. for the third-body efficiencies)
    scalar T = c[this->nSpecie()];
    scalar p = c[this->nSpecie() + 1];
    
    for(label i=0; i<this->nEqns(); i++)
    {
        for(label j=0; j<this->nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    //dcdt has the size of c (i.e. speciesNumber+2)
    scalarField tdcdt(omega(c, T, p));
    forAll(tdcdt, i)
    {
        dcdt[i] = tdcdt[i];
    }
    
    scalarField c2(completeC_.size(), 0.0);
    if(DAC_)
    {
        //when using DAC, the ODE solver submit a reduced set of species
        //the complete set is used and only the species in the simplified 
        //mechanism are updated
        c2 = completeC_;
        
        //update the concentration of the species in the simplified mechanism
        //the other species remain the same and are used only for third-body efficiencies
        for(label i=0; i<NsDAC_; i++)
        {
            c2[simplifiedToCompleteIndex(i)] = max(0.0, c[i]);
        }
    }
    else
    {
        for(label i=0; i<this->nSpecie(); i++)
        {
            c2[i] = max(0.0, c[i]);
        }
    }	
    
    for (label ri=0; ri<this->nReaction(); ri++)
    {
        if (!reactionsDisabled_[ri])
        {
            const Reaction<ThermoType>& R = this->reactions()[ri];
            
            scalar kf0 = R.kf(T, p, c2);
            scalar kr0 = R.kr(T, p, c2);
            
            forAll(R.lhs(), j)
            {
                label sj = R.lhs()[j].index;
                if (DAC_) sj = completeToSimplifiedIndex(sj);
                scalar kf = kf0;
                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    scalar el = R.lhs()[i].exponent;
                    if (i == j)
                    {
                        if (el < 1.0)
                        {
                            if (c2[si]>SMALL)
                            {
                                kf *= el*pow(c2[si]+VSMALL, el-1.0);
                            }
                            else
                            {
                                kf = 0.0;
                            }
                        }
                        else
                        {
                            kf *= el*pow(c2[si], el-1.0);
                        }
                    }
                    else
                    {
                        kf *= pow(c2[si], el);
                    }
                }
                
                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (DAC_) si = completeToSimplifiedIndex(si);
                    scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc[si][sj] -= sl*kf;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (DAC_) si = completeToSimplifiedIndex(si);
                    scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc[si][sj] += sr*kf;
                }
            }
            
            forAll(R.rhs(), j)
            {
                label sj = R.rhs()[j].index;
                if (DAC_) sj = completeToSimplifiedIndex(sj);
                scalar kr = kr0;
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    scalar er = R.rhs()[i].exponent;
                    if (i==j)
                    {
                        if (er<1.0)
                        {
                            if (c2[si]>SMALL)
                            {
                                kr *= er*pow(c2[si]+VSMALL, er-1.0);
                            }
                            else
                            {
                                kr = 0.0;
                            }
                        }
                        else
                        {
                            kr *= er*pow(c2[si], er-1.0);
                        }
                    }
                    else
                    {
                        kr *= pow(c2[si], er);
                    }
                }
                
                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (DAC_) si = completeToSimplifiedIndex(si);
                    scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc[si][sj] += sl*kr;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (DAC_) si = completeToSimplifiedIndex(si);
                    scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc[si][sj] -= sr*kr;
                }
            }
        }
    }

    // calculate the dcdT elements numerically
    scalar delta = 1.0e-8;
    scalarField dcdT0 = omega(c, T-delta, p);
    scalarField dcdT1 = omega(c, T+delta, p);

    for(label i=0; i<this->nEqns(); i++)
    {
        dfdc[i][this->nSpecie()] = 0.5*(dcdT1[i]-dcdT0[i])/delta;
    }

   /* // calculate the dcdp elements numerically
    scalar deltap = 1.0e-5;
    scalarField dcdp0 = omega(c2, T, p-deltap*p);
    scalarField dcdp1 = omega(c2, T, p+deltap*p);

    for(label i=0; i<nEqns(); i++)
    {
        dfdc[i][this->nSpecie()+1] = 0.5*(dcdp1[i]-dcdp0[i])/(deltap*p);
    }*/
} // end jacobian


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::setActive(label i)
{
    this->Y()[i].writeOpt()=IOobject::AUTO_WRITE;
    activeSpecies_[i]=true;
}

template<class CompType, class ThermoType>
Foam::label Foam::TDACChemistryModel<CompType, ThermoType>::tabSize()
{
    return tabPtr_->size();
}

template<class CompType, class ThermoType>
Foam::label Foam::TDACChemistryModel<CompType, ThermoType>::tabDepth()
{
    return tabPtr_->depth();
}

template<class CompType, class ThermoType>
bool Foam::TDACChemistryModel<CompType, ThermoType>::isActive(label i)
{
    return activeSpecies_[i];
}


// ************************************************************************* //
