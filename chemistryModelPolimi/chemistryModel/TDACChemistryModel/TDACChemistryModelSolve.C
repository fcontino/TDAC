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

Description
	This file implement the solve function of the chemistryModel class.
	This function uses computeA function which compute the mapping gradient
	matrix (implemented at the end of this file).
\*---------------------------------------------------------------------------*/

#include "TDACChemistryModel.H"
#include "chemistrySolverTDAC.H"
//#include "mechanismReduction.H"
//#include "tabulation.H"
#include "chemPointBase.H"
#include <sys/time.h>
#include "clockTime.H"
#include "Random.H"

/*---------------------------------------------------------------------------*\
	Solve function
 	Compute the Rates of Reaction (RR) of the species			  
	Input : initial time and time step							  
	Output: new timestep (not always used)						  
	The RR are stored in the RR_ variable returned by RR(i) funct.[kg/(m3.s)]
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<CompType, ThermoType>::solve(const scalar t0, const scalar deltaT)
{
		
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();

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
    label meshSize = rho.size();

    scalar deltaTMin = GREAT;
    scalarField Wi(this->nSpecie());
    scalarField invWi(this->nSpecie());
    for(label j=0; j<this->nSpecie(); j++)
    {
       Wi[j] = this->specieThermo()[j].W();
       invWi[j] = 1.0/this->specieThermo()[j].W();
    }

    tmp<volScalarField> thc = this->thermo().hc();
    const scalarField& hc = thc();
	
    //Update the mesh size inside chemistryModel
    label sizeOld = this->deltaTChem_.size();
    scalar minDeltaTChem(min(this->deltaTChem_));
    this->deltaTChem_.setSize(meshSize);
    
    if(sizeOld < meshSize)//new cells have been added
    {
		Info << "new size, sizeOld = " << sizeOld << ", meshSize = " << meshSize << endl;
		for(label i = sizeOld; i < meshSize; i++)
		{         
	    	this->deltaTChem_[i] = minDeltaTChem;    
		}
    }


        
    if (!this->chemistry())
    {
		return GREAT;
    }
	
	//in case of layering to avoid segmentation fault
    for(label i=0; i<this->nSpecie(); i++)
    {
		this->RR()[i].setSize(rho.size());
    }

    nFound_ = 0;
    nGrown_  = 0;

    //Random access to mesh cells to avoid problem using ISAT
    labelList cellIndexTmp = identity(meshSize);//cellIndexTmp[i]=i
    Random randGenerator(unsigned(time(NULL)));
    label j;
    for (label i=0; i<meshSize; i++)
    {
    	j=randGenerator.integer(i,meshSize-1);
        label tmp = cellIndexTmp[i];
        cellIndexTmp[i] = cellIndexTmp[j];
        cellIndexTmp[j] = tmp;
    }    


    //Start loop to solve chemistry in all cells
    reduceMechCpuTime_=0.0;
    addNewLeafCpuTime_=0.0;
    solveChemistryCpuTime_=0.0;
    searchISATCpuTime_=0.0;
    scalar computeRRTime = 0.0;        
    scalar loopInit = 0.0;
    nNsDAC_=0;
    meanNsDAC_=0;

    forAll(rho, ci)
    {	
        clockTime_.timeIncrement();
        label celli(cellIndexTmp[ci]);
        
        scalar rhoi = rho[celli];
        scalar Ti = this->thermo().T()[celli];
        scalar hi = this->thermo().hs()[celli] + hc[celli];
        scalar pi = this->thermo().p()[celli];

        scalarField phiq(this->nEqns());
        for(label i=0; i<this->nSpecie(); i++)
        {
            phiq[i] = this->Y()[i][celli];
        }

        //Species are stored in mass fraction in the cells
        //c arrays indicate the molar concentration of the species
        // c = (Y * rho)/W [kmol/m3]
        //phiq array store the mass fraction, the temperature and pressure
        //(i.e. the composition) of the query point
        scalarField c(this->nSpecie());
        scalarField c0(this->nSpecie());
        for(label i=0; i<this->nSpecie(); i++)
        {
            c[i] = rhoi*phiq[i]*invWi[i];
        }
        phiq[this->nSpecie()]=Ti;
        phiq[this->nSpecie()+1]=pi;
        
	//store the initial molar concentration to compute dc=c-c0
	c0 = c;
		
     	//time step and chemical time step
        scalar t = t0;
        scalar tauC = this->deltaTChem_[celli];
        scalar dt = min(deltaT, tauC);
        scalar timeLeft = deltaT;

        scalar cTot = 0.0;

	loopInit += clockTime_.timeIncrement();

	/*---------------------------------------------------------------------------*\
            Calculate the mapping of the query composition with the
            ISAT algorithm:
            1) phi0 the nearest point stored in the tree for a query point phiq
            according to the specific binary tree search (it is not 
            guaranted that it is the best point to approximate phiq)
            2) check if phiq lies in the region of accurate linear interpolation
            of phi0, this region is approximated by an hyper-ellipsoid 	
            (ellipsoid of accuracy EOA)
            3) then there are three possible methods:
            ==>RETRIEVE
                -phiq is in the EOA -> the mapping of phiq is obtained by linear 
                interp.: R(phiq)=phiq(t+dt)=R(phi0)+A(phi0)*(phiq-phi0),
                where A(phi0)=dR(phi0)/dphi is the mapping gradient matrix in phi0
            ==>GROW
                -phiq is not in the EOA but the mapping calculated by integrating
                the chemical equations is within the user-defined tolerance
                -> the EOA is grown in order to include the previous EOA and
                the query point phiq which ends up to be in the region of
                accuracy but not in its conservative approximation (i.e. EOA)
            ==>ADD
                -phiq is not in the region of accurate linear interpolation
                -> a new leaf is added to the tree to store the point phiq,
                the mapping R(phiq), the mapping gradient matrix A(phiq) and
                the specification of the ellipsoid of accuracy
        \*---------------------------------------------------------------------------*/
	if(isTabUsed_)
	{
	    clockTime_.timeIncrement();		
	    //phi0 will store the composition of the nearest stored point 
	    chemPointBase *phi0;
   
            //The tabulation algorithm try to retrieve the mapping
	    if (tabPtr_->retrieve(phiq,phi0))
	    {   
                nFound_ ++;
                //Rphiq array store the mapping of the query point
                scalarField Rphiq(this->nSpecie());					
                tabPtr_->calcNewC(phi0, phiq, Rphiq);
                searchISATCpuTime_ += clockTime_.timeIncrement();
                //Rphiq is in mass fraction, it is converted to molar 
                //concentration to obtain c (used to compute RR)
                for (label i=0; i<this->nSpecie(); i++) c[i] = rhoi*Rphiq[i]*invWi[i];
            }
            //Retrieve has failed. The mapping is computed
            else
            {
                searchISATCpuTime_ += clockTime_.timeIncrement();
                //When using mechanism reduction, the mechanism
                //is reduced before solving the ode including only
                //the active species
                if (DAC_) mechRed_->reduceMechanism(c, Ti, pi);
		reduceMechCpuTime_ += clockTime_.timeIncrement();

		while(timeLeft > SMALL)
		{
                    if (DAC_)
                    {
                        //The complete set of molar concentration is used even if only active species are updated

                        completeC_ = c;

                        tauC = this->solver().solve(simplifiedC_, Ti, pi, t, dt);

                        for (label i=0; i<NsDAC(); i++)  c[simplifiedToCompleteIndex(i)] = simplifiedC_[i];
                    }
                    else
                    {
			//Without dynamic reduction, the ode is directly solved
			//including all the species specified in the mechanism

			//the value of c is updated in the solve function of the chemistrySolverTDAC
                        tauC = this->solver().solve(c, Ti, pi, t, dt);
                    }
		
                    t += dt;
			
	            // update the temperature
                    cTot = sum(c);
                        
        	    ThermoType mixture(0.0*this->specieThermo()[0]);
                    for(label i=0; i<completeC_.size(); i++)
                    {
                        mixture += (c[i]/cTot)*this->specieThermo()[i];
                    }
                    Ti = mixture.TH(hi, Ti);

		    timeLeft -= dt;
                    this->deltaTChem()[celli] = tauC;
                    dt = min(timeLeft, tauC);
                    dt = max(dt, SMALL);
               }
                if (DAC_) 
		{
                    //after solving the number of species should be set back to the total number
                    nSpecie_ = mechRed_->nSpecie();
                    nNsDAC_++;
                    meanNsDAC_+=NsDAC();
                    //extend the array of active species to the full composition space
                    for (label i=0; i<NsDAC(); i++)  c[simplifiedToCompleteIndex(i)] = simplifiedC_[i];
                }
                
        	deltaTMin = min(tauC, deltaTMin);
        	
                //Rphiq array store the mapping of the query point
                scalarField Rphiq(this->nSpecie());
                //Transform c array containing the mapping in molar concentration [mol/m3]
                //to Rphiq array in mass fraction
                for(label i=0; i<this->nSpecie(); i++)
                {
                    Rphiq[i] = c[i]/rhoi*Wi[i];
                }
                solveChemistryCpuTime_ += clockTime_.timeIncrement();
				
                //check if the mapping is in the region of accurate linear interpolation
                //GROW (the grow operation is done in the checkSolution function)
                if(tabPtr_->grow(phi0, phiq, Rphiq))
                {
		    addNewLeafCpuTime_ += clockTime_.timeIncrement();
                    nGrown_ ++;
                }
                //ADD if the growth failed, a new leaf is created and added to the binary tree
                else
                {
                    //Compute the mapping gradient matrix
                    //Only computed with an add operation to avoid
                    //computing it each time
                    label Asize = this->nEqns();
                    if (DAC_) Asize = NsDAC_+2;
                    List<List<scalar> > A(Asize, List<scalar>(Asize,0.0));
                    scalarField Rcq(this->nEqns());
                    scalarField cq(this->nSpecie());					
                    for (label i=0; i<this->nSpecie(); i++)
                    {
                        Rcq[i] = rhoi*Rphiq[i]*invWi[i];
                        cq[i] = rhoi*phiq[i]*invWi[i];
                    }
                    Rcq[this->nSpecie()]=Ti;
                    Rcq[this->nSpecie()+1]=pi;
                    computeA(A, Rcq, cq, t0, deltaT, Wi, rhoi);
                    //add the new leaf which will contain phiq, R(phiq) and A(phiq)
                    //replace the leaf containing phi0 by a node splitting the
                    //composition space between phi0 and phiq (phi0 contains a reference to the node)
                    tabPtr_->add(phiq, Rphiq, A, phi0, this->nEqns());
                    addNewLeafCpuTime_ += clockTime_.timeIncrement();
                }
            }
        }//end if(isISATUsed_)
        //If ISAT is not used, direct integration is used for every cells
	else
        {
	    if (DAC_) mechRed_->reduceMechanism(c, Ti, pi);
	    while(timeLeft > SMALL)
	    {
		if (DAC_)
                {
                    //the value of c is updated in the solve function of the chemistrySolverTDAC
		    completeC_ = c;
                    tauC = this->solver().solve(simplifiedC_, Ti, pi, t, dt);
            	    for (label i=0; i<NsDAC(); i++)  c[simplifiedToCompleteIndex(i)] = simplifiedC_[i];                             
                }
                else
                {
                    //Without dynamic reduction, the ode is directly solved
                    //including all the species specified in the mechanism

                    //the value of c is updated in the solve function of the chemistrySolverTDAC
                    tauC = this->solver().solve(c, Ti, pi, t, dt);
                }
		t += dt;
		
		// update the temperature
		cTot = sum(c);
                       
       	    	ThermoType mixture(0.0*this->specieThermo()[0]);
         	for(label i=0; i<completeC_.size(); i++)
           	{
               	    mixture += (c[i]/cTot)*this->specieThermo()[i];
           	}        
           	Ti = mixture.TH(hi, Ti);
                timeLeft -= dt;
		this->deltaTChem_[celli] = tauC;
           	dt = min(timeLeft, tauC);
           	dt = max(dt, SMALL);
	    }
            if (DAC_)
            {
                //after solving the number of species should be set back to the total number
            	nSpecie_ = mechRed_->nSpecie();
                nNsDAC_++;
                meanNsDAC_+=NsDAC();
                //extend the array of active species to the full composition space
                for (label i=0; i<NsDAC(); i++)  c[simplifiedToCompleteIndex(i)] = simplifiedC_[i];
            }
	    deltaTMin = min(tauC, deltaTMin);        
        }//end of the part where it computes c after the time step
        
        
	clockTime_.timeIncrement();
	//Compute the rate of reaction according to dc=c-c0
	//In the CFD solver the following equation is solved:
	//d(Yi*rho)/dt +convection+diffusion = RR*turbulentCoeff(=1 if not used)
	//Therefore, the unit of RR should be [kg/(m3.s)]
        scalarField dc = c - c0;
 	for(label i=0; i<this->nSpecie(); i++)
	{
	     this->RR()[i][celli] = dc[i]*Wi[i]/deltaT;
	}
	computeRRTime += clockTime_.timeIncrement();
    }//End of loop over all cells
    //Display information about ISAT (if used)
    if(isTabUsed_)
    {
        scalar foundRatio =  (static_cast<scalar> (nFound_))/meshSize;
	/*Info << "Tabulation found " << foundRatio*100 << "% of the cells in the binary tree" << endl;
	Info << "Tolerance tabulation = " << tabPtr_->tolerance()<<endl;
	Info << "Chemistry library size = " ;
	Info << tabPtr_->size() << endl;
	Info << "Points Found = " << nFound_ << endl;
	Info << "Points Grown = " << nGrown_ << endl;
*/
	Pout << "Tabulation found " << foundRatio*100 << "% of the cells in the binary tree" << endl;
        Pout << "Tolerance tabulation = " << tabPtr_->tolerance()<<endl;
        Pout << "Chemistry library size = " ;
        Pout << tabPtr_->size() << endl;
        Pout << "Points Found = " << nFound_ << endl;
        Pout << "Points Grown = " << nGrown_ << endl;

    }
    if (DAC_ && nNsDAC_!=0)
        meanNsDAC_/=nNsDAC_;
    else
	meanNsDAC_=NsDAC();

    // Don't allow the time-step to change more than a factor of 2
    deltaTMin = min(deltaTMin, 2*deltaT);
    return deltaTMin;
} //end solve function





/*---------------------------------------------------------------------------*\
	Function to compute the mapping gradient matrix
	Input :	A the mapping gradient matrix (empty matrix which will contain it)
			cq the composition of the point where A is computed
				cq[0->nSpecie()-1] in [mol/m3], cq[nSpecie()]=T, cq[nSpecie()+1]=p
			t0 the initial time
			dt CFD time-step
	Output : void (the mapping gradient matrix is stored in A)
	
	The following paper describes the ode to integrate in order
	to calculate the mapping gradient matrix A:
	S.B. Pope,"Computationally efficient implementation of combustion chemistry
	using in situ adaptative tabulation", Combustion Theory Modelling, 1997
	(see p46 for the ode)
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::computeA
(
	List<List<scalar> >& A,
	const scalarField& Rcq,
		  scalarField& cq,
	const scalar& t0,
	const scalar& dt,
	const scalarField& Wi,
	const scalar& rhoi
)
{

	label speciesNumber=this->nSpecie();
	if (DAC_) speciesNumber = NsDAC_;
	//Matrix<scalar> J(speciesNumber+2, speciesNumber+2);

	jacobianForA(t0+dt, Rcq, A);
	//the jacobian is computed according to the molar concentration
	//the following conversion allow to use A with mass fraction

	for (register label i=0; i<speciesNumber; i++) 
	{	
		label si=i;
		if (DAC_) si = simplifiedToCompleteIndex(i);
		for (register label j=0; j<speciesNumber; j++)
		{	
			label sj=j;
			if (DAC_) sj = simplifiedToCompleteIndex(j);
			A[i][j] *= -dt*Wi[si]/Wi[sj];
		}
		A[i][i] += 1;
		//columns for pressure and temperature
		A[i][speciesNumber] *= -dt*Wi[si]/rhoi; 
		A[i][speciesNumber+1] *= -dt*Wi[si]/rhoi;
	}
	//Before inversion
	A[speciesNumber][speciesNumber] += 1;
	A[speciesNumber+1][speciesNumber+1] += 1;
	gaussj(A, speciesNumber+2);
	
//After inversion the last two lines of A are set to 0
// only A[this->nSpecie()][this->nSpecie()] and A[this->nSpecie()+1][this->nSpecie()+1] !=0
	A[speciesNumber][speciesNumber] = 0.0;
	A[speciesNumber+1][speciesNumber+1] = 0.0;
	
} //end computeA function


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::gaussj
(
List<List<scalar> >& A,
List<List<scalar> >& B,
label n
)
{
    //icol and irow initialized to 0 to make compiler happy (see below)
    label i, icol(0), irow(0), j, k, l, ll;
    scalar big, dum, pivinv;
    Field<label> indxc(n), indxr(n), ipiv(n);
    for (j=0; j<n; j++) ipiv[j]=0;
    for (i=0; i<n; i++)
    {
        big=0.0;
        for (j=0; j<n; j++)
        {
            if (ipiv[j] != 1)
            {
                for (k=0; k<n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(A[j][k]) >= big) //irow and icol are always initialized
                        {
                            big=fabs(A[j][k]);
                            irow=j;
                            icol=k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (l=0; l<n; l++) Swap(A[irow][l],A[icol][l]);
            for (l=0; l<n; l++) Swap(B[irow][l],B[icol][l]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (A[icol][icol] == 0.0) Info << "singular" << endl;
        pivinv = 1.0/A[icol][icol];
        A[icol][icol] = 1.0;
        for (l=0; l<n; l++) A[icol][l] *= pivinv;
        for (l=0; l<n; l++) B[icol][l] *= pivinv;
        for (ll=0; ll<n; ll++)
        {
            if (ll != icol)
            {
                dum = A[ll][icol];
                A[ll][icol] = 0.0;
                for (l=0; l<n; l++) A[ll][l] -= A[icol][l]*dum;
                for (l=0; l<n; l++) B[ll][l] -= B[icol][l]*dum;				
            }
        }
    }
    for (l=n-1; l>=0; l--)
    {
        if (indxr[l] != indxc[l])
        {
            for (k=0; k<n; k++) Swap (A[k][indxr[l]],A[k][indxc[l]]);
        }
    }
}//end gaussj(A,B)

template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::gaussj
(
List<List<scalar> >& A,
label n
)
{
	List<List<scalar> > B(n,List<scalar>(n, 0.0));
	gaussj(A,B, n);
}

template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::jacobianForA
(
    const scalar t,
    const scalarField& c2,
    List<List<scalar> >& dfdc
) const
{
	//if the DAC algorithm is used, the computed Jacobian
	//is compact (size of the reduced set of species)
	//but according to the informations of the complete set
	//(i.e. for the third-body efficiencies)
	label speciesNumber;
	if (DAC_) speciesNumber = NsDAC_;
	else speciesNumber = this->nSpecie();
	
    scalar T = c2[this->nSpecie()];
    scalar p = c2[this->nSpecie() + 1];
    for(label i=0; i<speciesNumber+2; i++)
    {
        for(label j=0; j<speciesNumber+2; j++)
        {
            dfdc[i][j] = 0.0;
        }
    }
	
	
    for (label ri=0; ri<this->reactions().size(); ri++)
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
	scalarField dcdT0;
	scalarField dcdT1;
	if (DAC_)
	{
		scalarField c1(speciesNumber,0.0);
		for (label i=0; i<speciesNumber; i++) c1[i] = c2[simplifiedToCompleteIndex(i)];
		dcdT0 = this->omega(c1, T-delta, p);
		dcdT1 = this->omega(c1, T+delta, p);
	}
	else
	{
		dcdT0 = this->omega(c2, T-delta, p);
		dcdT1 = this->omega(c2, T+delta, p);
	}

    for(label i=0; i<speciesNumber+2; i++)
    {
        dfdc[i][speciesNumber] = 0.5*(dcdT1[i]-dcdT0[i])/delta;
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

