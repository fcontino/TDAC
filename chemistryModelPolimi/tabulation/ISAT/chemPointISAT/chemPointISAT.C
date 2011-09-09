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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "chemPointISAT.H"
#include "binaryNode.H"
#include "TDACChemistryModel.H"
#include <limits>
#include "clockTime.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
Foam::scalar Foam::chemPointISAT<CompType, ThermoType>::epsTol_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
namespace Foam
{

template<class CompType, class ThermoType>
chemPointISAT<CompType, ThermoType>::chemPointISAT
(
TDACChemistryModel<CompType, ThermoType>& chemistry,
const scalarField& phi,
const scalarField& Rphi,
const List<List<scalar> >& A,
const scalarField& scaleFactor,
const scalar& epsTol,
const label& spaceSize,
binaryNode<CompType, ThermoType>* node
)
:
    chemistry_(&chemistry),
    phi_(phi),
    Rphi_(Rphi),
    A_(A),
    scaleFactor_(scaleFactor),
    node_(node),
    spaceSize_(spaceSize),
     nUsed_(0),
    nGrown_(0),    
    DAC_(chemistry.DAC()),
    NsDAC_(chemistry.NsDAC()),
    completeToSimplifiedIndex_(spaceSize-2),
    simplifiedToCompleteIndex_(NsDAC_),
    inertSpecie_(-1),
    timeTag_(chemistry_->time().timeOutputValue()),
    lastTimeUsed_(chemistry_->time().timeOutputValue()),
    lastError_(0.0)/*,
    failedSpeciesFile_(chemistry.thermo().T().mesh().time().path()+"/failedSpecies.out"),
    failedSpecies_(failedSpeciesFile_.c_str(), ofstream::app),
    refTime_(&chemistry.thermo().T().mesh().time())*/
{
    epsTol_=epsTol;
    clockTime cpuCP = clockTime();
    cpuCP.timeIncrement();
    
    if (DAC_)
    {
        for (label i=0; i<spaceSize-2; i++)
            completeToSimplifiedIndex_[i] = chemistry.completeToSimplifiedIndex(i);
        for (label i=0; i<NsDAC_; i++)
            simplifiedToCompleteIndex_[i] = chemistry.simplifiedToCompleteIndex(i);
    }
    
    label dim = spaceSize;
    if (DAC_) dim = NsDAC_+2;
    
    LT_ = List<List<scalar> >(dim,List<scalar>(dim,0.0));
    //QT_ = List<List<scalar> >(dim,List<scalar>(dim,0.0)); 	
    
    //SVD decomposition A= U*D*V^T 
    List<List<scalar> > Atilde(A);
    List<List<scalar> > B(dim,List<scalar>(dim,0.0));	
    scalarField diag(dim,0.0);
    svd(Atilde, dim-2, dim, diag, B);
    //replace the value of vector diag by max(diag, 1/2)
    for (register label i=0; i<dim; i++) diag[i]=max(diag[i], 0.5);
    
    //reconstruct Atilde = U*D'*V (ellipsoid in with length d'[i] and principal semi-axes in the direction of the column of )
    for (register label i=0; i<dim-2; i++)
    {
        scalarField AtildeI(dim);
        for (register label n=0; n<dim; n++) AtildeI[n] = Atilde[i][n];
        for (register label j=0; j<dim; j++)
        {	
            Atilde[i][j]=0.0;
            for (register label k=0; k<dim; k++) Atilde[i][j] += AtildeI[k]*diag[k]*B[j][k];
            //added to use qrDecompose on the reduced composition space
            label si=i;
            if(DAC_) si = simplifiedToCompleteIndex(i);
            Atilde[i][j] /= (epsTol*scaleFactor[si]); 
        }
    }
    qrDecompose(dim,Atilde);
    word inertSpecieName(chemistry.thermo().lookup("inertSpecie"));
    forAll(chemistry.Y(),Yi)
    {
        if(chemistry.Y()[Yi].name()==inertSpecieName)
        {
            inertSpecie_=Yi;
        }
    }
}


template<class CompType, class ThermoType>
chemPointISAT<CompType, ThermoType>::chemPointISAT
(
    chemPointISAT<CompType, ThermoType>& p
)
:
//	chemistry_(p.chemistry()),
    chemPointBase(),
    phi_(p.phi()),
    Rphi_(p.Rphi()),
    LT_(p.LT()),
    //QT_(p.QT()),
    A_(p.A()),
    scaleFactor_(p.scaleFactor()),
    node_(p.node()),
    spaceSize_(p.spaceSize()),
    nUsed_(p.nUsed()),
    nGrown_(p.nGrown()),
    //epsTol_(p.epsTol()),
    DAC_(p.DAC()),
    NsDAC_(p.NsDAC()),
    completeToSimplifiedIndex_(p.completeToSimplifiedIndex()),
    simplifiedToCompleteIndex_(p.simplifiedToCompleteIndex()),
    inertSpecie_(p.inertspecie()),
    timeTag_(p.timeTag()),
    lastTimeUsed_(p.lastTimeUsed())/*,
    failedSpeciesFile_(p.failedSpeciesFile()),
    failedSpecies_(failedSpeciesFile_.c_str(), ofstream::app)*/
{
   epsTol_ = p.epsTol();

}    

/*---------------------------------------------------------------------------*\
	To RETRIEVE the mapping from the chemPoint phi, the query point phiq has to 
    be in the EOA of phi. It follows that, dphi=phiq-phi and to test if phiq
    is in the ellipsoid there are two methods. First, compare r=||dphi|| with 
    rmin and rmax. If r < rmin, phiq is in the EOA. If r > rmax, phiq is out of
    the EOA. This operations is O(spaceSize) and is performed first.
    If rmin < r < rmax, then the second method is used:
    	||L^T.dphi|| <= 1 to be in the EOA.
    	
    Note : the use of rmin and rmax is not implemented yet	
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
bool chemPointISAT<CompType, ThermoType>::inEOA(const scalarField& phiq)
{
    lastError_=0.0;
    const List<List<scalar> >& LTvar = LT();
    scalarField dphi=phiq-phi();
    label dim = spaceSize()-2;
    if (DAC_) dim = NsDAC_;
    
    for (register label i=0; i<spaceSize()-2; i++)
    {
        //skip the inertSpecie
        if (i==inertSpecie_)
            continue;
        
        scalar epsTemp=0.0;
        scalarList curEps(spaceSize(),0.0);
      
        //without DAC OR with DAC and on an active species line
        //multiply L by dphi to get the distance in the active species direction
        //else (with DAC and inactive species), just multiply the diagonal element 
        //and dphi
        if (!(DAC_) || (DAC_ && completeToSimplifiedIndex(i)!=-1))
        {
            label si;
            if(DAC_) si=completeToSimplifiedIndex(i);
            else si=i;
            for(register label j=si; j<dim; j++)//LT is upper triangular
            {
                label sj;
                if(DAC_) sj =simplifiedToCompleteIndex(j);
                else	sj=j;

                curEps[sj]=LTvar[si][j]*dphi[sj];
                epsTemp += curEps[sj];
            }
            curEps[spaceSize()-2]=LTvar[si][NsDAC_]*dphi[spaceSize()-2];
            epsTemp += curEps[spaceSize()-2];
            curEps[spaceSize()-1]=LTvar[si][NsDAC_+1]*dphi[spaceSize()-1];
            epsTemp += curEps[spaceSize()-1];
        }
        else
        {
            epsTemp = dphi[i]/(epsTol_*scaleFactor_[i]);
            curEps[i] = epsTemp;
        }

        lastError_ += sqr(epsTemp);
        
        if(fabs(epsTemp) > 1.0)
        {
            if(chemistry_->analyzeTab())
            {
                chemistry_->addToSpeciesNotInEOA(i); //not in the EOA for the ith species direction in the composition space
                scalar maxEps = 1.0;
                label epsIndex = -1;
                forAll(curEps,cei)
                {
                    if(fabs(curEps[cei]) > maxEps)
                    {
                        maxEps = fabs(curEps[cei]);
                        epsIndex = cei;
                    }
                }
                if (epsIndex != -1)
                {
                    chemistry_->addToSpeciesImpact(epsIndex);
                }
            }   
            //should break when optimized but to analyze ISAT we need to go to the end
            break; 
        }
        else if(lastError_ > 1.0)
        {
            //the error can only grow, therefore, once it is above 1.0, we can stop the loop
            break;
        }
    }
    
    //sqrt(eps2) is not required since it is compared to 1	
    if(lastError_ > 1.0)
    {	
        return false;
    }
    else
    {   
    	if(nUsed_ < INT_MAX)
    	{
            nUsed_++;
    	}    

        return true;
    }
}

/*---------------------------------------------------------------------------*\
	If phiq is not in the EOA, then the mapping is computed. But as the EOA
    is a conservative approximation of the region of accuracy surrounding the
    point phi, we could expand it by comparing the computed results with the 
    one obtained by linear interpolation. The error epsGrow is calculated:
    	epsGrow = ||B.(dR - dRl)||,
    with dR = Rphiq - Rphi, dRl = A.dphi and B the diagonal scale factor
    matrix.
    If epsGrow <= epsTol, the EOA is too conservative and a GROW is perforned,
    otherwise, the newly computed mapping is associated to the initial 
    composition and added to the tree.
    
    A bool argument "isInEOA" is also used when a shrink process is required.
    If a query point phiq lies in this EOA but is not in the region of
    accuracy, then the grow algorithm is used to shrink the EOA with phiq on
    its boundary
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
bool chemPointISAT<CompType, ThermoType>::checkSolution(const scalarField& phiq, const scalarField& Rphiq)
{
    scalar eps2 = 0.0;
    scalarField dR = Rphiq - Rphi();
    scalarField dphi = phiq - phi();
    const scalarField& scaleFactorV = scaleFactor();
    const List<List<scalar> >& Avar = A();
    scalar dRl = 0.0;
    label dim = spaceSize()-2;
    if (DAC_) dim = NsDAC_;
    
    for (register label i=0; i<spaceSize()-2; i++)
    {
        dRl = 0.0;
        if (DAC_)
        {
            label si = completeToSimplifiedIndex_[i];
            if (si!=-1)
            {
                for (register label j=0; j<dim; j++)
                {
                    label sj;
                    if(DAC_) sj =simplifiedToCompleteIndex(j);
                    else    sj=j;
                    dRl += Avar[si][j]*dphi[sj];
                }
                dRl += Avar[si][NsDAC_]*dphi[spaceSize()-2];
                dRl += Avar[si][NsDAC_+1]*dphi[spaceSize()-1];
            }
            else
                dRl = dphi[i];
        }
        else
        {
            for (register label j=0; j<spaceSize(); j++)
            {
                dRl += Avar[i][j]*dphi[j];
            }
        }
        eps2 += sqr((dR[i]-dRl)/scaleFactorV[i]);
    }	
    
    eps2 = sqrt(eps2);
    
    
    if(eps2 > epsTol())
    {
    	return false;
    }
    else
    {
        // if the solution is in the ellipsoid of accuracy
        // GROW operation performed
	nGrown_++;
	grow(phiq); //phiq is on the boundary of the EOA
        return true;    
    }
}

/*---------------------------------------------------------------------------*\
	More details about the minimum-volume ellipsoid covering an ellispoid E and
  	a point p are found in [1]. Here is the main steps to obtain the modified 
  	matrix L' describind the new ellipsoid.
  		1) calculate the point p' in the transformed space :
  			p' = L^T.(p-phi)
  		2) compute the rank-one decomposition:
  			G = I + gamma.p'.p'^T,
  		   with gamma = (1/|p'|-1)*1/|p'|^2	
  		3) compute L':
  			L'L'^T = (L.G)(L.G)^T, 
  			L'^T is then obtained by QR decomposition of (L.G)^T = G^T.L^T
	[1] Stephen B. Pope, "Algorithms for ellipsoids", FDA 08-01,
		Cornell University, 2008  			
  			
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void chemPointISAT<CompType, ThermoType>::grow(const scalarField& phiq)
{
    List<List<scalar> >& LTvar = LT();
    List<List<scalar> >& Avar  = A();
    scalarField dphi = phiq - phi();
    label dim = spaceSize();
    label initNsDAC(NsDAC_);
    
    //first step when DAC is used: check if some species should be "activated"
    //this avoid growing the EOA and then fail inEOA because one inactive species
    //is sligthly different from the stored one
    if(DAC_)
    {
        label activeAdded(0);
        for (label i=0; i<spaceSize()-2; i++)
        {
            if(dphi[i]!=0.0 && completeToSimplifiedIndex(i)==-1 && i!=inertSpecie_)//if previously inactive but dphi!=0 and not inert
            {
                NsDAC_++;
                simplifiedToCompleteIndex_.setSize(NsDAC_,i); //add the new active species
                completeToSimplifiedIndex_[i]=NsDAC_-1;
                activeAdded++;
            }
        }
        
        //update LT and A :
        //-add new column and line for the new active species
        //-transfer last two lines of the previous matrix (p and T) to the end (change the diagonal position)
        //-set all element of the new lines and columns to zero except diagonal (=1/(epsTol*scaleFactor))
        if(NsDAC_ > initNsDAC)
        {
            LTvar.setSize(NsDAC_+2, List<scalar>(initNsDAC+2,0.0));
            Avar.setSize(NsDAC_+2, List<scalar>(initNsDAC+2,0.0));
            for (label i=0; i<NsDAC_+2; i++)
            {
                LTvar[i].setSize(NsDAC_+2, 0.0);
                Avar[i].setSize(NsDAC_+2, 0.0);
            }
            
            for (label i=0; i<NsDAC_-activeAdded; i++)
            {
                //star with last column, otherwise problems when activeAdded=1
                for (label j=1; j>=0; j--)
                {
                    LTvar[i][NsDAC_+j]=LTvar[i][NsDAC_+j-activeAdded];
                    Avar[i][NsDAC_+j]=Avar[i][NsDAC_+j-activeAdded];
                    Avar[NsDAC_+j][i]=Avar[NsDAC_+j-activeAdded][i];
                    LTvar[i][NsDAC_+j-activeAdded]=0.0;
                    Avar[i][NsDAC_+j-activeAdded]=0.0;
                    Avar[NsDAC_+j-activeAdded][i]=0.0;
                }
            }
            for (label i=NsDAC_+1; i>=NsDAC_; i--)
            {
                for (label j=NsDAC_+1; j>=NsDAC_; j--)
                {
                    LTvar[i][j]=LTvar[i-activeAdded][j-activeAdded];
                    Avar[i][j]=Avar[i-activeAdded][j-activeAdded];
                    LTvar[i-activeAdded][j-activeAdded]=0.0;
                    Avar[i-activeAdded][j-activeAdded]=0.0;
                }
            }
            for (label i=NsDAC_-activeAdded; i<NsDAC_;i++)
            {	
                LTvar[i][i]=1.0/(epsTol_*scaleFactor_[simplifiedToCompleteIndex(i)]);
                Avar[i][i]=1.0;

            }
        }//end if(NsDAC_>initNsDAC)
        dim = NsDAC_+2;
    }//end if(DAC_)
    //beginning of grow algorithm
    scalarField phiTilde(dim, 0.0);
    scalar	normPhiTilde = 0.0;	
    //p' = L^T.(p-phi)
    for (register label i=0; i<dim; i++)
    {
        for(label j=i; j<dim-2; j++)//LT is upper triangular
        {
            label sj = j;
            if(DAC_) sj=simplifiedToCompleteIndex(j);
            phiTilde[i] += LTvar[i][j]*dphi[sj];
        }
        phiTilde[i] += LTvar[i][dim-2]*dphi[spaceSize()-2];
        phiTilde[i] += LTvar[i][dim-1]*dphi[spaceSize()-1];
        normPhiTilde += sqr(phiTilde[i]);
    }
    scalar invSqrNormPhiTilde = 1.0/normPhiTilde;
    normPhiTilde = sqrt(normPhiTilde);
    //gamma = (1/|p'| - 1)/|p'|^2
    
    scalar gamma = (1/normPhiTilde - 1)*invSqrNormPhiTilde;
    
    
    //NEW VERSION : with updating QR decomposition
    scalarField u(gamma*phiTilde);
    scalarField v(dim,0.0);
    for (register label i=0; i<dim; i++)
    {
        for (register label j=0; j<=i;j++)
            v[i] += phiTilde[j]*LTvar[j][i];
    }
    
    qrUpdate(dim, u, v);
}



template<class CompType, class ThermoType>
void chemPointISAT<CompType, ThermoType>::setFree()
{
    node_ = NULL;
    nUsed_ = 0;
       
}

template<class CompType, class ThermoType>
void chemPointISAT<CompType, ThermoType>::clearData()
{
    nUsed_ = 0;
    phi_.clear();        
    Rphi_.clear();
    LT_.clear();        
    A_.clear();
    epsTol_ = 0;
}



/*---------------------------------------------------------------------------*\
	QR decomposition of the matrix A 
	Input : nCols cols number
			A the matrix to decompose A = Q.R
			R empty matrix in which the upper triangular matrix is stored
	Output: void
	Note : in this implementation we are only interested in the R matrix
		   which will be the transpose of L
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void chemPointISAT<CompType, ThermoType>::qrDecompose
(
    const label nCols,
    	  List<List<scalar> >& Q
)
{
    scalarField c(nCols);
    scalarField d(nCols);
    scalar scale, sigma, sum;
    
    for (label k=0; k<nCols-1; k++)
    {
    	scale = 0.0;
    	for (label i=k; i<nCols; i++) scale=max(scale, fabs(Q[i][k]));
    	if (scale == 0.0)
    	{
            c[k]=d[k]=0.0;
    	}
    	else
    	{
            for (label i=k; i<nCols; i++) Q[i][k] /= scale;
            sum = 0.0;
            for (label i=k; i<nCols; i++) sum += sqr(Q[i][k]);
            sigma = sign(Q[k][k])*sqrt(sum);
            Q[k][k] += sigma;
            c[k]=sigma*Q[k][k];
            d[k]=-scale*sigma;
            for (label j=k+1; j<nCols; j++)
            {
                sum=0.0;
                for ( label i=k; i<nCols; i++) sum += Q[i][k]*Q[i][j];
                scalar tau = sum/c[k];
                for ( label i=k; i<nCols; i++) Q[i][j] -= tau*Q[i][k];
            }
    	}
    }
    d[nCols-1] = Q[nCols-1][nCols-1];
    
    
    //form R
    List<List<scalar> >& R(LT());
    for (label i=0; i<nCols; i++)
    {
        R[i][i] = d[i];
    	for ( label j=0; j<i; j++) 
            R[i][j]=0.0;
        for (label j=i+1; j<nCols; j++) 
            R[i][j]=Q[i][j];
    }    
}//end qrDecompose

/*---------------------------------------------------------------------------*\
 QR update of the matrix A
 Input : nCols cols number
 A the matrix to decompose A = Q.R
 R empty matrix in which the upper triangular matrix is stored
 Output: void
 Note : in this implementation we are only interested in the R matrix
 which will be the transpose of L
 \*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void chemPointISAT<CompType, ThermoType>::qrUpdate
(
 const label n,
 const scalarField &u,
 const scalarField &v
)
{
    label k,i;
    scalarField w(u);
    for (k=n-1;k>=0;k--) 
        if (w[k] != 0.0) break; 
    if (k < 0) k=0; 
    for (i=k-1;i>=0;i--) { 
        rotate(i,w[i],-w[i+1], n); 
        if (w[i] == 0.0) 
            w[i]=fabs(w[i+1]);
        else if (fabs(w[i]) > fabs(w[i+1])) 
            w[i]=fabs(w[i])*sqrt(1.0+sqr(w[i+1]/w[i])); 
        else w[i]=fabs(w[i+1])*sqrt(1.0+sqr(w[i]/w[i+1])); 
    } 
    List<List<scalar> >&R(LT());
    for (i=0;i<n;i++) R[0][i] += w[0]*v[i]; 
    for (i=0;i<k;i++) 
        rotate(i,R[i][i],-R[i+1][i], n);
}		
    
//rotate function used by qrUpdate	
template<class CompType, class ThermoType>
void chemPointISAT<CompType, ThermoType>::rotate(const label i, const scalar a, const scalar b, label n)
{
    label j;
    scalar c, fact, s, w, y;
    if (a == 0.0)
    {
        c=0.0;
        s=(b >= 0.0 ? 1.0 : -1.0);
    }
    else if (fabs(a) > fabs(b))
    {
        fact = b/a;
        c=sign(a)/sqrt(1.0+(fact*fact));
        s=fact*c;
    }
    else
    {
        fact=a/b;
        s=sign(b)/sqrt(1.0+(fact*fact));
        c=fact*s;
    }
    List<List<scalar> >& R(LT());
    for (j=i;j<n;j++)
    {
        y=R[i][j];
        w=R[i+1][j];
        R[i][j]=c*y-s*w;
        R[i+1][j]=s*y+c*w;
    }
}	
	
	
/*---------------------------------------------------------------------------*\
	Singular Value Decomposition (SVD) for a square matrix
	needed to compute the length of the hyperellipsoid semi-axes
	SVD decompose a matrix A into:
		A = U * D * V^T , 
	with the singular value in the diagonal matrix D and U and V orthogonal 
	Input :  A (scalarMatrix) the matrix m x n to apply the decomposition on
		     m (label) the number of rows of the matrix A
		     n (label) the number of columns of the matrix A
	Output:  U (scalarMatrix) replace A
		     V (scalarMatrix) not the transpose V^T
		     d (scalarField) the diagonal element of matrix D
	The algorithm is based on "Numerical recipes in C", second edition (1992)	
	chapter 2, pp 67-70	  
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void chemPointISAT<CompType, ThermoType>::svd(List<List<scalar> >& A, label m, label n, scalarField& d, List<List<scalar> >& V)
{
    //UPDATED VERSION NR3
    bool flag;
    label i,its,j,jj,k,l,nm;
    scalar anorm,c,f,g,h,s,scale,x,y,z;
    scalarField rv1(n);
    scalar eps = std::numeric_limits<scalar>::epsilon();
    g = scale = anorm = 0.0;
    
    
    //Householder reduction to bidiagonal form 
    for( i = 0; i<n; i++)
    {
        l=i+2; //change from i+1 to i+2
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m)
        {
            for (k=i;k<m;k++)	scale += fabs(A[k][i]);
            if (scale != 0.0)
            {
                for ( k=i;k<m;k++)
                {
                    A[k][i] /= scale;
                    s += A[k][i]*A[k][i];
                }
                f = A[i][i];
                g = -sign(f)*sqrt(s);
                h = f*g-s;
                A[i][i]=f-g;
                for (j=l-1;j<n;j++)
                {
                    for (s=0.0,k=i;k<m;k++) s += A[k][i]*A[k][j];
                    f = s/h;
                    for (k=i; k<m;k++) A[k][j] += f*A[k][i];
                }
                for (k=i; k<m;k++) A[k][i] *= scale;
            }
        }
        d[i] = scale * g;
        g=s=scale=0.0;
        
        if (i+1 <= m && i+1 != n)
        {
            for (k=l-1; k<n; k++) scale += fabs(A[i][k]);
            if (scale != 0.0) 
            {
                for (k=l-1; k<n; k++)
                {
                    A[i][k] /= scale;
                    s += A[i][k]*A[i][k];
                }
                f = A[i][l-1];
                g = -sign(f)*sqrt(s);
                h = f*g-s;
                A[i][l-1] = f-g;
                for (k=l-1; k<n; k++) rv1[k]=A[i][k]/h;
                for (j=l-1; j<m; j++) 
                {
                    for (s=0.0,k=l-1; k<n; k++) s += A[j][k]*A[i][k];
                    for (k=l-1; k<n; k++) A[j][k] += s*rv1[k];
                }
                for (k=l-1; k<n; k++) A[i][k] *= scale;
            }
        }
        anorm = max(anorm, (fabs(d[i])+fabs(rv1[i])));	
    }//end Householder reduction
    
    //Accumulation of right-hand transformations
    for (i=n-1; i>=0; i--)
    {
        if (i < n-1)
        {
            if (g != 0.0)
            {
                for (j=l; j<n; j++) V[j][i] = (A[i][j]/A[i][l])/g; 
                for (j=l; j<n; j++)
                {
                    for (s=0.0,k=l; k<n; k++) s += A[i][k]*V[k][j];
                    for (k=l; k<n; k++) V[k][j] += s*V[k][i];			
                }
            }
            for (j=l; j<n; j++) V[i][j]=V[j][i]=0.0;
        }
        V[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    //Accumulation of left-hand transformations
    for (i = min(m,n)-1; i>=0; i--)
    {
        l=i+1;
        g=d[i];
        for (j=l; j<n; j++) A[i][j] = 0.0;
        if (g != 0.0) 
        {
            g = 1.0/g;
            for (j=l; j<n; j++)
            {
                for (s=0.0, k=l; k<m; k++) s+= A[k][i]*A[k][j];
                f = (s/A[i][i])*g;
                for (k=i; k<m; k++) A[k][j] += f*A[k][i];
            }	
            for (j=i; j<m; j++) A[j][i] *= g;
        }
        else for (j=i; j<m; j++) A[j][i]=0.0;
        ++A[i][i];
    }
    
    //Diagonalization of the bidiagonal form :
    //Loop over singular values, and over allowed iteration
    for (k=n-1; k>=0; k--)
    {
        for (its=0; its<30; its++)
        {
            flag=true;
            // Test for splitting (rv1[1] always zero)
            for (l=k; l>=0; l--)
            {
                nm = l-1;
                if (l == 0 || fabs(rv1[l]) <= eps*anorm)
                {
                    flag = false;
                    break;
                }
                if (fabs(d[nm]) <= eps*anorm) break;
            }
            //Cancellation of rv1[l], if l>1
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i=l; i<k+1; i++)
                {
                    f = s*rv1[i];
                    rv1[i] = c*rv1[i];
                    if (fabs(f) <= eps*anorm) break;
                    g = d[i];
                    h = pythag(f,g);
                    d[i] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;
                    for (j=0; j<m; j++)
                    {
                        y = A[j][nm];
                        z = A[j][i];
                        A[j][nm] = y*c + z*s;
                        A[j][i] = z*c - y*s;
                    }
                }
            }
            
            z = d[k];
            if (l == k) //Convergence
            {
                if (z < 0.0) //Singular value is made nonnegative
                {
                    d[k] = -z;
                    for (j=0; j<n; j++) V[j][k] = -V[j][k];
                }
                break;
            }
            if (its == 29)
            {
                //FatalErrorIn("void Foam::chemistryOnlineLibrary::svd(scalarMatrix& A, label n, scalarField& d, scalarMatrix& V)")
                //<< "No convergence in 30 iterations" << abort(FatalError);
                Info << "No convergence in 30 iterations" << endl;
            }
            x = d[l];
            nm = k-1;
            y = d[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f,1.0);
            f = ((x-z)*(x+z)+h*((y/(f+sign(f)*g))-h))/x;
            c=s=1.0; 
            //Next QR transformation
            for (j=l; j<=nm; j++)
            {
                i = j+1;
                g = rv1[i];
                y = d[i];
                h = s*g;
                g = c*g;
                z = pythag(f,h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y*s;
                y *= c;
                for (jj=0; jj<n; jj++)
                {
                    x = V[jj][j];
                    z = V[jj][i];
                    V[jj][j] = x*c + z*s;
                    V[jj][i] = z*c - x*s;
                }
                z = pythag(f,h);
                d[j] = z;
                if (z)
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g + s*y;
                x = c*y - s*g;
                for (jj=0; jj<m; jj++)
                {
                    y = A[jj][j];
                    z = A[jj][i];
                    A[jj][j] = y*c + z*s;
                    A[jj][i] = z*c - y*s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            d[k] = x;
        }
    }
    
}//end svd

// pythag function used in svd
// compute (a^2+b^2)^1/2 without descrutive underflow or overflow
template<class CompType, class ThermoType>
scalar chemPointISAT<CompType, ThermoType>::pythag(scalar a, scalar b)
{
    scalar absa, absb;
    //absa = mag(a);
    //absb = mag(b);
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb) return absa*sqrt(1.0+sqr(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+sqr(absa/absb)));
}


}
