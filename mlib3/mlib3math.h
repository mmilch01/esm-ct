/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#ifndef _MATH_H_INCLUDED_
#define _MATH_H_INCLUDED_

//#ifdef WIN32
#include <limits>
#define _inf	numeric_limits<double>::infinity()
#define _nan	numeric_limits<double>::quiet_NaN()

#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>

#include "newm/include.h"
#include "newm/newmat.h"
#include "newm/newmatap.h"
#include "newm/newmatio.h"

#include "mlib3common.h"

#define sign(a) (((a)>0)?1:(((a)==0)?0:-1))
class MMath;

typedef struct Range1
{
	Range1(){st=-1; en=-1;};
	Range1(int s, int e){Set(s,e);};
	void operator=(Range1& r){st=r.st; en=r.en;};
	void Set(int s, int e){st=(Real)s;en=(Real)e;};
	bool Valid(){return (st>-1 && en>-1);};
	Real st,en;
} Range1;

typedef struct Index3
{
	Index3(){d[0]=d[1]=d[2]=-1;};
	Index3(int d1,int d2=-1,int d3=-1) {Set(d1,d2,d3);};
	void Set(int d1,int d2=-1,int d3=-1){d[0]=d1;d[1]=d2;d[2]=d3;};	
	int d[3];
} Index3;

typedef struct Range3
{
	void Set(int x0, int x1, int y0, int y1, int z0, int z1)
	{
		r[0].st = x0; r[0].en = x1;
		r[1].st = y0; r[1].en = y1;
		r[2].st = z0; r[2].en = z1;
	}
	Range3(int x0,int x1,int y0, int y1, int z0, int z1)
	{
		Set(x0, x1, y0, y1, z0, z1);
	}
	Range3(){for (int i=0; i<3; i++) r[i].st=r[i].en=-1;};
	Range3(Range1 r1, Range1 r2=Range1(-1,-1), Range1 r3=Range1(-1,-1)) {r[0]=r1;r[1]=r2;r[2]=r3;}
	Range3(Range1 r1, Range1 r2, Range1 r3, int* permute)
	{
		r[permute[0]]=r1; r[permute[1]]=r2; r[permute[3]]=r3;
	}

	bool IsPtInside(ColumnVector& pt){
		return (pt(1)>=r[0].st && pt(1)<=r[0].en && pt(2)>=r[1].st && pt(2)<=r[1].en && pt(3)>=r[2].st && pt(3)<=r[2].en );
	}
	bool Contains2D(Range3& r1)
	{
		return (IsPtInside2D(r1[0].st,r1[1].st) && IsPtInside2D(r1[0].en,r1[1].en));
	}
	bool IsPtInside2D(Real x, Real y){
		return (x >= r[0].st && x <= r[0].en && y >= r[1].st && y <= r[1].en);
	}
	bool IsPtInside(Real x, Real y, Real z){
		return (x>=r[0].st && x<=r[0].en && y>=r[1].st && y<=r[1].en && z>=r[2].st && z<=r[2].en );
	}

	Range3(Range1* rn){for (int i=0; i<3; i++) r[i]=rn[i];};
	void GetBoundingBox(ColumnVector& v1, ColumnVector& v2)
	{
		r[0].st=min(v1(1),v2(1));r[0].en=max(v1(1),v2(1));
		r[1].st=min(v1(2),v2(2));r[1].en=max(v1(2),v2(2));
		r[2].st=min(v1(3),v2(3));r[2].en=max(v1(3),v2(3));
	}
	Real SX(){return r[0].en-r[0].st;};
	Real SY(){return r[1].en-r[1].st;};
	Real SZ(){return r[2].en-r[2].st;};

	void GetSize(Index3& ind){for (int i=0; i<3;i++) {ind.d[i]=(int)(r[i].en-r[i].st+.5); if(ind.d[i]>=0) ind.d[i]++;}};
	void Expand(int dx, int dy, int dz)
	{
		r[0].st-=dx; r[0].en+=dx;
		r[1].st-=dy; r[1].en+=dy;
		r[2].st-=dz; r[2].en+=dz;
	}
	Real Volume()
	{
		Real v=1;
		for(int i=0; i<3; i++) 
			v*=r[i].en-r[i].st;
		return fabs(v);
	}
#define _ima(x,y) ((x)>(y))?(x):(y)
#define _imi(x,y) ((x)<(y))?(x):(y)
	void Intersect (int minX, int maxX, int minY, int maxY, int minZ, int maxZ)
	{
		r[0].st=_ima(r[0].st,minX);  r[0].en=_imi(r[0].en,maxX);
		r[1].st=_ima(r[1].st,minY);  r[1].en=_imi(r[1].en,maxY);
		r[2].st=_ima(r[2].st,minZ);  r[2].en=_imi(r[2].en,maxZ);
	}
	void Intersect(Range3& r1)
	{
		Intersect(r1.r[0].st,r1.r[0].en,r1.r[1].st,r1.r[1].en,r1.r[2].st,r1.r[2].en);
	}
	Range3& operator =(Range3& rn){for (int i=0; i<3; i++) r[i]=rn.r[i]; return (*this);};
	void operator *= (Real f){ 
		for (int i = 0; i<3; i++) {
			r[i].st*=f; r[i].en*=f;
		} 
	};
	void operator += (Real f) {
		for (int i = 0; i<3; i++) {
			r[i].st += f; r[i].en += f;
		}
	};

	Range1 r[3];
	Range1& operator [](int i){ return (Range1&)(r[i]); }
} Range3;

class Volume;
class MMath
{
public:
	enum TRTYPE {ROTATE,TRANSLATE,SCALE};
	enum AXIS {X,Y,Z};
	static const long_Real UndefDouble, PI, RAD;
	static const int PLANE_XY,PLANE_YZ,PLANE_XZ,COMPARE_NMI,COMPARE_CORR,COMPARE_ABSDIFF;
	static void IdentityTransform(Matrix& T) { T.ReSize(4, 4); T << 1 << 0 << 0 << 0 << 0 << 1 << 0 << 0 << 0 << 0 << 1 << 0 << 0 << 0 << 0 << 1; };
	static bool	Gauss(Volume& a, Volume& b, int sz) {
		Real *krn = new Real[(sz * 2 + 1)*(sz * 2 + 1)];
		MMath::GaussCoeff(krn, sz * 2);
		MMath::ConvolveSep3D_sm(a, b, krn, sz);
		delete[] krn;
		return true;
	};
	static void ConvolveSep3D_sm(Volume& v0, Volume& v1, Real *krn, int krn_sz, Volume* mask = NULL);
	static void GaussCoeff(Real* buf, int n);
	static void VectorToE1(ColumnVector& A, ColumnVector& B, Matrix& T);
	static void ConvolveBinom1D_3(ColumnVector& v, ColumnVector* mask);
	static void ConvolveVolume1D(Volume& v0, Volume& v1, int dim, Real* krn, int krn_sz, Volume* mask = NULL);
	static void	Convolve1D(ColumnVector& v, ColumnVector& res, Real* krn, int krn_sz);
	static void TranslateMatrix(ColumnVector& v, Matrix& T);
};

//Press et al. Numerical Recipies, 3rd Edition 2007, p. 506
struct Amoeba
{
	//Multidimensional minimization by the downhill simplex method of Nelder and Mead.
	const Real ftol; //невязка
	bool bVerbose;
	int nfunc; //The number of function evaluations.
	int mpts;
	int ndim;
	int maxIt; //макс.число итераций
	Real fmin; //Function value at the minimum.
	ColumnVector y; //Function values at the vertices of the simplex.
	Matrix p; //Current simplex.
	Amoeba(const Real ftoll, const int iter, const bool ver=true) : ftol(ftoll),maxIt(iter),bVerbose(ver) {};
	//The constructor argument ftoll is the fractional convergence tolerance to be achieved in the function value (n.b.!).
	template <class T>
	ColumnVector minimize(const ColumnVector &point, const Real del, T &func)
	//Multidimensional minimization of the function or functor func(x), where x[0..ndim-1]
	//is a vector in ndim dimensions, by the downhill simplex method of Nelder and Mead.
	//The initial simplex is specified as in equation (10.5.1) by a point[0..ndim-1] and a
	//constant displacement del along each coordinate direction. Returned is the location of the
	//minimum.
	{
		ColumnVector dels(point.size()); dels=del;
		return minimize(point,dels,func);
	};
	template <class T>
	ColumnVector minimize(const ColumnVector &point, const ColumnVector &dels, T &func)
	//Alternative interface that takes different displacements dels[0..ndim-1] in different directions
	//for the initial simplex.
	{
		int ndim=point.size();
		Matrix pp(ndim+1,ndim);
		for (int i=1;i<=ndim+1;i++) {
			for (int j=1;j<=ndim;j++)
			pp(i,j)=point(j);
			if (i !=1 ) pp(i,i-1) += dels(i-1);
			}

		return minimize(pp,func);
	};
	template<class T>
	inline void SWAP(T &a, T &b) {T dum=a; a=b; b=dum;};

	template <class T>
	ColumnVector minimize(const Matrix &pp, T &func)
	//Most general interface: initial simplex specified by the matrix pp[0..ndim][0..ndim-1].
	//Its ndim+1 rows are ndim-dimensional vectors that are the vertices of the starting simplex.
	{
		const int NMAX=5000; //Maximum allowed number of function evaluations
		const Real TINY=1.0e-10;
		int ihi,ilo,inhi;
		mpts=pp.nrows();
		ndim=pp.ncols();
		ColumnVector psum(ndim),pmin(ndim),x(ndim);
		p=pp;
		y.resize(mpts);
		int nFlat=0;
		for (int i=1;i<=mpts;i++) 
		{
			for (int j=1;j<=ndim;j++)
				x(j)=p(i,j);
			y(i)=func(x);
			if (y(i)==y(1)) nFlat++;
		}
		//return if neighborhood is flat.
		if (nFlat==mpts) {
//			cout << "flat neighborhood, returning zero displacement" << endl;
			for (int i=1; i<=ndim; i++)
				pmin(i)=0;
			return pmin;
		}
		nfunc=0;
		get_psum(p,psum);
		int iter=0; 

		//Number of iterations not to exceed maxIt.
		for (;;) {
			ilo=1;
			//First we must determine which point is the highest (worst), next-highest, and
			//lowest (best), by looping over the points in the simplex.
			ihi = y(1)>y(2) ? (inhi=2,1) : (inhi=1,2);
			for (int i=1;i<=mpts;i++) {
				if (y(i) <= y(ilo)) ilo=i;
				if (y(i) > y(ihi)) {
					inhi=ihi;
					ihi=i;
				} 
				else if (y(i) > y(inhi) && i != ihi) inhi=i;
			}
			Real rtol=2.0*fabs(y(ihi)-y(ilo))/(fabs(y(ihi))+fabs(y(ilo))+TINY);
			//Compute the fractional range from highest to lowest and return if satisfactory.
			if ( rtol < ftol || iter++ >= maxIt ) { //If returning, put best point and value in slot 0.
				SWAP(y(1),y(ilo));
				for (int i=1;i<=ndim;i++) {
					SWAP(p(1,i),p(ilo,i));
					pmin(i)=p(1,i);
				}
				fmin=y(1);
				if (bVerbose){
#ifdef _WIN32
					cout << fixed << setprecision(4) << "D: " << pmin.t();
#endif
					cout << fixed << "eta: " << y(1) << " tol: " << setprecision(4) << rtol << " it: " << iter << " ev: " << nfunc << endl;
				}
				return pmin;
			}
			if (nfunc >= NMAX) throw("NMAX exceeded");
			nfunc += 2;
			//Begin a new iteration. First extrapolate by a factor 1 through the face of the
			//simplex across from the high point, i.e., reflect the simplex from the high point.
			Real ytry=amotry(p,y,psum,ihi,-1.0,func);
			if (ytry <= y(ilo))
				//Gives a result better than the best point, so try an additional extrapolation	by a factor 2.
				ytry=amotry(p,y,psum,ihi,2.0,func);
			else if (ytry >= y(inhi)) 
			{
				//The reflected point is worse than the second-highest, so look for an intermediate
				//lower point, i.e., do a one-dimensional contraction.
				Real ysave=y(ihi);
				ytry=amotry(p,y,psum,ihi,0.5,func);
				if (ytry >= ysave) { //Can’t seem to get rid of that high point.
					//Better contract around the lowest (best) point.
					for (int i=1;i<=mpts;i++) {
						if (i != ilo) {
							for (int j=1;j<=ndim;j++)
								p(i,j)=psum(j)=0.5*(p(i,j)+p(ilo,j));
							y(i)=func(psum);
						}
					}
					nfunc += ndim; //Keep track of function evaluations.
					get_psum(p,psum); //Recompute psum.
				}
			} else --nfunc; //Correct the evaluation count.
		} //Go back for the test of doneness and the next iteration.

	};
	inline void get_psum(const Matrix &p, ColumnVector &psum)
	//Utility function.
	{
		for (int j=1;j<=ndim;j++) {
			Real sum=0.0;
			for (int i=1;i<=mpts;i++)
				sum += p(i,j);
			psum(j)=sum;
		}
	};

	template <class T>
	Real amotry(Matrix &p, ColumnVector &y, ColumnVector &psum, const int ihi, const Real fac, T &func)
	//Helper function: Extrapolates by a factor fac through the face of the simplex across from
	//the high point, tries it, and replaces the high point if the new point is better.
	{
		ColumnVector ptry(ndim);
		Real fac1=(1.0-fac)/ndim;
		Real fac2=fac1-fac;
		for (int j=1;j<=ndim;j++)
			ptry(j)=psum(j)*fac1-p(ihi,j)*fac2;
		Real ytry=func(ptry); //Evaluate the function at the trial point.
		if (ytry < y(ihi)) { //If it’s better than the highest, then replace the highest.
			y(ihi)=ytry;
			for (int j=1;j<=ndim;j++) {
				psum(j) += ptry(j)-p(ihi,j);
				p(ihi,j)=ptry(j);
			}
		}
		return ytry;
	};
};

//Press et al. Numerical Recipies, 3rd Edition 2007, p. 815
struct Fitlin {
	/*Object for general linear least-squares fitting by solving the normal equations, also including
	the ability to hold specified parameters at fixed, specified values. Call constructor to bind data
	vectors and fitting functions. Then call any combination of hold, free, and fit as often as
	desired. fit sets the output quantities a, covar, and chisq.*/
	int ndat, ma;
	ColumnVector &x,&y,sig;
	ColumnVector (*funcs)(const Real);
	ColumnVector ia;
	ColumnVector a; //Output values. a is the vector of fitted coefficients, covar is its covariance matrix, and chisq is the value of X2 for the fit.
	Matrix covar;
	Real chisq;

	/*Constructor. Binds references to the data arrays xx, yy, and ssig, and to a user-supplied
	function funks(x) that returns a ColumnVector containing ma basis functions evaluated at x = X.
	Initializes all parameters as free (not held).*/
	Fitlin(ColumnVector &xx, ColumnVector &yy, ColumnVector funks(const Real)) : ndat(xx.size()), x(xx), y(yy), funcs(funks) 
	{
		sig.ReSize(xx.size()); sig=1;
		ma = funcs(x(1)).size();
		a.resize(ma);
		covar.resize(ma,ma);
		ia.resize(ma);
		for (int i=1;i<=ma;i++) ia(i) = true;		
	}
	Fitlin(ColumnVector &xx, ColumnVector &yy, ColumnVector &ssig, ColumnVector funks(const Real)) : ndat(xx.size()), x(xx), y(yy), sig(ssig), funcs(funks) {
		ma = funcs(x(1)).size();
		a.resize(ma);
		covar.resize(ma,ma);
		ia.resize(ma);
		for (int i=1;i<=ma;i++) ia(i) = true;
	}
	void hold(const int i, const Real val) {ia(i)=false; a(i)=val;}
	void free(const int i) {ia(i)=true;}
	template<class T>
	inline void SWAP(T &a, T &b) {T dum=a; a=b; b=dum;};
	template<class T>
	inline T SQR(T x){ return (x)*(x); };
	/*
	Optional functions for holding a parameter, identified by a value i in the range 0; : : : ; ma-1,
	fixed at the value val, or for freeing a parameter that was previously held fixed. hold and
	free may be called for any number of parameters before calling fit to calculate best-fit
	values for the remaining (not held) parameters, and the process may be repeated multiple
	times. Alternatively, you can set the boolean vector ia directly, before calling fit.
	*/
	void fit() {
		/*
		Solve the normal equations for chi2 minimization to fit for some or all of the coefficients
		a[0..ma-1] of a function that depends linearly on a, y =sum_i(a_i * funks_i(x)). Set answer
		values for a[0..ma-1], chi2 = chisq, and the covariance matrix covar[0..ma-1][0..ma-1].
		(Parameters held fixed by calls to hold will return zero covariances.)
		*/
		int i,j,k,l,m,mfit=0;
		Real ym,wt,sum,sig2i;
		ColumnVector afunc(ma);
		for (j=1;j<=ma;j++) if (ia(j)) mfit++;
		if (mfit == 0) throw("lfit: no parameters to be fitted");
		Matrix temp(mfit,mfit),beta(mfit,1);temp=0;beta=0;
		for (i=1;i<=ndat;i++) { //Loop over data to accumulate coefficients of
			afunc = funcs(x(i));   //the normal equations.
			ym=y(i);
			if (mfit < ma) { //Subtract off dependences on known pieces
				for (j=1;j<=ma;j++) //of the fitting function.
					if (!ia(j)) ym -= a(j)*afunc(j);
			}
			sig2i=1.0/SQR(sig(i));
			for (j=1,l=1;l<=ma;l++) { //Set up matrix and r.h.s. for matrix inversion.
				if (ia(l)) {
					wt=afunc(l)*sig2i;
					for (k=1,m=1;m<=l;m++)
						if (ia(m)) temp(j,k++) += wt*afunc(m);
					beta(j++,1) += ym*wt;
				}
			}
		}
		for (j=2;j<=mfit;j++) for (k=1;k<j;k++) temp(k,j)=temp(j,k);

		ColumnVector b=beta;
		beta=temp.i()*b; //Matrix solution.

		for (j=1,l=1;l<=ma;l++) if (ia(l)) a(l)=beta(j++,1);

		//Spread the solution to appropriate positions in a, and evaluate chi2 of the fit.
		chisq=0.0;
		for (i=1;i<=ndat;i++) {
			afunc = funcs(x(i));
			sum=0.0;
			for (j=1;j<=ma;j++) sum += a(j)*afunc(j);
			chisq += SQR((y(i)-sum)/sig(i));
		}
		for (j=1;j<=mfit;j++) for (k=1;k<=mfit;k++) covar(j,k)=temp(j,k);
		for (i=mfit+1;i<=ma;i++) //Rearrange covariance matrix into the correct order.
			for (j=1;j<=i+1;j++) covar(i,j)=covar(j,i)=0.0;
		k=mfit;
		for (j=ma;j>0;j--) {
			if (ia(j)) {
				for (i=1;i<=ma;i++) SWAP(covar(i,k),covar(i,j));
				for (i=1;i<=ma;i++) SWAP(covar(k,i),covar(j,i));
				k--;
			}
		}
	}
};
#endif
