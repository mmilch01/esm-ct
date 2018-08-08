/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#ifndef _VOLUME_H_INCLUDED_
#define _VOLUME_H_INCLUDED_
#include <iostream>
#include <memory.h>
#include <cmath>

#include "mlib3common.h"
#include "mlib3analyze.h"
#include "mlib3math.h"
#include "newm/newmatio.h"

#ifdef _4DFP
extern "C"
{
	#include <Getifh.h>
}
#endif


class Volume
{
public:
	enum COORD_FORMAT {FSL, FDFP};
	enum INTERP_METHOD {TRILINEAR,NNEIB};
	enum VOXEL_CENSORING_METHOD {NONE, NEAREST};

	static INTERP_METHOD s_interp_method;
	static VOXEL_CENSORING_METHOD s_vox_cens_method;

	int m_nElements, m_dims[3];
	int m_offset[3];
	COORD_FORMAT m_coord_format;
	Real m_voxel_dims[3];
#ifdef _4DFP
	IFH* m_pIFH;
#endif
	Volume(int* dims){InitVars();InitMemory(dims);};
	Volume(int dim1,int dim2,int dim3){InitVars();InitMemory(dim1,dim2,dim3);};
	Volume(){ InitVars();};
	~Volume(){
		if(m_pBuf) delete[] m_pBuf;
		if (m_pIFH)	delete m_pIFH;
	};
	void		InitVars();
	void		InitMemory(Analyze& an, void* pBuf=NULL,int nSlice=-1);
	void		InitMemory(int dim1,int dim2,int dim3);
	void		InitMemory(Index3& ind){InitMemory(ind.d[0],ind.d[1],ind.d[2]);};
	void		InitMemory(int* dim){InitMemory(dim[0],dim[1],dim[2]);};
	void		SetVoxelDims(Real* dims){memcpy(m_voxel_dims,dims,sizeof(Real)*3);};
	void		SetVoxelDims(Real dx,Real dy,Real dz)
				{
					m_voxel_dims[0]=dx;m_voxel_dims[1]=dy;m_voxel_dims[2]=dz;
					m_bVoxDimInited=true;
				};
	inline void ElCurSet(int ind){m_ElCur=ind;};
	inline void	SetX(Real x, int i){(*this)(i,m_CurPt.d[1],m_CurPt.d[2])=x;};
	inline void	SetY(Real y, int i){(*this)(m_CurPt.d[0],i,m_CurPt.d[2])=y;};
	inline void	SetZ(Real z, int i){(*this)(m_CurPt.d[0],m_CurPt.d[1],i)=z;};
	void		operator&= (Volume &v){for (int i=0; i<m_nElements; i++) if(v.m_pBuf[i]==0) m_pBuf[i]=0;};
	void		operator&= (Real d) { for (int i = 0; i<m_nElements; i++) if (m_pBuf[i] != d) m_pBuf[i] = 0; };
	void		operator/= (Real d){for (int i=0; i<m_nElements; i++) m_pBuf[i]/=d;};
	void		operator/= (Volume& v) {for (int i=0; i<m_nElements; i++) m_pBuf[i]/=v.m_pBuf[i];};
	void		operator*= (Real d){for (int i=0; i<m_nElements; i++) m_pBuf[i]*=d;};
	void		operator*= (Volume& v)
	{
		for (int i=0; i<m_nElements; i++) 
			m_pBuf[i]*=v.m_pBuf[i];
	};
	inline Real&
				operator () (ColumnVector& c){return (*this)((int)c(1),(int)c(2),(int)c(3)); };
	inline Real&
		operator () (int i){return m_pBuf[i];};
	void		Multiply(Real d){for (int i=0; i<m_nElements; i++) m_pBuf[i]*=d;};
	void		operator+= (Real d){for (int i=0; i<m_nElements; i++) m_pBuf[i]+=d;};
	void		operator+= (Volume& v){
		for (int i=0; i<m_nElements; i++) m_pBuf[i]+=v.m_pBuf[i];
	};	
	void		operator-= (Real d){for (int i=0; i<m_nElements; i++) m_pBuf[i]-=d;};
	void		operator-= (Volume& v){for (int i=0; i<m_nElements; i++) m_pBuf[i]-=v.m_pBuf[i];};
	void		operator = (Real d){for (int i=0; i<m_nElements; i++) m_pBuf[i]=d;};

	int			operator != (Real d) {int n=0; for (int i=0; i<m_nElements; i++) if (m_pBuf[i] != d ) n++; return n;};
	int			operator == (Real d) {int n=0; for (int i=0; i<m_nElements; i++) if (m_pBuf[i] == d ) n++; return n;};
	int			operator > (Real d) {int n=0; for (int i=0; i<m_nElements; i++) if (m_pBuf[i] > d ) n++; return n;};
	int			operator < (Real d) {int n=0; for (int i=0; i<m_nElements; i++) if (m_pBuf[i] < d ) n++; return n;};
	int			operator >= (Real d) {int n=0; for (int i=0; i<m_nElements; i++) if (m_pBuf[i] >= d ) n++; return n;};
	int			operator <= (Real d) {int n=0; for (int i=0; i<m_nElements; i++) if (m_pBuf[i] <= d ) n++; return n;};

	void		Log(){for (int i=0; i<m_nElements; i++) m_pBuf[i]=log(m_pBuf[i]);};
	void		Exp(){for (int i=0; i<m_nElements; i++) m_pBuf[i]=exp(m_pBuf[i]);};
	void		Abs(){for (int i=0; i<m_nElements; i++) if(m_pBuf[i]<0) m_pBuf[i]=-m_pBuf[i];};
	void		ReplaceNaNs(Real val)
	{
#ifdef _WIN32
			for(int n=0; n<m_nElements; n++) if( _isnan(m_pBuf[n])) m_pBuf[n]=val; 
#else
			for(int n=0; n<m_nElements; n++) if( std::isnan(m_pBuf[n])) m_pBuf[n]=val;
#endif
	}
	void		Angle(Volume& re, Volume& im)
	{
		static const double pi=3.141592653589793238462643383;
		InitMemory(re.m_dims);
		double x,y;
		for(int i=0; i<m_nElements; i++) 
		{
			x=re.m_pBuf[i]; y=im.m_pBuf[i];
			m_pBuf[i]=(Real)atan2(y,x);
		}
	};
	void		Binarize(){for (int i=0; i<m_nElements; i++) m_pBuf[i]=(m_pBuf[i]!=0)?1.0:0.0;};
	void		Binarize(Real r){for (int i=0; i<m_nElements; i++) m_pBuf[i]=(m_pBuf[i]==r)?1.0:0.0;};
	void		ZeroLt(Real r){for (int i=0; i<m_nElements; i++) m_pBuf[i]=(m_pBuf[i]<=r)?0:m_pBuf[i];};
	void		ZeroGt(Real r){for (int i=0; i<m_nElements; i++) m_pBuf[i]=(m_pBuf[i]>=r)?0:m_pBuf[i];};
	void		ZeroLtGt(Real l, Real h){for (int i=0; i<m_nElements; i++) m_pBuf[i]=(m_pBuf[i]<=l || m_pBuf[i]>=h)?0:m_pBuf[i];};
	void		InvertBinary(){for (int i=0; i<m_nElements; i++) m_pBuf[i]=(m_pBuf[i]>0)?0.0:1.0;};
	void		Stats(Matrix& stats,Volume* mask=NULL);
	void		GetMaxMin(Real& mx, Real& mn){Matrix stats; Stats(stats); mx=stats(2,1); mn=stats(1,1);};
	Real		GetMax() { Matrix stats; Stats(stats); return stats(2,1); };
	Real*		GetBuf(){return m_pBuf;};
	void		SetCurPt(int x, int y, int z){m_CurPt=Index3(x,y,z);};
	void		NearestValidVoxel(ColumnVector& X){
					X(1) = min((Real)(m_dims[0] - 1), (Real)(max(X(1), 0.0)));
					X(2) = min((Real)(m_dims[1] - 1), (Real)(max(X(2), 0.0)));
					X(3) = min((Real)(m_dims[2] - 1), (Real)(max(X(3), 0.0)));
	};
	void		SetVRange(double v0, double v1)
	{
		Matrix st;
		Stats(st);
		Real mn=st(1,1), mx=st(2,1);
		Real sl=(Real)((v1-v0)/(mx-mn)), inter=(Real)v0-sl*mn;
		for (int i=0; i<m_nElements; i++) m_pBuf[i]=m_pBuf[i]*sl+inter;
	}
	void		ReSize(int newW, int newH, int newD, int x0, int y0, int z0, Real fill, Volume &resized);
	bool		WriteAnalyze(char* fname, bool bSwap=false);
	bool		Read(char* root, COORD_FORMAT fmt)
	{
#ifdef _4DFP
		if (fmt==FDFP){ m_coord_format=FDFP; return Read4dfp(root); }
#else
		if (fmt == FSL) {	m_coord_format = FSL; return  ReadAnalyze(root); }
		return false;
#endif
	};

	bool		Read(char* root)
	{
#ifdef _4DFP
		m_coord_format=FSL;
		if (Read4dfp(root)) {m_coord_format=FDFP; return true;}
		return  ReadAnalyze(root);
#else
		m_coord_format = FSL;
		return ReadAnalyze(root);
#endif
	};
	bool		Write(char* root, COORD_FORMAT fmt=FSL, int argc=0, char* argv[] = NULL)
	{
		if (fmt==FSL) return WriteAnalyze(root);
		else if(fmt==FDFP)
		{
#ifdef _4DFP
			return Write4dfp(root,argc,argv);
#else
			return false;
#endif
		}
		return false;
	};
	bool		ReadAnalyze(char* fname);
	bool		ReadAnalyzeHeader(char *fname);
	bool		IsVoxDims(){return m_bVoxDimInited;};
#ifdef _4DFP
	bool		Read4dfp(char* root);
	bool		Write4dfp(char* root, int argc, char *argv[]);
	void		Get4dfpAxisFlips(Real* flips)
	{
		if (!m_pIFH) {flips[0]=flips[1]=flips[2]=0.0;return;}		
		flips[0]=flips[1]=flips[2]=1.0;
		switch (m_pIFH->orientation)
		{
			case 2: flips[0]=flips[1]=-1; break;
			case 3: break;
			case 4: flips[1]=flips[2]=-1; break;
		}
	};

#endif
	inline int	SX(){return m_dims[0];};
	inline int	SY(){return m_dims[1];};
	inline int  SZ(){return m_dims[2];};

	inline Real	DX(){return m_voxel_dims[0];};
	inline Real	DY(){return m_voxel_dims[1];};
	inline Real	DZ(){return m_voxel_dims[2];};
	void	Metric_Vox(DiagonalMatrix& M2V, DiagonalMatrix& V2M){
		V2M.ReSize(4); V2M(1)=DX(); V2M(2)=DY(); V2M(3)=DZ(); V2M(4)=1;	M2V=V2M.i();
	};
	Real 		Sum(Volume* mask=NULL)
	{
		Real s=0; 
		if (!mask) {for(int i=0; i<m_nElements;i++) s+=m_pBuf[i];}
		else {for(int i=0; i<m_nElements;i++) {if ((*mask)(i)) s+=m_pBuf[i];}}
		return s;
	};
	Real		SumAbs(Volume* mask=NULL)
	{
		Real s=0; 
		if (!mask) {for(int i=0; i<m_nElements;i++) s+=fabs(m_pBuf[i]);}
		else {for(int i=0; i<m_nElements;i++) {if ((*mask)(i)) s+=fabs(m_pBuf[i]);}}
		return s;
	};

#define _mn(x,y) ((x)<(y))?(x):(y)
#define _mx(x,y) ((x)>(y))?(x):(y)
#define _floor(x) (int)(x)
#define _ceil(x) (int)((Real)(x)+1)

	Real	Interp3(ColumnVector& X);
	inline Real Interp3(Real x, Real y, Real z)
	{
		static Real z_old, y_old, x_old, z_new,x_new,y_new,x0,x1,y0,y1,z0,z1,xd,yd,zd,xd1,yd1,zd1;
		static Real c00,c10,c01,c11,c0,c1;

		z_old=z;
		z0=floor(z_old); z1=_mn(SZ()-1,_ceil(z_old));
		zd=(z_old-z0); zd1=1-zd;
		y_old=y;
		y0=floor(y_old); y1=_mn(SY()-1,_ceil(y_old));
		yd=(y_old-y0); yd1=1-yd;
		x_old=x;
		x0=floor(x_old); x1=_mn(SX()-1,_ceil(x_old));
		xd=(x_old-x0);xd1=1-xd;
		c00=(*this)(x0,y0,z0)*xd1+(*this)(x1,y0,z0)*xd;
		c10=(*this)(x0,y1,z0)*xd1+(*this)(x1,y1,z0)*xd;
		c01=(*this)(x0,y0,z1)*xd1+(*this)(x1,y0,z1)*xd;
		c11=(*this)(x0,y1,z1)*xd1+(*this)(x1,y1,z1)*xd;
		c0=c00*yd1+c10*yd;
		c1=c01*yd1+c11*yd;
		return c0*zd1+c1*zd;
	};

	inline Real	VoxSz(){ if (!IsVoxDims()) return 1; return m_voxel_dims[0]*m_voxel_dims[1]*m_voxel_dims[2];};
	inline bool IsInsideBoundary(int x,int y,int z){
		return ( x>=0 && y>=0 && z>=0 && x<m_dims[0] && y<m_dims[1] && z<m_dims[2] );
	};
	inline bool IsInsideBoundary(ColumnVector& X){
		return ( X(1)>=0 && X(2)>=0 && X(3)>=0 && X(1)<m_dims[0] && X(2)<m_dims[1] && X(3)<m_dims[2] );
	}
	void		VoxDims(ColumnVector& dims)
	{
		if (!IsVoxDims()) { dims=1; return; }
		dims(1)=m_voxel_dims[0];	dims(2)=m_voxel_dims[1];	dims(3)=m_voxel_dims[2];
	}
	bool		Range_thr(Range3& subrange, Real thr)
	{		
		int minX=m_dims[0],minY=m_dims[1],minZ=m_dims[0],maxX=0,maxY=0,maxZ=0;
		for (int z=0; z<SZ(); z++)
			for(int y=0; y<SY(); y++)
				for(int x=0; x<SX(); x++)
				{
					if ( (*this)(x,y,z)>=thr )
					{
						minX=min(minX,x);minY=min(minY,y);minZ=min(minZ,z);
						maxX=max(maxX,x);maxY=max(maxY,y);maxZ=max(maxZ,z);
					}
				}
		if (maxX==0) return false;
		Range3 r3(minX,maxX,minY,maxY,minZ,maxZ);
		subrange=r3;
		return true;
	}
	operator	Real* (){return m_pBuf;};
	inline Real&
		operator[](int n){return m_pBuf[n];};
	Real& At(int dim1, int dim2, int dim3)
				{return m_pBuf[dim3*m_dims[1]*m_dims[0]+dim2*m_dims[0]+dim1];};
	inline Real& 	
				operator()(int dim1,int dim2,int dim3)
				{
					return m_pBuf[dim3*m_dims[1]*m_dims[0]+dim2*m_dims[0]+dim1];
				};
	inline Real& 
		GetOff(int x, int y, int z){return m_pBuf[(m_offset[2]+z)*m_dims[1]*m_dims[0]+(m_offset[1]+y)*m_dims[0]+m_offset[0]+x];};
	inline Real&
				El(int i){return m_pBuf[i];};
	inline Real& 
				ElCurd(){return m_pBuf[m_ElCur--];};
	inline Real& 
				ElCuri(){return m_pBuf[m_ElCur++];};
	inline Real& 
				ElCur(){return m_pBuf[m_ElCur];};
	inline Real& 
				GetX(int x){return (*this)(x,m_CurPt.d[1],m_CurPt.d[2]);};
	inline Real&
				GetY(int y){return (*this)(m_CurPt.d[0],y,m_CurPt.d[2]);};
	inline Real& 
				GetZ(int z){return (*this)(m_CurPt.d[0],m_CurPt.d[1],z);};
	inline bool FitPt(int x, int y,int z )
	{
				return (x>=0 && y>=0 && z>=0 && x<SX() && y<SY() && z<SZ());
	}
	void		GetSlice(Matrix& m, int ind)
				{
					m.ReSize(SY(),SX());
					memcpy((void*)m.Store(),(void*)(m_pBuf+m_dims[0]*m_dims[1]*ind), m_dims[0]*m_dims[1]*sizeof(Real));
				};
	void		SetSlice(Matrix& m, int ind)
				{
					memcpy((void*)(m_pBuf+m_dims[0]*m_dims[1]*ind), (void*)m.Store(), m_dims[0]*m_dims[1]*sizeof(Real));
				};
	void		Center(ColumnVector& c, COORD_FORMAT f)
	{
		c.ReSize(4);
		if (f==FSL)
			c << 0 << 0 << 0 << 1;
		else
			c << ((Real)(m_dims[0]-1))*.5 << ((Real)(m_dims[1]-1))*.5 << ((Real)(m_dims[2]-1))*.5 << 1;
	};
	void		Get1D(ColumnVector& v, int x, int y, int z, int dim){Get1D(v.Store(),x,y,z,dim);};
	void		Get1D(Real* buf, int x, int y, int z, int dim);
	void		Set1D(ColumnVector& v, int x, int y, int z, int dim){Set1D(v.Store(),x,y,z,dim);}
	void		Set1D(Real* buf, int x, int y, int z, int dim);
	Volume&		operator = (Volume& a);
	void		Subtract(Volume&a, Volume& res)
	{
		res.InitMemory(m_dims);
		for (int i=0; i<m_nElements; i++) res.m_pBuf[i]=m_pBuf[i]-a.m_pBuf[i];
	}
	void		Divide (Volume& a, Volume& res)
	{
		res.InitMemory(m_dims); 
		for (int i=0; i<m_nElements; i++) res.m_pBuf[i]=m_pBuf[i]/a.m_pBuf[i];
	};
	void		Print()	{for(int i=0; i<SZ();i++) PrintSlice(i);};
	void		AddPt(ColumnVector& pt, Real val, Real bg);
	void		AddCurve(ColumnVector& A, ColumnVector& B, Matrix& coefs, Real val, void (*func)(Matrix&,Real,ColumnVector&));
	void		AddLine	(ColumnVector& A, ColumnVector& B, Real val);
	void		PrintSlice(int s)
	{
		Matrix sl(SX(),SY());
		GetSlice(sl,s);
		cout << "Slice "<<s<<":"<<endl<<sl<<endl;
	};
	void		Square() {for (int i=0; i<m_nElements; i++) m_pBuf[i]*=m_pBuf[i];};
	void		Sqrt(){for (int i=0; i<m_nElements; i++) m_pBuf[i]=sqrt(m_pBuf[i]);};
	//complex modulus.
	static void	AbsC(Volume& re, Volume& im, Volume& res)
	{
		res.InitMemory(re.m_dims);
		res.SetVoxelDims(re.m_voxel_dims);
		for(int i=0; i<re.m_nElements; i++)
			res[i]=sqrt(re[i]*re[i]+im[i]*im[i]);
	};
	static void	ExpToReIm(Volume& abs, Volume& phi, Volume& re, Volume& im)
	{
		for (int i=0; i<abs.m_nElements; i++) 
		{
			re.m_pBuf[i]=abs.m_pBuf[i]*cos(phi.m_pBuf[i]);
			im.m_pBuf[i]=abs.m_pBuf[i]*sin(phi.m_pBuf[i]);
		}
	};
	static void	Swap4(char* a);
#ifdef _4DFP
	bool	Vrtflip();
	inline ColumnVector	X4dfp2ecat(ColumnVector& X){
		ColumnVector X1; X1=X;
		switch (m_pIFH->orientation)
		{
			case 2: X1(1)=SX()-1-X(1); X1(3)=SZ()-1-X(3); return X1; //transvers
			case 3: X1(1)=SX()-1-X(1); return X1; //coronal
			default: return X1; //sagittal
		}
	}
#endif

private:

	bool	m_bVoxDimInited;
	int m_ElCur;
	Real *m_pBuf;
	Index3 m_CurPt;
};
#endif //_VOLUME_H_INCLUDED_
