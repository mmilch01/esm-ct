/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#include "mlib3volume.h"
#include "mlib3consoleutil.h"
#ifdef _4DFP
extern "C"
{
	#include <Getifh.h>
	#include <rec.h>
}
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

Volume::INTERP_METHOD Volume::s_interp_method=Volume::TRILINEAR;
Volume::VOXEL_CENSORING_METHOD Volume::s_vox_cens_method=Volume::NONE;

void Volume::InitVars()
{
#ifdef _4DFP
	m_pIFH=NULL;
#endif
	m_bVoxDimInited=true;
	m_voxel_dims[0]=m_voxel_dims[1]=m_voxel_dims[2]=1;
	m_pBuf=NULL;
	m_coord_format=FSL;
	m_ElCur=0;
	m_nElements=0;
	m_offset[0]=m_offset[1]=m_offset[2]=0;
}
/******************************************************
* Create volume with variable number of dimensions.
*******************************************************/
void Volume::InitMemory(int dim1, int dim2, int dim3)
{
	int dims[]={dim1, dim2, dim3}, nEl=1, i;
	//memory already allocated
	if(m_dims[0]==dim1 && m_dims[1]==dim2 && m_dims[2]==dim3 && m_pBuf) return;
	memcpy(m_dims,dims,3*sizeof(int));
	for (i=0; i<3; i++) nEl*=dims[i];
	bool bChangeSize=true;
	if(m_nElements==nEl && m_pBuf!=NULL) bChangeSize=false;		
	m_nElements=nEl;
	if(bChangeSize)
	{
		if(m_pBuf) delete[] m_pBuf;
		m_pBuf=new Real[m_nElements];
	}
	m_ElCur=0;
	(*this)=0;
}
// newW/H/D - new dimensions
// x/y/z0 - translate vector from old origin. 
// fill - fill value
void	Volume::ReSize(int newW, int newH, int newD, int x0, int y0, int z0, Real fill, Volume &resized)
{
	resized.InitMemory(newW,newH,newD);
	resized.SetVoxelDims(m_voxel_dims);
	int cx,cy,cz;
	bool wby,wbz,wbyz;
//	Real val;
	for(int uz=0; uz<newD; uz++)
	{
		cz=uz+z0;
		wbz=(cz>=0 && cz<SZ());
		for(int uy=0; uy<newH; uy++)
		{
			cy=uy+y0;
			wby=(cy>=0 && cy<SY());
			wbyz=(wbz && wby);
			if(!wbyz)
				for(int ux=0; ux<newW; ux++)	resized(ux,uy,uz)=fill;
			else 
				for(int ux=0; ux<newW; ux++)	
				{
					cx=ux+x0;
					if(cx>=0 && cx<SX())
					{
//						val=(*this)(cx,cy,cz);
						resized(ux,uy,uz)=(*this)(cx,cy,cz);
					}
					else
						resized(ux,uy,uz)=fill;
				}
		}
	}
}

Real Volume::Interp3(ColumnVector& X)
{
		if (s_interp_method==NNEIB)
		{
			if (!IsInsideBoundary(round(X(1)), round(X(2)), round(X(3)))) return 0;
			return (*this)(round(X(1)),round(X(2)),round(X(3)));
		}
		static Real z_old, y_old, x_old, z_new,x_new,y_new,x0,x1,y0,y1,z0,z1,xd,yd,zd,xd1,yd1,zd1;
		static Real c00,c10,c01,c11,c0,c1;
		z_old=X(3);
		z0=floor(z_old); z1=_mn(SZ()-1,_ceil(z_old));
		zd=(z_old-z0); zd1=1-zd;
		y_old=X(2);
		y0=floor(y_old); y1=_mn(SY()-1,_ceil(y_old));
		yd=(y_old-y0); yd1=1-yd;
		x_old=X(1);
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

bool	Volume::ReadAnalyzeHeader(char *fname)
{
	Analyze an;
	if ( ! an.ReadHeaderOnly(fname) ) return false;
	InitMemory(an,NULL);
	return true;
}
bool	Volume::ReadAnalyze(char* fname)
{
	Analyze an;
	unsigned char *imbuf;
	if(!(imbuf=an.ReadAll(fname))) 
	{
		if(imbuf) delete[] imbuf;
		return false;
	}
	InitMemory(an,imbuf);
	delete[] imbuf;
	return true;
}
void Volume::Swap4(char* a)
{
	char t;
	t = a[0]; a[0] = a[3]; a[3] = t;
	t = a[1]; a[1] = a[2]; a[2] = t;
}

void	Volume::Get1D(Real* buf, int x, int y, int z, int dim)
{
	if(dim==1) //x
	{
		memcpy(buf,m_pBuf+m_dims[0]*m_dims[1]*z+m_dims[0]*y,m_dims[0]*sizeof(Real));						
	}
	else if(dim==2) //y
	{
		Real* st=m_pBuf+m_dims[0]*m_dims[1]*z+x, *p;
		int i;
		for(i=0, p=st; i<m_dims[1]; i++, p+=m_dims[0]) buf[i]=p[0];
	}
	else if(dim==3) //z
	{
		Real* st=m_pBuf+m_dims[0]*y+x, *p;
		int pl=m_dims[0]*m_dims[1], i;
		for(i=0, p=st; i<m_dims[2]; i++, p+=pl) buf[i]=p[0];
	}
}

void	Volume::Set1D(Real* buf, int x, int y, int z, int dim)
{
	if(dim==1) //x
	{
		memcpy(m_pBuf+m_dims[0]*m_dims[1]*z+m_dims[0]*y,buf,m_dims[0]*sizeof(Real));
	}
	else if(dim==2) //y
	{
		Real* st=m_pBuf+m_dims[0]*m_dims[1]*z+x, *p;
		int i;
		for(i=0, p=st; i<m_dims[1]; i++, p+=m_dims[0]) p[0]=buf[i];
	}
	else if(dim==3) //z
	{
		Real* st=m_pBuf+m_dims[0]*y+x, *p;
		int pl=m_dims[0]*m_dims[1],i;
		for(i=0, p=st; i<m_dims[2]; i++, p+=pl) p[0]=buf[i];
	}
}
#ifdef _4DFP
bool	Volume::Read4dfp(char* root0)
{
//	if(cc.m_bAHeader) return ReadAnalyze(cc.m_root);

	char name[MAXL];
	char root[MAXL];
	char* occ=strstr(root0,".4dfp");

	if ( occ && strlen(occ)==5 )
	{
		strncpy(root, root0, strlen(root0)-5);
		root[strlen(root0) - 5]=0;
	}
	else strcpy(root, root0);

	sprintf (name,"%s.4dfp.ifh",root);
	if (!m_pIFH) m_pIFH=new IFH;
	if(Getifh(name, m_pIFH) != 0) { delete m_pIFH; m_pIFH=NULL; return false; }
	bool isbig=(m_pIFH->imagedata_byte_order[0]=='b');
	InitMemory(m_pIFH->matrix_size);
	if(m_nElements<1) return false;
	bool bSwap=(!IsBigEndian() && isbig) || (IsBigEndian() && !isbig);

	sprintf(name,"%s.4dfp.img",root);
	FILE* fp;
	if(!(fp=fopen(name,"rb"))) return false;
	if(sizeof(Real)==4)
	{
		if(!(fread(m_pBuf,4,m_nElements,fp)==m_nElements)) {fclose(fp); return false;}
		if(bSwap) for(int i=0; i<m_nElements; i++) Swap4((char*)(&m_pBuf[i]));
	}
	else if(sizeof(Real)==8)
	{
		float *buf=new float[m_nElements];
		if(!(fread(buf,4,m_nElements,fp)==m_nElements)) {fclose(fp); return false;}
		if(bSwap) for(int i=0; i<m_nElements; i++) Swap4((char*)(&buf[i]));
		for(int i=0; i<m_nElements; i++)
		{
			m_pBuf[i]=buf[i];
			i=i;
		}
		delete[] buf;
	}
	fclose(fp);	
	for(int i=0; i<3; i++) m_voxel_dims[i]=fabs((Real)(m_pIFH->mmppix[i]));
	return true;
}
bool Volume::Write4dfp(char* root, int argc, char *argv[])
{
	char name[MAXL];
	bool bSwap=false, bBigEndian=IsBigEndian(), bRes=true;
	if (!m_pIFH) return false;

	if(((m_pIFH->imagedata_byte_order[0]=='b') && !bBigEndian) ||
		((m_pIFH->imagedata_byte_order[0]=='l') && bBigEndian)) bSwap=true;
	sprintf (name,"%s.4dfp",root);
	bRes=WriteAnalyze(name, bSwap);
	if (bRes) bRes &= (0==Writeifh(argv[0],name,m_pIFH,IsBigEndian()? 'b':'l'));
	char rcsid[]="$Id: " __FILE__ ", " __DATE__ ", " __TIME__;
	bRes=(0==startrece(name,argc,argv,rcsid,IsBigEndian()? 'b':'l'));
	if (bRes) bRes &= (0==endrec());
	return bRes;
}
#endif
bool	Volume::WriteAnalyze(char* fname, bool bSwap/*=false*/)
{
	Analyze an;
	an.SetDefaultHeader();

	an.m_dsr.dime.dim[0]=4;
	an.m_dsr.dime.dim[1]=m_dims[0];
	an.m_dsr.dime.dim[2]=m_dims[1];
	an.m_dsr.dime.dim[3]=m_dims[2];
	an.m_dsr.dime.dim[4]=1;
	an.m_dsr.dime.datatype=Analyze::DT_FLOAT;
	an.m_dsr.dime.bitpix=32;

	an.m_dsr.dime.pixdim[1]=(float)m_voxel_dims[0];
	an.m_dsr.dime.pixdim[2]=(float)m_voxel_dims[1];
	an.m_dsr.dime.pixdim[3]=(float)m_voxel_dims[2];

	if(sizeof(Real)==8) //double, 8-byte
	{
		float* tmp=new float[m_nElements];
		for(int i=0; i< m_nElements; i++) 
		{
			tmp[i]=(float)(m_pBuf[i]);
			if(bSwap) Swap4((char*)(&tmp[i]));
		}
		bool res=an.WriteAll(fname, (unsigned char*)tmp, bSwap);
		delete[] tmp;
		return res;
	}
	else //float, 4-byte 
	{
		if(bSwap) for (int i=0; i<m_nElements; i++) Swap4((char*)(&m_pBuf[i]));
		return an.WriteAll(fname, (unsigned char*)m_pBuf, bSwap);
	}
}

/******************************************************
* Create an array from Analyze header and pixel buffer.
*******************************************************/
void Volume::InitMemory(Analyze& an, void* pBuf/*=NULL*/, int nSlice/*=-1*/)
{
	int dims[3];
	int i;
	for (i=1; i<4; i++) dims[i-1]=an.m_dsr.dime.dim[i];
	bool bFull=(nSlice<0);
	if(bFull && pBuf) InitMemory(dims);
	else {m_dims[0]=dims[0]; m_dims[1]=dims[1]; m_dims[2]=dims[2];}

	int ind1,ind2,ind3;

	if(an.m_dsr.dime.pixdim[1]*an.m_dsr.dime.pixdim[2]*an.m_dsr.dime.pixdim[3]>0)
	{
		m_bVoxDimInited=true;
		m_voxel_dims[0]=an.m_dsr.dime.pixdim[1];
		m_voxel_dims[1]=an.m_dsr.dime.pixdim[2];
		m_voxel_dims[2]=an.m_dsr.dime.pixdim[3];
	}
	if(!pBuf) return;
	if(an.m_dsr.dime.datatype==Analyze::DT_FLOAT)
	{
		float* tmpB=(float*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}
	else if(an.m_dsr.dime.datatype==Analyze::DT_SIGNED_INT)
	{
		int* tmpB=(int*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}

	else if(an.m_dsr.dime.datatype==Analyze::DT_SIGNED_SHORT)
	{
		short* tmpB=(short*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}
	else if(an.m_dsr.dime.datatype==Analyze::DT_UNSIGNED_CHAR)
	{
		unsigned char* tmpB=(unsigned char*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=(Real)(tmpB[ind3+x]);
				}
			}
		}
	}
	else if(an.m_dsr.dime.datatype==Analyze::DT_DOUBLE)
	{
		double* tmpB=(double*) pBuf;
		for (int z=(bFull? 0:nSlice);z<(bFull?m_dims[2]:(nSlice+1)); z++)
		{
			ind2=m_dims[1]*m_dims[0]*z;
			for(int y=0; y<m_dims[1]; y++)
			{
				ind1=ind2+m_dims[0]*y;
				ind3=(bFull?ind2:0)+m_dims[0]*y;
				for(int x=0; x<m_dims[0]; x++)
				{
					m_pBuf[ind1+x]=double(tmpB[ind3+x]);
				}
			}
		}
	}
	m_ElCur=0;
}
Volume& Volume::operator = (Volume& a)
{
	InitMemory(a.m_dims);
	memcpy(m_pBuf,a.m_pBuf,sizeof(Real)*m_nElements);
	memcpy(m_dims,a.m_dims,sizeof(int)*3);
	m_ElCur=a.m_ElCur;
	m_bVoxDimInited=a.m_bVoxDimInited;
	memcpy(m_voxel_dims,a.m_voxel_dims,sizeof(double)*3);
	memcpy(m_offset,a.m_offset,sizeof(double)*3);
	m_nElements=a.m_nElements;
	m_coord_format=a.m_coord_format;
#ifdef _4DFP
	if (a.m_pIFH)
	{
		m_pIFH=new IFH;
		memcpy(m_pIFH,a.m_pIFH,sizeof(IFH));
	}
#endif
	return (*this);
}
#define _r(x) (int)(x+0.5)

void Volume::AddPt(ColumnVector& pt, Real val, Real bg)
{
#define _n(x,y,z) sqrt((x)*(x)+(y)*(y)+(z)*(z))
#define _rr(x) (int)((x)+.5)
	static ColumnVector X[3], DX[3];
	for(int i=0; i<3; i++){
		X[i].ReSize(3); DX[i].ReSize(3);
	}
	static int neib[27][3]={
		{0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},{0,2,0},{1,2,0},{2,2,0},
		{0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1},{0,2,1},{1,2,1},{2,2,1},
		{0,0,2},{1,0,2},{2,0,2},{0,1,2},{1,1,2},{2,1,2},{0,2,2},{1,2,2},{2,2,2}
	};

	for(int i=1;i<=3; i++){
		X[0](i)=_rr(pt(i))-1;
		X[1](i)=_rr(pt(i));
		X[2](i)=_rr(pt(i))+1;
		DX[0](i)=pt(i)-X[0](i);
		DX[1](i)=pt(i)-X[1](i);
		DX[2](i)=pt(i)-X[2](i);
	}
	Real D, dVal=val-bg;
	int x,y,z;
	for (int i=0; i<27; i++){
		x=neib[i][0];
		y=neib[i][1];
		z=neib[i][2];
		D=_n(DX[x](1),DX[y](2),DX[z](3));
		(*this)(X[x](1),X[y](2),X[z](3))+=exp(-D*D)*dVal+bg;
	}
}

void Volume::AddCurve(ColumnVector& A, ColumnVector& B, Matrix& coefs, Real val, void (*func)(Matrix&,Real,ColumnVector&))
{
	Real t=0;
	int i,j=0;
	ColumnVector Vab=B-A;
	Real L=_r(Vab.MaximumAbsoluteValue1(i));
	ColumnVector X(3),Xst(3);
	X=A;
	Real dVox=_mx(Vab(1),Vab(2)); dVox=_mx(Vab(3),dVox);
	Real st=1.0/(dVox*3);
	ColumnVector curPt(3);
	do{
		func(coefs,t,curPt);
		if (IsInsideBoundary(curPt))
			AddPt(curPt,val,0);
		t+=st;
	}while (t <= 1 + 1e-6);
}

void Volume::AddLine(ColumnVector& A, ColumnVector& B, Real val)
{
	ColumnVector Vab=B-A;
	int i,j=0;
	Real L=_r(Vab.MaximumAbsoluteValue1(i));
	ColumnVector X(3),Xst(3);
	X=A;
	Xst=Vab/L;
	do{
		AddPt(X,val,0);
		X+=Xst;
		j++;
	}while(j<=L);
}
#ifdef _4DFP
bool	Volume::Vrtflip()
{
	if (!m_pIFH) return false;
	int k=m_pIFH->orientation-1;
	if ( k < 0 ) return false;
	Matrix flips(3,3);
	float *mmpixt=(m_pIFH->mmppix), *centert=(m_pIFH->center);
	flips << -1 << 1 << -1 << -1 << 1 << 1 << 1 << 1 << 1;
	for (int i=1; i<=3; i++)
	{
		mmpixt[i-1]*=flips(k,i);
		centert[i-1]*=flips(k,i);
		if (flips(k,i)<0) 
//			centert[i-1]=mmpixt[i-1]*m_pIFH->matrix_size[i-1]+1-centert[i-1];
			centert[i-1]=mmpixt[i-1]*(m_pIFH->matrix_size[i-1]+1)-centert[i-1];
	}
	return true;
}
#endif
/////////////////////////////////////////////
// Minimum and maximum
void Volume::Stats(Matrix& stats,Volume* mask)
{
	Real mn, mx;
	mn=mx=m_pBuf[0];
	int indmn=0, indmx=0;
	for (int i=0; i<m_nElements; i++) 
	{
		if (mask && !(*mask)(i)) continue;
		if(m_pBuf[i]<mn){ mn=m_pBuf[i]; indmn=i;};
		if(m_pBuf[i]>mx){ mx=m_pBuf[i]; indmx=i;};
	}
	stats.resize(2,2);
	stats(1,1)=mn;
	stats(1,2)=(Real)indmn;
	stats(2,1)=mx;
	stats(2,2)=(Real)indmx;
}