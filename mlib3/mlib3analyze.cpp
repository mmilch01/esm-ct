/*******************************************************************************
Copyright (c) 2009
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "mlib3analyze.h"
const unsigned short Analyze::DT_NONE = 0;
const unsigned short Analyze::DT_UNKNOWN = 0;
const unsigned short Analyze::DT_BINARY = 1;
const unsigned short Analyze::DT_UNSIGNED_CHAR = 2;
const unsigned short Analyze::DT_SIGNED_SHORT = 4;
const unsigned short Analyze::DT_SIGNED_INT = 8;
const unsigned short Analyze::DT_FLOAT = 16;
const unsigned short Analyze::DT_COMPLEX = 32;
const unsigned short Analyze::DT_DOUBLE = 64;
const unsigned short Analyze::DT_RGB = 128;
const unsigned short Analyze::DT_ALL = 256;

/*************************************************************
*	Read Analyze header .hdr from file pointer
**************************************************************/
bool Analyze::ReadHeader(FILE* fp)
{
	int nRead=fread(&m_dsr,1,sizeof(struct dsr),fp);
	int sz=sizeof(struct dsr);
	int err=ferror(fp);
	int eof=feof(fp);
	if(err || eof || nRead <=0 ) return false;
//	if(nRead != sizeof(struct dsr)) return false;
	m_bSwapEndian=false;
	if(m_dsr.dime.dim[0]<0 || m_dsr.dime.dim[0]>15) 
	{
		m_bSwapEndian=true;
		swap_hdr(&m_dsr);
	}
	return true;
}

/*************************************************************
*	Read Analyze header .hdr from file name
**************************************************************/
bool Analyze::ReadHeader(char* file)
{
	FILE* fp;
	if(!(fp=fopen(file,"rb"))) return false;
	bool bRes=ReadHeader(fp);
	fclose(fp);
	return bRes;
}

/*************************************************************
*	Swap bytes in analyze header (little<->big endian)
**************************************************************/
void Analyze::swap_hdr(struct dsr *pntr)
{
	swap_long((unsigned char*)&pntr->hk.sizeof_hdr);
	swap_long((unsigned char*)&pntr->hk.extents);
	swap_short((unsigned char*)&pntr->hk.session_error);
	swap_short((unsigned char*)&pntr->dime.dim[0]);
	swap_short((unsigned char*)&pntr->dime.dim[1]);
	swap_short((unsigned char*)&pntr->dime.dim[2]);
	swap_short((unsigned char*)&pntr->dime.dim[3]);
	swap_short((unsigned char*)&pntr->dime.dim[4]);
	swap_short((unsigned char*)&pntr->dime.dim[5]);
	swap_short((unsigned char*)&pntr->dime.dim[6]);
	swap_short((unsigned char*)&pntr->dime.dim[7]);

	swap_short((unsigned char*)&pntr->dime.unused1);
	swap_short((unsigned char*)&pntr->dime.datatype);
	swap_short((unsigned char*)&pntr->dime.bitpix);
	swap_long((unsigned char*)&pntr->dime.pixdim[0]);
	swap_long((unsigned char*)&pntr->dime.pixdim[1]);
	swap_long((unsigned char*)&pntr->dime.pixdim[2]);
	swap_long((unsigned char*)&pntr->dime.pixdim[3]);
	swap_long((unsigned char*)&pntr->dime.pixdim[4]);
	swap_long((unsigned char*)&pntr->dime.pixdim[5]);
	swap_long((unsigned char*)&pntr->dime.pixdim[6]);
	swap_long((unsigned char*)&pntr->dime.pixdim[7]);
	swap_long((unsigned char*)&pntr->dime.vox_offset);
	swap_long((unsigned char*)&pntr->dime.funused1);
	swap_long((unsigned char*)&pntr->dime.funused2);
	swap_long((unsigned char*)&pntr->dime.cal_max);
	swap_long((unsigned char*)&pntr->dime.cal_min);
	swap_long((unsigned char*)&pntr->dime.compressed);
	swap_long((unsigned char*)&pntr->dime.verified);
	swap_short((unsigned char*)&pntr->dime.dim_un0);
	swap_long((unsigned char*)&pntr->dime.glmax);
	swap_long((unsigned char*)&pntr->dime.glmin);
}
/*************************************************************
*	Fill default header values
**************************************************************/
void Analyze::SetDefaultHeader()
{
	memset((void*)(&m_dsr),0,sizeof(m_dsr));
	m_dsr.hk.regular='r';
	m_dsr.hk.sizeof_hdr=sizeof(struct dsr);
	m_dsr.hk.extents=16384;

}

/*************************************************************
*	Endian swap utility functions.
**************************************************************/
void Analyze::swap_long(unsigned char* pntr)
{
	unsigned char b0, b1, b2, b3;
	b0 = *pntr;
	b1 = *(pntr+1);
	b2 = *(pntr+2);
	b3 = *(pntr+3);
	*pntr = b3;
	*(pntr+1) = b2;
	*(pntr+2) = b1;
	*(pntr+3) = b0;
}
void Analyze::swap_short(unsigned char* pntr)
{
	unsigned char b0, b1;
	b0 = *pntr;
	b1 = *(pntr+1);
	*pntr = b1;
	*(pntr+1) = b0;
}
/*************************************************************
*	Write analyze header .hdr file
**************************************************************/
bool Analyze::WriteHeader(const char* file, bool bSwap/*=false*/)
{
	FILE* fp;
	if(!(fp=fopen(file,"wb"))) return false;
	if(bSwap) swap_hdr(&m_dsr);
	bool bRes=(fwrite(&m_dsr,sizeof(m_dsr),1,fp)==1);
	fclose(fp);
	if(bSwap) swap_hdr(&m_dsr); //swap the header back
	return bRes;
}
bool Analyze::WritePixelBuf(FILE* fp, unsigned char* buf, int sz)
{
//	int sz=m_dsr.dime.dim[1]*m_dsr.dime.dim[2]*m_dsr.dime.dim[3];
	if(sz<=0) return false;
	int bytes_per_pixel=2;
	switch (m_dsr.dime.datatype)
	{
		case DT_SIGNED_SHORT: bytes_per_pixel=2; break;
		case DT_UNSIGNED_CHAR: bytes_per_pixel=1; break;
		case DT_FLOAT: bytes_per_pixel=4; break;
		case DT_DOUBLE: bytes_per_pixel=8; break;
		default: return false;
	}
	return fwrite(buf,bytes_per_pixel,sz,fp)==sz;
}

bool Analyze::WritePixels(char* file, unsigned char* buf)
{
	FILE* fp;
	int sz=m_dsr.dime.dim[1]*m_dsr.dime.dim[2]*m_dsr.dime.dim[3]*m_dsr.dime.dim[4];
	if(!(fp=fopen(file,"wb"))) return false;
	bool bRes=WritePixelBuf(fp,buf,sz);
	fclose(fp);
	return bRes;
}
bool Analyze::WriteAll(char* root, unsigned char* buf, bool bSwapHeader/*=false*/)
{
	char hdr[512],img[512];
	GetFileNames(root,hdr,img);
	return (WriteHeader(hdr, bSwapHeader) && WritePixels(img,buf));
}
unsigned char* Analyze::ReadPixels(char* filename, int& bits_per_pixel)
{
	FILE* fp;
	if(!(fp=fopen(filename,"rb"))) {fclose(fp); return NULL;}
	unsigned char* res=ReadPixelsFromFp(fp,bits_per_pixel);
	fclose(fp);
	return res;
}
unsigned char* Analyze::ReadPixelsFromFp(FILE* fp, int& bits_per_pixel)
{
	long nVoxels=m_dsr.dime.dim[1]*m_dsr.dime.dim[2]*m_dsr.dime.dim[3];
	if(m_dsr.dime.datatype!=Analyze::DT_SIGNED_SHORT && //only supported are 16 bit int and float
		m_dsr.dime.datatype!=Analyze::DT_FLOAT && 
		m_dsr.dime.datatype!=Analyze::DT_UNSIGNED_CHAR &&
		m_dsr.dime.datatype!=Analyze::DT_DOUBLE &&
		m_dsr.dime.datatype!=Analyze::DT_SIGNED_INT) return NULL;
	if (nVoxels<=0) return NULL;

	int nSlices=m_dsr.dime.dim[3];
	if(nSlices<1) return NULL;

	unsigned char* pd=NULL;
//	if(m_dsr.dime.dim[0]>3 && m_dsr.dime.dim[4]>1) //multiple volumes
//		nVoxels*=m_dsr.dime.dim[4];

	if(m_dsr.dime.datatype==Analyze::DT_SIGNED_INT)
	{
		pd=(unsigned char*)(new int[nVoxels]);
		bits_per_pixel=32;
		if(!ReadToBuffer(fp,pd,nVoxels,4,m_dsr.dime.datatype)){delete[] pd; return NULL;}
	}
	if(m_dsr.dime.datatype==Analyze::DT_SIGNED_SHORT)
	{
		pd=(unsigned char*)(new short[nVoxels]);
		bits_per_pixel=16;
		if(!ReadToBuffer(fp,pd,nVoxels,2,m_dsr.dime.datatype)){delete[] pd; return NULL;}
	}
	else if(m_dsr.dime.datatype==Analyze::DT_UNSIGNED_CHAR) //8 bit int
	{
		pd=(unsigned char*)(new unsigned char[nVoxels]);
		bits_per_pixel=8;
		if(!ReadToBuffer(fp,pd,nVoxels,1,m_dsr.dime.datatype)){delete[] pd; return NULL;}
	}
	else if(m_dsr.dime.datatype==Analyze::DT_FLOAT)
	{
		pd=(unsigned char*)(new float[nVoxels]);
		bits_per_pixel=32;
		if(!ReadToBuffer(fp,pd,nVoxels,4,m_dsr.dime.datatype)){delete[] pd; return NULL;}
	}
	else if(m_dsr.dime.datatype==Analyze::DT_DOUBLE)
	{
		pd=(unsigned char*)(new double[nVoxels]);
		bits_per_pixel=64;
		if(!ReadToBuffer(fp,pd,nVoxels,8,m_dsr.dime.datatype)){delete[] pd; return NULL;}
	}
	return pd;
}

void Analyze::GetFileNames(char* root, char* hdrfile, char* imfile)
{
	int len=strlen(root);
	strcpy(imfile,root);
	strcpy(hdrfile,root);
	if(!strcmp(root+len-4,".img") || !strcmp(root+len-4,".hdr"))
	{
		len-=4;
		imfile[len]=0;
		hdrfile[len]=0;
	}
	strcat(imfile,".img");
	strcat(hdrfile,".hdr");
}
unsigned char* Analyze::ReadAll(char* root)
{
	char imfile[512], hdrfile[512];
	GetFileNames(root,hdrfile,imfile);
	if(!ReadHeader(hdrfile)) return NULL;
	int bits_per_pixel;
	return ReadPixels(imfile,bits_per_pixel);
}
bool Analyze::ReadToBuffer(FILE* fp, unsigned char* buf, long nUnits, int unit_size, int dt)
{
	if(unit_size==1 || (!m_bSwapEndian && unit_size==2) )
		return (fread(buf,unit_size,nUnits,fp)==nUnits);
	if(unit_size==2)
	{
		short* ptr=(short*)buf;
		short wd;
		for(int i=0; i<nUnits; i++)
		{
			if(fread(&wd,2,1,fp)!=1) return false;
			swap_short((unsigned char*)&wd);
			ptr[i]=wd;
		}
		return true;
	}
	else if (unit_size==4)
	{
		if (dt==Analyze::DT_FLOAT)
		{
			float dwd; 
			for(int i=0; i<nUnits; i++)
			{
				if(fread(&dwd,4,1,fp)!=1) return false;
				if(m_bSwapEndian) swap_long((unsigned char*)&dwd);
				((float*)(buf))[i]=dwd;
			}
			return true;
		}
		else if(dt==Analyze::DT_SIGNED_INT)
		{
			int dwd;
			for(int i=0; i<nUnits; i++)
			{
				if(fread(&dwd,4,1,fp)!=1) return false;
				if(m_bSwapEndian) swap_long((unsigned char*)&dwd);
				((int*)(buf))[i]=dwd;
			}
			return true;
		}
		else return false;
	}
	else if (unit_size==8)
	{
		double qdwd; 
		for(int i=0; i<nUnits; i++)
		{
			if(fread(&qdwd,8,1,fp)!=1) return false;
//			if(m_bSwapEndian) swap_long((unsigned char*)&dwd);
			((double*)(buf))[i]=qdwd;
		}
		return true;
	}
	return false;
}
int Analyze::PixelCount()
{
	int cnt=1;
	for (int i=1; i<8; i++)
	{
		if(m_dsr.dime.pixdim[i]>0)
			cnt*=(int)(m_dsr.dime.pixdim[i]+0.5);
		else break;
	}
	return cnt;
}
