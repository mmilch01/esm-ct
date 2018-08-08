/*******************************************************************************
Copyright (c) 2009
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#if !defined(_ANALYZE_H__INCLUDED_)
#define _ANALYZE_H__INCLUDED_
#include <stdio.h>

class Analyze 
{
public:
	static const unsigned short DT_NONE, DT_UNKNOWN, DT_BINARY, DT_UNSIGNED_CHAR, DT_SIGNED_SHORT, DT_SIGNED_INT,
		DT_FLOAT, DT_COMPLEX, DT_DOUBLE, DT_RGB, DT_ALL;
	struct header_key				// header key 
	{								// off + size 
		int sizeof_hdr;				// 0 + 4 
		char data_type[10];			// 4 + 10 
		char db_name[18];			// 14 + 18 
		int extents;				// 32 + 4 
		short int session_error;	// 36 + 2 
		char regular;				// 38 + 1 
		char hkey_un0;				// 39 + 1 
	};								// total=40 bytes 
	struct image_dimension
	{								// off + size 
		short int dim[8];			// 0 + 16 

		unsigned char vox_units[4]; // 16 + 4
		unsigned char cal_units[8]; // 20 + 8
		short int unused1;			// 28 + 2

		short int datatype;			// 30 + 2 
		short int bitpix;			// 32 + 2 
		short int dim_un0;			// 34 + 2 
		float pixdim[8];			// 36 + 32 
		/*
		pixdim[] specifies the voxel dimensitons:
		pixdim[1] - voxel width
		pixdim[2] - voxel height
		pixdim[3] - interslice distance
		...etc
		*/
		float vox_offset;			// 68 + 4 
		float funused1;				// 72 + 4 
		float funused2;				// 76 + 4 
		float funused3;				// 80 + 4 
		float cal_max;				// 84 + 4 
		float cal_min;				// 88 + 4 
		float compressed;			// 92 + 4 
		float verified;				// 96 + 4 
		int glmax,glmin;			// 100 + 8 
	};								// total=108 bytes 
	struct data_history
	{								// off + size 
		char descrip[80];			// 0 + 80 
		char aux_file[24];			// 80 + 24 
		char orient;				// 104 + 1 
		char originator[10];		// 105 + 10 
		char generated[10];			// 115 + 10 
		char scannum[10];			// 125 + 10 
		char patient_id[10];		// 135 + 10 
		char exp_date[10];			// 145 + 10 
		char exp_time[10];			// 155 + 10 
		char hist_un0[3];			// 165 + 3 
		int views;					// 168 + 4 
		int vols_added;				// 172 + 4 
		int start_field;			// 176 + 4 
		int field_skip;				// 180 + 4 
		int omax, omin;				// 184 + 8 
		int smax, smin;				// 192 + 8 
	};
	struct dsr
	{
		struct header_key hk;		// 0 + 40 
		struct image_dimension dime;// 40 + 108 
		struct data_history hist;	// 148 + 200 
	} m_dsr;
	Analyze(){SetDefaultHeader();};
	~Analyze(){};
	static void swap_long(unsigned char* pntr);
	static void swap_short(unsigned char* pntr);
	void SetDefaultHeader();
	void	GetFileNames(char* root, char* imfile, char* hdrfile);
	bool ReadHeader(char* file);
	bool ReadHeaderOnly(char* root)
	{
		char imfile[512], hdrfile[512];
		GetFileNames(root,hdrfile,imfile);
		if(!ReadHeader(hdrfile)) return false;
		return true;
	}

	bool ReadHeader(FILE* fp);
	bool WriteHeader(const char* file, bool bSwap=false);
	bool WritePixelBuf(FILE* fp, unsigned char* buf, int len);
	bool WritePixels(char* file, unsigned char* buf);
	bool WriteAll(char* root, unsigned char* buf, bool bSwapHeader=false);
	unsigned char* ReadPixels(char* filename, int& bits_per_pixel);
	unsigned char* ReadPixelsFromFp(FILE* fp, int& bits_per_pixel);
	unsigned char* ReadAll(char* root);
	int  PixelCount();
private:
	bool	m_bSwapEndian;

	void	swap_hdr(struct dsr *pntr);
	bool	ReadToBuffer(FILE* fp, unsigned char* buf, long nUnits, int unit_size,int datatype);
};
#endif //#define _ANALYZE_INCLUDED
