/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#ifndef _CONSOLECONFIG_INCLUDED_
#define _CONSOLECONFIG_INCLUDED_
#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>

#include "mlib3volume.h"

#ifdef WIN32
#endif

#ifdef _4DFP
extern "C"
{
	#include <rec.h>
	#include <ifh.h>
}
#endif

#define MAXL 256 //max string length

class ConsoleConfig
{
public:
	bool m_bAHeader;
	char m_program[MAXL], m_root[MAXL], m_control, **m_argv, m_rcsid[MAXL];
	int m_argc;

	ConsoleConfig(char* rcsid)
	{
		m_control='\0';
		m_bAHeader=false;
		strcpy(m_rcsid,rcsid);
		cout << m_rcsid << endl;
	};
	~ConsoleConfig(){};
	static vector<string> Tokenize(const string& str, const string& delimiters)
	{
		vector<string> tokens;
	    	
		// skip delimiters at beginning.
    		string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	    	
		// find first "non-delimiter".
    		string::size_type pos = str.find_first_of(delimiters, lastPos);

    		while (string::npos != pos || string::npos != lastPos)
    		{
        		// found a token, add it to the vector.
        		tokens.push_back(str.substr(lastPos, pos - lastPos));
			
        		// skip delimiters.  Note the "not_of"
        		lastPos = str.find_first_not_of(delimiters, pos);
			
        		// find next "non-delimiter"
        		pos = str.find_first_of(delimiters, lastPos);
    		}

		return tokens;
	}

	static void GetCoords(const char* param, double* buf, int cnt)
	{			
		string ln(param);
		string delim(",-");
		vector<string> vals=Tokenize(ln,delim);
//		int x0=0,y0=0,z0=0,x1=0,y1=0,z1=0;
		if(cnt>=3)
		{
			buf[0]=atof(vals.at(0).data());
			buf[1]=atof(vals.at(1).data());
			buf[2]=atof(vals.at(2).data());
			if(cnt>=6)
			{
				buf[3]=atof(vals.at(3).data());
				buf[4]=atof(vals.at(4).data());
				buf[5]=atof(vals.at(5).data());
			}
		}
	}
	void WriteRec(char* suffix)
	{
#ifdef _4DFP
#ifdef __GNUG__
		char name[MAXL], imgfile[MAXL];
		if(suffix) sprintf(name,"%s_%s.4dfp.rec",m_root,suffix);
		else sprintf (name,"%s.4dfp.rec",m_root);
		startrece (name, m_argc, m_argv, m_rcsid,IsBigEndian()?'b':'l');
		sprintf(imgfile,"%s.4dfp.img",m_root);
		catrec(imgfile);
		endrec();
#endif
#endif
	}
	static bool ReadVolume4dfp(Volume &v, char* name)
	{
		char root[MAXL],fname[MAXL];
		GetRoot(name,root);
		if(v.Read4dfp(root)) return true;
		sprintf(fname,"%s.4dfp",root);
		return v.Read4dfp(fname);
	};

	bool		ReadVolume(Volume &v, char* fname=NULL)
	{
		char name[MAXL];
		GetRoot((fname)?fname:m_root,name);
#ifdef _4DFP
		if(!m_bAHeader)
		{
			if(v.Read4dfp(name)) return true;
			sprintf(name,"%s.4dfp",name);
			if(v.Read4dfp(name)) return true;
		}
#endif
		GetRoot((fname)?fname:m_root,name);
		if(!v.ReadAnalyze(name))
		{
			sprintf(name,"%s.4dint",name);
			if(!v.ReadAnalyze(name))
			{
				fprintf(stderr, "\nError reading %s\n",m_root);
				return false;
			}
		}
		return true;
	}
	void		SaveVolume(Volume& v, char* suf, bool bWriteRec=false)
	{
		char nm[MAXL];
		GetRootShort(nm);
		printf("Writing %s_%s\n",nm,suf);
#ifdef _4DFP
		if(!m_bAHeader)
		{
			if(bWriteRec) WriteRec(suf);
			v.Write4dfp(m_root,m_argc,m_argv);
			return;
		}
#endif
		sprintf(nm+strlen(nm),"_%s",suf);
		v.WriteAnalyze(nm);
	}
	static void SaveMatr(Matrix& M, ofstream& ofs, int prec=6){
		int acc; Real val;
		int nRow = M.Nrows(), nCol = M.Ncols();
		vector<int> wid(nCol);

		for (int i = 0; i<nCol; i++) {
			wid[i] = 1;
			for (int j = 1; j <= nRow; j++) {
				val = M(j, i + 1); acc = 1;
				while (fabs(val) >= 10) { val *= .1; acc++; }
				if (val<0) acc++;
				wid[i] = _mx(wid[i], acc);
			}
		}
		for (int i = 1; i <= nRow; i++) {
			for (int j = 1; j <= nCol; j++)
				ofs << fixed << setprecision(prec) << setw(prec + 2 + wid[j - 1]) << M(i, j);
			ofs << endl;
		}
	}
	static void PrintMatr(Matrix& M, int prec=6){
		int acc; Real val;
		int nRow=M.Nrows(), nCol=M.Ncols();
		int old_precision = cout.precision(prec);
		vector<int> wid(nCol);

		for (int i = 0; i<nCol; i++) {
			wid[i] = 1;
			for (int j = 1; j <= nRow; j++) {
				val = M(j, i + 1); acc = 1;
				while (fabs(val) >= 10) { val *= .1; acc++; }
				if (val<0) acc++;
				wid[i] = _mx(wid[i], acc);
			}
		}
		for (int i = 1; i <= nRow; i++) {
			for (int j = 1; j <= nCol; j++)
				cout << fixed << setw(prec+2+wid[j - 1]) << M(i, j);
			cout << endl;
		}
		cout.precision(old_precision);
	}
	static void PrintMatr(string msg, Matrix& m, int prec=6){
		cout << endl << msg.c_str() << endl;
		PrintMatr(m,prec);
	};
	static void PrintCV(string msg, ColumnVector& cv, int prec){
		cout << msg.c_str() << " ";
		for(int i=0; i<cv.nrows(); i++)
			cout << fixed << setprecision(prec) << cv(i+1) << " ";
		cout << endl;
	}
	static void MPrint(double* data, int nCols, int n, int prec=6)
	{
		for(int i=0; i<n; i++)
		{
			cout << setprecision(prec) << data[i] << " ";
			if((i>0 || nCols==1) && !(i%nCols) && i!=(n-1)) cout << "..." << endl;
		}
		cout << endl;
	}
	void		GetRootShort(char* root)
	{
		char* ptr1=strrchr((char*)m_root,'/'),
			*ptr2=strrchr((char*)m_root,'\\');
		if(ptr1<ptr2) ptr1=ptr2;
		if(!ptr1) ptr1=(char*)m_root;
		ptr2=strstr(ptr1,".4dfp");
		if(!ptr2) ptr2=strstr(ptr1,".4dint");
		if(!ptr2) sprintf (root, "%s", m_root);
		else
		{
			int len=(int)(ptr2-m_root);
			strncpy(root,m_root,len);
			root[len]=0;
		}
	}
	static void GetRoot(const char* name, char* root)
	{
		char	*str;
		strcpy (root, name);
		if (str = strrchr (root, '.')) 
			*str='\0';
	};

	static char* GetFile(char* name)
	{
		char* str;
		if(str=strrchr(name, '\\'))
			return ++str;
		else if (str=strrchr(name,'/'))
			return ++str;
		else 
			return name;
	}

	static void GetDir(const char* name, char* dir)
	{
		char* str;
		strcpy(dir,name);
		if(str=strrchr(dir, '\\'))
			*str='\0';
		else if (str=strrchr(dir,'/'))
			*str='\0';
		else
			dir[0]=0;
	}
	bool ProcessCommandLine(int argc, char* argv[])
	{
		m_argc=argc; m_argv=argv;
		char* ptr, str[MAXL],c;
		if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
		strcpy (m_program, ptr);

		int k,i;
		for (k = 0, i = 1; i < argc; i++) 
		{
			if (*(argv[i]) == '-')
			{
				strcpy(str,argv[i]); 
				ptr=str;
				while(c=*(ptr++))
					switch(c)
					{
						case '@': m_control= *ptr++; *ptr='\0'; break;
						case 'a': m_bAHeader=true; break;
					}

			}
			else if(k==0) 
			{
				GetRoot (argv[i], m_root); k++;
			}
		}
		if(k<1) return false;

		return true;
	};

};
#endif //CONSOLECONFIG_INCLUDED
