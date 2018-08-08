/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#ifndef _WIN32
#if !defined(HAVE_CONFIG_H)
#define HAVE_CONFIG_H
#endif
#endif

#include "objectives.h"

#define _VERBOSE
#define _DBG

static char g_vol_in[MAXL]; //input image
static char g_vol_out[MAXL]; //output image
static bool g_b4dfp=false;
static vector <Real> g_coords; //initial shaft coordinates
static Real g_DiskR=2; //radius of integration disk in mm.
static bool g_bOptimizeHUdiff=false; //flag to optimize for HU difference in the proximal end.

static Real g_ElecDistance=7.5; //distance covered by all electrodes. 10.5 for model 3387, 7.5 for model 3389.
static Real g_ElecSpacing=0.5;
static Real g_ElecLength=1.5;

static int g_iGaussKrnSz=2; //Gaussian kernel pre-blur size


enum COMMAND_TYPE {POS_ELECTRODE_CT};

static COMMAND_TYPE g_cmdType;
static int g_argc;
static char **g_argv;



bool GetArgs(int argc, char* argv[])
{
	AnyOption opt;

    /* 3. SET THE USAGE/HELP   */
	opt.addUsage( "Detect image feature" );
    opt.addUsage( "\nUsage: imfeat <command> [options] [arguments]" );

	opt.addUsage( "\nCommands:" );
	opt.addUsage( "\tOutput DBS electrode position for a CT image\r ep" );
	opt.addUsage( "\timfeat ep [options] <x0> <y0> <z0> <x1> <y1> <z1>" );
	opt.addUsage( "\tx0..z1: proximal (0) and distal (1) implant tips in image space" );
	opt.addUsage( "\trequired arguments: -in, -out" );
	opt.addUsage( "\toptions: -4, -h, -d, -g, -r" );


	opt.addUsage( "\nArguments and options:" );
	opt.addUsage( "\t\t\tinput image\r -in <image>" );
	opt.addUsage( "\t\t\toutput image\r -out <image>" );
	opt.addUsage( "\t\t\tassume 4dfp and t4 format ([Analyze 7.5 and FSL])\r -4" );
	opt.addUsage( "\t\t\toptimize for HU difference in distal end (clinical images)\r -h" );
	opt.addUsage( "\t\t\telectrode distance in mm [7.5]\r -d <float>" );
	opt.addUsage( "\t\t\telectrode spacing in mm [0.5]\r -s <float>" );

	opt.addUsage( "\t\t\tgaussian preblur kernel size in voxels (0 for no blur) [2]\r -g <int>" );
	opt.addUsage( "\t\t\tradius of integration disk in mm [2]\r -r <float>" );

	opt.noPOSIX();
	opt.setOption("in");
	opt.setOption("out");
	opt.setFlag("4");
	opt.setFlag("h");

	opt.setOption("d");
	opt.setOption("s");

	opt.setOption("g");
	opt.setOption("r");

	string arg;
	opt.processCommandArgs(argc,argv);
	g_b4dfp=opt.getFlag("4");
	g_bOptimizeHUdiff=opt.getFlag("h");

	if (argc<2) 
	{
		opt.printUsage(); return false;
	}
	//split a mask

	if (!strcmp(argv[1],"ep")) 
		g_cmdType=POS_ELECTRODE_CT;
	else{
		opt.printAutoUsage(); return false;
	}

	if ( g_cmdType==POS_ELECTRODE_CT )
	{
		if(argc<7){
			opt.printUsage(); return false;
		}
		Real l;
		for (int i=argc-6; i<argc; i++)
		{
			l=atof(argv[i]);
			if ( l==0 && argv[i][0] != 0 ){
				opt.printUsage(); return false;
			}			
			g_coords.push_back(l);
		}
	}
	else{
		opt.printUsage(); return false;
	}

	char* s=opt.getValue("in");
	if (s)	strcpy(g_vol_in,s);
	else {opt.printUsage(); return false;}
	s=opt.getValue("out");
	if (s)	strcpy(g_vol_out,s);
	else {opt.printUsage(); return false;}
	s=opt.getValue("d"); if (s)	g_ElecDistance=atof(s);
	s=opt.getValue("s"); if (s)	g_ElecSpacing=atof(s);

	s=opt.getValue("r"); if (s)	g_DiskR=atof(s);
	s=opt.getValue("g"); if (s)	g_iGaussKrnSz=atoi(s);

	return true;
}

//utility functions
template<class T>
void save_var1(ofstream& ofs, char* name, T val,int prec=4)
{
	ofs << fixed << setprecision(4) << name << "=" << val << endl;
}
void save_var(ofstream& ofs, char* name, ColumnVector& arr)
{
	ofs << fixed << setprecision(4) << name << "=(";
	for (int i=0; i<arr.size(); i++) ofs << arr(i+1) << " ";
	ofs << ")" << endl;
}

void cpcoord(ColumnVector& A, ColumnVector& B, vector<Real>& coords)
{
	for(int i=1; i<=6;i++){
		if(i<4)	A(i) = coords[i-1];
		else	B(i-3) = coords[i-1];
	}
}

void cpcoord(ColumnVector& A, ColumnVector& B, ColumnVector& coords)
{
	for(int i=1; i<=6;i++){
		if(i<4)	A(i) = coords(i);
		else	B(i-3) = coords(i);
	}
}


bool OptimizeShaftPosition(Volume& vin, Volume& vout, ColumnVector& A, ColumnVector& B)
{
	vout=vin; vout=0;
	Io io(&vin, g_coords,g_DiskR);
	Amoeba am(0.01,1000);
	ColumnVector C(4),Dels(4);
	C=0;
	Real d=4;
	Dels << d << d << d << d;
	C=am.minimize(C,Dels,io);
	A.ReSize(3); B.ReSize(3);
	io.GetTips(C,A,B);
	ConsoleConfig::PrintCV("optimized A: ",A, 2);
	ConsoleConfig::PrintCV("optimized B: ",B, 2);

	//сгенерировать конечную маску электрода. / generate final electrode mask
	vout.AddLine(A,B,1);
	return true;
}
bool OptimizeElecPos(ColumnVector& coords, string& root, string& suff, ColumnVector& A, ColumnVector& B)
{
	bool bRes=true;
	Volume v_in, v_out;
	bRes &= v_in.Read4dfp((char*)(root.c_str()));
	g_coords.resize(6);
	for(int i=1; i<=6;i++){
		if(i<4)	A(i) = coords(i);
		else	B(i-3) = coords(i);
		g_coords[i-1]=coords(i);
	}
	//debug
//	return true;
	v_out=v_in; v_out=0;
	ColumnVector A0=A, B0=B;
	bRes &= OptimizeShaftPosition(v_in,v_out,A,B);
	v_in=0;
	v_in.AddLine(A0,B0,1);
	v_in.Write4dfp((char*)(root+suff+"_init_el").c_str(),g_argc,g_argv);
	v_out.Write4dfp((char*)(root+suff+"_final_el").c_str(),g_argc,g_argv);
	return bRes;
}

bool OptimizeNonlinearTrajectory(ColumnVector& A, ColumnVector& B, Matrix& coefs, Volume& v, ofstream *ofs)
{
	DiagonalMatrix V2MM(4),MM2V(4); v.Metric_Vox(MM2V,V2MM);

//1. position the optimization points.
	ColumnVector ABmm=V2MM*(V4D(B)-V4D(A)), Amm=V2MM*V4D(A), Bmm=V2MM*V4D(B);
	Real L=norm3(ABmm);
	Real N0=36; //начальное число точек в разбиении / initial # of points in partition
	vector<Real> T; //точки сетки / point on grid
	Real t=0;//0.5/L;
	Real step_mul=4.0; //коэф-т изменения шага сетки / coefficient of grid step
	Real step=g_ElecDistance/(N0*L); //начальный шаг / initial step
	int nsteps=0;
	while (t < 1){
		T.push_back(1-t);
		if (nsteps>=N0){  
			step*=step_mul;
			N0=_mx(2,N0/2);
			nsteps=0;
		}
		t+=step; nsteps++;
	}
//2. Optimize in selected points.
	vector<ColumnVector> D; //векторы смещения на разбиении / shift vectors on partition
	ColumnVector d(2),dg(3); //смещение в перпендик. плоскости / shift in normal plane
	Amoeba am(0.01,1000,false);
	ColumnVector vT,vX,vY,vZ; //координаты оптим. точек на кривой / coordinates of optimized points on curve
	vX.ReSize(T.size()); vY.ReSize(T.size());vZ.ReSize(T.size());

	for (int i=0; i<T.size(); i++)
	{
		Tos tos(&v,A,B,g_DiskR,T[i]);
		d=0;
		d=am.minimize(d,2,tos); 
		dg=tos.DiskDisplacementToGlobalPt(d);
		vX(i+1)=dg(1); vY(i+1)=dg(2); vZ(i+1)=dg(3);
	}
//3. Update the trajectory.	
	V2CV(vT,T);
	Fitlin flx(vT,vX,quadfit); flx.fit();
	Fitlin fly(vT,vY,quadfit); fly.fit();
	Fitlin flz(vT,vZ,quadfit); flz.fit();
	coefs.resize(flx.a.nrows(),3);
	coefs.Columns(1,1)=flx.a;
	coefs.Columns(2,2)=fly.a;
	coefs.Columns(3,3)=flz.a;

	//cout << "LSfit coefs: " << endl << setprecision(6) << coefs << endl;
	bool bRes=true;
	Real maxC=10000;
	if (flx.chisq > maxC || fly.chisq > maxC || flz.chisq > maxC)  bRes=false;
	cout << "chisq: " << flx.chisq << " " << fly.chisq << " " << flz.chisq << endl;
	if (ofs){
		save_var1(*ofs,"chisq_x",flx.chisq);
		save_var1(*ofs,"chisq_y",fly.chisq);
		save_var1(*ofs,"chisq_z",flz.chisq);
	}

	//save output
	polyeval(coefs,0,A);
	polyeval(coefs,1,B);
	cout << "optimized nonlinear trajectory" << endl;
	ConsoleConfig::PrintCV("A:", A, 2);
	ConsoleConfig::PrintCV("B:", B, 2);

	return bRes;
}

//Optimize nonlinear trajectory in all directions.
bool OptimizeNonlinCombined(Volume& v, ColumnVector& A, ColumnVector& B, Matrix& coefs, ofstream* ofs=NULL)
{
	coefs.ReSize(3,3);
	DiagonalMatrix V2MM(4),MM2V(4); v.Metric_Vox(MM2V,V2MM);
	ColumnVector ABmm=V2MM*(V4D(B)-V4D(A));
	Real L=norm3(ABmm);

	//1. Fit quadratic curve.
	cout << "first pass of trajectory fit" << endl;
	ConsoleConfig::PrintCV("A0:",A,2);
	ConsoleConfig::PrintCV("B0:",B,2);
	if (ofs){
		save_var(*ofs,"A0",A);
		save_var(*ofs,"B0",B);
	}

	OptimizeNonlinearTrajectory(A,B,coefs,v,NULL);

	//2. Locate approximate beginning and end of electrodes along the shaft.
	ColumnVector Amm=V2MM*V4D(A), Bmm=V2MM*V4D(B);
	Tos to1(&v,A,B,g_DiskR,0);
	ColumnVector Atip(3), Btip(3), Atipmm, Btipmm, tmp;;
//	Real t=to1.ParamFromLength(g_ElecDistance,coefs);	//точно / exactly
		
	Real t=1-g_ElecDistance/L;                         //приблизительно / roughly
	polyeval(coefs,t,tmp);

	cout << "t=" << t << ", error=" << norm3(Bmm-V2MM*V4D(tmp))-g_ElecDistance << " mm" << endl;


	polyeval(coefs,t,Atip); Atip=V4D(Atip);
	polyeval(coefs,1,Btip); Btip=V4D(Btip);

	ConsoleConfig::PrintCV("A1:",A, 2);
	ConsoleConfig::PrintCV("B1:",B, 2);

	Atipmm=V2MM*Atip; Btipmm=V2MM*Btip;

	//3. Optimize along the axis parallel to shaft.
	cout << endl<<"optimizing contact end along the trajectory" << endl;
	Tol tol(&v, Atip, Btip, g_DiskR,g_bOptimizeHUdiff);
	ColumnVector LD(1),LDo; LD=0;
	Amoeba am(0.00001,1000);
	LDo=am.minimize(LD,0.1,tol);
	cout << "Amoeba dL=" << fixed << LDo(1)*g_ElecDistance << " mm." << endl;
	Real iV, mX=0, maxI=-1;
	for (Real t=-0.1; t<=0.1; t+=0.01){
		LD(1)=t;
		iV=tol(LD);
		if(mX>iV){
			mX=iV; maxI=t;
		}
//		cout << "t=" << t << ", mX=" << fixed << setprecision(3) << iV << endl;
	}
	cout  << "Brute force D=" << fixed << setprecision(4) << maxI << ", " << mX << ", offset=" << maxI*g_ElecDistance << " mm." << endl;	

	//един. вектор в напр. AtipBtip1 / unit vector aligned with
	ColumnVector AtipBtipmm=Btipmm-Atipmm;
	ColumnVector Atip1,Btip1;

	//4. Update the approximate B tip.
	Atip1=MM2V*(Atipmm-AtipBtipmm*LDo(1));
	Btip1=MM2V*(Btipmm-AtipBtipmm*LDo(1));

//	cout << setprecision(2) << "A2: " << Atip1.t();
//	cout << setprecision(2) << "B2: " << Btip1.t();

	//debug
	Matrix T;
	ColumnVector Pl(4);
	Pl<<0<<0<<0<<1;
	tol.GetTransform(LDo,T);
//	cout << setprecision(2) << "A2: " << (T*Pl).t();
	Pl<<1<<0<<0<<1;
//	cout << setprecision(2) << "B2: " << (T*Pl).t();

	Real temp=tol(LDo);
	ConsoleConfig::PrintCV("B2:",Btip1,2);
	To to2(&v,Atip1,Btip1,g_DiskR);
	cout << setprecision(2) << "v(contact 3 start)=" << to2(0) << ", v(contact 0 end)=" << to2(1) << endl;


	//5. Fit the curve one more time.
	cout << endl << "final pass of trajectory fit" << endl;
	B.rows(1,3)=Btip1.rows(1,3);
	bool bRes=OptimizeNonlinearTrajectory(A,B,coefs,v,ofs);

	//save optimized tips in A and B.
//	B=Btip1;
	Amm=V2MM*V4D(A); Bmm=V2MM*V4D(B);
	t=1-g_ElecDistance/norm3(Bmm-Amm); 
	polyeval(coefs,t,tmp);
	Real terr=norm3(Bmm-V2MM*V4D(tmp))-g_ElecDistance;
	cout << "Final t=" << t << ", error=" << terr << " mm" << endl;
	if (terr>0.1) bRes=false;
	if (ofs){
		save_var1(*ofs,"t_c3_start",t);
		save_var1(*ofs,"err_c3_start_mm",terr);
		ColumnVector elst(4),elstmm(4);		
		polyeval(coefs,t,elst);
		elstmm=V2MM*V4D(elst);
		save_var(*ofs,"contact3_start",elst);
		ColumnVector Cmm,C;
		Real e=(g_ElecDistance-1.5*g_ElecLength-g_ElecSpacing)/g_ElecDistance;
		C=elstmm*(1-e)+Bmm*e;
		ColumnVector Cvox;
		Cvox=(MM2V*C).rows(1,3);
		save_var(*ofs,"contact1_center",Cvox);
	}

	Tos to3(&v,A,B,g_DiskR,0);

	ConsoleConfig::PrintCV("proximal_end_final:",A,2);
	ConsoleConfig::PrintCV("contact0_end_final:",B,2);
	if (ofs){
		save_var(*ofs,"proximal_end",A);
		save_var(*ofs,"contact0_end",B);
	}
	Real vst=to3.ValFromNonlinearCurve(t,coefs),ven=to3.ValFromNonlinearCurve(1,coefs);
	cout << setprecision(2) << "v(contact3_start)=" << to3.ValFromNonlinearCurve(t,coefs) << ", v(contact0_end)=" << to3.ValFromNonlinearCurve(1,coefs) << endl;
	Real K=mean_curvature(coefs,to3.V2MM);
	Real l=norm3(to3.ABmm)*.5,maxdiv=0,rr=1e+10;
	cout << "mean curvature=" << setprecision(4) << scientific << K << endl;
	
	if(K!=0){
		rr=1.0/K;
		maxdiv=rr-sqrt(rr*rr-l*l);
		if (maxdiv > 10) bRes = false;
		cout << "mean r=" << setprecision(4) << fixed << rr << endl;
		cout << "maxdev=" << setprecision(4) << fixed << maxdiv << " mm" <<endl;
	}
	if (ofs){
		save_var1(*ofs,"v_contact3_start",vst);
		save_var1(*ofs,"v_contact0_end",ven);
		save_var1(*ofs,"mean_curv_mm",K,18);
		save_var1(*ofs,"maxdev_mm",maxdiv,4);
	}

	return bRes;
}

void write_profile_nonlin(Volume& v_in, string in, string suff, ColumnVector& A, ColumnVector& B, Matrix coefs)
{
#ifndef _i
#define _i(x) (int)(x+.5)
#endif
	Tos to1(&v_in,A,B,g_DiskR,0);
	string csvfile=in;  csvfile.append("_profile_"); csvfile.append(suff);
	string coordfile=in;  coordfile.append("_coords_"); coordfile.append(suff);

	ofstream ofs; 
	ofs.open(csvfile.c_str(),ios_base::out);
	cout << "writing " << csvfile << endl;
	ColumnVector out,out1;
	Real tmm0=norm3(to1.ABmm);
	Matrix Coords(3,1201);

	//найти t, соотв. началу электрода. / find t corresponding to electrode start
	int te=_i(1000.0*(1-g_ElecDistance/norm3(to1.ABmm)));

	//записать дисковый интеграл вдоль параметризованной кривой. / write disk integral along parameterized curve
	for (int j=0; j<=1200; j++)
	{
		Real t=(Real)j/1000.0;
		polyeval(coefs,t,out);
//		cout << norm3(out-(A+to1.AB*t)) << endl;
		Coords.columns(j+1,j+1)=out;

		if (j==te)
			ofs << t << ", " << 0 << endl;
		else if (j==1000)
			ofs << t << ", " << 0 << endl;
		else {
			ofs << t << ", " << to1.ValFromNonlinearCurve(t,coefs) << endl;
		}
	}
	ofs.close();
	ofs.open(coordfile.c_str(),ios_base::out);
	//записать смещения по x,y,z / write shifts by x, y, z
	for (int i=1; i<=3; i++){
		for (int j=1; j<=1200; j++){
			ofs << j-1 << ", " << Coords(i,j) << endl;
		}
	}
	ofs.close();

}

void write_profile(Volume& v_in, string in, string suff, ColumnVector& A, ColumnVector& B)
{
	To to1(&v_in,A,B,g_DiskR);
	string csvfile=in;
	csvfile.append("_profile_");
	csvfile.append(suff);
	ofstream ofs; 
	ofs.open(csvfile.c_str(),ios_base::out);
	cout << "writing " << csvfile << endl;
	Real tmm0=norm3(to1.ABmm);
	int tst=_i(1000.0*(1-g_ElecDistance/norm3(to1.ABmm)));
	for (int j=0; j<=1200; j++)
	{
		Real t=(Real)j/1000.0;
		if (j==tst)
			ofs << t << ", " << 0 << endl;
		//конец массива электродов. / end of electrode array
		else if (j==1000)
			ofs << t << ", " << 0 << endl;
		else {
			ofs << t << ", " << to1(t) << endl;
		}
	}
	ofs.close();
}

bool Test()
{
	return true;
}
int main (int argc, char *argv[])
{
	g_argc=argc; g_argv=argv;
	if (!Test()) exit (1);
	static char rcsid[]="$Id: " __FILE__ ", " __DATE__ ", " __TIME__;
	std::cout << rcsid << endl;	
	if(!GetArgs(argc, argv)) return 0;

	Volume::COORD_FORMAT fmt=(g_b4dfp)?Volume::FDFP:Volume::FSL;
	g_argc=argc; g_argv=argv;

	Volume v_orig, v_in, v_out;
	cout << "reading " << g_vol_in << endl;
	if (! v_orig.Read(g_vol_in)){
			std::cout << "imfeat ERROR: cannot read input image " << g_vol_in << endl; exit(1);
	}
	
	string outfile=g_vol_out;outfile.append("_imfeat.txt");
	ColumnVector A(3),B(3);
	cpcoord(A,B,g_coords);
	if (g_iGaussKrnSz>0){
		cout << "gaussian blur" << endl;
		MMath::Gauss(v_orig,v_in,g_iGaussKrnSz);
	}
	write_profile(v_in,g_vol_out,"init.csv",A,B);
	Matrix coefs;
	ofstream outfs; outfs.open(outfile.c_str(),ios_base::out);

	if (!OptimizeNonlinCombined(v_in,A,B,coefs,&outfs)){
		cout << "imfeat ERROR: cannot optimize shaft direction" << endl;
		outfs.close();
		exit (-1);
	}
	outfs.close();
	write_profile(v_in,g_vol_out,"final_lin.csv",A,B);
	write_profile_nonlin(v_in,g_vol_out,"final_nonlin.csv",A,B,coefs);
	
	Volume v_elec; v_elec=v_orig; v_elec=0;
	v_elec.AddCurve(A,B,coefs,1.0,&polyeval);

	string elecvol=g_vol_out; elecvol.append("_trajectory");
	cout << "writing " << elecvol << endl;
	v_elec.Write((char*)elecvol.c_str(),fmt,argc,argv);

	exit(0);
}