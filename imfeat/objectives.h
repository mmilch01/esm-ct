/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#include <mlib3.h>

Real norm3(ColumnVector c)
{
	return sqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3));
}

ColumnVector V4D(ColumnVector& v3d){ ColumnVector res(4); res.rows(1,3)=v3d.rows(1,3); res(4)=1; return res; }
ColumnVector V3D(ColumnVector& v){ ColumnVector res(3); res.rows(1,3)=v.rows(1,3); return res; }

void polyeval(Matrix& coefs, const Real t, ColumnVector& out)
{
	out.ReSize(3);
	Real acc;
	for(int i=1; i<=3; i++) {
		acc=1; out(i)=0;
		for (int j=1; j<=coefs.Nrows(); j++){
			out(i)+=acc*coefs(j,i); acc*=t;
		}
	}
}

void diff_curve(Matrix& c, Matrix& diff)
{
		diff.resize(c.Nrows(),c.Ncols());
		for(int i=1; i<=3;i++){
			for(int j=1; j<c.Nrows(); j++){
				diff(j,i)=c(j+1,i)*j;
			}
			diff(c.Nrows(),i)=0;
		}
//		cout << "diff_curve: " << setprecision(4)<<fixed<<endl << diff;
}

Real mean_curvature(Matrix& coefs,DiagonalMatrix& V2MM)
{
	Real t,acc=0;
	Matrix diff,diff2;
	diff_curve(coefs,diff);
	diff_curve(diff,diff2);
	ColumnVector dX(3),ddX(3);
	Real sqd,sqdd, dp;
	DiagonalMatrix V2MM1(3); 
	V2MM1(1)=V2MM(1);V2MM1(2)=V2MM(2);V2MM1(3)=V2MM(3);
	for(int i=0; i<=1000; i++){
		t=(Real)i/1000.0;
		polyeval(diff,t,dX); dX=V2MM1*dX;
		polyeval(diff2,t,ddX);ddX=V2MM1*ddX;
		sqd=dX.SumSquare(); sqdd=ddX.SumSquare();
		dp=DotProduct(dX,ddX);
		if(sqd>1e-6)
			acc+=fabs(sqd*sqdd-dp*dp)/(sqd*sqd*sqd);
	}
	return sqrt(acc/1001.0);
}

void V2CV(ColumnVector& out, vector<Real> v)
{
	out.ReSize(v.size());
	for(int i=0; i<v.size(); i++) out(i+1)=v[i];	
}

struct LengthIntegrand{
	Matrix coefs;
	LengthIntegrand(Matrix& c,DiagonalMatrix& V2MM)
	{	
		coefs.resize(c.Nrows()-1,c.Ncols());
		for(int i=1; i<=3;i++)
			for(int j=1; j<=coefs.Nrows(); j++){
//				coefs(i,j)=V2MM(i)*c(i,j+1)*j;
				coefs(j,i)=V2MM(i)*c(j+1,i)*j;
			}
		cout << "Length integrand: " << endl << setprecision(4) << coefs;
	};
	Real operator()(Real t){
		ColumnVector cv;
		polyeval(coefs,t,cv);
		return norm3(cv);
	};
};
/*
struct Length{
	LengthIntegrand *li;
	Real tlen;
	Length(LengthIntegrand &lint, Real target){li=&lint;tlen=target;};
	Real operator()(ColumnVector t)
	{
		Real res=qsimp(*li,t(1),1,1e-6);
		return -fabs(res-tlen);		
	};
};
*/
//значения базисных функций (квадратичных полиномов) в точке x. / values of basis quadratic polinomials in x
ColumnVector quadfit(const Real x){
	ColumnVector ans(4);
	ans(1)=1; ans(2)=x; ans(3)=x*x; ans(4)=x*x*x;
//	for (int i=2; i<=3; i++) ans(i)=ans(i-1)*x;		
	return ans;
}

//дисковый интеграл с центром смещённым относительно оси AB. / disk integral with center shifted relative to AB axis
struct Tos
{
	Real R; //радиус в мм / radius in mm
	Volume* pV;
	ColumnVector A, B, AB, Amm, Bmm, ABmm, Avox, Bvox;
	Real sigma; //параметр гауссиана / Gaussian parameter
	Real t; //параметр вдоль отрезка / parameter along segment
	DiagonalMatrix V2MM,MM2V;
	Tos(Volume* v, ColumnVector& PA, ColumnVector& PB, Real r, Real tt)
	{
		pV=v; A=PA.Rows(1,3); B=PB.Rows(1,3);
		AB=B-A;
		R=r;
		sigma=2*R; 
		v->Metric_Vox(MM2V,V2MM);
		Amm=V2MM*V4D(A);	Bmm=V2MM*V4D(B); ABmm=V2MM*V4D(AB);
		Avox=V4D(A); Bvox=V4D(B);
		t=tt;
	};
	void GetTransform(ColumnVector D, Matrix& T)
	{
		Matrix AB2E, TD;
		MMath::IdentityTransform(TD);
		TD(2,4)=D(1);
		TD(3,4)=D(2);
		//Т: AB -> (1,0,0). в вокс. координатах /in voxel coords
		MMath::VectorToE1(Amm,Bmm,T);		
		AB2E=TD*T*V2MM;
		//точка A теперь в (0,0,0), точка P- в (t,0,0), нормальная плоскость - плоскость YZ / now A is at (0,0,0), P in (t,0,0), normal plane is YZ
		T=AB2E.i();
	};
	ColumnVector DiskDisplacementToGlobalPt(ColumnVector d)
	{
		ColumnVector P(4),D(4);
		Matrix E2AB;
		GetTransform(d,E2AB);
		P(1)=t; P(2)=0; P(3)=0; P(4)=1;
		return E2AB*P;
	};
	//найти параметр для заданной длины в мм, отсчитываемой от конца отрезка. / find parameter for given length in mm, counted from the segment's end
	/*
	Real ParamFromLength(Real l, Matrix& params)
	{
		LengthIntegrand li(params,V2MM);
		Length len(li,l);
		Amoeba am(0.01,1000);
		ColumnVector t(1), res(1); t(1)=1-l/norm3(ABmm);
		res=am.minimize(t,1,len);
		return res(1);
	};
	*/
	//Смещение задаётся параметрической кривой, t от 0 до 1. / shift curve is parameterized by t from 0 to 1
	Real ValFromNonlinearCurve(Real tt, Matrix& params)
	{
		ColumnVector Pt(4),P0t(4);
		polyeval(params,tt,Pt);
		Pt=V4D(Pt);
		Real t0=DotProduct(Pt-Avox,Bvox-Avox)/(norm3(AB)*norm3(AB));
		P0t=V2MM*(Avox+(Bvox-Avox)*t0); P0t(4)=1;
		Matrix T,AB2E;
		MMath::VectorToE1(Amm,Bmm,T);
		ColumnVector D=T*(P0t-V2MM*Pt);
		t=t0;
		D(1)=D(2); D(2)=D(3);
		return -(*this)(D);
	};
	//D - малый двумерный вектор смещения в мм / small 2D shift vector in mm
	Real operator() (ColumnVector D)
	{

		Matrix E2AB;
		GetTransform(D,E2AB);


		Real dy=0.1,dz=0.1;
		Real R2=R*R;
//		cout << setprecision(4) << AB2E << endl;
//		cout << setprecision(4) << E2AB << endl;

		//точка на диске со смещённым отн. оси центром, в глобальных координатах / point on disk shifted relative to axis, in global coords
		ColumnVector Pg(4); 

		//та же точка в локальных координатах AB=[1,0,0] / same point in local coords
		ColumnVector Pl(4); Pl(4)=1; Pl(1)=t;
		Real acc=0,n=0;		
		Real s2=1.0/(sigma*sigma), l;

		for(Real z=-R; z<=R; z+=dz){
			Pl(3)=z;//+D(2);
			for(Real y=-R; y<=R; y+=dy){
				if ((l=z*z+y*y)>R2) continue;
				Pl(2)=y;//+D(1);
//				cout <<setprecision(4)<<Pl.t();
				Pg=E2AB*Pl;
//				cout << setprecision(4) << Pg.t();
				if (!pV->IsInsideBoundary(Pg)) continue;
				acc+=exp(-l*s2)*pV->Interp3(Pg);
				n++;
//debug
//				(*pV)(Pg)=0;
			}
		}
//		pV->Write4dfp("C:\\Temp\\elpos\\disk",g_argc,g_argv);

		return (n==0)? 0:-acc/n;
	}
};


//функтор, вычисляющий сферический/дисковый интеграл в точке. / functor that calculates spherical/disk integral at point
struct To
{
	Real R; //радиус в мм
	Volume* pV;
	ColumnVector A, B, AB, Amm, Bmm, ABmm, Avox, Bvox;
	Real sigma;
	DiagonalMatrix V2MM,MM2V;
	To(Volume *v, ColumnVector& PA, ColumnVector& PB, Real r)
	{
		pV=v; A=PA.rows(1,3); B=PB.rows(1,3);
		AB=B-A;
		R=r;
		sigma=2*R;
		v->Metric_Vox(MM2V,V2MM);
		Avox=V4D(A); Bvox=V4D(B);
		Amm=V2MM*Avox; Bmm=V2MM*Bvox;
		ABmm=V2MM*(Bvox-Avox);
	};
	//дисковый интеграл / disk integral
	Real operator() (Real t)
	{
		//перевод в мм диапазон / change to mm scale
		Real dy=0.1,dz=0.1;
		Real R2=R*R;
		Matrix AB2E, E2AB, T;

//		cout << setprecision (6) << "Avox: " << Avox.t() <<", Bvox: " << Bvox.t() << endl;
//		cout << setprecision (3) << V2MM << MM2V << endl;
		//Т: AB -> (1,0,0). в вокс. координатах / in voxel coords

		MMath::VectorToE1(Amm,Bmm,T);
		AB2E=T*V2MM;

		//точка A теперь в (0,0,0), точка P- в (t,0,0), нормальная плоскость - плоскость YZ / now A is at (0,0,0), P in (t,0,0), normal plane is YZ
		E2AB=AB2E.i();

//		cout << setprecision(4) << AB2E << endl;
//		cout << setprecision(4) << E2AB << endl;

		//точка на диске в глобальных координатах / point on disk in global coords
		ColumnVector Pg(4); 

		//та же точка в локальных координатах, параметризованная вдоль оси X / same point in local coords, parameterized along X axis
		ColumnVector Pl(4); Pl(4)=1; Pl(1)=t;
		Real acc=0,n=0;		
		Real s2=1.0/(sigma*sigma), l;

		for(Real z=-R; z<=R; z+=dz){
			Pl(3)=z;
			for(Real y=-R; y<=R; y+=dy){
				if ((l=z*z+y*y)>R2) continue;
				Pl(2)=y;
//				cout <<setprecision(4)<<Pl.t();
				Pg=E2AB*Pl;
//				cout << setprecision(4) << Pg.t();
				if (!pV->IsInsideBoundary(Pg)) continue;
				acc+=exp(-l*s2)*pV->Interp3(Pg);
				n++;
//debug
//				(*pV)(Pg)=0;
			}
		}
//		pV->Write4dfp("C:\\Temp\\elpos\\disk",g_argc,g_argv);

		return (n==0)? 0:acc/n;
	};
	//сферический интеграл / spherical integral
	Real operator_sphere(Real t)
	{
		ColumnVector P=A+AB*t;
		Real wx=R/pV->DX(),wy=R/pV->DY(),wz=R/pV->DZ(); //размеры в вокселях. / size in voxels
		Range3 r(-wx,+wx,-wy,+wy,-wz,+wz);
		Real dx=wx/ceil(wx),dy=wy/ceil(wy),dz=wz/ceil(wz); //шаг в вокс. / step in voxels
//		Real xr,yr,zr;
		Real sx,sz,sy;
		Real R2=R*R;
		Real n=0,sum=0;
		ColumnVector Poff(3);
		for (Real z=0,zr=r.r[2].st; z<=2*R; z++,zr+=dz){
			sz=zr*pV->DZ();sz*=sz;
			Poff(3)=zr+P(3);
			for(Real y=0,yr=r.r[1].st; y<=2*R; y++,yr+=dy){
				sy=yr*pV->DY();sy*=sy; sy+=sz;
				Poff(2)=yr+P(2);
				for(Real x=0,xr=r.r[0].st; x<=2*R; x++,xr+=dx)
				{
					sx=xr*pV->DX();sx*=sx;
					if(sy+sx>R2) continue;
					Poff(1)=xr+P(1);
					if (!pV->IsInsideBoundary(Poff)) continue;
					sum+=pV->Interp3(Poff);
					n++;
				}
			}
		}
		return sum/n;
	};
};


//функтор, вычисляющий интеграл вдоль оси (смещения в перпендикулярных оси плоскостях) 
//functor that calculates the integral along an axis (shifts are given in planes perpendicular to the axis)
struct Io
{
#ifndef _i
#define _i(x) (int)((x)+.5)
#endif
	ColumnVector Vab, A,B;
	Real Da,Db;
	Real R;
	int ax,ay,az; //оси.
	Volume* pV;
	Real sigma;
	Io(Volume* v, vector<Real>& coords,Real r){
		A.ReSize(3);B.ReSize(3);
		A<<coords[0]<<coords[1]<<coords[2];
		B<<coords[3]<<coords[4]<<coords[5];
		Vab=B-A;
		Da=-DotProduct(Vab,A);
		Db=-DotProduct(Vab,B);
		R=r;
		sigma=1*R;
		int amax, amin, amed;
		int ind[3];
		Vab.MaximumAbsoluteValue1(amax);
		Vab.MinimumAbsoluteValue1(amin);
		ind[0]=amax; ind[2]=amin;
		for (int i=0; i<3; i++){
			if (i!=amax && i!=amin){
				ind[1]=i; break;
			}
		}
		ax=ind[0]; ay=ind[1];
		pV=v;
	};
	void GetTips(const ColumnVector& disp, ColumnVector& Af, ColumnVector& Bf)
	{
		Displacement2Pos3(disp(1),disp(2),Af,A,Da);
		Displacement2Pos3(disp(3),disp(4),Bf,B,Db);
	};
	void Displacement2Pos3(Real x, Real y, ColumnVector& X, ColumnVector& X0, Real D)
	{
		X.ReSize(3);
		X(1)=x+X0(1);X(2)=y+X0(2);
		X(3)=-(Vab(1)*X(1)+Vab(2)*X(2)+D)/Vab(3);
	};

	Real operator() (const ColumnVector &vec)
	{
		ColumnVector Xa(3),Xb(3);
		GetTips(vec,Xa,Xb);

		ColumnVector Xab=Xb-Xa;
		Real L=norm3(Xb-Xa);
		Real iL=1.0/L;
		Real L2=L*L;
		Real iL2=1.0/L2;
//		Xab*=iL;
		Real acc=0;
		Real cnt=0;
		Range3 box;
		box.GetBoundingBox(Xa,Xb);
		box.Expand(_i(R),_i(R),_i(R));
		box.Intersect(0,pV->SX()-1,0,pV->SY()-1,0,pV->SZ()-1);
		Real t,syz=0;
		ColumnVector Xd(3),Dl(3);
		Real R2=R*R;
		Real dist;
		Real s2=1.0/(sigma*sigma);
		for (int z=box.r[2].st; z<=box.r[2].en; z++){
			Xd(3)=z-Xa(3);
			for(int y=box.r[1].st; y<=box.r[1].en; y++){
				Xd(2)=y-Xa(2);
				syz=Xd(2)*Xab(2)+Xd(3)*Xab(3);
				for(int x=box.r[0].st; x<=box.r[0].en; x++){
					Xd(1)=x-Xa(1);
					t=(Xd(1)*Xab(1)+syz)*iL2;
					Dl=Xd-Xab*t;
					dist=Dl(1)*Dl(1)+Dl(2)*Dl(2)+Dl(3)*Dl(3);
					if (dist <= R2){
						acc+=(*pV)(x,y,z)*exp(-dist*s2);
						cnt++;
					}
				}
			}
		}
		if (cnt==0) return 0;
//debug
		cout << "io=" << fixed << setprecision(4) <<  -acc/cnt << endl;
		ConsoleConfig::PrintCV("vec:", (ColumnVector&)vec, 2);
		
		return -acc/cnt;
	};
};

//функтор, вычисляющий интеграл вдоль оси (смещение вдоль оси)
//functor that calculates integral along axis (shifts along the axis).

struct Tol
{
	Real R; //радиус в мм / radius in mm
	Volume* pV;
	ColumnVector A, B, AB;
	Real sigma; //параметр гауссиана / Gaussian parameter
	Real t; //параметр вдоль отрезка / parameter along the segment
	Real thresh; //порог интенсивности / intensity threshold
	bool bOptimizeHUdiff;
	Real t_before, t_between, t_after, t_thresh;
	Tol(Volume* v, ColumnVector& PA, ColumnVector& PB, Real r, bool bOHUD):R(r)
	{
		pV=v; A=PA; B=PB;
		AB=B-A;
		sigma=2*R;
		thresh=0;

		//compute values before and after the distal tip.
		bOptimizeHUdiff=false;
		t_before=getval(1.5);
		t_after=getval(-1.5);
		t_between=getval(0);
		t_thresh=.5*(t_between+t_after);
		bOptimizeHUdiff=bOHUD;

		/*
		if (bAutoThresh){
			Real sh, elec;
			ColumnVector D(1); D(1)=0;
			elec=-(*this)(D);
			D(1)=-1.2;
			sh=-(*this)(D);
//			if (sh>0.9*elec) thresh=sh;
			thresh=sh;
			cout << fixed << setprecision(2) << "elec=" << elec << ", " << "sh=" << sh << ", thresh=" << thresh << endl;
		}
		*/
	};
	void GetTransform(ColumnVector D, Matrix& T)
	{
		//перевод в мм диапазон / translate to mm scale
		DiagonalMatrix V2MM(4),MM2V(4); pV->Metric_Vox(MM2V,V2MM);
		Matrix AB2E, E2AB, TD;
		ColumnVector Avox(4), Bvox(4),Dvox(4);
		Avox=V4D(A); Bvox=V4D(B); 
		Dvox(1)=D(1);Dvox(2)=0; Dvox(3)=0; Dvox(4)=1;

//		cout << setprecision (6) << "Avox: " << Avox.t() <<", Bvox: " << Bvox.t() << endl;
		ColumnVector Amm,Bmm,Dmm;
//		cout << setprecision (3) << V2MM << MM2V << endl;
		Amm=V2MM*Avox; Bmm=V2MM*Bvox; 
		//смещение вдоль оси X (эквив. смещению вдоль оси AB) / shift along X axis (equiv to shift along AB axis)
		MMath::IdentityTransform(TD);
		TD(1,4)=D(1);

		//Т: AB -> (1,0,0). в вокс. координатах / in voxel coords
		MMath::VectorToE1(Amm,Bmm,T);
		
		AB2E=TD*T*V2MM;

		//точка A теперь в (0,0,0), точка P- в (t,0,0), нормальная плоскость - плоскость YZ / now A is at (0,0,0), P in (t,0,0), normal plane is YZ
		E2AB=AB2E.i();
		T=E2AB;
	};
	ColumnVector DiskDisplacementToGlobalPt(ColumnVector d)
	{
		ColumnVector P(4),D(4);
		Matrix T;
		GetTransform(d,T);
		D(3)=d(2);D(2)=d(1);D(1)=t;D(4)=1;
		P=T*D;
	};
	Real getval(Real d)
	{
		ColumnVector D(1); D=d; return (*this)(D);
	}

	//D - одномерный вектор смещения в безразмерных единицах. /D is 1D shift vector, dimensionless
	Real operator() (ColumnVector D)
	{
#define _sq(x) (x)*(x)

		Matrix E2AB;

		GetTransform(D,E2AB);

		Real dy=0.1,dz=0.1, dt=0.1;
		Real R2=R*R;
//		cout << setprecision(4) << AB2E << endl;
//		cout << setprecision(4) << E2AB << endl;

		//точка на диске со смещённым отн. оси центром, в глобальных координатах / point on disk shifted relative to the axis, global coords
		ColumnVector Pg(4); 

		//та же точка в локальных координатах AB=[1,0,0] / the same point in local coords
		ColumnVector Pl(4); Pl(4)=1;
		Real acc=0,n=0,val;
		Real s2=1.0/(sigma*sigma), l;

		//интегрируем вдоль оси. / integrate along axis
		for (Real t=0; t<=1+dt*0.1; t+=dt){
			Pl(1)=t;
			for(Real z=-R; z<=R; z+=dz){
				Pl(3)=z;
				for(Real y=-R; y<=R; y+=dy){
					if ((l=z*z+y*y)>R2) continue;
					Pl(2)=y;
//					cout <<setprecision(4)<<Pl.t();
					Pg=E2AB*Pl;
//					cout << setprecision(4) << Pg.t();
					if (!pV->IsInsideBoundary(Pg)) continue;
					if ((val=exp(-l*s2)*pV->Interp3(Pg))<thresh) continue;
					acc+=val;
					n++;
				}
			}
	//		pV->Write4dfp("C:\\Temp\\elpos\\disk",g_argc,g_argv);
		}
		//находим значение в крайних точках. / find value in extreme points
		Real acc1=0,acc2=0,n1=0,n2=0;
		for (Real t=0; t<=1+dt*0.1; t+=1){
			Pl(1)=t;
			for(Real z=-R; z<=R; z+=dz){
				Pl(3)=z;
				for(Real y=-R; y<=R; y+=dy){
					if ((l=z*z+y*y)>R2) continue;
					Pl(2)=y;
	//				cout <<setprecision(4)<<Pl.t();
					Pg=E2AB*Pl;
	//				cout << setprecision(4) << Pg.t();
					if (!pV->IsInsideBoundary(Pg)) continue;
					if ((val=exp(-l*s2)*pV->Interp3(Pg))<thresh) continue;
					if (t==0){ 
						acc1+=val; n1++;
					}
					else	 {
						acc2+=val; n2++; 
					}
				}
			}
		}
		if (bOptimizeHUdiff){
			Real coef=_mx(0,(1-0.001*fabs(acc2/n2+t_thresh)));
	//		Real coef=_mx(0,(1-0.001*fabs(acc1/n1-acc2/n2)));
			//более гладкая функция / smoother function
			coef=-2*coef*coef*coef+3*coef*coef;
//			return -coef;
			return (n==0)? 0:-coef*acc/n;
		}
		else
			return (n==0)?0:-acc/n;
	}
};
