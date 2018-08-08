/*******************************************************************************
Copyright (c) 2016, 2017, 2018
Authored by: Mikhail Milchenko, mmilchenko@wustl.edu
This software is free for non-commercial use by acadimic, government, and non-profit/not-for-profit institutions. If this code is modified and/or redistributed, this note must appear in the beginning of each of re-distributed source code files. A commercial version of the software is available, and is licensed through the Office of Technology Management at Washington University School of Medicine. For more information regarding the commercial license, please email your request to otm@dom.wustl.edu.

BY DOWNLOADING THIS SOFTWARE YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED 'AS IS', WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
**********************************************************************************/

#include "mlib3math.h"
#include "mlib3volume.h"
#include <fstream>
#include <iostream>
#include <iomanip>


//#include <limits>
#include "newm/newmatio.h"
#include "mlib3consoleutil.h"

//////////////////////////////////////////////////////////////////////
//
//	Math class - includes basic math functions
//
//////////////////////////////////////////////////////////////////////
const long_Real MMath::UndefDouble = -1e-36;
const long_Real MMath::PI = 3.1415926535897932384626433832795;
const long_Real MMath::RAD = 180.0/MMath::PI;

const int MMath::PLANE_YZ=0;
const int MMath::PLANE_XZ=1;
const int MMath::PLANE_XY=2;

const int MMath::COMPARE_CORR=0;
const int MMath::COMPARE_NMI=1;
const int MMath::COMPARE_ABSDIFF=2;

void MMath::GaussCoeff(Real* buf, int n)
{
	int r = n + 1, c = 1;
	buf[0] = 1;
	for (int i = 1; i<n; i++)
	{
		buf[i] = (buf[i - 1] * (r - c)) / c; c++;
	}
	buf[n] = 1;
	Real cf = pow(2.0, n);
	for (int i = 0; i <= n; i++) buf[i] = buf[i] / cf;
}
void MMath::VectorToE1(ColumnVector& A, ColumnVector& B, Matrix& T)
{
	ColumnVector v = B - A;

	Real n = v.NormFrobenius();
	Real nxy = sqrt(v(1)*v(1) + v(2)*v(2));
	Real nxz = sqrt(v(1)*v(1) + v(3)*v(3));
	Real x = v(1), y = v(2), z = v(3);
	Real cz, sz;
	if (nxy>0) { cz = x / nxy; sz = y / nxy; }
	else { cz = 1; sz = 0; }
	Real cy = nxy / n;
	Real sy = z / n;
	Matrix M1, M2, M3;
	IdentityTransform(M1); M1(1, 1) = cz; M1(1, 2) = sz; M1(2, 1) = -sz; M1(2, 2) = cz;
	IdentityTransform(M2); M2(1, 1) = cy; M2(1, 3) = sy; M2(3, 1) = -sy; M2(3, 3) = cy;
	IdentityTransform(M3); M3(1, 1) = 1 / n;

	//Translate to origin first.
	ColumnVector Aneg = -A;
	Matrix M0; TranslateMatrix(Aneg, M0);

	T = M3*M2*M1*M0;
}

void MMath::TranslateMatrix(ColumnVector& v, Matrix& T) {
	IdentityTransform(T);
	T(1, 4) = v(1); T(2, 4) = v(2); T(3, 4) = v(3);
}

void MMath::ConvolveSep3D_sm(Volume& v0, Volume& v1, Real* krn, int half_krn_sz, Volume* mask)
{
	v1 = v0; v1 = 0;
	ConvolveVolume1D(v0, v1, 1, krn, half_krn_sz, mask);
	Volume v2; v2.InitMemory(v0.m_dims);
	ConvolveVolume1D(v1, v2, 2, krn, half_krn_sz, mask);
	if (v1.SZ()>1)
		ConvolveVolume1D(v2, v1, 3, krn, half_krn_sz, mask);
}
void MMath::Convolve1D(ColumnVector& v, ColumnVector& res, Real* krn, int half_krn_sz)
{
	int sz = v.size();
	res.resize(sz);

	Real sum = 0;
	//	bool bMaskStop;
	//main part
	for (int x = half_krn_sz; x<sz - half_krn_sz; x++)
	{
		sum = 0;
		for (int i = x - half_krn_sz, j = 0; i <= x + half_krn_sz; i++, j++) sum += v(i + 1)*krn[j];
		res(x + 1) = sum;
	}
	//left
	for (int x = 0; x<half_krn_sz; x++)
	{
		sum = 0;
		for (int i = x - half_krn_sz, j = 0; i <= x + half_krn_sz; i++, j++) sum += (i<0) ? v(1)*krn[j] : v(i + 1)*krn[j];
		res(x + 1) = sum;
	}
	//right
	for (int x = sz - half_krn_sz; x<sz; x++)
	{
		sum = 0;
		for (int i = x - half_krn_sz, j = 0; i <= x + half_krn_sz; i++, j++) sum += (i >= sz) ? v(sz)*krn[j] : v(i + 1)*krn[j];
		res(x + 1) = sum;
	}
}

void MMath::ConvolveVolume1D(Volume& v0, Volume& v1, int dim, Real* krn, int half_krn_sz, Volume* mask)
{
	int i1, i2, i3, d1, d2, d3;
	v1.InitMemory(v0.m_dims);
	if (dim == 1) { i1 = 2; i2 = 1; i3 = 0; }
	else if (dim == 2) { i1 = 2; i2 = 0; i3 = 1; }
	else if (dim == 3) { i1 = 1; i2 = 0; i3 = 2; }
	d1 = v0.m_dims[i1]; d2 = v0.m_dims[i2]; d3 = v0.m_dims[i3];
	int x[3];
	ColumnVector r0, r1, mask1D;
	r0.resize(d3); r1.resize(d3); mask1D.resize(d3);
	x[i3] = 0;
	for (x[i1] = 0; x[i1]<d1; x[i1]++)
		for (x[i2] = 0; x[i2]<d2; x[i2]++)
		{
			v0.Get1D(r0, x[0], x[1], x[2], dim);
			if (mask)
			{
				mask->Get1D(mask1D, x[0], x[1], x[2], dim);
				ConvolveBinom1D_3(r0, &mask1D);
			}
			else Convolve1D(r0, r1, krn, half_krn_sz);

			v1.Set1D((mask) ? r0 : r1, x[0], x[1], x[2], dim);
		}
}
void MMath::ConvolveBinom1D_3(ColumnVector& v, ColumnVector* mask)
{
	int sz = v.size();

	Real sum = 0, v1, v2, v3, m1, m2, m3;
	//main part
	v1 = 0; v2 = v(1); v3 = v2;
	m1 = 0; m2 = (*mask)(1); m3 = m2;
	for (int x = 1; x<sz - 1; x++)
	{
		if ((*mask)(x + 1) == 0) { m1 = 0; continue; }
		m2 = (*mask)(x + 1); m3 = (*mask)(x + 2);
		v2 = v(x + 1); v3 = v(x + 2);
		if (m1)
		{
			if (m3) v(x + 1) = (v1 + v2 + v2 + v3) / 4;
			else   v(x + 1) = (v1 + v2 + v2) / 3;
		}
		else if (m3) v(x + 1) = (v2 + v2 + v3) / 3;
		v1 = v2; m1 = m2;
	}
	if ((*mask)(sz) && m1) v(sz) = (v1 + 2 * v(sz)) / 3;
}