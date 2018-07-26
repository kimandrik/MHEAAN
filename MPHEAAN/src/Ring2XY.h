/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MHEAAN_Ring2XY_H_
#define MHEAAN_Ring2XY_H_

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <complex>
#include <map>
#include <math.h>

#include "MatrixContext.h"
#include "RingMultiplier.h"

using namespace std;
using namespace NTL;

static RR const Pi = ComputePi_RR();

class Ring2XY {
public:

	long logNx;
	long Nx;
	long Mx;
	long logNxh;
	long Nxh;

	long logN;
	long N;
	long Nh;

	long logNy;
	long Ny;
	long My;

	long logQ; ///< log of Q
	long logQQ; ///< log of PQ

	ZZ Q; ///< Q corresponds to the highest modulus
	ZZ QQ; ///< PQ = Q * Q
	ZZ* qvec;

	double sigma; ///< standard deviation for Gaussian distribution
	long h; ///< parameter for HWT distribution

	RingMultiplier multiplier;

	long* gxPows; ///< auxiliary information about rotation group indexes for batch encoding
	long* gyPows; ///< auxiliary information about rotation group indexes for batch encoding

	complex<double>* ksixPows; ///< storing ksi pows for fft calculation
	complex<double>* ksiyPows; ///< storing ksi pows for fft calculation

	complex<double>* dftomegaPows;
	complex<double>* omegaPows;

	map<pair<long, long>, MatrixContext> matrixContext;

	Ring2XY(long logNx, long logQ, double sigma = 3.2, long h = 64);


	//----------------------------------------------------------------------------------
	//   AUXILIARY CONTEXT
	//----------------------------------------------------------------------------------


	void addMatrixContext(long lognx);


	//----------------------------------------------------------------------------------
	//   ENCODING
	//----------------------------------------------------------------------------------


	void arrayBitReverse(complex<double>* vals, const long n);
	void DFTY(complex<double>* vals);
	void IDFTY(complex<double>* vals);

	void EMBX(complex<double>* vals, const long nx);
	void IEMBX(complex<double>* vals, const long nx);

	void EMBY(complex<double>* vals);
	void IEMBY(complex<double>* vals);

	void EMBXY(complex<double>* vals, const long nx);
	void IEMBXY(complex<double>* vals, const long nx);

	void encode(ZZ* mxy, complex<double>* vals, long nx, long logp);
	void encode(ZZ* mxy, double* vals, long nx, long logp);
	void decode(ZZ* mxy, complex<double>*vals, long nx, long logp, long logq);


	//----------------------------------------------------------------------------------
	//   MULTIPLICATION
	//----------------------------------------------------------------------------------


	void multXpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q);
	void multXpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q);

	void multYpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q);
	void multYpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q);

	void mult(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q);
	void multAndEqual(ZZ* a, const ZZ* b, const ZZ& q);

	void square(ZZ* x, const ZZ* a, const ZZ& q);
	void squareAndEqual(ZZ* a, const ZZ& q);


	//----------------------------------------------------------------------------------
	//   OTHER
	//----------------------------------------------------------------------------------


	void mod(ZZ* res, ZZ* p, ZZ& q);
	void modAndEqual(ZZ* p, ZZ& q);

	void negate(ZZ* res, ZZ* p);
	void negateAndEqual(ZZ* p);

	void add(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q);
	void addAndEqual(ZZ* p1, ZZ* p2, ZZ& q);

	void sub(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q);
	void subAndEqual(ZZ* p1, ZZ* p2, ZZ& q);
	void subAndEqual2(ZZ* p1, ZZ* p2, ZZ& q);

	void multByMonomial(ZZ* res, ZZ* p, const long degx, const long degy);
	void multByMonomialAndEqual(ZZ* p, const long degx, const long degy);

	void multByConst(ZZ* res, ZZ* p, const ZZ& cnst, ZZ& q);
	void multByConstAndEqual(ZZ* p, const ZZ& cnst, ZZ& q);


	//----------------------------------------------------------------------------------
	//   SHIFTING
	//----------------------------------------------------------------------------------


	void leftShift(ZZ* res, ZZ* p, const long bits, ZZ& q);
	void leftShiftAndEqual(ZZ* p, const long bits, ZZ& q);
	void doubleAndEqual(ZZ* p, ZZ& q);

	void rightShift(ZZ* res, ZZ* p, const long bits);
	void rightShiftAndEqual(ZZ* p, const long bits);


	//----------------------------------------------------------------------------------
	//   ROTATION & CONJUGATION & TRANSPOSITION
	//----------------------------------------------------------------------------------


	void leftRotate(ZZ* res, ZZ* p, const long rx, const long ry);
	void conjugate(ZZ* res, ZZ* p);
	void transpose(ZZ* res, ZZ* p);


	//----------------------------------------------------------------------------------
	//   SAMPLING
	//----------------------------------------------------------------------------------


	void sampleGauss(ZZ* res);
	void sampleHWT(ZZ* res);
	void sampleZO(ZZ* res);
	void sampleUniform(ZZ* res, const long bits);

};

#endif
