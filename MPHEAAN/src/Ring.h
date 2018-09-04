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
#include "BootContext.h"
#include "RingMultiplier.h"

using namespace std;
using namespace NTL;

static RR Pi = ComputePi_RR();

class Ring {
public:

	long logN0;
	long N0;
	long M0;
	long logN0h;
	long N0h;

	long logN;
	long N;
	long Nh;

	long logN1;
	long N1;
	long M1;

	long logQ; ///< log of Q
	long logQQ; ///< log of PQ

	ZZ Q; ///< Q corresponds to the highest modulus
	ZZ QQ; ///< PQ = Q * Q
	ZZ* qvec;

	double sigma; ///< standard deviation for Gaussian distribution
	long h; ///< parameter for HWT distribution

	RingMultiplier multiplier;

	long* gM0Pows; ///< auxiliary information about rotation group indexes for batch encoding
	long* gM1Pows; ///< auxiliary information about rotation group indexes for batch encoding

	complex<double>* ksiM0Pows; ///< storing ksi pows for fft calculation
	complex<double>* ksiM1Pows;
	complex<double>* ksiN1Pows; ///< storing ksi pows for fft calculation

	complex<double>* dftomegaPows;
	complex<double>* omegaPows;

	map<pair<long, long>, BootContext> bootContextMap;

	map<pair<long, long>, MatrixContext> matrixContext;

	Ring(long logN0, long logQ, double sigma = 3.2, long h = 64);


	//----------------------------------------------------------------------------------
	//   AUXILIARY CONTEXT
	//----------------------------------------------------------------------------------


	void addBootContext(long lognx, long logny, long logp);

	void addMatrixContext(long lognx);


	//----------------------------------------------------------------------------------
	//   ENCODING
	//----------------------------------------------------------------------------------


	void arrayBitReverse(complex<double>* vals, long n);
	void DFTX1(complex<double>* vals);
	void IDFTX1(complex<double>* vals);

	void EMBX0(complex<double>* vals, long nx);
	void IEMBX0(complex<double>* vals, long nx);

	void EMBX1(complex<double>* vals);
	void IEMBX1(complex<double>* vals);

	void EMB(complex<double>* vals, long nx);
	void IEMB(complex<double>* vals, long nx);

	void encode(ZZ* mxy, complex<double>* vals, long nx, long ny, long logp);
	void encode(ZZ* mxy, double* vals, long nx, long ny, long logp);
	void decode(ZZ* mxy, complex<double>*vals, long nx, long ny, long logp, long logq);


	//----------------------------------------------------------------------------------
	//   MULTIPLICATION
	//----------------------------------------------------------------------------------

	long MaxBits(const ZZ* f, long n);

	uint64_t* toNTTX0(ZZ* a, long np);

	uint64_t* toNTTX1(ZZ* a, long np);
	uint64_t* toNTTX1Lazy(ZZ* a, long np);

	uint64_t* toNTT(ZZ* a, long np);
	uint64_t* toNTTLazy(ZZ* a, long np);

	uint64_t* addNTT(uint64_t* ra, uint64_t* rb, long np);

	void multX0(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	void multX0AndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTTX0(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX0AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX0D(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q);

	void multX1(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	void multX1AndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTTX1(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX1AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multDNTTX1(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q);

	void mult(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	void multAndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTT(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTLazy(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTLazyAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multDNTT(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q);

	void square(ZZ* x, ZZ* a, long np, ZZ& q);
	void squareAndEqual(ZZ* a, long np, ZZ& q);
	void squareNTT(ZZ* x, uint64_t* ra, long np, ZZ& q);


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

	ZZ* multByMonomial(ZZ* p, long degx, long degy, ZZ& q);
	void multByMonomialAndEqual(ZZ* p, long degx, long degy, ZZ& q);

	void multByConst(ZZ* res, ZZ* p, ZZ& cnst, ZZ& q);
	void multByConstAndEqual(ZZ* p, ZZ& cnst, ZZ& q);


	//----------------------------------------------------------------------------------
	//   SHIFTING
	//----------------------------------------------------------------------------------


	void leftShift(ZZ* res, ZZ* p, long bits, ZZ& q);
	void leftShiftAndEqual(ZZ* p, long bits, ZZ& q);
	void doubleAndEqual(ZZ* p, ZZ& q);

	void rightShift(ZZ* res, ZZ* p, long bits);
	void rightShiftAndEqual(ZZ* p, long bits);


	//----------------------------------------------------------------------------------
	//   ROTATION & CONJUGATION & TRANSPOSITION
	//----------------------------------------------------------------------------------


	ZZ* leftRotate(ZZ* p, long rx, long ry);
	ZZ* conjugate(ZZ* p);


	//----------------------------------------------------------------------------------
	//   SAMPLING
	//----------------------------------------------------------------------------------


	void sampleGauss(ZZ* res);
	void sampleHWT(ZZ* res);
	void sampleZO(ZZ* res);
	void sampleUniform(ZZ* res, long bits);

};

#endif
