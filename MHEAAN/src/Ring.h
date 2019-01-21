/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MHEAAN_Ring2XY_H_
#define MHEAAN_Ring2XY_H_

#include "Params.h"

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <complex>
#include <map>
#include <math.h>
#include <vector>

#include "BootContext.h"
#include "RingMultiplier.h"
#include "SqrMatContext.h"


using namespace std;
using namespace NTL;

static RR Pi = ComputePi_RR();

class Ring {
public:
	RingMultiplier multiplier;

	ZZ qvec[logQQ + 1];

	uint64_t gM0Pows[N0h + 1]; ///< auxiliary information about rotation group indexes for batch encoding
	uint64_t gM1Pows[M1]; ///< auxiliary information about rotation group indexes for batch encoding

	complex<double> ksiM0Pows[M0 + 1]; ///< storing ksi pows for fft calculation
	complex<double> ksiM1Pows[M1 + 1];
	complex<double> ksiN1Pows[N1 + 1]; ///< storing ksi pows for fft calculation

	complex<double>* dftM1Pows[logN1 + 1];
	complex<double>* dftM1NTTPows[logN1 + 1];

	map<pair<long, long>, BootContext&> bootContextMap;
	map<long, SqrMatContext&> sqrMatContextMap;

	Ring();

	//----------------------------------------------------------------------------------
	//   ENCODING
	//----------------------------------------------------------------------------------


	void arrayBitReverse(complex<double>* vals, long n);
	void DFTX1(complex<double>* vals, long n1);
	void IDFTX1(complex<double>* vals, long n1);

	void EMBX0(complex<double>* vals, long n0);
	void IEMBX0(complex<double>* vals, long n0);

	void EMBX1(complex<double>* vals, long n1);
	void IEMBX1(complex<double>* vals, long n1);

	void EMB(complex<double>* vals, long n0, long n1);
	void IEMB(complex<double>* vals, long n0, long n1);

	void encode(ZZ* mx, complex<double>* vals, long n0, long n1, long logp);
	void encode(ZZ* mx, double* vals, long n0, long n1, long logp);

	complex<double>* decode(ZZ* mx, long n0, long n1, long logp);


	//----------------------------------------------------------------------------------
	//   MULTIPLICATION
	//----------------------------------------------------------------------------------

	long MaxBits(ZZ* f, long n);
	void addNTTAndEqual(uint64_t* ra, uint64_t* rb, long np);

	void toNTTX0(uint64_t* ra, ZZ* a, long np);
	void toNTTX1(uint64_t* ra, ZZ* a, long np);
	void toNTT(uint64_t* ra, ZZ* a, long np);

	void multX0(ZZ* x, ZZ* a, ZZ* b, long np, const ZZ& q);
	void multX0AndEqual(ZZ* a, ZZ* b, long np, const ZZ& q);
	void multNTTX0(ZZ* x, ZZ* a, uint64_t* rb, long np, const ZZ& q);
	void multNTTX0AndEqual(ZZ* a, uint64_t* rb, long np, const ZZ& q);
	void multDNTTX0(ZZ* x, uint64_t* ra, uint64_t* rb, long np, const ZZ& q);

	void multX1(ZZ* x, ZZ* a, ZZ* b, long np, const ZZ& q);
	void multX1AndEqual(ZZ* a, ZZ* b, long np, const ZZ& q);
	void multNTTX1(ZZ* x, ZZ* a, uint64_t* rb, long np, const ZZ& q);
	void multNTTX1AndEqual(ZZ* a, uint64_t* rb, long np, const ZZ& q);
	void multDNTTX1(ZZ* x, uint64_t* ra, uint64_t* rb, long np, const ZZ& q);

	void mult(ZZ* x, ZZ* a, ZZ* b, long np, const ZZ& q);
	void mult(ZZ* x, ZZ* a, long* b, long np, const ZZ& q);
	void multAndEqual(ZZ* a, ZZ* b, long np, const ZZ& q);
	void multNTT(ZZ* x, ZZ* a, uint64_t* rb, long np, const ZZ& q);
	void multNTT(ZZ* x, long* a, uint64_t* rb, long np, const ZZ& q);
	void multNTTAndEqual(ZZ* a, uint64_t* rb, long np, const ZZ& q);
	void multDNTT(ZZ* x, uint64_t* ra, uint64_t* rb, long np, const ZZ& q);

	void square(ZZ* x, ZZ* a, long np, const ZZ& q);
	void square(ZZ* x, long* a, const ZZ& q);
	void squareAndEqual(ZZ* a, long np, const ZZ& q);
	void squareNTT(ZZ* x, uint64_t* ra, long np, const ZZ& q);


	//----------------------------------------------------------------------------------
	//   OTHER
	//----------------------------------------------------------------------------------


	void normalizeAndEqual(ZZ* p, const ZZ& q);

	void mod(ZZ* res, ZZ* p, const ZZ& q);
	void modAndEqual(ZZ* p, const ZZ& q);

	void negate(ZZ* res, ZZ* p);
	void negateAndEqual(ZZ* p);

	void add(ZZ* res, ZZ* p1, ZZ* p2, const ZZ& q);
	void addAndEqual(ZZ* p1, ZZ* p2, const ZZ& q);

	void sub(ZZ* res, ZZ* p1, ZZ* p2, const ZZ& q);
	void subAndEqual(ZZ* p1, ZZ* p2, const ZZ& q);
	void subAndEqual2(ZZ* p1, ZZ* p2, const ZZ& q);

	void multByMonomial(ZZ* res, ZZ* p, long deg0, long deg1, const ZZ& q);
	void multByMonomialAndEqual(ZZ* p, long deg0, long deg1, const ZZ& q);

	void multByConst(ZZ* res, ZZ* p, ZZ& cnst, const ZZ& q);
	void multByConstAndEqual(ZZ* p, ZZ& cnst, const ZZ& q);


	//----------------------------------------------------------------------------------
	//   SHIFTING
	//----------------------------------------------------------------------------------


	void leftShift(ZZ* res, ZZ* p, long bits, const ZZ& q);
	void leftShiftAndEqual(ZZ* p, long bits, const ZZ& q);

	void rightShift(ZZ* res, ZZ* p, long bits);
	void rightShiftAndEqual(ZZ* p, long bits);


	//----------------------------------------------------------------------------------
	//   ROTATION & CONJUGATION & TRANSPOSITION
	//----------------------------------------------------------------------------------


	void leftRotate(ZZ* res, ZZ* p, long r0, long r1);
	void leftRotate(ZZ* res, long* p, long r0, long r1);

	void conjugate(ZZ* res, ZZ* p);
	void conjugate(ZZ* res, long* p);


	//----------------------------------------------------------------------------------
	//   SAMPLING
	//----------------------------------------------------------------------------------


	void sampleRLWE(ZZ* ax, ZZ* bx, ZZ* sx, long logq);

	void addGauss(ZZ* ax, const ZZ& q);
	void sampleHWT(ZZ* res);
	void sampleZO(long* res);
	void sampleUniform(ZZ* res, long logq);

};

#endif
