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
#include <vector>

#include "BootContext.h"
#include "RingMultiplier.h"
#include "SqrMatContext.h"

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

	long pbnd;

	long logQ; ///< log of Q
	long logQQ; ///< log of PQ

	ZZ Q; ///< Q corresponds to the highest modulus
	ZZ QQ; ///< PQ = Q * Q
	ZZ* qvec;

	double sigma; ///< standard deviation for Gaussian distribution
	long h; ///< parameter for HWT distribution

	RingMultiplier* multiplier;

	uint64_t* gM0Pows; ///< auxiliary information about rotation group indexes for batch encoding
	uint64_t* gM1Pows; ///< auxiliary information about rotation group indexes for batch encoding

	complex<double>* ksiM0Pows; ///< storing ksi pows for fft calculation
	complex<double>* ksiM1Pows;
	complex<double>* ksiN1Pows; ///< storing ksi pows for fft calculation

	complex<double>** dftM1Pows;
	complex<double>** dftM1NTTPows;

	map<pair<long, long>, BootContext*> bootContextMap;

	map<long, SqrMatContext*> sqrMatContextMap;

	Ring(long logN0, long logN1, long logQ, double sigma = 3.2, long h = 64);


	//----------------------------------------------------------------------------------
	//   AUXILIARY CONTEXT
	//----------------------------------------------------------------------------------


	void addBootContext(long logn0, long logn1, long logp);

	void addSqrMatContext(long logn, long logp);


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

	ZZ* encode(complex<double>* vals, long n0, long n1, long logp);
	ZZ* encode(double* vals, long n0, long n1, long logp);
	complex<double>* decode(ZZ* mxy, long n0, long n1, long logp, long logq);


	//----------------------------------------------------------------------------------
	//   MULTIPLICATION
	//----------------------------------------------------------------------------------

	long MaxBits(const ZZ* f, long n);
	void addNTTAndEqual(uint64_t* ra, uint64_t* rb, long np);

	void toNTTX0(uint64_t* ra, ZZ* a, long np);
	uint64_t* toNTTX0(ZZ* a, long np);
	void toNTTX1(uint64_t* ra, ZZ* a, long np);
	uint64_t* toNTTX1(ZZ* a, long np);
	void toNTT(uint64_t* ra, ZZ* a, long np);
	uint64_t* toNTT(ZZ* a, long np);

	void multX0(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	ZZ* multX0(ZZ* a, ZZ* b, long np, ZZ& q);
	void multX0AndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTTX0(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	ZZ* multNTTX0(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX0AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multDNTTX0(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q);
	ZZ* multDNTTX0(uint64_t* ra, uint64_t* rb, long np, ZZ& q);

	void multX1(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	ZZ* multX1(ZZ* a, ZZ* b, long np, ZZ& q);
	void multX1AndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTTX1(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	ZZ* multNTTX1(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX1AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multDNTTX1(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q);
	ZZ* multDNTTX1(uint64_t* ra, uint64_t* rb, long np, ZZ& q);

	void mult(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	ZZ* mult(ZZ* a, ZZ* b, long np, ZZ& q);
	void multAndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTT(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	ZZ* multNTT(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multDNTT(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q);
	ZZ* multDNTT(uint64_t* ra, uint64_t* rb, long np, ZZ& q);

	void square(ZZ* x, ZZ* a, long np, ZZ& q);
	ZZ* square(ZZ* a, long np, ZZ& q);
	void squareAndEqual(ZZ* a, long np, ZZ& q);
	void squareNTT(ZZ* x, uint64_t* ra, long np, ZZ& q);
	ZZ* squareNTT(uint64_t* ra, long np, ZZ& q);


	//----------------------------------------------------------------------------------
	//   OTHER
	//----------------------------------------------------------------------------------


	void mod(ZZ* res, ZZ* p, ZZ& q);
	ZZ* mod(ZZ* p, ZZ& q);
	void modAndEqual(ZZ* p, ZZ& q);

	void negate(ZZ* res, ZZ* p);
	ZZ* negate(ZZ* p);
	void negateAndEqual(ZZ* p);

	void add(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q);
	ZZ* add(ZZ* p1, ZZ* p2, ZZ& q);
	void addAndEqual(ZZ* p1, ZZ* p2, ZZ& q);

	void sub(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q);
	ZZ* sub(ZZ* p1, ZZ* p2, ZZ& q);
	void subAndEqual(ZZ* p1, ZZ* p2, ZZ& q);
	void subAndEqual2(ZZ* p1, ZZ* p2, ZZ& q);

	void multByMonomial(ZZ* res, ZZ* p, long deg0, long deg1, ZZ& q);
	ZZ* multByMonomial(ZZ* p, long deg0, long deg1, ZZ& q);
	void multByMonomialAndEqual(ZZ* p, long deg0, long deg1, ZZ& q);

	void multByConst(ZZ* res, ZZ* p, ZZ& cnst, ZZ& q);
	ZZ* multByConst(ZZ* p, ZZ& cnst, ZZ& q);
	void multByConstAndEqual(ZZ* p, ZZ& cnst, ZZ& q);


	//----------------------------------------------------------------------------------
	//   SHIFTING
	//----------------------------------------------------------------------------------


	void leftShift(ZZ* res, ZZ* p, long bits, ZZ& q);
	ZZ* leftShift(ZZ* p, long bits, ZZ& q);
	void leftShiftAndEqual(ZZ* p, long bits, ZZ& q);

	void rightShift(ZZ* res, ZZ* p, long bits);
	ZZ* rightShift(ZZ* p, long bits);
	void rightShiftAndEqual(ZZ* p, long bits);


	//----------------------------------------------------------------------------------
	//   ROTATION & CONJUGATION & TRANSPOSITION
	//----------------------------------------------------------------------------------


	void leftRotate(ZZ* res, ZZ* p, long r0, long r1);
	ZZ* leftRotate(ZZ* p, long r0, long r1);

	void conjugate(ZZ* res, ZZ* p);
	ZZ* conjugate(ZZ* p);


	//----------------------------------------------------------------------------------
	//   SAMPLING
	//----------------------------------------------------------------------------------


	void sampleGauss(ZZ* res);
	ZZ* sampleGauss();
	void sampleHWT(ZZ* res);
	ZZ* sampleHWT();
	void sampleZO(ZZ* res);
	ZZ* sampleZO();
	void sampleUniform(ZZ* res, long bits);
	ZZ* sampleUniform(long bits);

};

#endif
