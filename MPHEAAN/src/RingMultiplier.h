/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MHEAAN_RINGMULTIPLIER_H_
#define MHEAAN_RINGMULTIPLIER_H_

#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class RingMultiplier {
public:

	long logN0;
	long N0;

	long logN1;
	long N1;

	long logN;
	long N;

	long* g1Pows;
	uint64_t** rootM1DFTPows;
	uint64_t** rootM1DFTPowsInv;
	uint64_t** rootM1Pows;
	uint64_t** rootM1PowsInv;

	uint64_t* pVec;
	uint64_t* prVec;
	long* pTwok;

	uint64_t* pInvVec;

	uint64_t** scaledRootM0Pows;
	uint64_t** scaledRootN1Pows;

	uint64_t** scaledRootM0PowsInv;
	uint64_t** scaledRootN1PowsInv;

	uint64_t* scaledN0Inv;
	uint64_t* scaledN1Inv;

	_ntl_general_rem_one_struct** red_ss_array;
	mulmod_precon_t** coeffpinv_array;
	ZZ* pProd;
	ZZ* pProdh;
	ZZ** pHat;
	uint64_t** pHatInvModp;

	RingMultiplier(long logN0 = 0, long nprimes = 0);

	void arrayBitReverse(uint64_t* a, long n);

	void butt1(uint64_t& a1, uint64_t& a2, uint64_t& p, uint64_t& pInv, uint64_t& W);
	void butt2(uint64_t& a1, uint64_t& a2, uint64_t& p, uint64_t& pInv, uint64_t& W);

	void divByN(uint64_t& a, uint64_t& p, uint64_t& pInv, uint64_t& NScaleInv);

	void NTTX0(uint64_t* a, long index);
	void INTTX0(uint64_t* a, long index);
	void NTTPO2X1(uint64_t* a, long index);
	void INTTPO2X1(uint64_t* a, long index);
	void NTTX1(uint64_t* a, long index);
	void INTTX1(uint64_t* a, long index);
	void NTTX1Lazy(uint64_t* a, long index);
	void INTTX1Lazy(uint64_t* a, long index);
	void NTT(uint64_t* a, long index);
	void INTT(uint64_t* a, long index);
	void NTTLazy(uint64_t* a, long index);
	void INTTLazy(uint64_t* a, long index);

	uint64_t* toNTTX0(ZZ* a, long np);
	uint64_t* toNTTX1(ZZ* a, long np);
	uint64_t* toNTTX1Lazy(ZZ* a, long np);
	uint64_t* toNTT(ZZ* a, long np);
	uint64_t* toNTTLazy(ZZ* a, long np);

	uint64_t* addNTT(uint64_t* ra, uint64_t* rb, long np);

	void reconstruct(ZZ* x, uint64_t* rx, long np, ZZ& q);

	void multX0(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	void multX0AndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTTX0(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX0AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multDNTTX0(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q);

	void multX1(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q);
	void multX1AndEqual(ZZ* a, ZZ* b, long np, ZZ& q);
	void multNTTX1(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX1AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX1Lazy(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q);
	void multNTTX1LazyAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q);
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

	void mulMod(uint64_t& r, uint64_t a, uint64_t b, uint64_t p);

	void mulModBarrett(uint64_t& r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr, long twok);

	uint64_t powMod(uint64_t x, uint64_t y, uint64_t p);

	uint64_t inv(uint64_t x);

	uint32_t bitReverse(uint32_t x);

	void findPrimeFactors(vector<uint64_t> &s, uint64_t number);

	uint64_t findPrimitiveRoot(uint64_t m);

	uint64_t findMthRootOfUnity(uint64_t M, uint64_t p);

};

#endif /* RINGMULTIPLIER_H_ */
