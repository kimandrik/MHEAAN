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

	long logNx;
	long Nx;

	long logNy;
	long Ny;

	long logN;
	long N;

	long* gyPows;
	uint64_t** dftomegaPows;
	uint64_t** dftomegaPowsInv;
	uint64_t** omegaPows;
	uint64_t** omegaPowsInv;

	uint64_t* pVec;
	uint64_t* prVec;
	long* pTwok;

	uint64_t* pInvVec;

	uint64_t** scaledRootxPows;
	uint64_t** scaledRootyPows;

	uint64_t** scaledRootxPowsInv;
	uint64_t** scaledRootyPowsInv;

	uint64_t* scaledNxInv;
	uint64_t* scaledNyInv;

	_ntl_general_rem_one_struct** red_ss_array;
	mulmod_precon_t** coeffpinv_array;
	ZZ* pProd;
	ZZ* pProdh;
	ZZ** pHat;
	uint64_t** pHatInvModp;

	RingMultiplier(long logNx = 0, long logQ = 0);

	void arrayBitReverse(uint64_t* a, long n);

	void butt1(uint64_t& a1, uint64_t& a2, uint64_t& p, uint64_t& pInv, uint64_t& W);
	void butt2(uint64_t& a1, uint64_t& a2, uint64_t& p, uint64_t& pInv, uint64_t& W);

	void divByN(uint64_t& a, uint64_t& p, uint64_t& pInv, uint64_t& NScaleInv);

	void NTTX(uint64_t* a, long index);
	void INTTX(uint64_t* a, long index);

	void NTTPO2Y(uint64_t* a, long index);
	void INTTPO2Y(uint64_t* a, long index);

	void NTTY(uint64_t* a, long index);
	void INTTY(uint64_t* a, long index);

	void NTTY1(uint64_t* a, long index);
	void INTTY1(uint64_t* a, long index);

	void NTTXY(uint64_t* a, long index);
	void INTTXY(uint64_t* a, long index);

	void NTTXY1(uint64_t* a, long index);
	void INTTXY1(uint64_t* a, long index);

	long MaxBits(const ZZ* f, long n);

	void multXpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q);
	void multXpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q);

	void multYpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q);
	void multYpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q);

	void mult(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q);
	void multAndEqual(ZZ* a, const ZZ* b, const ZZ& q);

	void square(ZZ* x, const ZZ* a, const ZZ& q);
	void squareAndEqual(ZZ* a, const ZZ& q);

};

#endif /* RINGMULTIPLIER_H_ */
