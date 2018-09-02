/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "RingMultiplier.h"

#include <NTL/lip.h>
#include <NTL/sp_arith.h>
#include <NTL/tools.h>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "Primes.h"

RingMultiplier::RingMultiplier(long logNx, long logQ) : logNx(logNx) {
	Nx = 1 << logNx;

	logNy = 8;
	Ny = 1 << logNy;

	logN = logNx + logNy;
	N = 1 << logN;

	long Mx = 1 << (logNx + 1);

	long bound = 2 + logN + 4 * logQ;
	long nprimes = ceil(bound / 59.0);

	pVec = new uint64_t[nprimes];
	prVec = new uint64_t[nprimes];
	pTwok = new long[nprimes];
	pInvVec = new uint64_t[nprimes];

	scaledRootxPows = new uint64_t*[nprimes];
	scaledRootxPowsInv = new uint64_t*[nprimes];

	scaledRootyPows = new uint64_t*[nprimes];
	scaledRootyPowsInv = new uint64_t*[nprimes];

	scaledNxInv = new uint64_t[nprimes];
	scaledNyInv = new uint64_t[nprimes];

	dftomegaPows = new uint64_t*[nprimes]();
	dftomegaPowsInv = new uint64_t*[nprimes]();
	omegaPows = new uint64_t*[nprimes]();
	omegaPowsInv = new uint64_t*[nprimes]();

	red_ss_array = new _ntl_general_rem_one_struct*[nprimes];

	long My = Ny + 1;
	gyPows = new long[My];
	long gy = 3;
	long gyPow = 1;
	for (long i = 0; i < Ny; ++i) {
		gyPows[i] = gyPow;
		gyPow *= gy;
		gyPow %= My;
	}
	gyPows[Ny] = 1;

	for (long i = 0; i < nprimes; ++i) {
		pVec[i] = pPrimesVec[i];
		red_ss_array[i] = _ntl_general_rem_one_struct_build(pVec[i]);
		pInvVec[i] = inv(pVec[i]);
		pTwok[i] = (2 * ((long) log2(pVec[i]) + 1));
		prVec[i] = (static_cast<unsigned __int128>(1) << pTwok[i]) / pVec[i];
		uint64_t root = findMthRootOfUnity(Mx, pVec[i]);
		uint64_t rootinv = powMod(root, pVec[i] - 2, pVec[i]);

		uint64_t NxInvModp = powMod(Nx, pVec[i] - 2, pVec[i]);
		mulMod(scaledNxInv[i], NxInvModp, (1ULL << 32), pVec[i]);
		mulMod(scaledNxInv[i], scaledNxInv[i],(1ULL << 32), pVec[i]);

		uint64_t NyInvModp = powMod(Ny, pVec[i] - 2, pVec[i]);
		mulMod(scaledNyInv[i], NyInvModp, (1ULL << 32), pVec[i]);
		mulMod(scaledNyInv[i], scaledNyInv[i],(1ULL << 32), pVec[i]);

		scaledRootxPows[i] = new uint64_t[Nx]();
		scaledRootxPowsInv[i] = new uint64_t[Nx]();
		uint64_t power = 1;
		uint64_t powerInv = 1;
		for (long j = 0; j < Nx; ++j) {
			uint32_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logNx);
			uint64_t rootpow = power;
			mulMod(scaledRootxPows[i][jprime], rootpow, (1ULL << 32), pVec[i]);
			mulMod(scaledRootxPows[i][jprime], scaledRootxPows[i][jprime], (1ULL << 32), pVec[i]);
			uint64_t rootpowinv = powerInv;
			mulMod(scaledRootxPowsInv[i][jprime], rootpowinv, (1ULL << 32), pVec[i]);
			mulMod(scaledRootxPowsInv[i][jprime], scaledRootxPowsInv[i][jprime], (1ULL << 32), pVec[i]);
			mulMod(power, power, root, pVec[i]);
			mulMod(powerInv, powerInv, rootinv, pVec[i]);
		}

		root = findMthRootOfUnity(Ny, pVec[i]);
		rootinv = powMod(root, pVec[i] - 2, pVec[i]);
		scaledRootyPows[i] = new uint64_t[Ny]();
		scaledRootyPowsInv[i] = new uint64_t[Ny]();

		power = 1;
		powerInv = 1;
		for (long j = 0; j < Ny; ++j) {
			uint64_t rootpow = power;
			mulMod(scaledRootyPows[i][j], rootpow, (1ULL << 32), pVec[i]);
			mulMod(scaledRootyPows[i][j], scaledRootyPows[i][j], (1ULL << 32), pVec[i]);
			uint64_t rootpowinv = powerInv;
			mulMod(scaledRootyPowsInv[i][j], rootpowinv, (1ULL << 32), pVec[i]);
			mulMod(scaledRootyPowsInv[i][j], scaledRootyPowsInv[i][j], (1ULL << 32), pVec[i]);
			mulMod(power, power, root, pVec[i]);
			mulMod(powerInv, powerInv, rootinv, pVec[i]);
		}

		root = findMthRootOfUnity(My, pVec[i]);
		dftomegaPows[i] = new uint64_t[Ny]();
		dftomegaPowsInv[i] = new uint64_t[Ny]();
		omegaPows[i] = new uint64_t[Ny]();
		omegaPowsInv[i] = new uint64_t[Ny]();
		for (long j = 0; j < Ny; ++j) {
			dftomegaPows[i][j] = powMod(root, gyPows[j], pVec[i]);
			omegaPows[i][j] = dftomegaPows[i][j];
			omegaPowsInv[i][j] = powMod(root, My - gyPows[j], pVec[i]);
		}
		NTTPO2Y(dftomegaPows[i], i);
		for (long j = 0; j < Ny; ++j) {
			dftomegaPowsInv[i][j] = powMod(dftomegaPows[i][j], pVec[i] - 2, pVec[i]);
		}
	}

	coeffpinv_array = new mulmod_precon_t*[nprimes];
	pProd = new ZZ[nprimes];
	pProdh = new ZZ[nprimes];
	pHat = new ZZ*[nprimes];
	pHatInvModp = new uint64_t*[nprimes];
	for (long i = 0; i < nprimes; ++i) {
		pProd[i] = (i == 0) ? to_ZZ((long) pVec[i]) : pProd[i - 1] * (long) pVec[i];
		pProdh[i] = pProd[i] / 2;
		pHat[i] = new ZZ[i + 1];
		pHatInvModp[i] = new uint64_t[i + 1];
		coeffpinv_array[i] = new mulmod_precon_t[i + 1];
		for (long j = 0; j < i + 1; ++j) {
			pHat[i][j] = ZZ(1);
			for (long k = 0; k < j; ++k) {
				pHat[i][j] *= (long) pVec[k];
			}
			for (long k = j + 1; k < i + 1; ++k) {
				pHat[i][j] *= (long) pVec[k];
			}
			pHatInvModp[i][j] = to_long(pHat[i][j] % (long) pVec[j]);
			pHatInvModp[i][j] = powMod(pHatInvModp[i][j], pVec[j] - 2, pVec[j]);
			coeffpinv_array[i][j] = PrepMulModPrecon(pHatInvModp[i][j], pVec[j]);
		}
	}
}

void RingMultiplier::arrayBitReverse(uint64_t* vals, long n) {
	for (long i = 1, j = 0; i < n; ++i) {
		long bit = n >> 1;
		for (; j >= bit; bit>>=1) {
			j -= bit;
		}
		j += bit;
		if(i < j) {
			swap(vals[i], vals[j]);
		}
	}
}

void RingMultiplier::butt1(uint64_t& a1, uint64_t& a2, uint64_t& p, uint64_t& pInv, uint64_t& W) {
	uint64_t U = a1 + a2;
	if (U > p) U -= p;
	uint64_t T = a1 < a2 ? a1 + p - a2 : a1 - a2;
	unsigned __int128 UU = static_cast<unsigned __int128>(T) * W;
	uint64_t U0 = static_cast<uint64_t>(UU);
	uint64_t U1 = UU >> 64;
	uint64_t Q = U0 * pInv;
	unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
	uint64_t H = Hx >> 64;
	a1 = U;
	a2 = (U1 < H) ? U1 + p - H : U1 - H;
}

void RingMultiplier::butt2(uint64_t& a1, uint64_t& a2, uint64_t& p, uint64_t& pInv, uint64_t& W) {
	uint64_t T = a2;
	unsigned __int128 U = static_cast<unsigned __int128>(T) * W;
	uint64_t U0 = static_cast<uint64_t>(U);
	uint64_t U1 = U >> 64;
	uint64_t Q = U0 * pInv;
	unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
	uint64_t H = Hx >> 64;
	uint64_t V = U1 < H ? U1 + p - H : U1 - H;
	a2 = a1 < V ? a1 + p - V : a1 - V;
	a1 += V;
	if (a1 > p) a1 -= p;
}

void RingMultiplier::divByN(uint64_t& a, uint64_t& p, uint64_t& pInv, uint64_t& NScaleInv) {
	uint64_t T = a;
	unsigned __int128 U = static_cast<unsigned __int128>(T) * NScaleInv;
	uint64_t U0 = static_cast<uint64_t>(U);
	uint64_t U1 = U >> 64;
	uint64_t Q = U0 * pInv;
	unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
	uint64_t H = Hx >> 64;
	a = (U1 < H) ? U1 + p - H : U1 - H;
}

void RingMultiplier::NTTX(uint64_t* a, long index) {
	long t = Nx;
	long logt1 = logNx + 1;
	uint64_t p = pVec[index];
	uint64_t pInv = pInvVec[index];
	for (long m = 1; m < Nx; m <<= 1) {
		t >>= 1;
		logt1 -= 1;
		for (long i = 0; i < m; i++) {
			long j1 = i << logt1;
			long j2 = j1 + t - 1;
			uint64_t W = scaledRootxPows[index][m + i];
			for (long j = j1; j <= j2; j++) {
				butt2(a[j], a[j + t], p, pInv, W);
			}
		}
	}
}

void RingMultiplier::INTTX(uint64_t* a, long index) {
	uint64_t p = pVec[index];
	uint64_t pInv = pInvVec[index];
	long t = 1;
	for (long m = Nx; m > 1; m >>= 1) {
		long j1 = 0;
		long h = m >> 1;
		for (long i = 0; i < h; i++) {
			long j2 = j1 + t - 1;
			uint64_t W = scaledRootxPowsInv[index][h + i];
			for (long j = j1; j <= j2; j++) {
				butt1(a[j], a[j+t], p, pInv, W);
			}
			j1 += (t << 1);
		}
		t <<= 1;
	}

	uint64_t NxScale = scaledNxInv[index];
	for (long i = 0; i < Nx; i++) {
		divByN(a[i], p, pInv, NxScale);
	}
}

void RingMultiplier::NTTPO2Y(uint64_t* a, long index) {
	uint64_t p = pVec[index];
	uint64_t pInv = pInvVec[index];
	arrayBitReverse(a, Ny);
	for (long i = 0; i < logNy; ++i) {
		long ihpow = 1 << i;
		long ipow = 1 << (i + 1);
		for (long j = 0; j < Ny; j += ipow) {
			for (long k = 0; k < ihpow; ++k) {
				long idx = k << (logNy - i - 1);
				uint64_t W = scaledRootyPows[index][idx];
				butt2(a[j + k], a[j + k + ihpow], p, pInv, W);
			}
		}
	}
}

void RingMultiplier::INTTPO2Y(uint64_t* a, long index) {
	uint64_t p = pVec[index];
	uint64_t pInv = pInvVec[index];
	arrayBitReverse(a, Ny);
	for (long i = 0; i < logNy; ++i) {
		long ihpow = 1 << i;
		long ipow = 1 << (i + 1);
		for (long j = 0; j < Ny; j += ipow) {
			for (long k = 0; k < ihpow; ++k) {
				long idx = k << (logNy - i - 1);
				uint64_t W = scaledRootyPowsInv[index][idx];
				butt2(a[j + k], a[j + k + ihpow], p, pInv, W);
			}
		}
	}
	uint64_t NyScale = scaledNyInv[index];
	for (long i = 0; i < Ny; i++) {
		divByN(a[i], p, pInv, NyScale);
	}
}

void RingMultiplier::NTTY(uint64_t* a, long index) {
	uint64_t pi = pVec[index];
	uint64_t pri = prVec[index];
	uint64_t pti = pTwok[index];
	uint64_t* dftomegaPowsi = dftomegaPows[index];
	uint64_t* omegaPowsInvi = omegaPowsInv[index];

	uint64_t* tmp = new uint64_t[Ny];

	for (long i = 0; i < Ny; ++i) {
		tmp[i] = a[gyPows[Ny - i] - 1];
	}
	for (long i = 0; i < Ny; ++i) {
		a[i] = tmp[i];
	}
	delete[] tmp;
	NTTPO2Y(a, index);
	for (long i = 0; i < Ny; ++i) {
		mulModBarrett(a[i], a[i], dftomegaPowsi[i], pi, pri, pti);
	}
	INTTPO2Y(a, index);
	for (long i = 0; i < Ny; ++i) {
		mulModBarrett(a[i], a[i], omegaPowsInvi[i], pi, pri, pti);
	}
}

void RingMultiplier::INTTY(uint64_t* a, long index) {
	uint64_t pi = pVec[index];
	uint64_t pri = prVec[index];
	uint64_t pti = pTwok[index];
	uint64_t* omegaPowsi = omegaPows[index];
	uint64_t* dftomegaPowsInvi = dftomegaPowsInv[index];

	for (long i = 0; i < Ny; ++i) {
		mulModBarrett(a[i], a[i], omegaPowsi[i], pi, pri, pti);
	}

	NTTPO2Y(a, index);

	for (long i = 0; i < Ny; ++i) {
		mulModBarrett(a[i], a[i], dftomegaPowsInvi[i], pi, pri, pti);
	}

	INTTPO2Y(a, index);

	uint64_t* tmp = new uint64_t[Ny]();
	for (long i = 0; i < Ny; ++i) {
		tmp[gyPows[Ny - i] - 1] = a[i];
	}
	for (long i = 0; i < Ny; ++i) {
		a[i] = tmp[i];
	}
	delete[] tmp;
}

void RingMultiplier::NTTY1(uint64_t* a, long index) {
	uint64_t pi = pVec[index];
	uint64_t pri = prVec[index];
	uint64_t pti = pTwok[index];
	uint64_t* dftomegaPowsi = dftomegaPows[index];

	uint64_t* tmp = new uint64_t[Ny];

	for (long i = 0; i < Ny; ++i) {
		tmp[i] = a[gyPows[Ny - i] - 1];
	}
	for (long i = 0; i < Ny; ++i) {
		a[i] = tmp[i];
	}
	delete[] tmp;
	NTTPO2Y(a, index);
	for (long i = 0; i < Ny; ++i) {
		mulModBarrett(a[i], a[i], dftomegaPowsi[i], pi, pri, pti);
	}
	INTTPO2Y(a, index);
}

void RingMultiplier::INTTY1(uint64_t* a, long index) {
	uint64_t pi = pVec[index];
	uint64_t pri = prVec[index];
	uint64_t pti = pTwok[index];
	uint64_t* dftomegaPowsInvi = dftomegaPowsInv[index];

	NTTPO2Y(a, index);

	for (long i = 0; i < Ny; ++i) {
		mulModBarrett(a[i], a[i], dftomegaPowsInvi[i], pi, pri, pti);
	}

	INTTPO2Y(a, index);

	uint64_t* tmp = new uint64_t[Ny]();
	for (long i = 0; i < Ny; ++i) {
		tmp[gyPows[Ny - i] - 1] = a[i];
	}
	for (long i = 0; i < Ny; ++i) {
		a[i] = tmp[i];
	}
	delete[] tmp;
}

void RingMultiplier::NTTXY(uint64_t* a, long index) {
	for (long j = 0; j < Ny; ++j) {
		uint64_t* aj = a + (j << logNx);
		NTTX(aj, index);
	}
	uint64_t* tmp = new uint64_t[Ny];
	for (long j = 0; j < Nx; ++j) {
		for (long k = 0; k < Ny; ++k) {
			tmp[k] = a[j + (k << logNx)];
		}
		NTTY(tmp, index);
		for (long k = 0; k < Ny; ++k) {
			a[j + (k << logNx)] = tmp[k];
		}
	}
	delete[] tmp;
}

void RingMultiplier::INTTXY(uint64_t* a, long index) {
	uint64_t* tmp = new uint64_t[Ny];
	for (long j = 0; j < Nx; ++j) {
		for (long k = 0; k < Ny; ++k) {
			tmp[k] = a[j + (k << logNx)];
		}
		INTTY(tmp, index);
		for (long k = 0; k < Ny; ++k) {
			a[j + (k << logNx)] = tmp[k];
		}
	}
	delete[] tmp;

	for (long j = 0; j < Ny; ++j) {
		uint64_t* aj = a + (j << logNx);
		INTTX(aj, index);
	}
}

void RingMultiplier::NTTXY1(uint64_t* a, long index) {
	for (long j = 0; j < Ny; ++j) {
		uint64_t* aj = a + (j << logNx);
		NTTX(aj, index);
	}
	uint64_t* tmp = new uint64_t[Ny];
	for (long j = 0; j < Nx; ++j) {
		for (long k = 0; k < Ny; ++k) {
			tmp[k] = a[j + (k << logNx)];
		}
		NTTY1(tmp, index);
		for (long k = 0; k < Ny; ++k) {
			a[j + (k << logNx)] = tmp[k];
		}
	}
	delete[] tmp;
}

void RingMultiplier::INTTXY1(uint64_t* a, long index) {
	uint64_t* tmp = new uint64_t[Ny];
	for (long j = 0; j < Nx; ++j) {
		for (long k = 0; k < Ny; ++k) {
			tmp[k] = a[j + (k << logNx)];
		}
		INTTY1(tmp, index);
		for (long k = 0; k < Ny; ++k) {
			a[j + (k << logNx)] = tmp[k];
		}
	}
	delete[] tmp;

	for (long j = 0; j < Ny; ++j) {
		uint64_t* aj = a + (j << logNx);
		INTTX(aj, index);
	}
}

uint64_t* RingMultiplier::toNTTX(ZZ* a, long maxBnd) {
	long np = ceil(maxBnd / 59.0);

	uint64_t* ra = new uint64_t[np << logNx];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logNx);
		for (long n = 0; n < Nx; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTX(rai, i);
	}
	NTL_EXEC_RANGE_END;
	return ra;
}

uint64_t* RingMultiplier::toNTTY(ZZ* a, long maxBnd) {
	long np = ceil(maxBnd / 59.0);
	uint64_t* ra = new uint64_t[np << logNy];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logNy);
		for (long n = 0; n < Ny; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTY(rai, i);
	}
	NTL_EXEC_RANGE_END;
	return ra;
}

uint64_t* RingMultiplier::toNTTY1(ZZ* a, long maxBnd) {
	long np = ceil(maxBnd / 59.0);
	uint64_t* ra = new uint64_t[np << logNy];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logNy);
		for (long n = 0; n < Ny; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTY1(rai, i);
	}
	NTL_EXEC_RANGE_END;
	return ra;
}

uint64_t* RingMultiplier::toNTTXY(ZZ* a, long maxBnd) {
	long np = ceil(maxBnd / 59.0);
	uint64_t* ra = new uint64_t[np << logN];
	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY(rai, i);
	}
	NTL_EXEC_RANGE_END;
	return ra;
}

uint64_t* RingMultiplier::toNTTXY1(ZZ* a, long maxBnd) {
	long np = ceil(maxBnd / 59.0);
	uint64_t* ra = new uint64_t[np << logN];
	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t pi = pVec[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY1(rai, i);
	}
	NTL_EXEC_RANGE_END;
	return ra;
}

long RingMultiplier::MaxBits(const ZZ* f, long n) {
   long i, m;
   m = 0;

   for (i = 0; i < n; i++) {
      m = max(m, NumBits(f[i]));
   }
   return m;
}

void RingMultiplier::reconstruct(ZZ* x, uint64_t* rx, long np, ZZ& q) {
	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	NTL_EXEC_RANGE(N, first, last);
	for (long n = first; n < last; ++n) {
		ZZ& acc = x[n];
		QuickAccumBegin(acc, pProd[np - 1].size());
		for (long i = 0; i < np; i++) {
			long p = pVec[i];
			long tt = pHatInvModpnp[i];
			mulmod_precon_t ttpinv = coeffpinv_array[np-1][i];
			long s = MulModPrecon(rx[n + (i << logN)], tt, p, ttpinv);
			QuickAccumMulAdd(acc, pHatnp[i], s);
		}
		QuickAccumEnd(acc);
		rem(x[n], x[n], pProd[np - 1]);
		if (x[n] > pProdh[np - 1]) x[n] -= pProd[np - 1];
		x[n] %= q;
	}
	NTL_EXEC_RANGE_END;
}

void RingMultiplier::multXpoly(ZZ* x, ZZ* a, ZZ* b, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + MaxBits(b, Nx) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logNx];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNx);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		for (long iy = 0; iy < N; iy += Nx) {
			uint64_t* raij = rai + iy;
			NTTX(raij, i);
		}

		for (int ix = 0; ix < Nx; ++ix) {
			rbi[ix] = _ntl_general_rem_one_struct_apply(b[ix].rep, pi, red_ss);
		}
		NTTX(rbi, i);

		uint64_t* rxi = rx + (i << logN);
		for (long ix = 0; ix < Nx; ++ix) {
			for (long iy = 0; iy < N; iy += Nx) {
				long n = ix + iy;
				mulModBarrett(rxi[n], rai[n], rbi[ix], pi, pri, pTwoki);
			}
		}

		for (long iy = 0; iy < N; iy += Nx) {
			uint64_t* rxij = rxi + iy;
			INTTX(rxij, i);
		}
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] ra;
	delete[] rb;
	delete[] rx;
}

void RingMultiplier::multXpolyAndEqual(ZZ* a, ZZ* b, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + MaxBits(b, Nx) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logNx];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNx);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];

		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		for (long j = 0; j < Ny; ++j) {
			uint64_t* raij = rai + (j << logNx);
			NTTX(raij, i);
		}

		for (long ix = 0; ix < Nx; ++ix) {
			rbi[ix] = _ntl_general_rem_one_struct_apply(b[ix].rep, pi, red_ss);
		}
		NTTX(rbi, i);

		for (long ix = 0; ix < Nx; ++ix) {
			for (long iy = 0; iy < N; iy += Nx) {
				long n = ix + iy;
				mulModBarrett(rai[n], rai[n], rbi[ix], pi, pri, pTwoki);
			}
		}

		for (long j = 0; j < Ny; ++j) {
			uint64_t* raij = rai + (j << logNx);
			INTTX(raij, i);
		}

	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, q);

	delete[] ra;
	delete[] rb;
}

void RingMultiplier::multXpolyNTT(ZZ* x, ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {

}

void RingMultiplier::multXpolyNTTAndEqual(ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {

}

void RingMultiplier::multXpolyNTT2(ZZ* x, uint64_t* ra, long raBnd, uint64_t* rb, long rbBnd, ZZ& q) {

}


void RingMultiplier::multYpoly(ZZ* x, ZZ* a, ZZ* b, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + MaxBits(b, Ny) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logNy];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}

		uint64_t* rbi = rb + (i << logNy);
		for (int iy = 0; iy < Ny; ++iy) {
			rbi[iy] = _ntl_general_rem_one_struct_apply(b[iy].rep, pi, red_ss);
		}
		NTTY(rbi, i);

		uint64_t* rxi = rx + (i << logN);
		for (long iy = 0; iy < Ny; ++iy) {
			uint64_t rbiy = rbi[iy];
			uint64_t* rxiy = rxi + (iy << logNx);
			uint64_t* raiy = rai + (iy << logNx);
			for (long ix = 0; ix < Nx; ++ix) {
				mulModBarrett(rxiy[ix], raiy[ix], rbiy, pi, pri, pTwoki);
			}
		}

		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rxi[j + (k << logNx)];
			}
			INTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rxi[j + (k << logNx)] = tmp[k];
			}
		}
		delete[] tmp;
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] ra;
	delete[] rb;
	delete[] rx;
}

void RingMultiplier::multYpolyAndEqual(ZZ* a, ZZ* b, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + MaxBits(b, Ny) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logNy];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		uint64_t* tmp = new uint64_t[1 << logNy];
		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}

		uint64_t* rbi = rb + (i << logNy);
		for (long iy = 0; iy < Ny; ++iy) {
			rbi[iy] = _ntl_general_rem_one_struct_apply(b[iy].rep, pi, red_ss);
		}
		NTTY(rbi, i);

		for (long iy = 0; iy < Ny; ++iy) {
			uint64_t rbiy = rbi[iy];
			uint64_t* raiy = rai + (iy << logNx);
			for (long ix = 0; ix < Nx; ++ix) {
				mulModBarrett(raiy[ix], raiy[ix], rbiy, pi, pri, pTwoki);
			}
		}

		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			INTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}
		delete[] tmp;
	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, q);

	delete[] ra;
	delete[] rb;
}

void RingMultiplier::multYpolyNTT(ZZ* x, ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNy);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}

		uint64_t* rxi = rx + (i << logN);
		for (long iy = 0; iy < Ny; ++iy) {
			uint64_t rbiy = rbi[iy];
			uint64_t* rxiy = rxi + (iy << logNx);
			uint64_t* raiy = rai + (iy << logNx);
			for (long ix = 0; ix < Nx; ++ix) {
				mulModBarrett(rxiy[ix], raiy[ix], rbiy, pi, pri, pTwoki);
			}
		}

		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rxi[j + (k << logNx)];
			}
			INTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rxi[j + (k << logNx)] = tmp[k];
			}
		}
		delete[] tmp;
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] ra;
	delete[] rx;
}

void RingMultiplier::multYpolyNTTAndEqual(ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNy);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		uint64_t* tmp = new uint64_t[1 << logNy];
		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}

		for (long iy = 0; iy < Ny; ++iy) {
			uint64_t rbiy = rbi[iy];
			uint64_t* raiy = rai + (iy << logNx);
			for (long ix = 0; ix < Nx; ++ix) {
				mulModBarrett(raiy[ix], raiy[ix], rbiy, pi, pri, pTwoki);
			}
		}

		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			INTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}
		delete[] tmp;
	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, q);

	delete[] ra;
}

void RingMultiplier::multYpolyNTTD(ZZ* x, uint64_t* ra, long raBnd, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + raBnd + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNy);
		uint64_t* rxi = rx + (i << logN);
		for (long iy = 0; iy < Ny; ++iy) {
			uint64_t rbiy = rbi[iy];
			uint64_t* rxiy = rxi + (iy << logNx);
			uint64_t* raiy = rai + (iy << logNx);
			for (long ix = 0; ix < Nx; ++ix) {
				mulModBarrett(rxiy[ix], raiy[ix], rbiy, pi, pri, pTwoki);
			}
		}

		uint64_t* tmp = new uint64_t[1 << logNy];
		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rxi[j + (k << logNx)];
			}
			INTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rxi[j + (k << logNx)] = tmp[k];
			}
		}
		delete[] tmp;
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] rx;
}

void RingMultiplier::mult(ZZ* x, ZZ* a, ZZ* b, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + MaxBits(b, N) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logN];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];

		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
			rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss);
		}
		NTTXY(rai, i);
		NTTXY1(rbi, i);

		uint64_t* rxi = rx + (i << logN);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		INTTXY1(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] ra;
	delete[] rb;
	delete[] rx;
}

void RingMultiplier::multAndEqual(ZZ* a, ZZ* b, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + MaxBits(b, N) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
			rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss);
		}
		NTTXY(rai, i);
		NTTXY1(rbi, i);

		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		INTTXY1(rai, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, q);

	delete[] ra;
	delete[] rb;
}

void RingMultiplier::multNTTXY(ZZ* x, ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];

		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY1(rai, i);

		uint64_t* rxi = rx + (i << logN);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		INTTXY1(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] ra;
	delete[] rx;
}

void RingMultiplier::multNTTXYAndEqual(ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY1(rai, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		INTTXY1(rai, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, q);

	delete[] ra;
}

void RingMultiplier::multNTTXY1(ZZ* x, ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];

		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY(rai, i);

		uint64_t* rxi = rx + (i << logN);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		INTTXY1(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] ra;
	delete[] rx;
}

void RingMultiplier::multNTTXY1AndEqual(ZZ* a, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + MaxBits(a, N) + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY(rai, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		INTTXY1(rai, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, q);

	delete[] ra;
}

void RingMultiplier::multNTTXYD(ZZ* x, uint64_t* ra, long raBnd, uint64_t* rb, long rbBnd, ZZ& q) {
	long bound = 2 * logN + raBnd + rbBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];

		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		INTTXY1(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] rx;
}

void RingMultiplier::square(ZZ* x, ZZ* a, ZZ& q) {
	long bound = 2 * logN + 2 * MaxBits(a, N) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY(rai, i);

		uint64_t* rxi = rx + (i << logN);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rai[n], pi, pri, pTwoki);
		}

		INTTXY(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] ra;
	delete[] rx;
}

void RingMultiplier::squareAndEqual(ZZ* a, ZZ& q) {
	long bound = 2 * logN + 2 * MaxBits(a, N) + 2;
	long np = ceil(bound / 59.0);
	uint64_t* ra = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = 0; i < np; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		uint64_t* rai = ra + (i << logN);
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTTXY(rai, i);

		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rai[n], pi, pri, pTwoki);
		}

		INTTXY(rai, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, q);

	delete[] ra;
}

void RingMultiplier::squareNTT(ZZ* x, uint64_t* ra, long raBnd, ZZ& q) {
	long bound = 2 * logN + 2 * raBnd + 2;
	long np = ceil(bound / 59.0);

	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rai[n], pi, pri, pTwoki);
		}

		INTTXY(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, q);

	delete[] rx;
}

void RingMultiplier::mulMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t m) {
	unsigned __int128
	mul = static_cast<unsigned __int128>(a) * b;
	mul %= static_cast<unsigned __int128>(m);
	r = static_cast<uint64_t>(mul);
}

void RingMultiplier::mulModBarrett(uint64_t& r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr, long twok) {
	unsigned __int128
	mul = static_cast<unsigned __int128>(a) * b;
	uint64_t atop, abot;
	abot = static_cast<uint64_t>(mul);
	atop = static_cast<uint64_t>(mul >> 64);
	unsigned __int128
	tmp = static_cast<unsigned __int128>(abot) * pr;
	tmp >>= 64;
	tmp += static_cast<unsigned __int128>(atop) * pr;
	tmp >>= twok - 64;
	tmp *= p;
	tmp = mul - tmp;
	r = static_cast<uint64_t>(tmp);
	if (r >= p) {
		r -= p;
	}
}

uint64_t RingMultiplier::powMod(uint64_t x, uint64_t y, uint64_t modulus) {
	uint64_t res = 1;
	while (y > 0) {
		if (y & 1) {
			mulMod(res, res, x, modulus);
		}
		y = y >> 1;
		mulMod(x, x, x, modulus);
	}
	return res;
}

uint64_t RingMultiplier::inv(uint64_t x) {
	uint64_t res = 1;
	for (long i = 0; i < 62; ++i) {
		res *= x;
		x *= x;
	}
	return res;
}

uint32_t RingMultiplier::bitReverse(uint32_t x) {
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return ((x >> 16) | (x << 16));
}

void RingMultiplier::findPrimeFactors(vector<uint64_t> &s, uint64_t number) {
	while (number % 2 == 0) {
		s.push_back(2);
		number /= 2;
	}
	for (uint64_t i = 3; i < sqrt(number); i++) {
		while (number % i == 0) {
			s.push_back(i);
			number /= i;
		}
	}
	if (number > 2) {
		s.push_back(number);
	}
}

uint64_t RingMultiplier::findPrimitiveRoot(uint64_t modulus) {
	vector<uint64_t> s;
	uint64_t phi = modulus - 1;
	findPrimeFactors(s, phi);
	for (uint64_t r = 2; r <= phi; r++) {
		bool flag = false;
		for (uint64_t i : s) {
			if(powMod(r, phi / i, modulus) == 1) {
				flag = true;
				break;
			}
		}
		if (flag == false) {
			return r;
		}
	}
	return -1;
}

// Algorithm to find m-th primitive root in Z_mod
uint64_t RingMultiplier::findMthRootOfUnity(uint64_t M, uint64_t mod) {
	uint64_t res = findPrimitiveRoot(mod);
	uint64_t factor = (mod - 1) / M;
	res = powMod(res, factor, mod);
	return res;
}
