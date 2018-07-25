/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "RingMultiplier.h"
#include "Numb.h"
#include "Primes.h"

RingMultiplier::RingMultiplier(long logNx, long logNy, long logQ) : logNx(logNx), logNy(logNy) {
	logN = logNx + logNy;

	Nx = 1 << logNx;
	Ny = 1 << logNy;
	N = 1 << logN;

	long Mx = 1 << (logNx + 1);
	long My = 1 << (logNy + 1);

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

	red_ss_array = new _ntl_general_rem_one_struct*[nprimes];

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

		root = findMthRootOfUnity(My, pVec[i]);
		rootinv = powMod(root, pVec[i] - 2, pVec[i]);
		scaledRootyPows[i] = new uint64_t[Ny]();
		scaledRootyPowsInv[i] = new uint64_t[Ny]();

		power = 1;
		powerInv = 1;
		for (long j = 0; j < Ny; ++j) {
			uint32_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logNy);
			uint64_t rootpow = power;
			mulMod(scaledRootyPows[i][jprime], rootpow, (1ULL << 32), pVec[i]);
			mulMod(scaledRootyPows[i][jprime], scaledRootyPows[i][jprime], (1ULL << 32), pVec[i]);
			uint64_t rootpowinv = powerInv;
			mulMod(scaledRootyPowsInv[i][jprime], rootpowinv, (1ULL << 32), pVec[i]);
			mulMod(scaledRootyPowsInv[i][jprime], scaledRootyPowsInv[i][jprime], (1ULL << 32), pVec[i]);
			mulMod(power, power, root, pVec[i]);
			mulMod(powerInv, powerInv, rootinv, pVec[i]);
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
				uint64_t T = a[j + t];
				unsigned __int128
				U = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(U);
				uint64_t U1 = U >> 64;
				uint64_t Q = U0 * pInv;
				unsigned __int128
				Hx = static_cast<unsigned __int128>(Q) * p;
				uint64_t H = Hx >> 64;
				uint64_t V = U1 < H ? U1 + p - H : U1 - H;
				a[j + t] = a[j] < V ? a[j] + p - V : a[j] - V;
				a[j] += V;
				if (a[j] > p) a[j] -= p;
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
				uint64_t U = a[j] + a[j + t];
				if (U > p) U -= p;
				uint64_t T = a[j] < a[j + t] ? a[j] + p - a[j + t] : a[j] - a[j + t];
				unsigned __int128
				UU = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(UU);
				uint64_t U1 = UU >> 64;
				uint64_t Q = U0 * pInv;
				unsigned __int128
				Hx = static_cast<unsigned __int128>(Q) * p;
				uint64_t H = Hx >> 64;
				a[j] = U;
				a[j + t] = (U1 < H) ? U1 + p - H : U1 - H;
			}
			j1 += (t << 1);
		}
		t <<= 1;
	}

	uint64_t NxScale = scaledNxInv[index];
	for (long i = 0; i < Nx; i++) {
		uint64_t T = a[i];
		unsigned __int128
		U = static_cast<unsigned __int128>(T) * NxScale;
		uint64_t U0 = static_cast<uint64_t>(U);
		uint64_t U1 = U >> 64;
		uint64_t Q = U0 * pInv;
		unsigned __int128
		Hx = static_cast<unsigned __int128>(Q) * p;
		uint64_t H = Hx >> 64;
		a[i] = (U1 < H) ? U1 + p - H : U1 - H;
	}
}

void RingMultiplier::NTTY(uint64_t* a, long index) {
	long t = Ny;
	long logt1 = logNy + 1;
	uint64_t p = pVec[index];
	uint64_t pInv = pInvVec[index];
	for (long m = 1; m < Ny; m <<= 1) {
		t >>= 1;
		logt1 -= 1;
		for (long i = 0; i < m; i++) {
			long j1 = i << logt1;
			long j2 = j1 + t - 1;
			uint64_t W = scaledRootyPows[index][m + i];
			for (long j = j1; j <= j2; j++) {
				uint64_t T = a[j + t];
				unsigned __int128
				U = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(U);
				uint64_t U1 = U >> 64;
				uint64_t Q = U0 * pInv;
				unsigned __int128
				Hx = static_cast<unsigned __int128>(Q) * p;
				uint64_t H = Hx >> 64;
				uint64_t V = U1 < H ? U1 + p - H : U1 - H;
				a[j + t] = a[j] < V ? a[j] + p - V : a[j] - V;
				a[j] += V;
				if (a[j] > p) a[j] -= p;
			}
		}
	}
}

void RingMultiplier::INTTY(uint64_t* a, long index) {
	uint64_t p = pVec[index];
	uint64_t pInv = pInvVec[index];
	long t = 1;
	for (long m = Ny; m > 1; m >>= 1) {
		long j1 = 0;
		long h = m >> 1;
		for (long i = 0; i < h; i++) {
			long j2 = j1 + t - 1;
			uint64_t W = scaledRootyPowsInv[index][h + i];
			for (long j = j1; j <= j2; j++) {
				uint64_t U = a[j] + a[j + t];
				if (U > p) U -= p;
				uint64_t T = a[j] < a[j + t] ? a[j] + p - a[j + t] : a[j] - a[j + t];
				unsigned __int128
				UU = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(UU);
				uint64_t U1 = UU >> 64;
				uint64_t Q = U0 * pInv;
				unsigned __int128
				Hx = static_cast<unsigned __int128>(Q) * p;
				uint64_t H = Hx >> 64;
				a[j] = U;
				a[j + t] = (U1 < H) ? U1 + p - H : U1 - H;
			}
			j1 += (t << 1);
		}
		t <<= 1;
	}

	uint64_t NyScale = scaledNyInv[index];
	for (long i = 0; i < Ny; i++) {
		uint64_t T = a[i];
		unsigned __int128
		U = static_cast<unsigned __int128>(T) * NyScale;
		uint64_t U0 = static_cast<uint64_t>(U);
		uint64_t U1 = U >> 64;
		uint64_t Q = U0 * pInv;
		unsigned __int128
		Hx = static_cast<unsigned __int128>(Q) * p;
		uint64_t H = Hx >> 64;
		a[i] = (U1 < H) ? U1 + p - H : U1 - H;
	}
}

long RingMultiplier::MaxBits(const ZZ* f, long n) {
   long i, m;
   m = 0;

   for (i = 0; i < n; i++) {
      m = max(m, NumBits(f[i]));
   }
   return m;
}

void RingMultiplier::multXpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q) {
	long bound = logNx + MaxBits(a, N) + MaxBits(b, Nx) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logNx];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNx);
		uint64_t* rxi = rx + (i << logN);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		for (int ix = 0; ix < Nx; ++ix) {
			rbi[ix] = _ntl_general_rem_one_struct_apply(b[ix].rep, pi, red_ss);
		}
		NTTX(rbi, i);
		for (long iy = 0; iy < N; iy += Nx) {
			uint64_t* raij = rai + iy;
			NTTX(raij, i);
		}

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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
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

	delete[] ra;
	delete[] rb;
	delete[] rx;
}

void RingMultiplier::multYpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q) {
	long bound = logNy + MaxBits(a, N) + MaxBits(b, Ny) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logNy];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNy);
		uint64_t* rxi = rx + (i << logN);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}

		for (int iy = 0; iy < Ny; ++iy) {
			rbi[iy] = _ntl_general_rem_one_struct_apply(b[iy].rep, pi, red_ss);
		}

		NTTY(rbi, i);

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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
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

	delete[] ra;
	delete[] rb;
	delete[] rx;
}

void RingMultiplier::mult(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q) {
	long bound = logN + MaxBits(a, N) + MaxBits(b, N) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logN];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
			rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss);
		}
		for (long j = 0; j < Ny; ++j) {
			uint64_t* raij = rai + (j << logNx);
			uint64_t* rbij = rbi + (j << logNx);
			NTTX(raij, i);
			NTTX(rbij, i);
		}

		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rbi[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rbi[j + (k << logNx)] = tmp[k];
			}
		}

		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		for (long j = 0; j < Ny; ++j) {
			uint64_t* rxij = rxi + (j << logNx);
			INTTX(rxij, i);
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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
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

	delete[] ra;
	delete[] rb;
	delete[] rx;
}

void RingMultiplier::multXpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q) {
	long bound = logNx + MaxBits(a, N) + MaxBits(b, Nx) + 2;
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

		for (long ix = 0; ix < Nx; ++ix) {
			rbi[ix] = _ntl_general_rem_one_struct_apply(b[ix].rep, pi, red_ss);
		}

		for (long j = 0; j < Ny; ++j) {
			uint64_t* raij = rai + (j << logNx);
			NTTX(raij, i);
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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
		ZZ& acc = a[n];
		QuickAccumBegin(acc, pProd[np - 1].size());
		for (long i = 0; i < np; i++) {
			long p = pVec[i];
			long tt = pHatInvModpnp[i];
			mulmod_precon_t ttpinv = coeffpinv_array[np-1][i];
			long s = MulModPrecon(ra[n + (i << logN)], tt, p, ttpinv);
			QuickAccumMulAdd(acc, pHatnp[i], s);
		}
		QuickAccumEnd(acc);
		rem(a[n], a[n], pProd[np - 1]);
		if (a[n] > pProdh[np - 1]) a[n] -= pProd[np - 1];
		a[n] %= q;
	}


	delete[] ra;
	delete[] rb;
}

void RingMultiplier::multYpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q) {
	long bound = logNy + MaxBits(a, N) + MaxBits(b, Ny) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logNy];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logNy);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}

		for (long iy = 0; iy < Ny; ++iy) {
			rbi[iy] = _ntl_general_rem_one_struct_apply(b[iy].rep, pi, red_ss);
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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
		ZZ& acc = a[n];
		QuickAccumBegin(acc, pProd[np - 1].size());
		for (long i = 0; i < np; i++) {
			long p = pVec[i];
			long tt = pHatInvModpnp[i];
			mulmod_precon_t ttpinv = coeffpinv_array[np-1][i];
			long s = MulModPrecon(ra[n + (i << logN)], tt, p, ttpinv);
			QuickAccumMulAdd(acc, pHatnp[i], s);
		}
		QuickAccumEnd(acc);
		rem(a[n], a[n], pProd[np - 1]);
		if (a[n] > pProdh[np - 1]) a[n] -= pProd[np - 1];
		a[n] %= q;
	}


	delete[] ra;
	delete[] rb;
}

void RingMultiplier::multAndEqual(ZZ* a, const ZZ* b, const ZZ& q) {
	long bound = logN + MaxBits(a, N) + MaxBits(b, N) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rb = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];

		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
			rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss);
		}

		for (long j = 0; j < Ny; ++j) {
			uint64_t* raij = rai + (j << logNx);
			uint64_t* rbij = rbi + (j << logNx);
			NTTX(raij, i);
			NTTX(rbij, i);
		}

		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rbi[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rbi[j + (k << logNx)] = tmp[k];
			}
		}

		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rbi[n], pi, pri, pTwoki);
		}

		for (long j = 0; j < Ny; ++j) {
			uint64_t* raij = rai + (j << logNx);
			INTTX(raij, i);
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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
		ZZ& acc = a[n];
		QuickAccumBegin(acc, pProd[np - 1].size());
		for (long i = 0; i < np; i++) {
			long p = pVec[i];
			long tt = pHatInvModpnp[i];
			mulmod_precon_t ttpinv = coeffpinv_array[np-1][i];
			long s = MulModPrecon(ra[n + (i << logN)], tt, p, ttpinv);
			QuickAccumMulAdd(acc, pHatnp[i], s);
		}
		QuickAccumEnd(acc);
		rem(a[n], a[n], pProd[np - 1]);
		if (a[n] > pProdh[np - 1]) a[n] -= pProd[np - 1];
		a[n] %= q;
	}


	delete[] ra;
	delete[] rb;
}

void RingMultiplier::square(ZZ* x, const ZZ* a, const ZZ& q) {
	long bound = logN + 2 * MaxBits(a, N) + 2;
	long np = ceil(bound / 59.0);

	uint64_t* ra = new uint64_t[np << logN];
	uint64_t* rx = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
		uint64_t* rxi = rx + (i << logN);
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
		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}

		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rai[n], pi, pri, pTwoki);
		}

		for (long j = 0; j < Ny; ++j) {
			uint64_t* rxij = rxi + (j << logNx);
			INTTX(rxij, i);
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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
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

	delete[] ra;
	delete[] rx;
}

void RingMultiplier::squareAndEqual(ZZ* a, const ZZ& q) {
	long bound = logN + 2 * MaxBits(a, N) + 2;
	long np = ceil(bound / 59.0);
	uint64_t* ra = new uint64_t[np << logN];

	NTL_EXEC_RANGE(np, first, last);
	for (long i = 0; i < np; ++i) {
		uint64_t* tmp = new uint64_t[1 << logNy];
		uint64_t* rai = ra + (i << logN);
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
		for (long j = 0; j < Nx; ++j) {
			for (long k = 0; k < Ny; ++k) {
				tmp[k] = rai[j + (k << logNx)];
			}
			NTTY(tmp, i);
			for (long k = 0; k < Ny; ++k) {
				rai[j + (k << logNx)] = tmp[k];
			}
		}

		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rai[n], pi, pri, pTwoki);
		}

		for (long j = 0; j < Ny; ++j) {
			uint64_t* raij = rai + (j << logNx);
			INTTX(raij, i);
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

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	for (long n = 0; n < N; ++n) {
		ZZ& acc = a[n];
		QuickAccumBegin(acc, pProd[np - 1].size());
		for (long i = 0; i < np; i++) {
			long p = pVec[i];
			long tt = pHatInvModpnp[i];
			mulmod_precon_t ttpinv = coeffpinv_array[np-1][i];
			long s = MulModPrecon(ra[n + (i << logN)], tt, p, ttpinv);
			QuickAccumMulAdd(acc, pHatnp[i], s);
		}
		QuickAccumEnd(acc);
		rem(a[n], a[n], pProd[np - 1]);
		if (a[n] > pProdh[np - 1]) a[n] -= pProd[np - 1];
		a[n] %= q;
	}

	delete[] ra;
}

