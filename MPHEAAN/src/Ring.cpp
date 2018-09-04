/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ring.h"

#include <NTL/BasicThreadPool.h>
#include "EvaluatorUtils.h"
#include "StringUtils.h"


Ring::Ring(long logN0, long logQ, double sigma, long h) :
		logN0(logN0), logQ(logQ), sigma(sigma), h(h) {

	N0 = (1 << logN0);
	M0 = 1 << (logN0 + 1);
	N0h = N0 >> 1;
	logN0h = logN0 - 1;

	logN1 = 8;
	N1 = 1 << logN1;
	M1 = N1 + 1;

	logN = logN0 + logN1;
	N = (1 << logN);
	Nh = N >> 1;

	long nprimes = ceil((2 + 2 * logN + 4 * logQ) / 59.0);
	multiplier = RingMultiplier(logN0, nprimes);

	logQQ = 2 * logQ;

	Q = power2_ZZ(logQ);
	QQ = power2_ZZ(logQQ);

	gM0Pows = new long[N0h + 1];
	long gx = 5;
	long gxPow = 1;
	for (long i = 0; i < N0h; ++i) {
		gM0Pows[i] = gxPow;
		gxPow *= gx;
		gxPow %= M0;
	}
	gM0Pows[N0h] = gM0Pows[0];

	gM1Pows = new long[M1];
	dftomegaPows = new complex<double>[N1]();
	omegaPows = new complex<double>[N1]();
	long gy = 3;
	long gyPow = 1;
	for (long i = 0; i < N1; ++i) {
		double angle = 2.0 * M_PI * gyPow / M1;
		omegaPows[i].real(cos(angle));
		omegaPows[i].imag(sin(angle));
		dftomegaPows[i].real(cos(angle));
		dftomegaPows[i].imag(sin(angle));
		gM1Pows[i] = gyPow;
		gyPow *= gy;
		gyPow %= M1;
	}
	gM1Pows[N1] = gM1Pows[0];

	ksiM0Pows = new complex<double>[M0 + 1]();
	for (long j = 0; j < M0; ++j) {
		double angle = 2.0 * M_PI * j / M0;
		ksiM0Pows[j].real(cos(angle));
		ksiM0Pows[j].imag(sin(angle));
	}
	ksiM0Pows[M0] = ksiM0Pows[0];

	ksiN1Pows = new complex<double>[M1]();
	for (long j = 0; j < N1; ++j) {
		double angle = 2.0 * M_PI * j / N1;
		ksiN1Pows[j].real(cos(angle));
		ksiN1Pows[j].imag(sin(angle));
	}
	ksiN1Pows[N1] = ksiN1Pows[0];

	ksiM1Pows = new complex<double>[M1]();
	for (long j = 0; j < M1; ++j) {
		double angle = 2.0 * M_PI * j / M1;
		ksiM1Pows[j].real(cos(angle));
		ksiM1Pows[j].imag(sin(angle));
	}

	arrayBitReverse(dftomegaPows, N1);
	for (long len = 2; len <= N1; len <<= 1) {
		long lenh = len >> 1;
		long gap = N1 / len;
		for (long i = 0; i < N1; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = j * gap;
				complex<double> u = dftomegaPows[i + j];
				complex<double> v = dftomegaPows[i + j + lenh];
				v *= ksiN1Pows[idx];
				dftomegaPows[i + j] = u + v;
				dftomegaPows[i + j + lenh] = u - v;
			}
		}
	}

	qvec = new ZZ[logQQ + 1];
	qvec[0] = ZZ(1);
	for (long i = 1; i < logQQ + 1; ++i) {
		qvec[i] = qvec[i-1] << 1;
	}
}

void Ring::addBootContext(long lognx, long logny, long logp) {
	if (bootContextMap.find({lognx, logny}) == bootContextMap.end()) {
		long nx = 1 << lognx;
		long logkx = lognx >> 1;
		long kx = 1 << logkx;

		uint64_t** rpxVec = new uint64_t*[nx];
		uint64_t** rpxInvVec = new uint64_t*[nx];

		complex<double>* pxvals = new complex<double>[nx];
		uint64_t* rp1 = NULL;
		uint64_t* rp2 = NULL;

		long np0 = ceil((logQ + logp + logN0 + 3)/59.0);

		ZZ* pxVec = new ZZ[N0];

		long gapx = N0h >> lognx;
		long deg;
		for (long kxi = 0; kxi < nx; kxi += kx) {
			for (long pos = kxi; pos < kxi + kx; ++pos) {
				for (long i = 0; i < nx - pos; ++i) {
					deg = ((M0 - gM0Pows[i + pos]) * i * gapx) % M0;
					pxvals[i] = ksiM0Pows[deg];
				}
				for (long i = nx - pos; i < nx; ++i) {
					deg = ((M0 - gM0Pows[i + pos - nx]) * i * gapx) % M0;
					pxvals[i] = ksiM0Pows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pxvals, nx, 1, kxi, 0);
				IEMBX0(pxvals, nx);
				for (long ix = 0, jdx = N0h, idx = 0; ix < nx; ++ix, jdx += gapx, idx += gapx) {
					pxVec[idx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].real(), logp);
					pxVec[jdx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].imag(), logp);
				}
				rpxVec[pos] = toNTTX0(pxVec, np0);
			}
		}

		for (long kxi = 0; kxi < nx; kxi += kx) {
			for (long pos = kxi; pos < kxi + kx; ++pos) {
				for (long i = 0; i < nx - pos; ++i) {
					deg = (gM0Pows[i] * (i + pos) * gapx) % M0;
					pxvals[i] = ksiM0Pows[deg];
				}
				for (long i = nx - pos; i < nx; ++i) {
					deg = (gM0Pows[i] * (i + pos - nx) * gapx) % M0;
					pxvals[i] = ksiM0Pows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pxvals, nx, 1, kxi, 0);
				IEMBX0(pxvals, nx);
				for (long ix = 0, jdx = N0h, idx = 0; ix < nx; ++ix, jdx += gapx, idx += gapx) {
					pxVec[idx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].real(), logp);
					pxVec[jdx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].imag(), logp);
				}
				rpxInvVec[pos] = toNTTX0(pxVec, np0);
			}
		}

		delete[] pxvals;
		delete[] pxVec;

		bootContextMap.insert(pair<pair<long, long>, BootContext>({lognx, logny}, BootContext(rpxVec, rpxInvVec, rp1, rp2, logp)));
	}
}

void Ring::addMatrixContext(long lognx) {
	if (matrixContext.find({lognx, lognx}) == matrixContext.end()) {
		long nx = 1 << lognx;
		ZZ** pvec = new ZZ*[nx];

		long gap = (N0h >> lognx);
		long powsum = N0h - gap;
		pvec[0] = new ZZ[N];
		for (long ix = 0; ix < N0h; ix += gap) {
			pvec[0][ix + (powsum - ix) * N0] = ZZ(1);
		}

		multByMonomialAndEqual(pvec[0], M0 - powsum, 0, Q);

		for (long ix = 1; ix < nx; ++ix) {
			pvec[ix] = leftRotate(pvec[0], 0, ix);
		}
		matrixContext.insert(pair<pair<long, long>, MatrixContext>({lognx, lognx}, MatrixContext(pvec)));
	}
}

void Ring::arrayBitReverse(complex<double>* vals, long n) {
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

void Ring::DFTX1(complex<double>* vals) {
	arrayBitReverse(vals, N1);
	for (long len = 2; len <= N1; len <<= 1) {
		long lenh = len >> 1;
		long gap = N1 / len;
		for (long i = 0; i < N1; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = j * gap;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiN1Pows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Ring::IDFTX1(complex<double>* vals) {
	arrayBitReverse(vals, N1);
	for (long len = 2; len <= N1; len <<= 1) {
		long lenh = len >> 1;
		long gap = N1 / len;
		for (long i = 0; i < N1; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = N1 - (j * gap);
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiN1Pows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
	for (long iy = 0; iy < N1; ++iy) {
		vals[iy] /= N1;
	}
}

void Ring::EMBX0(complex<double>* vals, long nx) {
	arrayBitReverse(vals, nx);
	for (long len = 2; len <= nx; len <<= 1) {
		for (long i = 0; i < nx; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			long gap = M0 / lenq;
			for (long j = 0; j < lenh; ++j) {
				long idx = (gM0Pows[j] % lenq) * gap;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiM0Pows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Ring::IEMBX0(complex<double>* vals, long nx) {
	for (long len = nx; len >= 1; len >>= 1) {
		for (long i = 0; i < nx; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			long gap = M0 / lenq;
			for (long j = 0; j < lenh; ++j) {
				long idx = (lenq - (gM0Pows[j] % lenq)) * gap;
				complex<double> u = vals[i + j] + vals[i + j + lenh];
				complex<double> v = vals[i + j] - vals[i + j + lenh];
				v *= ksiM0Pows[idx];
				vals[i + j] = u;
				vals[i + j + lenh] = v;
			}
		}
	}
	arrayBitReverse(vals, nx);
	for (long i = 0; i < nx; ++i) {
		vals[i] /= nx;
	}
}

void Ring::EMBX1(complex<double>* vals) {
	complex<double>* tmp = new complex<double>[N1];
	for (long i = 0; i < N1; ++i) {
		tmp[i] = vals[gM1Pows[N1 - i] - 1];
	}
	for (long i = 0; i < N1; ++i) {
		vals[i] = tmp[i];
	}
	delete[] tmp;

	DFTX1(vals);
	for (long i = 0; i < N1; ++i) {
		vals[i] *= dftomegaPows[i];
	}
	IDFTX1(vals);
	for (long i = 0; i < N1; ++i) {
		vals[i] /= omegaPows[i];
	}
}

void Ring::IEMBX1(complex<double>* vals) {
	for (long i = 0; i < N1; ++i) {
		vals[i] *= omegaPows[i];
	}
	DFTX1(vals);
	for (long i = 0; i < N1; ++i) {
		vals[i] /= dftomegaPows[i];
	}
	IDFTX1(vals);

	complex<double>* tmp = new complex<double>[N1]();
	for (long i = 0; i < N1; ++i) {
		tmp[gM1Pows[N1 - i] - 1] = vals[i];
	}
	for (long i = 0; i < N1; ++i) {
		vals[i] = tmp[i];
	}
	delete[] tmp;
}

void Ring::EMB(complex<double>* vals, long nx) {
	complex<double>* tmp = new complex<double>[N1]();
	for (long iy = 0; iy < N1; ++iy) {
		EMBX0(vals + (iy * nx), nx);
	}
	for (long ix = 0; ix < nx; ++ix) {
		for (long iy = 0; iy < N1; ++iy) {
			tmp[iy] = vals[ix + (iy * nx)];
		}
		EMBX1(tmp);
		for (long iy = 0; iy < N1; ++iy) {
			vals[ix + (iy * nx)] = tmp[iy];
		}
	}
	delete[] tmp;
}

void Ring::IEMB(complex<double>* vals, long nx) {
	complex<double>* tmp = new complex<double>[N1]();

	for (long iy = 0; iy < N1; ++iy) {
		IEMBX0(vals + (iy * nx), nx);
	}
	for (long ix = 0; ix < nx; ++ix) {
		for (long iy = 0; iy < N1; ++iy) {
			tmp[iy] = vals[ix + (iy * nx)];
		}
		IEMBX1(tmp);
		for (long iy = 0; iy < N1; ++iy) {
			vals[ix + (iy * nx)] = tmp[iy];
		}
	}

	delete[] tmp;
}

void Ring::encode(ZZ* mx, complex<double>* vals, long nx, long ny, long logp) {
	long gapx = N0h / nx;
	long gapy = N1 / ny;

	complex<double>* uvals = new complex<double>[nx * N1];
	for (long j = 0; j < gapy; ++j) {
		copy(vals, vals + nx * ny, uvals + j * nx * ny);
	}

	IEMB(uvals, nx);
	for (long ix = 0, iix = N0h, irx = 0; ix < nx; ++ix, iix += gapx, irx += gapx) {
		for (long iy = 0; iy < N1; ++iy) {
			mx[irx + N0 * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].real(), logp);
			mx[iix + N0 * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].imag(), logp);
		}
	}
	delete[] uvals;
}

void Ring::encode(ZZ* mx, double* vals, long nx, long ny, long logp) {
	long gapx = N0h / nx;
	long gapy = N1 / ny;

	complex<double>* uvals = new complex<double>[nx * N1];
	for (long i = 0; i < nx * ny; ++i) {
		for (long j = 0; j < gapy; ++j) {
			uvals[i + j * nx * ny].real(vals[i]);
		}
	}

	IEMB(uvals, nx);
	for (long ix = 0, iix = N0h, irx = 0; ix < nx; ++ix, iix += gapx, irx += gapx) {
		for (long iy = 0; iy < N1; ++iy) {
			mx[irx + N0 * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].real(), logp);
			mx[iix + N0 * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].imag(), logp);
		}
	}
	delete[] uvals;
}

void Ring::decode(ZZ* mxy, complex<double>* vals, long nx, long ny, long logp, long logq) {
	ZZ q = qvec[logq];
	ZZ qh = qvec[logq - 1];
	long gapx = N0h / nx;
	ZZ tmp;

	complex<double>* fvals = new complex<double>[nx * N1];
	for (long ix = 0, iix = N0h, irx = 0; ix < nx; ++ix, iix += gapx, irx += gapx) {
		for (long iy = 0; iy < N1; ++iy) {
			rem(tmp, mxy[irx + N0 * iy], q);
			while (tmp < 0) tmp += q;
			while (tmp > qh) tmp -= q;
			fvals[ix + nx * iy].real(EvaluatorUtils::scaleDownToReal(tmp, logp));

			rem(tmp, mxy[iix + N0 * iy], q);
			while(tmp < 0) tmp += q;
			while (tmp > qh) tmp -= q;
			fvals[ix + nx * iy].imag(EvaluatorUtils::scaleDownToReal(tmp, logp));
		}
	}

	EMB(fvals, nx);
	copy(fvals, fvals + nx * ny, vals);
	delete[] fvals;
}


//----------------------------------------------------------------------------------
//   MULTIPLICATION
//----------------------------------------------------------------------------------


long Ring::MaxBits(const ZZ* f, long n) {
   long i, m;
   m = 0;

   for (i = 0; i < n; i++) {
      m = max(m, NumBits(f[i]));
   }
   return m;
}

uint64_t* Ring::toNTTX0(ZZ* a, long np) {
	return multiplier.toNTTX0(a, np);
}

uint64_t* Ring::toNTTX1(ZZ* a, long np) {
	return multiplier.toNTTX1(a, np);
}

uint64_t* Ring::toNTTX1Lazy(ZZ* a, long np) {
	return multiplier.toNTTX1Lazy(a, np);
}

uint64_t* Ring::toNTT(ZZ* a, long np) {
	return multiplier.toNTT(a, np);
}

uint64_t* Ring::toNTTLazy(ZZ* a, long np) {
	return multiplier.toNTTLazy(a, np);
}

uint64_t* Ring::addNTT(uint64_t* ra, uint64_t* rb, long np) {
	return multiplier.addNTT(ra, rb, np);
}

void Ring::multX0(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier.multX0(x, a, b, np, q);
}

void Ring::multX0AndEqual(ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier.multX0AndEqual(a, b, np, q);
}

void Ring::multNTTX0(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier.multNTTX0(x, a, rb, np, q);
}

void Ring::multNTTX0AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier.multNTTX0AndEqual(a, rb, np, q);
}

void Ring::multNTTX0D(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q) {
	multiplier.multDNTTX0(x, ra, rb, np, q);
}

void Ring::multX1(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier.multX1(x, a, b, np, q);
}

void Ring::multX1AndEqual(ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier.multX1AndEqual(a, b, np, q);
}

void Ring::multNTTX1(ZZ* x, ZZ* a, uint64_t* b, long np, ZZ& q) {
	multiplier.multNTTX1(x, a, b, np, q);
}

void Ring::multNTTX1AndEqual(ZZ* a, uint64_t* b, long np, ZZ& q) {
	multiplier.multNTTX1AndEqual(a, b, np, q);
}

void Ring::mult(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier.mult(x, a, b, np, q);
}

void Ring::multAndEqual(ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier.multAndEqual(a, b, np, q);
}

void Ring::multNTT(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier.multNTT(x, a, rb, np, q);
}

void Ring::multNTTAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier.multNTTAndEqual(a, rb, np, q);
}

void Ring::multNTTLazy(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier.multNTTLazy(x, a, rb, np, q);
}

void Ring::multNTTLazyAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier.multNTTLazyAndEqual(a, rb, np, q);
}

void Ring::multDNTT(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q) {
	multiplier.multDNTT(x, ra, rb, np, q);
}

void Ring::square(ZZ* x, ZZ* a, long np, ZZ& q) {
	multiplier.square(x, a, np, q);
}

void Ring::squareAndEqual(ZZ* a, long np, ZZ& q) {
	multiplier.squareAndEqual(a, np, q);
}

void Ring::squareNTT(ZZ* x, uint64_t* ra, long np, ZZ& q) {
	multiplier.squareNTT(x, ra, np, q);
}

//----------------------------------------------------------------------------------
//   OTHER
//----------------------------------------------------------------------------------


void Ring::mod(ZZ* res, ZZ* p, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		rem(res[i], p[i], q);
	}
}

void Ring::modAndEqual(ZZ* p, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		rem(p[i], p[i], q);
	}
}

void Ring::negate(ZZ* res, ZZ* p) {
	for (long i = 0; i < N; ++i) {
		res[i] = -p[i];
	}
}

void Ring::negateAndEqual(ZZ* p) {
	for (long i = 0; i < N; ++i) {
		p[i] = -p[i];
	}
}

void Ring::add(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(res[i], p1[i], p2[i], q);
	}
}

void Ring::addAndEqual(ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(p1[i], p1[i], p2[i], q);
	}
}

void Ring::sub(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(res[i], p1[i], -p2[i], q);
	}
}

void Ring::subAndEqual(ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(p1[i], p1[i], -p2[i], q);
	}
}

void Ring::subAndEqual2(ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(p2[i], p1[i], -p2[i], q);
	}
}

ZZ* Ring::multByMonomial(ZZ* p, long degx, long degy, ZZ& q) {
	ZZ* res = new ZZ[N];

	for (long ix = 0; ix < N0; ++ix) {
		for (long iy = 0; iy < N1; ++iy) {
			long ixres = (degx + ix) % M0;
			long iyres = (degy + iy) % M1;

			if(ixres < N0) {
				if(iyres != N1) {
					AddMod(res[ixres + (iyres << logN0)], res[ixres + (iyres << logN0)], p[ix + (iy << logN0)], q);
				} else {
					for (long j = 0; j < N1; ++j) {
						AddMod(res[ixres + (j << logN0)], res[ixres + (j << logN0)], -p[ix + (iy << logN0)], q);
					}
				}
			} else {
				if(iyres != N1) {
					AddMod(res[(ixres - N0) + (iyres << logN0)], res[(ixres - N0) + (iyres << logN0)], -p[ix + (iy << logN0)], q);
				} else {
					for (long j = 0; j < N1; ++j) {
						AddMod(res[(ixres - N0) + (j << logN0)], res[(ixres - N0) + (j << logN0)], p[ix + (iy << logN0)], q);
					}
				}
			}
		}
	}
	return res;
}

void Ring::multByMonomialAndEqual(ZZ* p, long degx, long degy, ZZ& q) {
	ZZ* res = multByMonomial(p, degx, degy, q);
	for (long i = 0; i < N; ++i) {
		p[i] = res[i];
	}
	delete[] res;
}

void Ring::multByConst(ZZ* res, ZZ* p, ZZ& cnst, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		MulMod(res[i], p[i], cnst, q);
	}
}

void Ring::multByConstAndEqual(ZZ* p, ZZ& cnst, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		MulMod(p[i], p[i], cnst, q);
	}
}


//----------------------------------------------------------------------------------
//   SHIFTING
//----------------------------------------------------------------------------------


void Ring::leftShift(ZZ* res, ZZ* p, long bits, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		res[i] = p[i] << bits;
		res[i] %= q;
	}
}

void Ring::leftShiftAndEqual(ZZ* p, long bits, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		p[i] <<= bits;
		p[i] %= q;
	}
}

void Ring::doubleAndEqual(ZZ* p, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		p[i] <<= 1;
		p[i] %= q;
	}
}

void Ring::rightShift(ZZ* res, ZZ* p, long bits) {
	for (long i = 0; i < N; ++i) {
		res[i] = p[i] >> bits;
	}
}

void Ring::rightShiftAndEqual(ZZ* p, long bits) {
	for (long i = 0; i < N; ++i) {
		p[i] >>= bits;
	}
}


//----------------------------------------------------------------------------------
//   ROTATION & CONJUGATION & TRANSPOSITION
//----------------------------------------------------------------------------------


ZZ* Ring::leftRotate(ZZ* p, long rx, long ry) {

	ZZ* xxx = new ZZ[N];
	long degx = gM0Pows[rx];
	for (long j = 0; j < N; j += N0) {
		for (long ix = 0; ix < N0; ++ix) {
			long ipow = ix * degx;
			long shift = ipow % M0;
			if (shift < N0) {
				xxx[shift + j] = p[ix + j];
			} else {
				xxx[shift - N0 + j] = -p[ix + j];
			}
		}
	}

	ZZ* res = new ZZ[N];

	long degy = gM1Pows[ry];
	for (long ix = 0; ix < N0; ++ix) {
		for (long iy = 0; iy < N1; ++iy) {
			long ipow = iy * degy;
			long shift = ipow % M1;
			if (shift < N1) {
				res[ix + (shift << logN0)] += xxx[ix + (iy << logN0)];
			} else {
				for (long t = 0; t < N1; ++t) {
					res[ix + (t << logN0)] -= xxx[ix + (iy << logN0)];
				}
			}
		}
	}
	delete[] xxx;
	return res;
}

ZZ* Ring::conjugate(ZZ* p) {
	ZZ* xxx = new ZZ[N];

	for (long j = 0; j < N; j += N0) {
		xxx[j] = p[j];
		for (long ix = 1; ix < N0; ++ix) {
			xxx[N0 - ix + j] = -p[ix + j];
		}
	}
	ZZ* res = new ZZ[N];

	for (long ix = 0; ix < N0; ++ix) {
		res[ix] += xxx[ix];
		for (long j = 2 * N0; j < N; j += N0) {
			res[ix + (N - j + N0)] += xxx[ix + j];
		}
		for (long j = 0; j < N; j += N0) {
			res[ix + j] -= xxx[ix + N0];
		}
	}

	delete[] xxx;
	return res;
}


//----------------------------------------------------------------------------------
//   SAMPLING
//----------------------------------------------------------------------------------


void Ring::sampleGauss(ZZ* res) {
	static double Pi = 4.0 * atan(1.0);
	static long const bignum = 0xfffffff;

	for (long i = 0; i < N; i+=2) {
		double r1 = (1 + RandomBnd(bignum)) / ((double)bignum + 1);
		double r2 = (1 + RandomBnd(bignum)) / ((double)bignum + 1);
		double theta=2 * Pi * r1;
		double rr= sqrt(-2.0 * log(r2)) * sigma;

		res[i] = (long) floor(rr * cos(theta) + 0.5);
		res[i + 1] = (long) floor(rr * sin(theta) + 0.5);
	}
}

void Ring::sampleHWT(ZZ* res) {
	long idx = 0;
	ZZ tmp = RandomBits_ZZ(h);
	while(idx < h) {
		long i = RandomBnd(N);
		if(res[i] == 0) {
			res[i] = (bit(tmp, idx) == 0) ? ZZ(1) : ZZ(-1);
			idx++;
		}
	}
}

void Ring::sampleZO(ZZ* res) {
	ZZ tmp = RandomBits_ZZ((N << 1));
	for (long i = 0; i < N; ++i) {
		res[i] = (bit(tmp, 2 * i) == 0) ? ZZ(0) : (bit(tmp, 2 * i + 1) == 0) ? ZZ(1) : ZZ(-1);
	}
}

void Ring::sampleUniform(ZZ* res, long bits) {
	for (long i = 0; i < N; i++) {
		res[i] = RandomBits_ZZ(bits);
	}
}
