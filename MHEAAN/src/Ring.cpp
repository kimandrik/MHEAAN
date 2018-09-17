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


Ring::Ring(long logN0, long logN1, long logQ, double sigma, long h) :
		logN0(logN0), logN1(logN1), logQ(logQ), sigma(sigma), h(h) {

	N0 = (1 << logN0);
	M0 = 1 << (logN0 + 1);
	N0h = N0 >> 1;
	logN0h = logN0 - 1;

	N1 = 1 << logN1;
	M1 = N1 + 1;

	logN = logN0 + logN1;
	N = (1 << logN);
	Nh = N >> 1;

	long nprimes = ceil((2 + 2 * logN + 4 * logQ) / 59.0);
	multiplier = new RingMultiplier(logN0, logN1, nprimes);

	logQQ = 2 * logQ;

	Q = power2_ZZ(logQ);
	QQ = power2_ZZ(logQQ);

	gM0Pows = new uint64_t[N0h + 1];
	uint64_t g0 = 5;
	uint64_t g0Pow = 1;
	for (long i = 0; i < N0h; ++i) {
		gM0Pows[i] = g0Pow;
		g0Pow *= g0;
		g0Pow %= M0;
	}
	gM0Pows[N0h] = gM0Pows[0];

	gM1Pows = new uint64_t[M1];
	dftomegaPows = new complex<double>[N1]();
	omegaPows = new complex<double>[N1]();
	uint64_t g1 = multiplier->findPrimitiveRoot(M1);
	uint64_t g1Pow = 1;
	for (long i = 0; i < N1; ++i) {
		double angle = 2.0 * M_PI * g1Pow / M1;
		omegaPows[i].real(cos(angle));
		omegaPows[i].imag(sin(angle));
		dftomegaPows[i].real(cos(angle));
		dftomegaPows[i].imag(sin(angle));
		gM1Pows[i] = g1Pow;
		g1Pow *= g1;
		g1Pow %= M1;
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

void Ring::addBootContext(long logn0, long logn1, long logp) {
	if (bootContextMap.find({logn0, logn1}) == bootContextMap.end()) {
		long n0 = 1 << logn0;
		long logk0 = logn0 >> 1;
		long k0 = 1 << logk0;

		uint64_t** rpVec = new uint64_t*[n0];
		uint64_t** rpInvVec = new uint64_t*[n0];
		uint64_t* rp1 = NULL;
		uint64_t* rp2 = NULL;

		long* bndVec = new long[n0];
		long* bndInvVec = new long[n0];
		long bnd1 = 0;
		long bnd2 = 0;

		long np;
		complex<double>* pvals = new complex<double>[n0];
		ZZ* pVec = new ZZ[N0];

		long gap0 = N0h >> logn0;
		long deg;
		for (long ki = 0; ki < n0; ki += k0) {
			for (long pos = ki; pos < ki + k0; ++pos) {
				for (long i = 0; i < n0 - pos; ++i) {
					deg = ((M0 - gM0Pows[i + pos]) * i * gap0) % M0;
					pvals[i] = ksiM0Pows[deg];
				}
				for (long i = n0 - pos; i < n0; ++i) {
					deg = ((M0 - gM0Pows[i + pos - n0]) * i * gap0) % M0;
					pvals[i] = ksiM0Pows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pvals, n0, 1, ki, 0);
				IEMBX0(pvals, n0);
				for (long i = 0, jd = N0h, id = 0; i < n0; ++i, jd += gap0, id += gap0) {
					pVec[id] = EvaluatorUtils::scaleUpToZZ(pvals[i].real(), logp);
					pVec[jd] = EvaluatorUtils::scaleUpToZZ(pvals[i].imag(), logp);
				}
				bndVec[pos] = MaxBits(pVec, N0);
				np = ceil((logQ + bndVec[pos] + logN0 + 3)/59.0);
				rpVec[pos] = toNTTX0(pVec, np);
			}
		}

		for (long ki = 0; ki < n0; ki += k0) {
			for (long pos = ki; pos < ki + k0; ++pos) {
				for (long i = 0; i < n0 - pos; ++i) {
					deg = (gM0Pows[i] * (i + pos) * gap0) % M0;
					pvals[i] = ksiM0Pows[deg];
				}
				for (long i = n0 - pos; i < n0; ++i) {
					deg = (gM0Pows[i] * (i + pos - n0) * gap0) % M0;
					pvals[i] = ksiM0Pows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pvals, n0, 1, ki, 0);
				IEMBX0(pvals, n0);
				for (long i = 0, jd = N0h, id = 0; i < n0; ++i, jd += gap0, id += gap0) {
					pVec[id] = EvaluatorUtils::scaleUpToZZ(pvals[i].real(), logp);
					pVec[jd] = EvaluatorUtils::scaleUpToZZ(pvals[i].imag(), logp);
				}
				bndInvVec[pos] = MaxBits(pVec, N0);
				np = ceil((logQ + bndInvVec[pos] + logN0 + 3)/59.0);
				rpInvVec[pos] = toNTTX0(pVec, np);
			}
		}

		delete[] pvals;
		delete[] pVec;

		bootContextMap.insert(pair<pair<long, long>, BootContext*>({logn0, logn1}, new BootContext(rpVec, rpInvVec, rp1, rp2,
				bndVec, bndInvVec, bnd1, bnd2, logp)));
	}
}

void Ring::addSqrMatContext(long logn, long logp) {
	if (sqrMatContextMap.find(logn) == sqrMatContextMap.end()) {
		long n = 1 << logn;
		ZZ** pvec = new ZZ*[n];

		long gap = (N0h >> logn);
		long powsum = N0h - gap;
		pvec[0] = new ZZ[N];
		for (long i = 0; i < N0h; i += gap) {
			pvec[0][i + (powsum - i) * N0] = ZZ(1);
		}

		multByMonomialAndEqual(pvec[0], M0 - powsum, 0, Q);

		for (long i = 1; i < n; ++i) {
			pvec[i] = leftRotate(pvec[0], 0, i);
		}
		sqrMatContextMap.insert(pair<long, SqrMatContext*>(logn, new SqrMatContext(pvec, logp)));
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
	for (long i = 0; i < N1; ++i) {
		vals[i] /= N1;
	}
}

void Ring::EMBX0(complex<double>* vals, long n0) {
	arrayBitReverse(vals, n0);
	for (long len = 2; len <= n0; len <<= 1) {
		for (long i = 0; i < n0; i += len) {
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

void Ring::IEMBX0(complex<double>* vals, long n0) {
	for (long len = n0; len >= 1; len >>= 1) {
		for (long i = 0; i < n0; i += len) {
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
	arrayBitReverse(vals, n0);
	for (long i = 0; i < n0; ++i) {
		vals[i] /= n0;
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

void Ring::EMB(complex<double>* vals, long n0) {
	complex<double>* tmp = new complex<double>[N1]();
	for (long i = 0; i < N1; ++i) {
		EMBX0(vals + (i * n0), n0);
	}
	for (long i = 0; i < n0; ++i) {
		for (long j = 0; j < N1; ++j) {
			tmp[j] = vals[i + (j * n0)];
		}
		EMBX1(tmp);
		for (long j = 0; j < N1; ++j) {
			vals[i + (j * n0)] = tmp[j];
		}
	}
	delete[] tmp;
}

void Ring::IEMB(complex<double>* vals, long n0) {
	complex<double>* tmp = new complex<double>[N1]();

	for (long iy = 0; iy < N1; ++iy) {
		IEMBX0(vals + (iy * n0), n0);
	}
	for (long ix = 0; ix < n0; ++ix) {
		for (long iy = 0; iy < N1; ++iy) {
			tmp[iy] = vals[ix + (iy * n0)];
		}
		IEMBX1(tmp);
		for (long iy = 0; iy < N1; ++iy) {
			vals[ix + (iy * n0)] = tmp[iy];
		}
	}

	delete[] tmp;
}

ZZ* Ring::encode(complex<double>* vals, long n0, long n1, long logp) {
	long gap0 = N0h / n0;
	long gap1 = N1 / n1;
	ZZ* mx = new ZZ[N];
	complex<double>* uvals = new complex<double>[n0 * N1];
	for (long j = 0; j < gap1; ++j) {
		copy(vals, vals + n0 * n1, uvals + j * n0 * n1);
	}

	IEMB(uvals, n0);
	for (long i = 0, ii = N0h, ir = 0; i < n0; ++i, ii += gap0, ir += gap0) {
		for (long j = 0; j < N1; ++j) {
			mx[ir + N0 * j] = EvaluatorUtils::scaleUpToZZ(uvals[i + n0 * j].real(), logp);
			mx[ii + N0 * j] = EvaluatorUtils::scaleUpToZZ(uvals[i + n0 * j].imag(), logp);
		}
	}
	delete[] uvals;
	return mx;
}

ZZ* Ring::encode(double* vals, long n0, long n1, long logp) {
	long gap0 = N0h / n0;
	long gap1 = N1 / n1;
	ZZ* mx = new ZZ[N];
	complex<double>* uvals = new complex<double>[n0 * N1];
	for (long i = 0; i < n0 * n1; ++i) {
		for (long j = 0; j < gap1; ++j) {
			uvals[i + j * n0 * n1].real(vals[i]);
		}
	}

	IEMB(uvals, n0);
	for (long i = 0, ii = N0h, ir = 0; i < n0; ++i, ii += gap0, ir += gap0) {
		for (long j = 0; j < N1; ++j) {
			mx[ir + N0 * j] = EvaluatorUtils::scaleUpToZZ(uvals[i + n0 * j].real(), logp);
			mx[ii + N0 * j] = EvaluatorUtils::scaleUpToZZ(uvals[i + n0 * j].imag(), logp);
		}
	}
	delete[] uvals;
	return mx;
}

complex<double>* Ring::decode(ZZ* mxy, long n0, long n1, long logp, long logq) {
	ZZ q = qvec[logq];
	ZZ qh = qvec[logq - 1];
	long gap0 = N0h / n0;
	ZZ tmp;
	complex<double>* vals = new complex<double>[n0 * n1];
	complex<double>* fvals = new complex<double>[n0 * N1];
	for (long i0 = 0, ii0 = N0h, ir0 = 0; i0 < n0; ++i0, ii0 += gap0, ir0 += gap0) {
		for (long i1 = 0; i1 < N1; ++i1) {
			rem(tmp, mxy[ir0 + N0 * i1], q);
			while (tmp < 0) tmp += q;
			while (tmp > qh) tmp -= q;
			fvals[i0 + n0 * i1].real(EvaluatorUtils::scaleDownToReal(tmp, logp));

			rem(tmp, mxy[ii0 + N0 * i1], q);
			while(tmp < 0) tmp += q;
			while (tmp > qh) tmp -= q;
			fvals[i0 + n0 * i1].imag(EvaluatorUtils::scaleDownToReal(tmp, logp));
		}
	}

	EMB(fvals, n0);
	copy(fvals, fvals + n0 * n1, vals);
	delete[] fvals;
	return vals;
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
	return multiplier->toNTTX0(a, np);
}

uint64_t* Ring::toNTTX1(ZZ* a, long np) {
	return multiplier->toNTTX1(a, np);
}

uint64_t* Ring::toNTTX1Lazy(ZZ* a, long np) {
	return multiplier->toNTTX1Lazy(a, np);
}

uint64_t* Ring::toNTT(ZZ* a, long np) {
	return multiplier->toNTT(a, np);
}

uint64_t* Ring::toNTTLazy(ZZ* a, long np) {
	return multiplier->toNTTLazy(a, np);
}

void Ring::addNTTAndEqual(uint64_t* ra, uint64_t* rb, long np) {
	multiplier->addNTTAndEqual(ra, rb, np);
}

void Ring::multX0(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier->multX0(x, a, b, np, q);
}

void Ring::multX0AndEqual(ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier->multX0AndEqual(a, b, np, q);
}

void Ring::multNTTX0(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier->multNTTX0(x, a, rb, np, q);
}

void Ring::multNTTX0AndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier->multNTTX0AndEqual(a, rb, np, q);
}

void Ring::multDNTTX0(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q) {
	multiplier->multDNTTX0(x, ra, rb, np, q);
}

void Ring::multX1(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier->multX1(x, a, b, np, q);
}

void Ring::multX1AndEqual(ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier->multX1AndEqual(a, b, np, q);
}

void Ring::multNTTX1(ZZ* x, ZZ* a, uint64_t* b, long np, ZZ& q) {
	multiplier->multNTTX1(x, a, b, np, q);
}

void Ring::multNTTX1AndEqual(ZZ* a, uint64_t* b, long np, ZZ& q) {
	multiplier->multNTTX1AndEqual(a, b, np, q);
}

void Ring::mult(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier->mult(x, a, b, np, q);
}

void Ring::multAndEqual(ZZ* a, ZZ* b, long np, ZZ& q) {
	multiplier->multAndEqual(a, b, np, q);
}

void Ring::multNTT(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier->multNTT(x, a, rb, np, q);
}

void Ring::multNTTAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier->multNTTAndEqual(a, rb, np, q);
}

void Ring::multNTTLazy(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier->multNTTLazy(x, a, rb, np, q);
}

void Ring::multNTTLazyAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& q) {
	multiplier->multNTTLazyAndEqual(a, rb, np, q);
}

void Ring::multDNTT(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& q) {
	multiplier->multDNTT(x, ra, rb, np, q);
}

void Ring::square(ZZ* x, ZZ* a, long np, ZZ& q) {
	multiplier->square(x, a, np, q);
}

void Ring::squareAndEqual(ZZ* a, long np, ZZ& q) {
	multiplier->squareAndEqual(a, np, q);
}

void Ring::squareNTT(ZZ* x, uint64_t* ra, long np, ZZ& q) {
	multiplier->squareNTT(x, ra, np, q);
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

ZZ* Ring::multByMonomial(ZZ* p, long deg0, long deg1, ZZ& q) {
	ZZ* res = new ZZ[N];

	for (long i = 0; i < N0; ++i) {
		for (long j = 0; j < N1; ++j) {
			long resdeg0 = (deg0 + i) % M0;
			long resdeg1 = (deg1 + j) % M1;

			if(resdeg0 < N0) {
				if(resdeg1 != N1) {
					AddMod(res[resdeg0 + (resdeg1 << logN0)], res[resdeg0 + (resdeg1 << logN0)], p[i + (j << logN0)], q);
				} else {
					for (long k = 0; k < N1; ++k) {
						AddMod(res[resdeg0 + (k << logN0)], res[resdeg0 + (k << logN0)], -p[i + (j << logN0)], q);
					}
				}
			} else {
				if(resdeg1 != N1) {
					AddMod(res[(resdeg0 - N0) + (resdeg1 << logN0)], res[(resdeg0 - N0) + (resdeg1 << logN0)], -p[i + (j << logN0)], q);
				} else {
					for (long k = 0; k < N1; ++k) {
						AddMod(res[(resdeg0 - N0) + (k << logN0)], res[(resdeg0 - N0) + (k << logN0)], p[i + (j << logN0)], q);
					}
				}
			}
		}
	}
	return res;
}

void Ring::multByMonomialAndEqual(ZZ* p, long deg0, long deg1, ZZ& q) {
	ZZ* res = multByMonomial(p, deg0, deg1, q);
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


ZZ* Ring::leftRotate(ZZ* p, long r0, long r1) {
	ZZ* tmp = new ZZ[N];
	long deg0 = gM0Pows[r0];
	for (long j = 0; j < N; j += N0) {
		for (long i = 0; i < N0; ++i) {
			long ipow = i * deg0;
			long shift = ipow % M0;
			if (shift < N0) {
				tmp[shift + j] = p[i + j];
			} else {
				tmp[shift - N0 + j] = -p[i + j];
			}
		}
	}

	ZZ* res = new ZZ[N];
	long degy = gM1Pows[r1];
	for (long i = 0; i < N0; ++i) {
		for (long j = 0; j < N1; ++j) {
			long jpow = j * degy;
			long shift = jpow % M1;
			if (shift < N1) {
				res[i + (shift << logN0)] += tmp[i + (j << logN0)];
			} else {
				for (long k = 0; k < N1; ++k) {
					res[i + (k << logN0)] -= tmp[i + (j << logN0)];
				}
			}
		}
	}
	delete[] tmp;
	return res;
}

ZZ* Ring::conjugate(ZZ* p) {
	ZZ* tmp = new ZZ[N];
	for (long j = 0; j < N; j += N0) {
		tmp[j] = p[j];
		for (long i = 1; i < N0; ++i) {
			tmp[N0 - i + j] = -p[i + j];
		}
	}

	ZZ* res = new ZZ[N];
	for (long i = 0; i < N0; ++i) {
		res[i] += tmp[i];
		for (long j = 2 * N0; j < N; j += N0) {
			res[i + (N - j + N0)] += tmp[i + j];
		}
		for (long j = 0; j < N; j += N0) {
			res[i + j] -= tmp[i + N0];
		}
	}

	delete[] tmp;
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
