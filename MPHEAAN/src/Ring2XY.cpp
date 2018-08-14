/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include <NTL/BasicThreadPool.h>
#include "Ring2XY.h"
#include "EvaluatorUtils.h"
#include "StringUtils.h"


Ring2XY::Ring2XY(long logNx, long logQ, double sigma, long h) :
		logNx(logNx), logQ(logQ), sigma(sigma), h(h), multiplier(logNx, logQ) {

	Nx = (1 << logNx);
	Mx = 1 << (logNx + 1);
	Nxh = Nx >> 1;
	logNxh = logNx - 1;

	logNy = 8;
	Ny = 1 << logNy;
	My = Ny + 1;

	logN = logNx + logNy;
	N = (1 << logN);
	Nh = N >> 1;

	logQQ = 2 * logQ;

	Q = power2_ZZ(logQ);
	QQ = power2_ZZ(logQQ);

	gxPows = new long[Nxh + 1];
	long gx = 5;
	long gxPow = 1;
	for (long i = 0; i < Nxh; ++i) {
		gxPows[i] = gxPow;
		gxPow *= gx;
		gxPow %= Mx;
	}
	gxPows[Nxh] = gxPows[0];

	gyPows = new long[My];
	dftomegaPows = new complex<double>[Ny]();
	omegaPows = new complex<double>[Ny]();
	long gy = 3;
	long gyPow = 1;
	for (long i = 0; i < Ny; ++i) {
		double angle = 2.0 * M_PI * gyPow / My;
		omegaPows[i].real(cos(angle));
		omegaPows[i].imag(sin(angle));
		dftomegaPows[i].real(cos(angle));
		dftomegaPows[i].imag(sin(angle));
		gyPows[i] = gyPow;
		gyPow *= gy;
		gyPow %= My;
	}
	gyPows[Ny] = 1;

	ksixPows = new complex<double>[Mx + 1]();
	for (long j = 0; j < Mx; ++j) {
		double angle = 2.0 * M_PI * j / Mx;
		ksixPows[j].real(cos(angle));
		ksixPows[j].imag(sin(angle));
	}
	ksixPows[Mx] = ksixPows[0];

	ksiyPows = new complex<double>[Ny + 1]();
	for (long j = 0; j < Ny; ++j) {
		double angle = 2.0 * M_PI * j / Ny;
		ksiyPows[j].real(cos(angle));
		ksiyPows[j].imag(sin(angle));
	}
	ksiyPows[Ny] = ksiyPows[0];

	arrayBitReverse(dftomegaPows, Ny);
	for (long len = 2; len <= Ny; len <<= 1) {
		long lenh = len >> 1;
		long gap = Ny / len;
		for (long i = 0; i < Ny; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = j * gap;
				complex<double> u = dftomegaPows[i + j];
				complex<double> v = dftomegaPows[i + j + lenh];
				v *= ksiyPows[idx];
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

void Ring2XY::addBootContext(long lognx, long logny, long logp) {
	if (bootContextMap.find({lognx, logny}) == bootContextMap.end()) {
		long nx = 1 << lognx;
		long ny = 1 << logny;

		long logkx = lognx >> 1;
		long logky = logny >> 1;

		long kx = 1 << logkx;
		long ky = 1 << logky;

		ZZ** pxVec = new ZZ*[nx];
		ZZ** pyrVec = new ZZ*[ny];
		ZZ** pyiVec = new ZZ*[ny];

		ZZ** pxInvVec = new ZZ*[nx];
		ZZ** pyrInvVec = new ZZ*[ny];
		ZZ** pyiInvVec = new ZZ*[ny];

		for (long ix = 0; ix < nx; ++ix) {
			pxVec[ix] = new ZZ[Nx];
			pxInvVec[ix] = new ZZ[Nx];
		}
		for (long iy = 0; iy < ny; ++iy) {
			pyrVec[iy] = new ZZ[Ny];
			pyiVec[iy] = new ZZ[Ny];
			pyrInvVec[iy] = new ZZ[Ny];
			pyiInvVec[iy] = new ZZ[Ny];
		}

		complex<double>* pxvals = new complex<double>[nx];
		complex<double>* pyvals = new complex<double>[ny];

		ZZ* p1 = new ZZ[N];
		ZZ* p2 = new ZZ[N];

		long gapx = Nxh >> lognx;
		long deg;
		for (long kxi = 0; kxi < nx; kxi += kx) {
			for (long pos = kxi; pos < kxi + kx; ++pos) {
				for (long i = 0; i < nx - pos; ++i) {
					deg = ((Mx - gxPows[i + pos]) * i * gapx) % Mx;
					pxvals[i] = ksixPows[deg];
				}
				for (long i = nx - pos; i < nx; ++i) {
					deg = ((Mx - gxPows[i + pos - nx]) * i * gapx) % Mx;
					pxvals[i] = ksixPows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pxvals, nx, 1, kxi, 0);
				IEMBX(pxvals, nx);
				for (long ix = 0, jdx = Nxh, idx = 0; ix < nx; ++ix, jdx += gapx, idx += gapx) {
					pxVec[pos][idx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].real(), logp);
					pxVec[pos][jdx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].imag(), logp);
				}
			}
		}

		for (long kxi = 0; kxi < nx; kxi += kx) {
			for (long pos = kxi; pos < kxi + kx; ++pos) {
				for (long i = 0; i < nx - pos; ++i) {
					deg = (gxPows[i] * (i + pos) * gapx) % Mx;
					pxvals[i] = ksixPows[deg];
				}
				for (long i = nx - pos; i < nx; ++i) {
					deg = (gxPows[i] * (i + pos - nx) * gapx) % Mx;
					pxvals[i] = ksixPows[deg];
				}
				EvaluatorUtils::rightRotateAndEqual(pxvals, nx, 1, kxi, 0);
				IEMBX(pxvals, nx);
				for (long ix = 0, jdx = Nxh, idx = 0; ix < nx; ++ix, jdx += gapx, idx += gapx) {
					pxInvVec[pos][idx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].real(), logp);
					pxInvVec[pos][jdx] = EvaluatorUtils::scaleUpToZZ(pxvals[ix].imag(), logp);
				}
			}
		}

		for (long kyi = 0; kyi < ny; kyi += ky) {
			for (long pos = kyi; pos < kyi + ky; ++pos) {
				for (long i = 0; i < ny - pos; ++i) {
					deg = ((My - gyPows[i + pos]) * i) % My;
					pyvals[i] = (ksiyPows[deg] - ksiyPows[gyPows[i + pos + 1]]) * (double)Ny / (double)My;
				}
				for (long i = ny - pos; i < ny; ++i) {
					deg = ((My - gyPows[i + pos - ny]) * i) % My;
					pyvals[i] = (ksiyPows[deg] - ksiyPows[gyPows[i + pos + 1 - ny]]) * (double)Ny / (double)My;
				}

				EvaluatorUtils::rightRotateAndEqual(pyvals, 1, ny, 0, kyi);

				IEMBY(pyvals);
				for (long iy = 0; iy < ny; ++iy) {
					pyrVec[pos][iy] = EvaluatorUtils::scaleUpToZZ(pyvals[iy].real(), logp);
					pyiVec[pos][iy] = EvaluatorUtils::scaleUpToZZ(pyvals[iy].imag(), logp);
				}
			}
		}

		for (long kyi = 0; kyi < ny; kyi += ky) {
			for (long pos = kyi; pos < kyi + ky; ++pos) {
				for (long iy = 0; iy < ny - pos; ++iy) {
					deg = (gyPows[iy] * (iy + pos)) % My;
					pyvals[iy] = ksiyPows[deg];
				}
				for (long iy = ny - pos; iy < ny; ++iy) {
					deg = (gyPows[iy] * (iy + pos - ny)) % My;
					pyvals[iy] = ksiyPows[deg];
				}

				EvaluatorUtils::rightRotateAndEqual(pyvals, 1, ny, 0, kyi);

				IEMBY(pyvals);
				for (long iy = 0; iy < ny; ++iy) {
					pyrInvVec[pos][iy] = EvaluatorUtils::scaleUpToZZ(pyvals[iy].real(), logp);
					pyiInvVec[pos][iy] = EvaluatorUtils::scaleUpToZZ(pyvals[iy].imag(), logp);
				}
			}
		}

		delete[] pxvals;
		delete[] pyvals;

		bootContextMap.insert(pair<pair<long, long>, BootContext>({lognx, logny}, BootContext(pxVec, pyrVec, pyiVec, pxInvVec, pyrInvVec, pyiInvVec, p1, p2, logp)));
	}
}

void Ring2XY::addMatrixContext(long lognx) {
	if (matrixContext.find({lognx, lognx}) == matrixContext.end()) {
		long nx = 1 << lognx;
		ZZ** pvec = new ZZ*[nx];

		long gap = (Nxh >> lognx);
		long powsum = Nxh - gap;
		pvec[0] = new ZZ[N];
		for (long ix = 0; ix < Nxh; ix += gap) {
			pvec[0][ix + (powsum - ix) * Nx] = ZZ(1);
		}

		multByMonomialAndEqual(pvec[0], Mx - powsum, 0);

		for (long ix = 1; ix < nx; ++ix) {
			pvec[ix] = new ZZ[N];
			leftRotate(pvec[ix], pvec[0], 0, ix);
		}
		matrixContext.insert(pair<pair<long, long>, MatrixContext>({lognx, lognx}, MatrixContext(pvec)));
	}
}

void Ring2XY::arrayBitReverse(complex<double>* vals, const long n) {
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

void Ring2XY::DFTY(complex<double>* vals) {
	arrayBitReverse(vals, Ny);
	for (long len = 2; len <= Ny; len <<= 1) {
		long lenh = len >> 1;
		long gap = Ny / len;
		for (long i = 0; i < Ny; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = j * gap;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiyPows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Ring2XY::IDFTY(complex<double>* vals) {
	arrayBitReverse(vals, Ny);
	for (long len = 2; len <= Ny; len <<= 1) {
		long lenh = len >> 1;
		long gap = Ny / len;
		for (long i = 0; i < Ny; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = Ny - (j * gap);
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiyPows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
	for (long iy = 0; iy < Ny; ++iy) {
		vals[iy] /= Ny;
	}
}

void Ring2XY::EMBX(complex<double>* vals, const long nx) {
	arrayBitReverse(vals, nx);
	for (long len = 2; len <= nx; len <<= 1) {
		for (long i = 0; i < nx; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			long gap = Mx / lenq;
			for (long j = 0; j < lenh; ++j) {
				long idx = (gxPows[j] % lenq) * gap;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksixPows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Ring2XY::IEMBX(complex<double>* vals, const long nx) {
	for (long len = nx; len >= 1; len >>= 1) {
		for (long i = 0; i < nx; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			long gap = Mx / lenq;
			for (long j = 0; j < lenh; ++j) {
				long idx = (lenq - (gxPows[j] % lenq)) * gap;
				complex<double> u = vals[i + j] + vals[i + j + lenh];
				complex<double> v = vals[i + j] - vals[i + j + lenh];
				v *= ksixPows[idx];
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

void Ring2XY::EMBY(complex<double>* vals) {
	complex<double>* tmp = new complex<double>[Ny];
	for (long i = 0; i < Ny; ++i) {
		tmp[i] = vals[gyPows[Ny - i] - 1];
	}
	for (long i = 0; i < Ny; ++i) {
		vals[i] = tmp[i];
	}
	delete[] tmp;

	DFTY(vals);
	for (long i = 0; i < Ny; ++i) {
		vals[i] *= dftomegaPows[i];
	}
	IDFTY(vals);
	for (long i = 0; i < Ny; ++i) {
		vals[i] /= omegaPows[i];
	}
}

void Ring2XY::IEMBY(complex<double>* vals) {
	for (long i = 0; i < Ny; ++i) {
		vals[i] *= omegaPows[i];
	}
	DFTY(vals);
	for (long i = 0; i < Ny; ++i) {
		vals[i] /= dftomegaPows[i];
	}
	IDFTY(vals);

	complex<double>* tmp = new complex<double>[Ny]();
	for (long i = 0; i < Ny; ++i) {
		tmp[gyPows[Ny - i] - 1] = vals[i];
	}
	for (long i = 0; i < Ny; ++i) {
		vals[i] = tmp[i];
	}
	delete[] tmp;
}

void Ring2XY::EMBXY(complex<double>* vals, const long nx) {
	complex<double>* tmp = new complex<double>[Ny]();
	for (long iy = 0; iy < Ny; ++iy) {
		EMBX(vals + (iy * nx), nx);
	}
	for (long ix = 0; ix < nx; ++ix) {
		for (long iy = 0; iy < Ny; ++iy) {
			tmp[iy] = vals[ix + (iy * nx)];
		}
		EMBY(tmp);
		for (long iy = 0; iy < Ny; ++iy) {
			vals[ix + (iy * nx)] = tmp[iy];
		}
	}
	delete[] tmp;
}

void Ring2XY::IEMBXY(complex<double>* vals, const long nx) {
	complex<double>* tmp = new complex<double>[Ny]();

	for (long iy = 0; iy < Ny; ++iy) {
		IEMBX(vals + (iy * nx), nx);
	}
	for (long ix = 0; ix < nx; ++ix) {
		for (long iy = 0; iy < Ny; ++iy) {
			tmp[iy] = vals[ix + (iy * nx)];
		}
		IEMBY(tmp);
		for (long iy = 0; iy < Ny; ++iy) {
			vals[ix + (iy * nx)] = tmp[iy];
		}
	}

	delete[] tmp;
}

void Ring2XY::encode(ZZ* mx, complex<double>* vals, long nx, long ny, long logp) {
	long gapx = Nxh / nx;
	long gapy = Ny / ny;

	complex<double>* uvals = new complex<double>[nx * Ny];
	for (long j = 0; j < gapy; ++j) {
		copy(vals, vals + nx * ny, uvals + j * nx * ny);
	}

	IEMBXY(uvals, nx);
	for (long ix = 0, iix = Nxh, irx = 0; ix < nx; ++ix, iix += gapx, irx += gapx) {
		for (long iy = 0; iy < Ny; ++iy) {
			mx[irx + Nx * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].real(), logp);
			mx[iix + Nx * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].imag(), logp);
		}
	}
	delete[] uvals;
}

void Ring2XY::encode(ZZ* mx, double* vals, long nx, long ny, long logp) {
	long gapx = Nxh / nx;
	long gapy = Ny / ny;

	complex<double>* uvals = new complex<double>[nx * Ny];
	for (long i = 0; i < nx * ny; ++i) {
		for (long j = 0; j < gapy; ++j) {
			uvals[i + j * nx * ny].real(vals[i]);
		}
	}

	IEMBXY(uvals, nx);
	for (long ix = 0, iix = Nxh, irx = 0; ix < nx; ++ix, iix += gapx, irx += gapx) {
		for (long iy = 0; iy < Ny; ++iy) {
			mx[irx + Nx * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].real(), logp);
			mx[iix + Nx * iy] = EvaluatorUtils::scaleUpToZZ(uvals[ix + nx * iy].imag(), logp);
		}
	}
	delete[] uvals;
}

void Ring2XY::decode(ZZ* mxy, complex<double>* vals, long nx, long ny, long logp, long logq) {
	ZZ q = qvec[logq];
	ZZ qh = qvec[logq - 1];
	long gapx = Nxh / nx;
	ZZ tmp;

	complex<double>* fvals = new complex<double>[nx * Ny];
	for (long ix = 0, iix = Nxh, irx = 0; ix < nx; ++ix, iix += gapx, irx += gapx) {
		for (long iy = 0; iy < Ny; ++iy) {
			rem(tmp, mxy[irx + Nx * iy], q);
			while (tmp < 0) tmp += q;
			while (tmp > qh) tmp -= q;
			fvals[ix + nx * iy].real(EvaluatorUtils::scaleDownToReal(tmp, logp));

			rem(tmp, mxy[iix + Nx * iy], q);
			while(tmp < 0) tmp += q;
			while (tmp > qh) tmp -= q;
			fvals[ix + nx * iy].imag(EvaluatorUtils::scaleDownToReal(tmp, logp));
		}
	}

	EMBXY(fvals, nx);
	copy(fvals, fvals + nx * ny, vals);
	delete[] fvals;
}


//----------------------------------------------------------------------------------
//   MULTIPLICATION
//----------------------------------------------------------------------------------


void Ring2XY::multXpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q) {
	multiplier.multXpoly(x, a, b, q);
}

void Ring2XY::multXpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q) {
	multiplier.multXpolyAndEqual(a, b, q);
}

void Ring2XY::multYpoly(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q) {
	multiplier.multYpoly(x, a, b, q);
}

void Ring2XY::multYpolyAndEqual(ZZ* a, const ZZ* b, const ZZ& q) {
	multiplier.multYpolyAndEqual(a, b, q);
}

void Ring2XY::mult(ZZ* x, const ZZ* a, const ZZ* b, const ZZ& q) {
	multiplier.mult(x, a, b, q);
}

void Ring2XY::multAndEqual(ZZ* a, const ZZ* b, const ZZ& q) {
	multiplier.multAndEqual(a, b, q);
}

void Ring2XY::square(ZZ* x, const ZZ* a, const ZZ& q) {
	multiplier.square(x, a, q);
}

void Ring2XY::squareAndEqual(ZZ* a, const ZZ& q) {
	multiplier.squareAndEqual(a, q);
}


//----------------------------------------------------------------------------------
//   OTHER
//----------------------------------------------------------------------------------


void Ring2XY::mod(ZZ* res, ZZ* p, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		rem(res[i], p[i], q);
	}
}

void Ring2XY::modAndEqual(ZZ* p, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		rem(p[i], p[i], q);
	}
}

void Ring2XY::negate(ZZ* res, ZZ* p) {
	for (long i = 0; i < N; ++i) {
		res[i] = -p[i];
	}
}

void Ring2XY::negateAndEqual(ZZ* p) {
	for (long i = 0; i < N; ++i) {
		p[i] = -p[i];
	}
}

void Ring2XY::add(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(res[i], p1[i], p2[i], q);
	}
}

void Ring2XY::addAndEqual(ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(p1[i], p1[i], p2[i], q);
	}
}

void Ring2XY::sub(ZZ* res, ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(res[i], p1[i], -p2[i], q);
	}
}

void Ring2XY::subAndEqual(ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(p1[i], p1[i], -p2[i], q);
	}
}

void Ring2XY::subAndEqual2(ZZ* p1, ZZ* p2, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		AddMod(p2[i], p1[i], -p2[i], q);
	}
}

void Ring2XY::multByMonomial(ZZ* res, ZZ* p, const long degx, const long degy) {
	long degxRem = degx % Mx;
	long degyRem = degy % My;

	bool isdegxSmall = degxRem < Nx;
	bool isdegySmall = degyRem < Ny;

	degxRem %= Nx;
	degyRem %= Ny;

	long degxRembar = Nx - degxRem;
	long degyRembar = Ny - degyRem;

	if ((isdegxSmall && isdegySmall) || (!isdegxSmall && !isdegySmall)) {
		for (long ix = 0; ix < degxRem; ++ix) {
			for (long iy = 0; iy < degyRem; ++iy) {
				res[ix + (iy << logNx)] = p[degxRembar + ix + ((degyRembar + iy) << logNx)];
			}
			for (long iy = degyRem; iy < Ny; ++iy) {
				res[ix + (iy << logNx)] = -p[degxRembar + ix + ((iy - degyRem) << logNx)];
			}
		}
		for (long ix = degxRem; ix < Nx; ++ix) {
			for (long iy = 0; iy < degyRem; ++iy) {
				res[ix + (iy << logNx)] = -p[(ix - degxRem) + ((degyRembar + iy) << logNx)];
			}
			for (long iy = degyRem; iy < Ny; ++iy) {
				res[ix + (iy << logNx)] = p[(ix - degxRem) + ((iy - degyRem) << logNx)];
			}
		}
	} else {
		for (long ix = 0; ix < degxRem; ++ix) {
			for (long iy = 0; iy < degyRem; ++iy) {
				res[ix + (iy << logNx)] = -p[degxRembar + ix + ((degyRembar + iy) << logNx)];
			}
			for (long iy = degyRem; iy < Ny; ++iy) {
				res[ix + (iy << logNx)] = p[degxRembar + ix + ((iy - degyRem) << logNx)];
			}
		}
		for (long ix = degxRem; ix < Nx; ++ix) {
			for (long iy = 0; iy < degyRem; ++iy) {
				res[ix + (iy << logNx)] = p[(ix - degxRem) + ((degyRembar + iy) << logNx)];
			}
			for (long iy = degyRem; iy < Ny; ++iy) {
				res[ix + (iy << logNx)] = -p[(ix - degxRem) + ((iy - degyRem) << logNx)];
			}
		}
	}
}

void Ring2XY::multByMonomialAndEqual(ZZ* p, const long degx, const long degy) {
	ZZ* res = new ZZ[N];
	multByMonomial(res, p, degx, degy);
	for (long i = 0; i < N; ++i) {
		p[i] = res[i];
	}
	delete[] res;
}

void Ring2XY::multByConst(ZZ* res, ZZ* p, const ZZ& cnst, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		MulMod(res[i], p[i], cnst, q);
	}
}

void Ring2XY::multByConstAndEqual(ZZ* p, const ZZ& cnst, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		MulMod(p[i], p[i], cnst, q);
	}
}


//----------------------------------------------------------------------------------
//   SHIFTING
//----------------------------------------------------------------------------------


void Ring2XY::leftShift(ZZ* res, ZZ* p, const long bits, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		res[i] = p[i] << bits;
		res[i] %= q;
	}
}

void Ring2XY::leftShiftAndEqual(ZZ* p, const long bits, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		p[i] <<= bits;
		p[i] %= q;
	}
}

void Ring2XY::doubleAndEqual(ZZ* p, ZZ& q) {
	for (long i = 0; i < N; ++i) {
		p[i] <<= 1;
		p[i] %= q;
	}
}

void Ring2XY::rightShift(ZZ* res, ZZ* p, const long bits) {
	for (long i = 0; i < N; ++i) {
		res[i] = p[i] >> bits;
	}
}

void Ring2XY::rightShiftAndEqual(ZZ* p, const long bits) {
	for (long i = 0; i < N; ++i) {
		p[i] >>= bits;
	}
}


//----------------------------------------------------------------------------------
//   ROTATION & CONJUGATION & TRANSPOSITION
//----------------------------------------------------------------------------------


void Ring2XY::leftRotate(ZZ* res, ZZ* p, const long rx, const long ry) {

	ZZ* xxx = new ZZ[N];
	long degx = gxPows[rx];
	for (long j = 0; j < N; j += Nx) {
		for (long ix = 0; ix < Nx; ++ix) {
			long ipow = ix * degx;
			long shift = ipow % Mx;
			if (shift < Nx) {
				xxx[shift + j] = p[ix + j];
			} else {
				xxx[shift - Nx + j] = -p[ix + j];
			}
		}
	}

	for (long i = 0; i < N; ++i) {
		res[i] = ZZ::zero();
	}

	long degy = gyPows[ry];
	for (long ix = 0; ix < Nx; ++ix) {
		for (long iy = 0; iy < Ny; ++iy) {
			long ipow = iy * degy;
			long shift = ipow % My;
			if (shift < Ny) {
				res[ix + (shift << logNx)] += xxx[ix + (iy << logNx)];
			} else {
				for (long t = 0; t < Ny; ++t) {
					res[ix + (t << logNx)] -= xxx[ix + (iy << logNx)];
				}
			}
		}
	}
	delete[] xxx;
}

void Ring2XY::conjugate(ZZ* res, ZZ* p) {
	ZZ* xxx = new ZZ[N];

	for (long j = 0; j < N; j += Nx) {
		xxx[j] = p[j];
		for (long ix = 1; ix < Nx; ++ix) {
			xxx[Nx - ix + j] = -p[ix + j];
		}
	}
	for (long i = 0; i < N; ++i) {
		res[i] = ZZ::zero();
	}

	for (long ix = 0; ix < Nx; ++ix) {
		res[ix] += xxx[ix];
		for (long j = 2 * Nx; j < N; j += Nx) {
			res[ix + (N - j + Nx)] += xxx[ix + j];
		}
		for (long j = 0; j < N; j += Nx) {
			res[ix + j] -= xxx[ix + Nx];
		}
	}

	delete[] xxx;
}


//----------------------------------------------------------------------------------
//   SAMPLING
//----------------------------------------------------------------------------------


void Ring2XY::sampleGauss(ZZ* res) {
	static double const Pi = 4.0 * atan(1.0);
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

void Ring2XY::sampleHWT(ZZ* res) {
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

void Ring2XY::sampleZO(ZZ* res) {
	ZZ tmp = RandomBits_ZZ((N << 1));
	for (long i = 0; i < N; ++i) {
		res[i] = (bit(tmp, 2 * i) == 0) ? ZZ(0) : (bit(tmp, 2 * i + 1) == 0) ? ZZ(1) : ZZ(-1);
	}
}

void Ring2XY::sampleUniform(ZZ* res, long bits) {
	for (long i = 0; i < N; i++) {
		res[i] = RandomBits_ZZ(bits);
	}
}
