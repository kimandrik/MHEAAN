/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "SchemeAlgo.h"

Ciphertext* SchemeAlgo::encryptSingleArray(complex<double>* vals, long size, long logp) {
	Ciphertext* res = new Ciphertext[size];
	for (long i = 0; i < size; ++i) {
		res[i] = scheme.encryptSingle(vals[i], logp, scheme.ring.logQ);
	}
	return res;
}

complex<double>* SchemeAlgo::decryptSingleArray(SecretKey& secretKey, Ciphertext* ciphers, long size) {
	complex<double>* res = new complex<double>[size];
	for (int i = 0; i < size; ++i) {
		res[i] = scheme.decryptSingle(secretKey, ciphers[i]);
	}
	return res;
}

Ciphertext SchemeAlgo::powerOf2(Ciphertext& cipher, const long logp, const long logDegree) {
	Ciphertext res = cipher;
	for (long i = 0; i < logDegree; ++i) {
		scheme.squareAndEqual(res);
		scheme.reScaleByAndEqual(res, logp);
	}
	return res;
}

Ciphertext* SchemeAlgo::powerOf2Extended(Ciphertext& cipher, const long logp,
		const long logDegree) {
	Ciphertext* res = new Ciphertext[logDegree + 1];
	res[0] = cipher;
	for (long i = 1; i < logDegree + 1; ++i) {
		res[i] = scheme.square(res[i - 1]);
		scheme.reScaleByAndEqual(res[i], logp);
	}
	return res;
}

//-----------------------------------------

Ciphertext SchemeAlgo::power(Ciphertext& cipher, const long logp, const long degree) {
	long logDegree = log2((double) degree);
	long po2Degree = 1 << logDegree;

	Ciphertext res = powerOf2(cipher, logp, logDegree);
	long remDegree = degree - po2Degree;
	if (remDegree > 0) {
		Ciphertext tmp = power(cipher, logp, remDegree);
		long bitsDown = tmp.logq - res.logq;
		scheme.modDownByAndEqual(tmp, bitsDown);
		scheme.multAndEqual(res, tmp);
		scheme.reScaleByAndEqual(res, logp);
	}
	return res;
}

Ciphertext* SchemeAlgo::powerExtended(Ciphertext& cipher, const long logp, const long degree) {
	Ciphertext* res = new Ciphertext[degree];
	long logDegree = log2((double) degree);
	Ciphertext* cpows = powerOf2Extended(cipher, logp, logDegree);
	long idx = 0;
	for (long i = 0; i < logDegree; ++i) {
		long powi = (1 << i);
		res[idx++] = cpows[i];
		for (int j = 0; j < powi - 1; ++j) {
			long bitsDown = res[j].logq - cpows[i].logq;
			res[idx] = scheme.modDownBy(res[j], bitsDown);
			scheme.multAndEqual(res[idx], cpows[i]);
			scheme.reScaleByAndEqual(res[idx++], logp);
		}
	}
	res[idx++] = cpows[logDegree];
	long degree2 = (1 << logDegree);
	for (int i = 0; i < (degree - degree2); ++i) {
		long bitsDown = res[i].logq - cpows[logDegree].logq;
		res[idx] = scheme.modDownBy(res[i], bitsDown);
		scheme.multAndEqual(res[idx], cpows[logDegree]);
		scheme.reScaleByAndEqual(res[idx++], logp);
	}
	return res;
}

//-----------------------------------------

Ciphertext SchemeAlgo::prodOfPo2(Ciphertext* ciphers, const long logp, const long logDegree) {
	Ciphertext* res = ciphers;
	for (long i = logDegree - 1; i >= 0; --i) {
		long powih = (1 << i);
		Ciphertext* tmp = new Ciphertext[powih];
		NTL_EXEC_RANGE(powih, first, last);
		for (long j = first; j < last; ++j) {
			tmp[j] = scheme.mult(res[2 * j], res[2 * j + 1]);
			scheme.reScaleByAndEqual(tmp[j], logp);
		}
		NTL_EXEC_RANGE_END;
		res = tmp;
	}
	return res[0];
}

Ciphertext SchemeAlgo::prod(Ciphertext* ciphers, const long logp, const long degree) {
	long logDegree = log2((double) degree) + 1;
	long idx = 0;
	bool isinit = false;
	Ciphertext res;
	for (long i = 0; i < logDegree; ++i) {
		if (bit(degree, i)) {
			long powi = (1 << i);
			Ciphertext* tmp = new Ciphertext[powi];
			for (long j = 0; j < powi; ++j) {
				tmp[j] = ciphers[idx + j];
			}
			Ciphertext iprod = prodOfPo2(tmp, logp, i);
			if (isinit) {
				long bitsDown = res.logq - iprod.logq;
				scheme.modDownByAndEqual(res, bitsDown);
				scheme.multAndEqual(res, iprod);
				scheme.reScaleByAndEqual(res, logp);
			} else {
				res = iprod;
				isinit = true;
			}
			idx += powi;
		}
	}
	return res;
}

Ciphertext SchemeAlgo::sum(Ciphertext* ciphers, const long size) {
	Ciphertext res = ciphers[0];
	for (long i = 1; i < size; ++i) {
		scheme.addAndEqual(res, ciphers[i]);
	}
	return res;
}

Ciphertext* SchemeAlgo::multVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size) {
	Ciphertext* res = new Ciphertext[size];
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		res[i] = scheme.mult(ciphers1[i], ciphers2[i]);
	}
	NTL_EXEC_RANGE_END;
	return res;
}

void SchemeAlgo::multAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size) {
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multAndEqual(ciphers1[i], ciphers2[i]);
	}
	NTL_EXEC_RANGE_END;
}

Ciphertext* SchemeAlgo::multAndModSwitchVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long precisionBits, const long size) {
	Ciphertext* res = new Ciphertext[size];
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		res[i] = scheme.mult(ciphers1[i], ciphers2[i]);
		scheme.reScaleByAndEqual(res[i], precisionBits);
	}
	NTL_EXEC_RANGE_END;
	return res;
}

void SchemeAlgo::multModSwitchAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long precisionBits, const long size) {
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multAndEqual(ciphers1[i], ciphers2[i]);
		scheme.reScaleByAndEqual(ciphers1[i], precisionBits);
	}
	NTL_EXEC_RANGE_END;
}

Ciphertext SchemeAlgo::innerProd(Ciphertext* ciphers1, Ciphertext* ciphers2, const long logp, const long size) {
	Ciphertext cip = scheme.mult(ciphers1[size - 1], ciphers2[size - 1]);

	NTL_EXEC_RANGE(size-1, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext cprodi = scheme.mult(ciphers1[i], ciphers2[i]);
		scheme.addAndEqual(cip, cprodi);
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(cip, logp);
	return cip;
}

Ciphertext SchemeAlgo::squareMatMult(Ciphertext& cipher1, Ciphertext& cipher2, long logp, long size) {
	long logSize = log2(size);

	MatrixContext ptmpm = scheme.ring.matrixContext.at({logSize,logSize});

	Ciphertext* cipherP = new Ciphertext[size];
	Ciphertext* cipherR = new Ciphertext[size];

	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		cipherR[i] = (i == 0) ? cipher1 : scheme.rightRotateFast(cipher1, i, 0);
		scheme.modDownByAndEqual(cipherR[i], logSize);

		cipherP[i] = scheme.multByPoly(cipher2, ptmpm.mvec[i], logSize);
		scheme.reScaleByAndEqual(cipherP[i], logSize);
		for (long j = 0; j < logSize; ++j) {
			Ciphertext rot = scheme.leftRotateFast(cipherP[i], 0, (1 << j));
			scheme.addAndEqual(cipherP[i], rot);
		}
		scheme.multAndEqual(cipherR[i], cipherP[i]);
	}
	NTL_EXEC_RANGE_END;

	for (long i = 1; i < size; ++i) {
		scheme.addAndEqual(cipherR[0], cipherR[i]);
	}

	Ciphertext res = scheme.reScaleBy(cipherR[0], logp);
	delete[] cipherR;
	delete[] cipherP;
	return res;
}

void SchemeAlgo::squareMatMultAndEqual(Ciphertext& cipher, long logp, long size) {
	long logSize = log2(size);

	MatrixContext ptmpm = scheme.ring.matrixContext.at({logSize, logSize});

	Ciphertext* cipherP = new Ciphertext[size];
	Ciphertext* cipherR = new Ciphertext[size];

	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		cipherR[i] = (i == 0) ? cipher : scheme.rightRotateFast(cipher, i, 0);
		scheme.modDownByAndEqual(cipherR[i], logp);

		cipherP[i] = scheme.multByPoly(cipher, ptmpm.mvec[i], logp);
		scheme.reScaleByAndEqual(cipherP[i], logp);
		for (long j = 0; j < logSize; ++j) {
			Ciphertext rot = scheme.leftRotateFast(cipherP[i], 0, (1 << j));
			scheme.addAndEqual(cipherP[i], rot);
		}
		scheme.multAndEqual(cipherR[i], cipherP[i]);
	}
	NTL_EXEC_RANGE_END;

	for (long i = 1; i < size; ++i) {
		scheme.addAndEqual(cipherR[0], cipherR[i]);
	}


	cipher = scheme.reScaleBy(cipherR[0], logp);
	delete[] cipherR;
	delete[] cipherP;
}

Ciphertext SchemeAlgo::matMult(Ciphertext& cipher1, Ciphertext& cipher2, const long logp, const long sizex, const long sizey, const long sizez) {
	return NULL;
}

Ciphertext SchemeAlgo::matInv(Ciphertext& cipher, long logp, long size, long r) {
	//TODO change method to right one
	long logSize = log2(size);
	Ciphertext cbar = scheme.negate(cipher);
	MatrixContext mpmpt = scheme.ring.matrixContext.at({logSize, logSize});
	scheme.addPolyAndEqual(cbar, mpmpt.mvec[0], logp);
	Ciphertext cpow = cbar;
	Ciphertext tmp = scheme.addPoly(cbar, mpmpt.mvec[0], logp);
	Ciphertext res = tmp;

	for (long i = 1; i < r; ++i) {
		squareMatMultAndEqual(cpow, logp, size);
		tmp = cpow;
		scheme.addPolyAndEqual(tmp, mpmpt.mvec[0], logp);
		scheme.modDownToAndEqual(res, tmp.logq);
		tmp = squareMatMult(tmp, res, logp, size);
		res = tmp;
	}
	return res;
}

//-----------------------------------------

Ciphertext SchemeAlgo::inverse(Ciphertext& cipher, long logp, long steps) {
	Ciphertext cbar = scheme.negate(cipher);
	scheme.addConstAndEqual(cbar, 1.0, logp);
	Ciphertext cpow = cbar;
	Ciphertext tmp = scheme.addConst(cbar, 1.0, logp);
	scheme.modDownByAndEqual(tmp, logp);
	Ciphertext res = tmp;

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, logp);
		tmp = cpow;
		scheme.addConstAndEqual(tmp, 1.0, logp);
		scheme.multAndEqual(tmp, res);
		scheme.reScaleByAndEqual(tmp, logp);
		res = tmp;
	}
	return res;
}

Ciphertext* SchemeAlgo::inverseExtended(Ciphertext& cipher, const long logp, const long steps) {
	Ciphertext* res = new Ciphertext[steps];
	Ciphertext cpow = cipher;
	Ciphertext tmp = scheme.addConst(cipher, 1.0, logp);
	scheme.modDownByAndEqual(tmp, logp);
	res[0] = tmp;

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, logp);
		tmp = cpow;
		scheme.addConstAndEqual(tmp, 1.0, logp);
		scheme.multAndEqual(tmp, res[i - 1]);
		scheme.reScaleByAndEqual(tmp, logp);
		res[i] = tmp;
	}
	return res;
}

//-----------------------------------------

Ciphertext SchemeAlgo::function(Ciphertext& cipher, string& funcName, long logp, long degree) {
	Ciphertext* cpows = powerExtended(cipher, logp, degree);

	long dlogp = 2 * logp;

	double* coeffs = taylorCoeffsMap.at(funcName);
	Ciphertext res = scheme.multByConst(cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(res, coeffs[0], dlogp);

	for (int i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext aixi = scheme.multByConst(cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(res, aixi.logq);
			scheme.addAndEqual(res, aixi);
		}
	}
	scheme.reScaleByAndEqual(res, logp);
	return res;
}

Ciphertext SchemeAlgo::functionLazy(Ciphertext& cipher, string& funcName,
		const long logp, const long degree) {
	Ciphertext* cpows = powerExtended(cipher, logp, degree);

	long dlogp = 2 * logp;

	double* coeffs = taylorCoeffsMap.at(funcName);

	Ciphertext res = scheme.multByConst(cpows[0], coeffs[1], logp);

	scheme.addConstAndEqual(res, coeffs[0], dlogp);

	for (int i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext aixi = scheme.multByConst(cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(res, aixi.logq);
			scheme.addAndEqual(res, aixi);
		}
	}
	return res;
}

Ciphertext* SchemeAlgo::functionExtended(Ciphertext& cipher, string& funcName, long logp, long degree) {
	Ciphertext* cpows = powerExtended(cipher, logp, degree);

	long dlogp = 2 * logp;
	double* coeffs = taylorCoeffsMap.at(funcName);
	Ciphertext aixi = scheme.multByConst(cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(aixi, coeffs[0], dlogp);

	Ciphertext* res = new Ciphertext[degree];
	res[0] = aixi;
	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			aixi = scheme.multByConst(cpows[i], coeffs[i + 1], logp);
			Ciphertext ctmp = scheme.modDownTo(res[i - 1], aixi.logq);
			scheme.addAndEqual(aixi, ctmp);
			res[i] = aixi;
		} else {
			res[i] = res[i - 1];
		}
	}
	NTL_EXEC_RANGE(degree, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(res[i], logp);
	}
	NTL_EXEC_RANGE_END;
	return res;
}

void SchemeAlgo::bitReverse(Ciphertext* ciphers, const long n) {
	for (long i = 1, j = 0; i < n; ++i) {
		long bit = n >> 1;
		for (; j >= bit; bit >>= 1) {
			j -= bit;
		}
		j += bit;
		if (i < j) {
			swap(ciphers[i], ciphers[j]);
		}
	}
}

void SchemeAlgo::DFTX(Ciphertext* ciphers, const long nx) {
	bitReverse(ciphers, nx);
	for (long len = 2; len <= nx; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.Mx / len;
		for (long i = 0; i < nx; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext u = ciphers[i + j];
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], shift * j, 0);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTXLazy(Ciphertext* ciphers, const long nx) {
	bitReverse(ciphers, nx);
	for (long len = 2; len <= nx; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.Mx - scheme.ring.Mx / len;
		for (long i = 0; i < nx; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext u = ciphers[i + j];
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], shift * j, 0);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTX(Ciphertext* ciphers, const long nx) {
	IDFTXLazy(ciphers, nx);

	long lognx = log2((double) nx);
	NTL_EXEC_RANGE(nx, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(ciphers[i], lognx);
	}
	NTL_EXEC_RANGE_END;
}

void SchemeAlgo::DFTY(Ciphertext* ciphers, const long ny) {
	bitReverse(ciphers, ny);
	for (long len = 2; len <= ny; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.My / len;
		for (long i = 0; i < ny; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext u = ciphers[i + j];
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], 0, shift * j);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTYLazy(Ciphertext* ciphers, const long ny) {
	bitReverse(ciphers, ny);
	for (long len = 2; len <= ny; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.My - scheme.ring.My / len;
		for (long i = 0; i < ny; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext u = ciphers[i + j];
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], 0, shift * j);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTY(Ciphertext* ciphers, const long ny) {
	IDFTYLazy(ciphers, ny);

	long logny = log2((double) ny);
	NTL_EXEC_RANGE(ny, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(ciphers[i], logny);
	}
	NTL_EXEC_RANGE_END;
}

void SchemeAlgo::DFTXY(Ciphertext* ciphers, const long nx, const long ny) {
	for (long iy = 0; iy < ny; ++iy) {
		DFTX(ciphers + iy * nx, nx);
	}

	Ciphertext* tmp = new Ciphertext[ny];
	for (long ix = 0; ix < nx; ++ix) {
		for (long iy = 0; iy < ny; ++iy) {
			tmp[iy] = ciphers[ix + iy * nx];
		}

		DFTY(tmp, ny);

		for (long iy = 0; iy < ny; ++iy) {
			ciphers[ix + iy * nx] = tmp[iy];
		}
	}
	delete[] tmp;
}

void SchemeAlgo::IDFTXYLazy(Ciphertext* ciphers, const long nx, const long ny) {
	for (long iy = 0; iy < ny; ++iy) {
		IDFTXLazy(ciphers + iy * nx, nx);
	}

	Ciphertext* tmp = new Ciphertext[ny];
	for (long ix = 0; ix < nx; ++ix) {
		for (long iy = 0; iy < ny; ++iy) {
			tmp[iy] = ciphers[ix + iy * nx];
		}

		IDFTYLazy(tmp, ny);

		for (long iy = 0; iy < ny; ++iy) {
			ciphers[ix + iy * nx] = tmp[iy];
		}
	}
	delete[] tmp;
}

void SchemeAlgo::IDFTXY(Ciphertext* ciphers, const long nx, const long ny) {
	IDFTXYLazy(ciphers, nx, ny);

	long n = nx * ny;
	long logn = log2((double) n);
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		scheme.divByPo2AndEqual(ciphers[i], logn);
	}
	NTL_EXEC_RANGE_END;
}
