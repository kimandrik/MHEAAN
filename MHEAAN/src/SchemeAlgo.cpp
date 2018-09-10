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

Ciphertext SchemeAlgo::sqrMatMult(Ciphertext& cipher1, Ciphertext& cipher2, long logp, long n) {
	long logn = log2(n);

	SqrMatContext sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);

	Ciphertext* cipherP = new Ciphertext[n];
	Ciphertext* cipherR = new Ciphertext[n];

	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		cipherR[i] = (i == 0) ? cipher1 : scheme.rightRotateFast(cipher1, i, 0);
		scheme.modDownByAndEqual(cipherR[i], logn);

		cipherP[i] = scheme.multPoly(cipher2, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(cipherP[i], logn);
		for (long j = 0; j < logn; ++j) {
			Ciphertext rot = scheme.leftRotateFast(cipherP[i], 0, (1 << j));
			scheme.addAndEqual(cipherP[i], rot);
		}
		scheme.multAndEqual(cipherR[i], cipherP[i]);
	}
	NTL_EXEC_RANGE_END;

	for (long i = 1; i < n; ++i) {
		scheme.addAndEqual(cipherR[0], cipherR[i]);
	}

	Ciphertext res = scheme.reScaleBy(cipherR[0], logp);
	delete[] cipherR;
	delete[] cipherP;
	return res;
}

void SchemeAlgo::sqrMatMultAndEqual(Ciphertext& cipher, long logp, long n) {
	long logn = log2(n);

	SqrMatContext sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);

	Ciphertext* cipherP = new Ciphertext[n];
	Ciphertext* cipherR = new Ciphertext[n];

	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		cipherR[i] = (i == 0) ? cipher : scheme.rightRotateFast(cipher, i, 0);
		scheme.modDownByAndEqual(cipherR[i], logp);

		cipherP[i] = scheme.multPoly(cipher, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(cipherP[i], logp);
		for (long j = 0; j < logn; ++j) {
			Ciphertext rot = scheme.leftRotateFast(cipherP[i], 0, (1 << j));
			scheme.addAndEqual(cipherP[i], rot);
		}
		scheme.multAndEqual(cipherR[i], cipherP[i]);
	}
	NTL_EXEC_RANGE_END;

	for (long i = 1; i < n; ++i) {
		scheme.addAndEqual(cipherR[0], cipherR[i]);
	}


	cipher = scheme.reScaleBy(cipherR[0], logp);
	delete[] cipherR;
	delete[] cipherP;
}


Ciphertext SchemeAlgo::matInv(Ciphertext& cipher, long logp, long n, long r) {
	//TODO change method to right one
	long logn = log2(n);
	Ciphertext cbar = scheme.negate(cipher);
	SqrMatContext sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	scheme.addPolyAndEqual(cbar, sqrMatContext.mvec[0], sqrMatContext.logp);
	Ciphertext cpow = cbar;
	Ciphertext tmp = scheme.addPoly(cbar, sqrMatContext.mvec[0], sqrMatContext.logp);
	Ciphertext res = tmp;

	for (long i = 1; i < r; ++i) {
		sqrMatMultAndEqual(cpow, logp, n);
		tmp = cpow;
		scheme.addPolyAndEqual(tmp, sqrMatContext.mvec[0], sqrMatContext.logp);
		scheme.modDownToAndEqual(res, tmp.logq);
		tmp = sqrMatMult(tmp, res, logp, n);
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
	Ciphertext res = scheme.multConst(cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(res, coeffs[0], dlogp);

	for (int i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext aixi = scheme.multConst(cpows[i], coeffs[i + 1], logp);
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

	Ciphertext res = scheme.multConst(cpows[0], coeffs[1], logp);

	scheme.addConstAndEqual(res, coeffs[0], dlogp);

	for (int i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext aixi = scheme.multConst(cpows[i], coeffs[i + 1], logp);
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
	Ciphertext aixi = scheme.multConst(cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(aixi, coeffs[0], dlogp);

	Ciphertext* res = new Ciphertext[degree];
	res[0] = aixi;
	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			aixi = scheme.multConst(cpows[i], coeffs[i + 1], logp);
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

void SchemeAlgo::DFTX0(Ciphertext* ciphers, long n0) {
	bitReverse(ciphers, n0);
	for (long len = 2; len <= n0; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.M0 / len;
		for (long i = 0; i < n0; i += len) {
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

void SchemeAlgo::IDFTX0Lazy(Ciphertext* ciphers, long n0) {
	bitReverse(ciphers, n0);
	for (long len = 2; len <= n0; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.M0 - scheme.ring.M0 / len;
		for (long i = 0; i < n0; i += len) {
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

void SchemeAlgo::IDFTX0(Ciphertext* ciphers, long n0) {
	IDFTX0Lazy(ciphers, n0);

	long lognx = log2((double) n0);
	NTL_EXEC_RANGE(n0, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(ciphers[i], lognx);
	}
	NTL_EXEC_RANGE_END;
}

void SchemeAlgo::DFTX1(Ciphertext* ciphers, long n1) {
	bitReverse(ciphers, n1);
	for (long len = 2; len <= n1; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.M1 / len;
		for (long i = 0; i < n1; i += len) {
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

void SchemeAlgo::IDFTX1Lazy(Ciphertext* ciphers, long n1) {
	bitReverse(ciphers, n1);
	for (long len = 2; len <= n1; len <<= 1) {
		long lenh = len >> 1;
		long shift = scheme.ring.M1 - scheme.ring.M1 / len;
		for (long i = 0; i < n1; i += len) {
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

void SchemeAlgo::IDFTX1(Ciphertext* ciphers, long n1) {
	IDFTX1Lazy(ciphers, n1);

	long logny = log2((double) n1);
	NTL_EXEC_RANGE(n1, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(ciphers[i], logny);
	}
	NTL_EXEC_RANGE_END;
}

void SchemeAlgo::DFT(Ciphertext* ciphers, long n0, long n1) {
	for (long iy = 0; iy < n1; ++iy) {
		DFTX0(ciphers + iy * n0, n0);
	}

	Ciphertext* tmp = new Ciphertext[n1];
	for (long ix = 0; ix < n0; ++ix) {
		for (long iy = 0; iy < n1; ++iy) {
			tmp[iy] = ciphers[ix + iy * n0];
		}

		DFTX1(tmp, n1);

		for (long iy = 0; iy < n1; ++iy) {
			ciphers[ix + iy * n0] = tmp[iy];
		}
	}
	delete[] tmp;
}

void SchemeAlgo::IDFTLazy(Ciphertext* ciphers, long n0, long n1) {
	for (long iy = 0; iy < n1; ++iy) {
		IDFTX0Lazy(ciphers + iy * n0, n0);
	}

	Ciphertext* tmp = new Ciphertext[n1];
	for (long ix = 0; ix < n0; ++ix) {
		for (long iy = 0; iy < n1; ++iy) {
			tmp[iy] = ciphers[ix + iy * n0];
		}

		IDFTX1Lazy(tmp, n1);

		for (long iy = 0; iy < n1; ++iy) {
			ciphers[ix + iy * n0] = tmp[iy];
		}
	}
	delete[] tmp;
}

void SchemeAlgo::IDFT(Ciphertext* ciphers, long n0, long n1) {
	IDFTLazy(ciphers, n0, n1);

	long n = n0 * n1;
	long logn = log2((double) n);
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		scheme.divPo2AndEqual(ciphers[i], logn);
	}
	NTL_EXEC_RANGE_END;
}
