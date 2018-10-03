/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "SchemeAlgo.h"

Ciphertext** SchemeAlgo::encryptSingleArray(complex<double>* vals, long size, long logp) {
	Ciphertext** res = new Ciphertext*[size];
	for (long i = 0; i < size; ++i) {
		res[i] = scheme.encryptSingle(vals[i], logp, logQ);
	}
	return res;
}

complex<double>* SchemeAlgo::decryptSingleArray(SecretKey& secretKey, Ciphertext** ciphers, long size) {
	complex<double>* res = new complex<double>[size];
	for (int i = 0; i < size; ++i) {
		res[i] = scheme.decryptSingle(secretKey, ciphers[i]);
	}
	return res;
}

Ciphertext* SchemeAlgo::powerOf2(Ciphertext* cipher, const long logp, const long logDegree) {
	Ciphertext* res = new Ciphertext(cipher);
	for (long i = 0; i < logDegree; ++i) {
		scheme.squareAndEqual(res);
		scheme.reScaleByAndEqual(res, logp);
	}
	return res;
}

Ciphertext** SchemeAlgo::powerOf2Extended(Ciphertext* cipher, const long logp, const long logDegree) {
	Ciphertext** res = new Ciphertext*[logDegree + 1];
	res[0] = new Ciphertext(cipher);
	for (long i = 1; i < logDegree + 1; ++i) {
		res[i] = scheme.square(res[i - 1]);
		scheme.reScaleByAndEqual(res[i], logp);
	}
	return res;
}

//-----------------------------------------

Ciphertext* SchemeAlgo::power(Ciphertext* cipher, const long logp, const long degree) {
	long logDegree = log2((double) degree);
	long po2Degree = 1 << logDegree;

	Ciphertext* res = powerOf2(cipher, logp, logDegree);
	long remDegree = degree - po2Degree;
	if (remDegree > 0) {
		Ciphertext* tmp = power(cipher, logp, remDegree);
		long bitsDown = tmp->logq - res->logq;
		scheme.modDownByAndEqual(tmp, bitsDown);
		scheme.multAndEqual(res, tmp);
		scheme.reScaleByAndEqual(res, logp);
		delete tmp;
	}
	return res;
}

Ciphertext** SchemeAlgo::powerExtended(Ciphertext* cipher, const long logp, const long degree) {
	Ciphertext** res = new Ciphertext*[degree];
	long logDegree = log2((double) degree);
	Ciphertext** cpows = powerOf2Extended(cipher, logp, logDegree);
	long idx = 0;
	for (long i = 0; i < logDegree; ++i) {
		long powi = (1 << i);
		res[idx++] = new Ciphertext(cpows[i]);
		for (long j = 0; j < powi - 1; ++j) {
			res[idx] = scheme.modDownTo(res[j], cpows[i]->logq);
			scheme.multAndEqual(res[idx], cpows[i]);
			scheme.reScaleByAndEqual(res[idx++], logp);
		}
	}
	res[idx++] = new Ciphertext(cpows[logDegree]);
	long degree2 = (1 << logDegree);
	for (int i = 0; i < (degree - degree2); ++i) {
		res[idx] = scheme.modDownTo(res[i], cpows[logDegree]->logq);
		scheme.multAndEqual(res[idx], cpows[logDegree]);
		scheme.reScaleByAndEqual(res[idx++], logp);
	}

	for (long i = 0; i < logDegree; ++i) {
		delete cpows[i];
	}
	delete[] cpows;
	return res;
}

//-----------------------------------------

Ciphertext* SchemeAlgo::prodOfPo2(Ciphertext** ciphers, const long logp, const long logDegree) {
	Ciphertext** res = ciphers;
	for (long i = logDegree - 1; i >= 0; --i) {
		long powih = (1 << i);
		Ciphertext** tmp = new Ciphertext*[powih];
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

Ciphertext* SchemeAlgo::prod(Ciphertext** ciphers, const long logp, const long degree) {
	long logDegree = log2((double) degree) + 1;
	long idx = 0;
	bool isinit = false;
	Ciphertext* res;
	for (long i = 0; i < logDegree; ++i) {
		if (bit(degree, i)) {
			long powi = (1 << i);
			Ciphertext** tmp = new Ciphertext*[powi];
			for (long j = 0; j < powi; ++j) {
				tmp[j] = new Ciphertext(ciphers[idx + j]);
			}
			Ciphertext* iprod = prodOfPo2(tmp, logp, i);
			if (isinit) {
				long bitsDown = res->logq - iprod->logq;
				scheme.modDownByAndEqual(res, bitsDown);
				scheme.multAndEqual(res, iprod);
				scheme.reScaleByAndEqual(res, logp);
			} else {
				res = new Ciphertext(iprod);
				isinit = true;
			}
			idx += powi;
		}
	}
	return res;
}

Ciphertext* SchemeAlgo::sum(Ciphertext** ciphers, const long size) {
	Ciphertext* res = new Ciphertext(ciphers[0]);
	for (long i = 1; i < size; ++i) {
		scheme.addAndEqual(res, ciphers[i]);
	}
	return res;
}

Ciphertext** SchemeAlgo::multVec(Ciphertext** ciphers1, Ciphertext** ciphers2, const long size) {
	Ciphertext** res = new Ciphertext*[size];
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		res[i] = scheme.mult(ciphers1[i], ciphers2[i]);
	}
	NTL_EXEC_RANGE_END;
	return res;
}

void SchemeAlgo::multAndEqualVec(Ciphertext** ciphers1, Ciphertext** ciphers2, const long size) {
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multAndEqual(ciphers1[i], ciphers2[i]);
	}
	NTL_EXEC_RANGE_END;
}

Ciphertext** SchemeAlgo::multAndModSwitchVec(Ciphertext** ciphers1, Ciphertext** ciphers2, const long precisionBits, const long size) {
	Ciphertext** res = new Ciphertext*[size];
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		res[i] = scheme.mult(ciphers1[i], ciphers2[i]);
		scheme.reScaleByAndEqual(res[i], precisionBits);
	}
	NTL_EXEC_RANGE_END;
	return res;
}

void SchemeAlgo::multModSwitchAndEqualVec(Ciphertext** ciphers1, Ciphertext** ciphers2, const long precisionBits, const long size) {
	NTL_EXEC_RANGE(size, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multAndEqual(ciphers1[i], ciphers2[i]);
		scheme.reScaleByAndEqual(ciphers1[i], precisionBits);
	}
	NTL_EXEC_RANGE_END;
}

Ciphertext* SchemeAlgo::innerProd(Ciphertext** ciphers1, Ciphertext** ciphers2, const long logp, const long size) {
	Ciphertext* cip = scheme.mult(ciphers1[size - 1], ciphers2[size - 1]);

	NTL_EXEC_RANGE(size-1, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext* cprodi = scheme.mult(ciphers1[i], ciphers2[i]);
		scheme.addAndEqual(cip, cprodi);
		delete cprodi;
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(cip, logp);
	return cip;
}

Ciphertext* SchemeAlgo::transpose(Ciphertext* cipher, long logp, long n) {
	long logn = log2(n);
	mutex m;
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	Ciphertext* res = new Ciphertext(cipher->logp + sqrMatContext.logp, cipher->logq, cipher->n0, cipher->n1);
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext* tmp = scheme.multPoly(cipher, sqrMatContext.mvec[i], sqrMatContext.logp);
		if(i > 0) scheme.leftRotateAndEqual(tmp,i, N1 - i);
		m.lock();
		scheme.addAndEqual(res, tmp);
		m.unlock();
		delete tmp;
	}
	NTL_EXEC_RANGE_END;
	scheme.reScaleByAndEqual(res, sqrMatContext.logp);
	return res;
}

Ciphertext* SchemeAlgo::sqrMatMult(Ciphertext* cipher1, Ciphertext* cipher2, long logp, long n) {
	long logn = log2(n);

	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	mutex m;
	Ciphertext* res = new Ciphertext(cipher1->logp + cipher2->logp, cipher1->logq - sqrMatContext.logp, cipher1->n0, cipher1->n1);
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext* tmp2 = scheme.multPoly(cipher2, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(tmp2, sqrMatContext.logp);

		for (long j = 0; j < logn; ++j) {
			Ciphertext* rot = scheme.leftRotate(tmp2, 0, (1 << j));
			scheme.addAndEqual(tmp2, rot);
			delete rot;
		}

		Ciphertext* tmp1 = (i == 0) ? new Ciphertext(cipher1) : scheme.rightRotate(cipher1, i, 0);
		scheme.modDownByAndEqual(tmp1, sqrMatContext.logp);
		scheme.multAndEqual(tmp1, tmp2);
		delete tmp2;
		m.lock();
		scheme.addAndEqual(res, tmp1);
		m.unlock();
		delete tmp1;
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(res, logp);

	return res;
}

Ciphertext* SchemeAlgo::sqrMatSqr(Ciphertext* cipher, long logp, long n) {
	long logn = log2(n);

	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	mutex m;
	Ciphertext* res = new Ciphertext(2 * cipher->logp, cipher->logq - sqrMatContext.logp, cipher->n0, cipher->n1);
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext* tmp = scheme.multPoly(cipher, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(tmp, sqrMatContext.logp);
		for (long j = 0; j < logn; ++j) {
			Ciphertext* rot = scheme.leftRotate(tmp, 0, (1 << j));
			scheme.addAndEqual(tmp, rot);
			delete rot;
		}
		Ciphertext* rtmp = (i == 0) ? new Ciphertext(cipher) : scheme.rightRotate(cipher, i, 0);
		scheme.modDownByAndEqual(rtmp, sqrMatContext.logp);
		scheme.multAndEqual(rtmp, tmp);
		delete tmp;
		m.lock();
		scheme.addAndEqual(res, rtmp);
		m.unlock();
		delete rtmp;
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(res, logp);

	return res;
}


Ciphertext* SchemeAlgo::matInv(Ciphertext* cipher, long logp, long n, long r) {
	long logn = log2(n);
	Ciphertext* cbar = scheme.negate(cipher);
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	scheme.addPolyAndEqual(cbar, sqrMatContext.mvec[0], sqrMatContext.logp);
	Ciphertext* cpow = new Ciphertext(cbar);
	Ciphertext* res = scheme.addPoly(cbar, sqrMatContext.mvec[0], sqrMatContext.logp);

	for (long i = 1; i < r; ++i) {
		Ciphertext* tmp = sqrMatSqr(cpow, logp, n);
		delete cpow;
		cpow = new Ciphertext(tmp);
		scheme.addPolyAndEqual(tmp, sqrMatContext.mvec[0], sqrMatContext.logp);
		scheme.modDownToAndEqual(res, tmp->logq);
		Ciphertext* x = sqrMatMult(tmp, res, logp, n);
		delete res;
		delete tmp;
		res = x;
	}
	return res;
}

//-----------------------------------------

Ciphertext* SchemeAlgo::inverse(Ciphertext* cipher, long logp, long steps) {
	Ciphertext* cbar = scheme.negate(cipher);
	scheme.addConstAndEqual(cbar, 1.0, logp);
	Ciphertext* cpow = new Ciphertext(cbar);
	Ciphertext* tmp = scheme.addConst(cbar, 1.0, logp);
	scheme.modDownByAndEqual(tmp, logp);
	Ciphertext* res = new Ciphertext(tmp);

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, logp);
		delete tmp;
		tmp = new Ciphertext(cpow);
		scheme.addConstAndEqual(tmp, 1.0, logp);
		scheme.multAndEqual(tmp, res);
		scheme.reScaleByAndEqual(tmp, logp);
		res = new Ciphertext(tmp);
	}
	delete tmp; delete cpow; delete cbar;

	return res;
}

Ciphertext** SchemeAlgo::inverseExtended(Ciphertext* cipher, const long logp, const long steps) {
	Ciphertext** res = new Ciphertext*[steps];
	Ciphertext* cpow = cipher;
	Ciphertext* tmp = scheme.addConst(cipher, 1.0, logp);
	scheme.modDownByAndEqual(tmp, logp);
	res[0] = new Ciphertext(tmp);

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, logp);
		delete tmp;
		tmp = new Ciphertext(cpow);
		scheme.addConstAndEqual(tmp, 1.0, logp);
		scheme.multAndEqual(tmp, res[i - 1]);
		scheme.reScaleByAndEqual(tmp, logp);
		res[i] = new Ciphertext(tmp);
	}
	delete tmp; delete cpow;
	return res;
}

//-----------------------------------------

Ciphertext* SchemeAlgo::function(Ciphertext* cipher, string& funcName, long logp, long degree) {
	Ciphertext** cpows = powerExtended(cipher, logp, degree);
	long dlogp = 2 * logp;
	double* coeffs = taylorCoeffsMap.at(funcName);
	Ciphertext* res = scheme.multConst(cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(res, coeffs[0], dlogp);

	for (int i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext* aixi = scheme.multConst(cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(res, aixi->logq);
			scheme.addAndEqual(res, aixi);
			delete aixi;
		}
	}
	scheme.reScaleByAndEqual(res, logp);
	for (long i = 0; i < degree; ++i) {
		delete cpows[i];
	}
	delete[] cpows;
	return res;
}

Ciphertext* SchemeAlgo::functionLazy(Ciphertext* cipher, string& funcName, long logp, long degree) {
	Ciphertext** cpows = powerExtended(cipher, logp, degree);

	long dlogp = 2 * logp;

	double* coeffs = taylorCoeffsMap.at(funcName);

	Ciphertext* res = scheme.multConst(cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(res, coeffs[0], dlogp);

	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext* aixi = scheme.multConst(cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(res, aixi->logq);
			scheme.addAndEqual(res, aixi);
			delete aixi;
		}
	}
	for (long i = 0; i < degree; ++i) {
		delete cpows[i];
	}
	delete[] cpows;

	return res;
}

Ciphertext** SchemeAlgo::functionExtended(Ciphertext* cipher, string& funcName, long logp, long degree) {
	Ciphertext** cpows = powerExtended(cipher, logp, degree);

	long dlogp = 2 * logp;
	double* coeffs = taylorCoeffsMap.at(funcName);
	Ciphertext** res = new Ciphertext*[degree];

	res[0] = scheme.multConst(cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(res[0], coeffs[0], dlogp);
	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			res[i] = scheme.multConst(cpows[i], coeffs[i + 1], logp);
			Ciphertext* ctmp = scheme.modDownTo(res[i - 1], res[i]->logq);
			scheme.addAndEqual(res[i], ctmp);
			delete ctmp;
		} else {
			res[i] = new Ciphertext(res[i - 1]);
		}
	}
	NTL_EXEC_RANGE(degree, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(res[i], logp);
	}
	NTL_EXEC_RANGE_END;
	return res;
}

void SchemeAlgo::bitReverse(Ciphertext** ciphers, const long n) {
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

void SchemeAlgo::DFTX0(Ciphertext** ciphers, long n0) {
	bitReverse(ciphers, n0);
	for (long len = 2; len <= n0; len <<= 1) {
		long lenh = len >> 1;
		long shift = M0 / len;
		for (long i = 0; i < n0; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext* u = new Ciphertext(ciphers[i + j]);
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], shift * j, 0);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
				delete u;
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTX0Lazy(Ciphertext** ciphers, long n0) {
	bitReverse(ciphers, n0);
	for (long len = 2; len <= n0; len <<= 1) {
		long lenh = len >> 1;
		long shift = M0 - M0 / len;
		for (long i = 0; i < n0; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext* u = new Ciphertext(ciphers[i + j]);
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], shift * j, 0);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
				delete u;
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTX0(Ciphertext** ciphers, long n0) {
	IDFTX0Lazy(ciphers, n0);

	long lognx = log2((double) n0);
	NTL_EXEC_RANGE(n0, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(ciphers[i], lognx);
	}
	NTL_EXEC_RANGE_END;
}

void SchemeAlgo::DFTX1(Ciphertext** ciphers, long n1) {
	bitReverse(ciphers, n1);
	for (long len = 2; len <= n1; len <<= 1) {
		long lenh = len >> 1;
		long shift = M1 / len;
		for (long i = 0; i < n1; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext* u = new Ciphertext(ciphers[i + j]);
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], 0, shift * j);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
				delete u;
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTX1Lazy(Ciphertext** ciphers, long n1) {
	bitReverse(ciphers, n1);
	for (long len = 2; len <= n1; len <<= 1) {
		long lenh = len >> 1;
		long shift = M1 - M1 / len;
		for (long i = 0; i < n1; i += len) {
			NTL_EXEC_RANGE(lenh, first, last);
			for (long j = first; j < last; ++j) {
				Ciphertext* u = new Ciphertext(ciphers[i + j]);
				scheme.multByMonomialAndEqual(ciphers[i + j + lenh], 0, shift * j);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + lenh]);
				scheme.subAndEqual2(u, ciphers[i + j + lenh]);
				delete u;
			}
			NTL_EXEC_RANGE_END;
		}
	}
}

void SchemeAlgo::IDFTX1(Ciphertext** ciphers, long n1) {
	IDFTX1Lazy(ciphers, n1);

	long logny = log2((double) n1);
	NTL_EXEC_RANGE(n1, first, last);
	for (long i = first; i < last; ++i) {
		scheme.reScaleByAndEqual(ciphers[i], logny);
	}
	NTL_EXEC_RANGE_END;
}

void SchemeAlgo::DFT(Ciphertext** ciphers, long n0, long n1) {
	for (long iy = 0; iy < n1; ++iy) {
		DFTX0(ciphers + iy * n0, n0);
	}

	Ciphertext** tmp = new Ciphertext*[n1];
	for (long ix = 0; ix < n0; ++ix) {
		for (long iy = 0; iy < n1; ++iy) {
			tmp[iy] = ciphers[ix + iy * n0];
		}

		DFTX1(tmp, n1);

		for (long iy = 0; iy < n1; ++iy) {
			ciphers[ix + iy * n0] = tmp[iy];
			delete tmp[iy];
		}
	}
	delete[] tmp;
}

void SchemeAlgo::IDFTLazy(Ciphertext** ciphers, long n0, long n1) {
	for (long iy = 0; iy < n1; ++iy) {
		IDFTX0Lazy(ciphers + iy * n0, n0);
	}

	Ciphertext** tmp = new Ciphertext*[n1];
	for (long ix = 0; ix < n0; ++ix) {
		for (long iy = 0; iy < n1; ++iy) {
			tmp[iy] = ciphers[ix + iy * n0];
		}

		IDFTX1Lazy(tmp, n1);

		for (long iy = 0; iy < n1; ++iy) {
			ciphers[ix + iy * n0] = tmp[iy];
			delete tmp[iy];
		}
	}
	delete[] tmp;
}

void SchemeAlgo::IDFT(Ciphertext** ciphers, long n0, long n1) {
	IDFTLazy(ciphers, n0, n1);

	long n = n0 * n1;
	long logn = log2((double) n);
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		scheme.divPo2AndEqual(ciphers[i], logn);
	}
	NTL_EXEC_RANGE_END;
}
