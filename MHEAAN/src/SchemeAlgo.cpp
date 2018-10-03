/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "SchemeAlgo.h"


void SchemeAlgo::powerOf2(Ciphertext* res, Ciphertext* cipher, const long logp, const long logDegree) {
	res->copy(cipher);
	for (long i = 0; i < logDegree; ++i) {
		scheme.squareAndEqual(res);
		scheme.reScaleByAndEqual(res, logp);
	}
}

void SchemeAlgo::powerOf2Extended(Ciphertext** res, Ciphertext* cipher, const long logp, const long logDegree) {
	res[0] = new Ciphertext(cipher);
	for (long i = 1; i < logDegree + 1; ++i) {
		res[i] = new Ciphertext(res[i - 1]);
		scheme.squareAndEqual(res[i]);
		scheme.reScaleByAndEqual(res[i], logp);
	}
}

//-----------------------------------------

void SchemeAlgo::power(Ciphertext* res, Ciphertext* cipher, const long logp, const long degree) {
	long logDegree = log2((double) degree);
	long po2Degree = 1 << logDegree;
	powerOf2(res, cipher, logp, logDegree);
	long remDegree = degree - po2Degree;
	if (remDegree > 0) {
		Ciphertext* tmp = new Ciphertext();
		power(tmp, cipher, logp, remDegree);
		scheme.modDownToAndEqual(tmp, res->logq);
		scheme.multAndEqual(res, tmp);
		scheme.reScaleByAndEqual(res, logp);
		delete tmp;
	}
}

void SchemeAlgo::powerExtended(Ciphertext** res, Ciphertext* cipher, const long logp, const long degree) {
	long logDegree = log2((double) degree);
	Ciphertext** cpows = new Ciphertext*[logDegree + 1];
	powerOf2Extended(cpows, cipher, logp, logDegree);
	long idx = 0;
	for (long i = 0; i < logDegree; ++i) {
		long powi = (1 << i);
		res[idx++] = new Ciphertext(cpows[i]);
		for (long j = 0; j < powi - 1; ++j) {
			res[idx] = new Ciphertext(res[j]);
			scheme.modDownToAndEqual(res[idx], cpows[i]->logq);
			scheme.multAndEqual(res[idx], cpows[i]);
			scheme.reScaleByAndEqual(res[idx++], logp);
		}
	}
	res[idx++] = new Ciphertext(cpows[logDegree]);
	long degree2 = (1 << logDegree);
	for (int i = 0; i < (degree - degree2); ++i) {
		res[idx] = new Ciphertext(res[i]);
		scheme.modDownToAndEqual(res[idx], cpows[logDegree]->logq);
		scheme.multAndEqual(res[idx], cpows[logDegree]);
		scheme.reScaleByAndEqual(res[idx++], logp);
	}

	for (long i = 0; i < logDegree; ++i) {
		delete cpows[i];
	}
	delete[] cpows;
}


//----------------------------------------------------------------------------------
//   METHODS ON ARRAYS OF CIPHERTEXTS
//----------------------------------------------------------------------------------


void SchemeAlgo::transpose(Ciphertext* res, Ciphertext* cipher, long logp, long n) {
	long logn = log2(n);
	mutex m;
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	res->free();
	res->copyParams(cipher);
	res->logp += sqrMatContext.logp;
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext* tmp = new Ciphertext(cipher);
		scheme.multPolyAndEqual(tmp, sqrMatContext.mvec[i], sqrMatContext.logp);
		if(i > 0) scheme.leftRotateAndEqual(tmp,i, N1 - i);
		m.lock();
		scheme.addAndEqual(res, tmp);
		m.unlock();
		delete tmp;
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(res, sqrMatContext.logp);
}

void SchemeAlgo::sqrMatMult(Ciphertext* res, Ciphertext* cipher1, Ciphertext* cipher2, long logp, long n) {
	long logn = log2(n);
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	mutex m;
	res->free();
	res->copyParams(cipher1);
	res->logp += cipher2->logp;
	res->logq -= sqrMatContext.logp;

	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext* tmp2 = new Ciphertext(cipher2);
		scheme.multPolyAndEqual(tmp2, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(tmp2, sqrMatContext.logp);
		Ciphertext* aux = new Ciphertext();
		for (long j = 0; j < logn; ++j) {
			scheme.leftRotate(aux, tmp2, 0, (1 << j));
			scheme.addAndEqual(tmp2, aux);
		}
		aux->copy(cipher1);
		if(i > 0) scheme.rightRotateAndEqual(aux, i, 0);
		scheme.modDownByAndEqual(aux, sqrMatContext.logp);
		scheme.multAndEqual(aux, tmp2);
		delete tmp2;
		m.lock();
		scheme.addAndEqual(res, aux);
		m.unlock();
		delete aux;
	}
	NTL_EXEC_RANGE_END;
	scheme.reScaleByAndEqual(res, logp);
}

void SchemeAlgo::sqrMatSqr(Ciphertext* res, Ciphertext* cipher, long logp, long n) {
	long logn = log2(n);
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	mutex m;
	res->free();
	res->copyParams(cipher);
	res->logp *= 2;
	res->logq -= sqrMatContext.logp;
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext* tmp = new Ciphertext(cipher);
		scheme.multPolyAndEqual(tmp, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(tmp, sqrMatContext.logp);

		Ciphertext* aux = new Ciphertext();
		for (long j = 0; j < logn; ++j) {
			scheme.leftRotate(aux, tmp, 0, (1 << j));
			scheme.addAndEqual(tmp, aux);
		}
		aux->copy(cipher);
		if (i > 0) scheme.rightRotateAndEqual(aux, i, 0);
		scheme.modDownByAndEqual(aux, sqrMatContext.logp);
		scheme.multAndEqual(aux, tmp);
		delete tmp;
		m.lock();
		scheme.addAndEqual(res, aux);
		m.unlock();
		delete aux;
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(res, logp);
}


void SchemeAlgo::matInv(Ciphertext* res, Ciphertext* cipher, long logp, long n, long r) {
	long logn = log2(n);
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);

	Ciphertext* cbar = new Ciphertext();
	scheme.negate(cbar, cipher);
	scheme.addPolyAndEqual(cbar, sqrMatContext.mvec[0], sqrMatContext.logp);

	Ciphertext* cpow = new Ciphertext(cbar);
	scheme.addPoly(res, cbar, sqrMatContext.mvec[0], sqrMatContext.logp);
	Ciphertext* x = new Ciphertext();
	Ciphertext* tmp = new Ciphertext();
	for (long i = 1; i < r; ++i) {
		sqrMatSqr(tmp, cpow, logp, n);
		cpow->copy(tmp);
		scheme.addPolyAndEqual(tmp, sqrMatContext.mvec[0], sqrMatContext.logp);
		scheme.modDownToAndEqual(res, tmp->logq);
		sqrMatMult(x, tmp, res, logp, n);
		res->copy(x);
	}
	delete cbar; delete cpow; delete tmp; delete x;
}

//-----------------------------------------

void SchemeAlgo::inverse(Ciphertext* res, Ciphertext* cipher, long logp, long steps) {
	Ciphertext* cbar = new Ciphertext();
	scheme.negate(cbar, cipher);
	scheme.addConstAndEqual(cbar, 1.0, logp);
	Ciphertext* cpow = new Ciphertext(cbar);
	Ciphertext* tmp = new Ciphertext(cbar);
	scheme.addConstAndEqual(tmp, 1.0, logp);
	scheme.modDownByAndEqual(tmp, logp);
	res->copy(tmp);

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, logp);
		tmp->copy(cpow);
		scheme.addConstAndEqual(tmp, 1.0, logp);
		scheme.multAndEqual(tmp, res);
		scheme.reScaleByAndEqual(tmp, logp);
		res->copy(tmp);
	}
	delete tmp; delete cpow; delete cbar;
}

//-----------------------------------------

void SchemeAlgo::function(Ciphertext* res, Ciphertext* cipher, string& funcName, long logp, long degree) {
	Ciphertext** cpows = new Ciphertext*[degree];
	powerExtended(cpows, cipher, logp, degree);
	long dlogp = 2 * logp;
	double* coeffs = taylorCoeffsMap.at(funcName);
	scheme.multConst(res, cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(res, coeffs[0], dlogp);

	Ciphertext* aixi = new Ciphertext();
	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			scheme.multConst(aixi, cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(res, aixi->logq);
			scheme.addAndEqual(res, aixi);
		}
	}
	delete aixi;

	scheme.reScaleByAndEqual(res, logp);
	for (long i = 0; i < degree; ++i) {
		delete cpows[i];
	}
	delete[] cpows;
}

void SchemeAlgo::functionLazy(Ciphertext* res, Ciphertext* cipher, string& funcName, long logp, long degree) {
	Ciphertext** cpows = new Ciphertext*[degree];
	powerExtended(cpows, cipher, logp, degree);
	long dlogp = 2 * logp;
	double* coeffs = taylorCoeffsMap.at(funcName);
	scheme.multConst(res, cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(res, coeffs[0], dlogp);
	Ciphertext* aixi = new Ciphertext();
	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			scheme.multConst(aixi, cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(res, aixi->logq);
			scheme.addAndEqual(res, aixi);
		}
	}
	delete aixi;
	for (long i = 0; i < degree; ++i) {
		delete cpows[i];
	}
	delete[] cpows;
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
