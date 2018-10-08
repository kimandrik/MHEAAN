/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "SchemeAlgo.h"


void SchemeAlgo::powerOf2AndEqual(Ciphertext& cipher, const long logp, const long logDegree) {
	for (long i = 0; i < logDegree; ++i) {
		scheme.squareAndEqual(cipher);
		scheme.reScaleByAndEqual(cipher, logp);
	}
}

void SchemeAlgo::powerOf2Extended(Ciphertext* res, Ciphertext& cipher, const long logp, const long logDegree) {
	res[0].copy(cipher);
	for (long i = 1; i < logDegree + 1; ++i) {
		res[i].copy(res[i - 1]);
		scheme.squareAndEqual(res[i]);
		scheme.reScaleByAndEqual(res[i], logp);
	}
}

//-----------------------------------------

void SchemeAlgo::powerAndEqual(Ciphertext& cipher, const long logp, const long degree) {
	long logDegree = log2((double) degree);
	long po2Degree = 1 << logDegree;
	Ciphertext res(cipher);
	powerOf2AndEqual(res, logp, logDegree);
	long remDegree = degree - po2Degree;
	if (remDegree > 0) {
		Ciphertext tmp(cipher);
		powerAndEqual(tmp, logp, remDegree);
		scheme.modDownToAndEqual(tmp, res.logq);
		scheme.multAndEqual(res, tmp);
		scheme.reScaleByAndEqual(res, logp);
	}
}

void SchemeAlgo::powerExtended(Ciphertext* res, Ciphertext& cipher, const long logp, const long degree) {
	long logDegree = log2((double) degree);
	Ciphertext* cpows = new Ciphertext[logDegree + 1];
	powerOf2Extended(cpows, cipher, logp, logDegree);
	long idx = 0;
	for (long i = 0; i < logDegree; ++i) {
		long powi = (1 << i);
		res[idx++].copy(cpows[i]);
		for (long j = 0; j < powi - 1; ++j) {
			res[idx].copy(res[j]);
			scheme.modDownToAndEqual(res[idx], cpows[i].logq);
			scheme.multAndEqual(res[idx], cpows[i]);
			scheme.reScaleByAndEqual(res[idx++], logp);
		}
	}
	res[idx++].copy(cpows[logDegree]);
	long degree2 = (1 << logDegree);
	for (int i = 0; i < (degree - degree2); ++i) {
		res[idx].copy(res[i]);
		scheme.modDownToAndEqual(res[idx], cpows[logDegree].logq);
		scheme.multAndEqual(res[idx], cpows[logDegree]);
		scheme.reScaleByAndEqual(res[idx++], logp);
	}

	delete[] cpows;
}

//-----------------------------------------

void SchemeAlgo::inverseAndEqual(Ciphertext& cipher, long logp, long steps) {
	Ciphertext cbar;
	scheme.negate(cbar, cipher);
	scheme.addConstAndEqual(cbar, 1.0, logp);
	Ciphertext cpow(cbar);
	Ciphertext tmp(cbar);
	scheme.addConstAndEqual(tmp, 1.0, logp);
	scheme.modDownByAndEqual(tmp, logp);
	cipher.copy(tmp);

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, logp);
		tmp.copy(cpow);
		scheme.addConstAndEqual(tmp, 1.0, logp);
		scheme.multAndEqual(tmp, cipher);
		scheme.reScaleByAndEqual(tmp, logp);
		cipher.copy(tmp);
	}
}

void SchemeAlgo::functionAndEqual(Ciphertext& cipher, string& funcName, long logp, long degree) {
	Ciphertext* cpows = new Ciphertext[degree];
	powerExtended(cpows, cipher, logp, degree);
	long dlogp = 2 * logp;
	double* coeffs = taylorCoeffsMap.at(funcName);
	scheme.multConst(cipher, cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(cipher, coeffs[0], dlogp);

	Ciphertext aixi;
	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			scheme.multConst(aixi, cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(cipher, aixi.logq);
			scheme.addAndEqual(cipher, aixi);
		}
	}
	scheme.reScaleByAndEqual(cipher, logp);

	delete[] cpows;
}

void SchemeAlgo::functionLazyAndEqual(Ciphertext& cipher, string& funcName, long logp, long degree) {
	Ciphertext* cpows = new Ciphertext[degree];
	powerExtended(cpows, cipher, logp, degree);
	long dlogp = 2 * logp;
	double* coeffs = taylorCoeffsMap.at(funcName);
	scheme.multConst(cipher, cpows[0], coeffs[1], logp);
	scheme.addConstAndEqual(cipher, coeffs[0], dlogp);
	Ciphertext aixi;
	for (long i = 1; i < degree; ++i) {
		if (abs(coeffs[i + 1]) > 1e-27) {
			scheme.multConst(aixi, cpows[i], coeffs[i + 1], logp);
			scheme.modDownToAndEqual(cipher, aixi.logq);
			scheme.addAndEqual(cipher, aixi);
		}
	}
	delete[] cpows;
}

void SchemeAlgo::transpose(Ciphertext& res, Ciphertext& cipher, long logp, long n) {
	long logn = log2(n);
	mutex m;
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	res.free();
	res.copyParams(cipher);
	res.logp += sqrMatContext.logp;
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext tmp(cipher);
		scheme.multPolyAndEqual(tmp, sqrMatContext.mvec[i], sqrMatContext.logp);
		if(i > 0) scheme.leftRotateAndEqual(tmp,i, N1 - i);
		m.lock();
		scheme.addAndEqual(res, tmp);
		m.unlock();
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(res, sqrMatContext.logp);
}

void SchemeAlgo::sqrMatMult(Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2, long logp, long n) {
	long logn = log2(n);
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	mutex m;
	res.free();
	res.copyParams(cipher1);
	res.logp += cipher2.logp;
	res.logq -= sqrMatContext.logp;

	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext tmp2(cipher2);
		scheme.multPolyAndEqual(tmp2, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(tmp2, sqrMatContext.logp);
		Ciphertext aux;
		for (long j = 0; j < logn; ++j) {
			scheme.leftRotate(aux, tmp2, 0, (1 << j));
			scheme.addAndEqual(tmp2, aux);
		}
		aux.copy(cipher1);
		if(i > 0) scheme.rightRotateAndEqual(aux, i, 0);
		scheme.modDownByAndEqual(aux, sqrMatContext.logp);
		scheme.multAndEqual(aux, tmp2);
		m.lock();
		scheme.addAndEqual(res, aux);
		m.unlock();
	}
	NTL_EXEC_RANGE_END;
	scheme.reScaleByAndEqual(res, logp);
}

void SchemeAlgo::sqrMatSqr(Ciphertext& res, Ciphertext& cipher, long logp, long n) {
	long logn = log2(n);
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);
	mutex m;
	res.free();
	res.copyParams(cipher);
	res.logp *= 2;
	res.logq -= sqrMatContext.logp;
	NTL_EXEC_RANGE(n, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext tmp(cipher);
		scheme.multPolyAndEqual(tmp, sqrMatContext.mvec[i], sqrMatContext.logp);
		scheme.reScaleByAndEqual(tmp, sqrMatContext.logp);

		Ciphertext aux;
		for (long j = 0; j < logn; ++j) {
			scheme.leftRotate(aux, tmp, 0, (1 << j));
			scheme.addAndEqual(tmp, aux);
		}
		aux.copy(cipher);
		if (i > 0) scheme.rightRotateAndEqual(aux, i, 0);
		scheme.modDownByAndEqual(aux, sqrMatContext.logp);
		scheme.multAndEqual(aux, tmp);
		m.lock();
		scheme.addAndEqual(res, aux);
		m.unlock();
	}
	NTL_EXEC_RANGE_END;

	scheme.reScaleByAndEqual(res, logp);
}


void SchemeAlgo::matInv(Ciphertext& res, Ciphertext& cipher, long logp, long n, long r) {
	long logn = log2(n);
	SqrMatContext& sqrMatContext = scheme.ring.sqrMatContextMap.at(logn);

	Ciphertext cbar;
	scheme.negate(cbar, cipher);
	scheme.addPolyAndEqual(cbar, sqrMatContext.mvec[0], sqrMatContext.logp);

	Ciphertext cpow(cbar);
	scheme.addPoly(res, cbar, sqrMatContext.mvec[0], sqrMatContext.logp);
	Ciphertext x;
	Ciphertext tmp;
	for (long i = 1; i < r; ++i) {
		sqrMatSqr(tmp, cpow, logp, n);
		cpow.copy(tmp);
		scheme.addPolyAndEqual(tmp, sqrMatContext.mvec[0], sqrMatContext.logp);
		scheme.modDownToAndEqual(res, tmp.logq);
		sqrMatMult(x, tmp, res, logp, n);
		res.copy(x);
	}
}

