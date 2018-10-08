/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MHEAAN_SCHEMEALGO_H_
#define MHEAAN_SCHEMEALGO_H_

#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <complex>

#include "EvaluatorUtils.h"
#include "SecretKey.h"
#include "Ciphertext.h"
#include "Scheme.h"

static string LOGARITHM = "Logarithm"; ///< log(x)
static string EXPONENT  = "Exponent"; ///< exp(x)
static string SIGMOID   = "Sigmoid"; ///< sigmoid(x) = exp(x) / (1 + exp(x))

class SchemeAlgo {
public:
	Scheme& scheme;

	map<string, double*> taylorCoeffsMap; ///< storing taylor coefficients for function calculation

	SchemeAlgo(Scheme& scheme) : scheme(scheme) {
		taylorCoeffsMap.insert(pair<string, double*>(LOGARITHM, new double[11]{0,1,-0.5,1./3,-1./4,1./5,-1./6,1./7,-1./8,1./9,-1./10}));
		taylorCoeffsMap.insert(pair<string, double*>(EXPONENT, new double[11]{1,1,0.5,1./6,1./24,1./120,1./720,1./5040, 1./40320,1./362880,1./3628800}));
		taylorCoeffsMap.insert(pair<string, double*>(SIGMOID, new double[11]{1./2,1./4,0,-1./48,0,1./480,0,-17./80640,0,31./1451520,0}));
	};


	//----------------------------------------------------------------------------------
	//   POWERS & PRODUCTS
	//----------------------------------------------------------------------------------


	void powerOf2AndEqual(Ciphertext& cipher, long precisionBits, long logDegree);

	void powerOf2Extended(Ciphertext* res, Ciphertext& cipher, long logp, long logDegree);

	void powerAndEqual(Ciphertext& cipher, long logp, long degree);

	void powerExtended(Ciphertext* res, Ciphertext& cipher, long logp, long degree);


	//----------------------------------------------------------------------------------
	//   FUNCTIONS
	//----------------------------------------------------------------------------------


	void inverseAndEqual(Ciphertext& cipher, long logp, long steps);

	void functionAndEqual(Ciphertext& cipher, string& funcName, long logp, long degree);

	void functionLazyAndEqual(Ciphertext& cipher, string& funcName, long logp, long degree);


	//----------------------------------------------------------------------------------
	//   METHODS ON ARRAYS OF CIPHERTEXTS
	//----------------------------------------------------------------------------------


	void transpose(Ciphertext& res, Ciphertext& cipher, long logp, long size);

	void sqrMatMult(Ciphertext& res, Ciphertext& cipher1, Ciphertext& cipher2, long logp, long size);

	void sqrMatSqr(Ciphertext& res, Ciphertext& cipher, long logp, long size);

	void matInv(Ciphertext& res, Ciphertext& cipher, long logp, long n, long r);

};

#endif
