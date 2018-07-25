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

#include "Common.h"
#include "EvaluatorUtils.h"
#include "Plaintext.h"
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
	//   ARRAY ENCRYPTION & DECRYPTION
	//----------------------------------------------------------------------------------


	Ciphertext* encryptSingleArray(complex<double>* vals, long size, long logp);

	complex<double>* decryptSingleArray(SecretKey& secretKey, Ciphertext* ciphers, long size);


	//----------------------------------------------------------------------------------
	//   POWERS & PRODUCTS
	//----------------------------------------------------------------------------------


	Ciphertext powerOf2(Ciphertext& cipher, long precisionBits, long logDegree);

	Ciphertext* powerOf2Extended(Ciphertext& cipher, long logp, long logDegree);

	Ciphertext power(Ciphertext& cipher, long logp, long degree);

	Ciphertext* powerExtended(Ciphertext& cipher, long logp, long degree);

	Ciphertext prodOfPo2(Ciphertext* ciphers, long logp, long logDegree);

	Ciphertext prod(Ciphertext* ciphers, long logp, long degree);


	//----------------------------------------------------------------------------------
	//   FUNCTIONS
	//----------------------------------------------------------------------------------


	Ciphertext inverse(Ciphertext& cipher, long logp, long steps);

	Ciphertext* inverseExtended(Ciphertext& cipher, long logp, long steps);

	Ciphertext function(Ciphertext& cipher, string& funcName, long logp, long degree);

	Ciphertext functionLazy(Ciphertext& cipher, string& funcName, long logp, long degree);

	Ciphertext* functionExtended(Ciphertext& cipher, string& funcName, long logp, long degree);


	//----------------------------------------------------------------------------------
	//   METHODS ON ARRAYS OF CIPHERTEXTS
	//----------------------------------------------------------------------------------


	Ciphertext sum(Ciphertext* ciphers, long size);

	Ciphertext* multVec(Ciphertext* ciphers1, Ciphertext* ciphers2, long size);

	void multAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, long size);

	Ciphertext* multAndModSwitchVec(Ciphertext* ciphers1, Ciphertext* ciphers2, long logp, long size);

	void multModSwitchAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, long logp, long size);

	Ciphertext innerProd(Ciphertext* ciphers1, Ciphertext* ciphers2, long logp, long size);


	//----------------------------------------------------------------------------------
	//   MATRIX
	//----------------------------------------------------------------------------------


	Ciphertext squareMatMult(Ciphertext& cipher1, Ciphertext& cipher2, long logp, long size);

	void squareMatMultAndEqual(Ciphertext& cipher, long logp, long size);

	Ciphertext matMult(Ciphertext& cipher1, Ciphertext& cipher2, long logp, long sizex, long sizey, long sizez);

	Ciphertext matInv(Ciphertext& cipher, long logp, long size, long r);


	//----------------------------------------------------------------------------------
	//   FFT & IFFT
	//----------------------------------------------------------------------------------

	void bitReverse(Ciphertext* ciphers, long n);

	void DFTX(Ciphertext* ciphers, long nx);
	void IDFTX(Ciphertext* ciphers, long nx);
	void IDFTXLazy(Ciphertext* ciphers, long nx);

	void DFTY(Ciphertext* ciphers, long ny);
	void IDFTY(Ciphertext* ciphers, long ny);
	void IDFTYLazy(Ciphertext* ciphers, long ny);

	void DFTXY(Ciphertext* ciphers, long nx, long ny);
	void IDFTXY(Ciphertext* ciphers, long nx, long ny);
	void IDFTXYLazy(Ciphertext* ciphers, long nx, long ny);

};

#endif
