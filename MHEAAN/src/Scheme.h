/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MHEAAN_SCHEME_H_
#define MHEAAN_SCHEME_H_

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include "SecretKey.h"
#include "Ciphertext.h"
#include "Plaintext.h"
#include "Key.h"
#include "EvaluatorUtils.h"
#include "Ring.h"

using namespace std;
using namespace NTL;

static long ENCRYPTION = 0;
static long MULTIPLICATION  = 1;
static long CONJUGATION = 2;

class Scheme {
private:
public:
	Ring& ring;

	map<long, Key> keyMap; ///< contain Encryption, Multiplication and Conjugation keys, if generated

	map<pair<long, long>, Key> leftRotKeyMap;

	Scheme(SecretKey& secretKey, Ring& ring);


	//----------------------------------------------------------------------------------
	//   KEYS GENERATION
	//----------------------------------------------------------------------------------


	void addEncKey(SecretKey& secretKey);

	void addMultKey(SecretKey& secretKey);

	void addLeftRotKey(SecretKey& secretKey, long rx, long ry);

	void addLeftXRotKeys(SecretKey& secretKey);

	void addLeftYRotKeys(SecretKey& secretKey);

	void addRightXRotKeys(SecretKey& secretKey);

	void addRightYRotKeys(SecretKey& secretKey);

	void addConjKey(SecretKey& secretKey);

	void addBootKey(SecretKey& secretKey, long logn0, long logn1, long logp);

	void addSqrMatKeys(SecretKey& secretKey, long n, long logp);


	//----------------------------------------------------------------------------------
	//   ENCODING & DECODING
	//----------------------------------------------------------------------------------


	Plaintext encode(complex<double>* vals, long n0, long n1, long logp, long logq);

	Plaintext encode(double* vals, long n0, long n1, long logp, long logq);

	Plaintext encodeSingle(complex<double> val, long logp, long logq);

	Plaintext encodeSingle(double val, long logp, long logq);

	complex<double>* decode(Plaintext& msg);

	complex<double> decodeSingle(Plaintext& msg);


	//----------------------------------------------------------------------------------
	//   ENCRYPTION & DECRYPTION
	//----------------------------------------------------------------------------------


	Ciphertext encryptMsg(Plaintext& msg);

	Ciphertext encrypt(complex<double>* vals, long n0, long n1, long logp, long logq);

	Ciphertext encrypt(double* vals, long n0, long n1, long logp, long logq);

	Ciphertext encryptSingle(complex<double> val, long logp, long logq);

	Ciphertext encryptSingle(double val, long logp, long logq);

	Ciphertext encryptZeros(long n0, long n1, long logp, long logq);

	Plaintext decryptMsg(SecretKey& secretKey, Ciphertext& cipher);

	complex<double>* decrypt(SecretKey& secretKey, Ciphertext& cipher);

	complex<double> decryptSingle(SecretKey& secretKey, Ciphertext& cipher);


	//----------------------------------------------------------------------------------
	//   HOMOMORPHIC OPERATIONS
	//----------------------------------------------------------------------------------


	Ciphertext negate(Ciphertext& cipher);
	void negateAndEqual(Ciphertext& cipher);

	Ciphertext add(Ciphertext& cipher1, Ciphertext& cipher2);
	void addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	Ciphertext addConst(Ciphertext& cipher, double cnst, long logp = -1);
	Ciphertext addConst(Ciphertext& cipher, RR& cnst, long logp = -1);

	void addConstAndEqual(Ciphertext& cipher, RR& cnst, long logp = -1);
	void addConstAndEqual(Ciphertext& cipher, double cnst, long logp = -1);

	Ciphertext addPoly(Ciphertext& cipher, ZZ* poly, long logp);
	void addPolyAndEqual(Ciphertext& cipher, ZZ* poly, long logp);

	Ciphertext sub(Ciphertext& cipher1, Ciphertext& cipher2);
	void subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);
	void subAndEqual2(Ciphertext& cipher1, Ciphertext& cipher2);

	Ciphertext imult(Ciphertext& cipher);
	void imultAndEqual(Ciphertext& cipher);

	Ciphertext idiv(Ciphertext& cipher);
	void idivAndEqual(Ciphertext& cipher);

	Ciphertext mult(Ciphertext& cipher1, Ciphertext& cipher2);
	void multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	Ciphertext square(Ciphertext& cipher);
	void squareAndEqual(Ciphertext& cipher);

	Ciphertext multConst(Ciphertext& cipher, RR& cnst, long logp);
	Ciphertext multConst(Ciphertext& cipher, double cnst, long logp);
	Ciphertext multConst(Ciphertext& cipher, complex<double> cnst, long logp);

	void multConstAndEqual(Ciphertext& cipher, RR& cnst, long logp);
	void multConstAndEqual(Ciphertext& cipher, double cnst, long logp);
	void multConstAndEqual(Ciphertext& cipher, complex<double> cnst, long logp);

	Ciphertext multPolyX0(Ciphertext& cipher, ZZ* poly, long logp);
	void multPolyX0AndEqual(Ciphertext& cipher, ZZ* poly, long logp);
	Ciphertext multPolyNTTX0(Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp);
	void multPolyNTTX0AndEqual(Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp);

	Ciphertext multPolyX1(Ciphertext& cipher, ZZ* rpoly, ZZ* ipoly, long logp);
	void multPolyX1AndEqual(Ciphertext& cipher, ZZ* rpoly, ZZ* ipoly, long logp);

	Ciphertext multPoly(Ciphertext& cipher, ZZ* poly, long logp);
	void multPolyAndEqual(Ciphertext& cipher, ZZ* poly, long logp);
	Ciphertext multPolyNTT(Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp);
	void multPolyNTTAndEqual(Ciphertext& cipher, uint64_t* rpoly, long bnd, long logp);

	Ciphertext multByMonomial(Ciphertext& cipher, const long d0, const long d1);
	void multByMonomialAndEqual(Ciphertext& cipher, const long d0, const long d1);

	Ciphertext multPo2(Ciphertext& cipher, long bits);
	void multPo2AndEqual(Ciphertext& cipher, long bits);

	void doubleAndEqual(Ciphertext& cipher);

	Ciphertext divPo2(Ciphertext& cipher, long logd);
	void divPo2AndEqual(Ciphertext& cipher, long logd);


	//----------------------------------------------------------------------------------
	//   RESCALING & MODULUS DOWN
	//----------------------------------------------------------------------------------


	Ciphertext reScaleBy(Ciphertext& cipher, long dlogq);
	Ciphertext reScaleTo(Ciphertext& cipher, long logq);

	void reScaleByAndEqual(Ciphertext& cipher, long dlogq);
	void reScaleToAndEqual(Ciphertext& cipher, long logq);

	Ciphertext modDownBy(Ciphertext& cipher, long dlogq);
	Ciphertext modDownTo(Ciphertext& cipher, long logq);

	void modDownByAndEqual(Ciphertext& cipher, long dlogq);
	void modDownToAndEqual(Ciphertext& cipher, long logq);


	//----------------------------------------------------------------------------------
	//   ROTATIONS & CONJUGATIONS
	//----------------------------------------------------------------------------------

	Ciphertext leftRotateFast(Ciphertext& cipher, long r0, long r1);
	Ciphertext rightRotateFast(Ciphertext& cipher, long r0, long r1);

	void leftRotateFastAndEqual(Ciphertext& cipher, long r0, long r1);
	void rightRotateFastAndEqual(Ciphertext& cipher, long r0, long r1);

	Ciphertext leftRotate(Ciphertext& cipher, long r0, long r1);
	Ciphertext rightRotate(Ciphertext& cipher, long r0, long r1);

	void leftRotateAndEqual(Ciphertext& cipher, long r0, long r1);
	void rightRotateAndEqual(Ciphertext& cipher, long r0, long r1);

	Ciphertext conjugate(Ciphertext& cipher);
	void conjugateAndEqual(Ciphertext& cipher);


	//----------------------------------------------------------------------------------
	//   BOOTSTRAPPING
	//----------------------------------------------------------------------------------


	void normalizeAndEqual(Ciphertext& cipher);

	void coeffToSlotX0AndEqual(Ciphertext& cipher);
	void coeffToSlotX1AndEqual(Ciphertext& cipher);
	void coeffToSlotAndEqual(Ciphertext& cipher);

	void slotToCoeffX0AndEqual(Ciphertext& cipher);
	void slotToCoeffX1AndEqual(Ciphertext& cipher);
	void slotToCoeffAndEqual(Ciphertext& cipher);

	void exp2piAndEqual(Ciphertext& cipher, long logp);

	void evalExpAndEqual(Ciphertext& cipher, long logT, long logI = 4);

	void bootstrapX0AndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI = 4);
	void bootstrapX1AndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI = 4);
	void bootstrapAndEqual(Ciphertext& cipher, long logq, long logQ, long logT, long logI = 4);

};

#endif
