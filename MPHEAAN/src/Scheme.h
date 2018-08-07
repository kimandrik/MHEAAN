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
#include "Ring2XY.h"

using namespace std;
using namespace NTL;

static long ENCRYPTION = 0;
static long MULTIPLICATION  = 1;
static long CONJUGATION = 2;

class Scheme {
private:
public:
	Ring2XY& ring;

	map<long, Key> keyMap; ///< contain Encryption, Multiplication and Conjugation keys, if generated

	map<pair<long, long>, Key> leftRotKeyMap;

	Scheme(SecretKey& secretKey, Ring2XY& ring);


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

	void addSquareMatrixKeys(SecretKey& secretKey, long nx);


	//----------------------------------------------------------------------------------
	//   ENCODING & DECODING
	//----------------------------------------------------------------------------------


	Plaintext encode(complex<double>* vals, long nx, long ny, long logp, long logq);

	Plaintext encode(double* vals, long nx, long ny, long logp, long logq);

	Plaintext encodeSingle(complex<double> val, long logp, long logq);

	Plaintext encodeSingle(double val, long logp, long logq);

	complex<double>* decode(Plaintext& msg);

	complex<double> decodeSingle(Plaintext& msg);


	//----------------------------------------------------------------------------------
	//   ENCRYPTION & DECRYPTION
	//----------------------------------------------------------------------------------


	Ciphertext encryptMsg(Plaintext& msg);

	Ciphertext encrypt(complex<double>* vals, long nx, long ny, long logp, long logq);

	Ciphertext encrypt(double* vals, long nx, long ny, long logp, long logq);

	Ciphertext encryptSingle(complex<double> val, long logp, long logq);

	Ciphertext encryptSingle(double val, long logp, long logq);

	Ciphertext encryptZeros(long nx, long ny, long logp, long logq);

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

	Ciphertext multByConst(Ciphertext& cipher, RR& cnst, long logp);
	Ciphertext multByConst(Ciphertext& cipher, double cnst, long logp);

	void multByConstAndEqual(Ciphertext& cipher, RR& cnst, long logp);
	void multByConstAndEqual(Ciphertext& cipher, double cnst, long logp);

	Ciphertext multByMonomial(Ciphertext& cipher, const long dx, const long dy);
	void multByMonomialAndEqual(Ciphertext& cipher, const long dx, const long dy);

	Ciphertext multByXPoly(Ciphertext& cipher, ZZ* xpoly, long logp);
	void multByXPolyAndEqual(Ciphertext& cipher, ZZ* xpoly, long logp);

	Ciphertext multByYPoly(Ciphertext& cipher, ZZ* ypoly, long logp);
	void multByYPolyAndEqual(Ciphertext& cipher, ZZ* ypoly, long logp);

	Ciphertext multByPoly(Ciphertext& cipher, ZZ* poly, long logp);
	void multByPolyAndEqual(Ciphertext& cipher, ZZ* poly, long logp);

	Ciphertext multByPo2(Ciphertext& cipher, long bits);
	void multByPo2AndEqual(Ciphertext& cipher, long bits);

	void doubleAndEqual(Ciphertext& cipher);

	Ciphertext divByPo2(Ciphertext& cipher, long logd);
	void divByPo2AndEqual(Ciphertext& cipher, long logd);


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
	//   ROTATIONS & CONJUGATIONS & TRANSPOSITION
	//----------------------------------------------------------------------------------

	Ciphertext leftRotateFast(Ciphertext& cipher, long rx, long ry);
	Ciphertext rightRotateFast(Ciphertext& cipher, long rx, long ry);

	void leftRotateFastAndEqual(Ciphertext& cipher, long rx, long ry);
	void rightRotateFastAndEqual(Ciphertext& cipher, long rx, long ry);

	Ciphertext leftRotate(Ciphertext& cipher, long rx, long ry);
	Ciphertext rightRotate(Ciphertext& cipher, long rx, long ry);

	void leftRotateAndEqual(Ciphertext& cipher, long rx, long ry);
	void rightRotateAndEqual(Ciphertext& cipher, long rx, long ry);

	Ciphertext conjugate(Ciphertext& cipher);
	void conjugateAndEqual(Ciphertext& cipher);

};

#endif
