/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MHEAAN_TESTSCHEME_H_
#define MHEAAN_TESTSCHEME_H_

class TestScheme {
public:


	//----------------------------------------------------------------------------------
	//   STANDARD TESTS
	//----------------------------------------------------------------------------------


	static void testEncrypt(long logN0, long logN1, long logQ, long logp, long logn0, long logn1);

	static void testEncryptSingle(long logN0, long logN1, long logQ, long logp);

	static void testStandard(long logN0, long logN1, long logQ, long logp, long logn0, long logn1);

	static void testimult(long logN0, long logN1, long logQ, long logp, long logn0, long logn1);


	//----------------------------------------------------------------------------------
	//   ROTATION & CONJUGATION & TRANSPOSITION TESTS
	//----------------------------------------------------------------------------------


	static void testRotateFast(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long r0, long r1);

	static void testRotate(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long r0, long r1);

	static void testConjugate(long logN0, long logN1, long logQ, long logp, long logn0, long logn1);


	//----------------------------------------------------------------------------------
	//   POWER & PRODUCT TESTS
	//----------------------------------------------------------------------------------


	static void testPowerOf2(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long logDegree);

	static void testPower(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long degree);

	static void testProdOfPo2(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long logDegree);

	static void testProd(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long degree);


	//----------------------------------------------------------------------------------
	//   FUNCTION TESTS
	//----------------------------------------------------------------------------------


	static void testInverse(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long invSteps);

	static void testLogarithm(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long degree);

	static void testExponent(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long degree);

	static void testExponentLazy(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long degree);

	static void testSigmoid(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long degree);

	static void testSigmoidLazy(long logN0, long logN1, long logQ, long logp, long logn0, long logn1, long degree);


	//----------------------------------------------------------------------------------
	//   MATRIX TESTS
	//----------------------------------------------------------------------------------


	static void testSqrMatMult(long logN0, long logN1, long logQ, long logp, long logn);

	static void testSqrMatPow(long logN0, long logN1, long logQ, long logp, long logn, long logDegree);

	static void testMatInv(long logN0, long logN1, long logQ, long logp, long logn, long steps);


	//----------------------------------------------------------------------------------
	//   OTHER TESTS
	//----------------------------------------------------------------------------------

	static void testBootstrap(long logN0, long logN1, long logq, long logQ, long logp, long logn0, long logn1, long logT, long logI);

	static void testCiphertextWriteAndRead(long logN0, long logN1, long logQ, long logp, long logn0, long logn1);

	static void test();

	static void test2();

};

#endif
