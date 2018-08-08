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


	static void testEncrypt(long logNx, long logNy, long logQ, long logp, long lognx, long logny);

	static void testEncryptSingle(long logNx, long logNy, long logQ, long logp);

	static void testStandard(long logNx, long logNy, long logQ, long logp, long lognx, long logny);

	static void testimult(long logNx, long logNy, long logQ, long logp, long lognx, long logny);


	//----------------------------------------------------------------------------------
	//   ROTATION & CONJUGATION & TRANSPOSITION TESTS
	//----------------------------------------------------------------------------------


	static void testRotateFast(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long rx, long ry);

	static void testRotate(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long rx, long ry);

	static void testConjugate(long logNx, long logNy, long logQ, long logp, long lognx, long logny);


	//----------------------------------------------------------------------------------
	//   POWER & PRODUCT TESTS
	//----------------------------------------------------------------------------------


	static void testPowerOf2(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long logDegree);

	static void testPower(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree);

	static void testProdOfPo2(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long logDegree);

	static void testProd(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree);


	//----------------------------------------------------------------------------------
	//   FUNCTION TESTS
	//----------------------------------------------------------------------------------


	static void testInverse(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long invSteps);

	static void testLogarithm(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree);

	static void testExponent(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree);

	static void testExponentLazy(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree);

	static void testSigmoid(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree);

	static void testSigmoidLazy(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree);


	//----------------------------------------------------------------------------------
	//   MATRIX TESTS
	//----------------------------------------------------------------------------------


	static void testSquareMatMult(long logNx, long logNy, long logQ, long logp, long lognx);

	static void testSquareMatPow(long logNx, long logNy, long logQ, long logp, long lognx, long logDegree);

	static void testSquareMatInv(long logNx, long logNy, long logQ, long logp, long lognx, long steps);

	static void testMatMult(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long lognz);


	//----------------------------------------------------------------------------------
	//   OTHER TESTS
	//----------------------------------------------------------------------------------

	static void testBootstrap(long logNx, long logNy, long logq, long logQ, long logp, long lognx, long logny, long logT, long logI = 4);

	static void testCiphertextWriteAndRead(long logNx, long logNy, long logQ, long logp, long lognx, long logny);

	static void test();
};

#endif
