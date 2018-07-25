/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "TestScheme.h"

int main() {


//----------------------------------------------------------------------------------
//   STANDARD TESTS
//----------------------------------------------------------------------------------

	TestScheme::testEncode(8, 1200, 30, 2);
//	TestScheme::testEncrypt(8, 8, 1200, 30, 2, 2);
//	TestScheme::testEncryptSingle(7, 7, 300, 30);
//	TestScheme::testStandard(8, 8, 1200, 30, 2, 2);
//	TestScheme::testimult(7, 7, 300, 30, 2, 2);


//----------------------------------------------------------------------------------
//   ROTATION & CONJUGATION & TRANSPOSITION TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testRotateFast(8, 8, 1200, 30, 4, 2, 1, 1);
//	TestScheme::testRotate(7, 7, 300, 30, 2, 2, 1, 1);
//	TestScheme::testConjugate(7, 7, 300, 30, 2, 2);
//	TestScheme::testTranspose(7, 7, 300, 30, 2, 2);


//----------------------------------------------------------------------------------
//   POWER & PRODUCT TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testPowerOf2(7, 7, 300, 30, 2, 2, 4);
//	TestScheme::testPower(7, 7, 300, 30, 2, 2, 13);
//	TestScheme::testProdOfPo2(7, 7, 300, 30, 2, 2, 4);
//	TestScheme::testProd(7, 7, 300, 30, 2, 2, 13);


//----------------------------------------------------------------------------------
//   FUNCTION TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testInverse(7, 7, 300, 25, 2, 2, 5);
//	TestScheme::testLogarithm(7, 7, 300, 30, 2, 2, 7);
//	TestScheme::testExponent(7, 7, 300, 30, 2, 2, 7);
//	TestScheme::testSigmoid(7, 7, 300, 30, 2, 2, 7);
//	TestScheme::testSigmoidLazy(7, 7, 300, 30, 2, 2, 7);


//----------------------------------------------------------------------------------
//   MATRIX TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testSquareMatMult(7, 7, 100, 30, 2);
//	TestScheme::testSquareMatPow(7, 7, 300, 30, 2, 4);
//	TestScheme::testSquareMatInv(7, 7, 300, 25, 2, 4);
//	TestScheme::testMatMult(7, 7, 300, 30, 2, 3, 4);


//----------------------------------------------------------------------------------
//   FFT TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testFFTBatch(7, 7, 300, 30, 2, 2, 2, 2);
//	TestScheme::testFFTBatchLazy(7, 7, 300, 30, 2, 2, 2, 2);
//	TestScheme::testFFTBatchLazyMultipleHadamard(7, 7, 300, 30, 2, 2, 2, 2, 2);


//----------------------------------------------------------------------------------
//   OTHER TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testCiphertextWriteAndRead(10, 65, 30, 2);

	return 0;
}
