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

//	TestScheme::testEncrypt(8, 8, 1200, 50, 2, 2);
//	TestScheme::testEncryptSingle(6, 8, 300, 30);
//	TestScheme::testStandard(10, 4, 1200, 50, 1, 3);
//	TestScheme::testimult(8, 8, 300, 30, 2, 2);


//----------------------------------------------------------------------------------
//   ROTATION & CONJUGATION & TRANSPOSITION TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testRotateFast(8, 8, 1200, 50, 2, 3, 1, 2);
//	TestScheme::testRotate(8, 8, 300, 30, 2, 2, 1, 1);
//	TestScheme::testConjugate(8, 8, 300, 30, 2, 2);


//----------------------------------------------------------------------------------
//   POWER & PRODUCT TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testPowerOf2(8, 8, 300, 30, 2, 2, 4);
//	TestScheme::testPower(8, 8, 300, 30, 2, 2, 13);
//	TestScheme::testProdOfPo2(8, 8, 300, 30, 2, 2, 4);
//	TestScheme::testProd(8, 8, 300, 30, 2, 2, 13);


//----------------------------------------------------------------------------------
//   FUNCTION TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testInverse(8, 8, 300, 25, 2, 2, 5);
//	TestScheme::testLogarithm(8, 8, 300, 30, 2, 2, 7);
//	TestScheme::testExponent(8, 8, 300, 30, 2, 2, 7);
//	TestScheme::testSigmoid(8, 8, 300, 30, 2, 2, 7);
//	TestScheme::testSigmoidLazy(8, 8, 300, 30, 2, 2, 7);


//----------------------------------------------------------------------------------
//   MATRIX TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testSqrMatMult(8, 8, 100, 30, 2);
//	TestScheme::testSqrMatPow(8, 8, 300, 30, 2, 4);
//	TestScheme::testMatInv(8, 8, 300, 25, 2, 4);


//----------------------------------------------------------------------------------
//   OTHER TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testBootstrap(8, 8, 41, 1240, 31, 7, 8, 7, 4);

//	TestScheme::testCiphertextWriteAndRead(10, 65, 30, 2);

	TestScheme::test();

	return 0;
}
