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

//	TestScheme::testEncrypt(300, 30, 7, 8);
//	TestScheme::testEncryptSingle(300, 30);
//	TestScheme::testStandard(300, 30, 2, 2);
//	TestScheme::testimult(300, 30, 2, 2);


//----------------------------------------------------------------------------------
//   ROTATION & CONJUGATION & TRANSPOSITION TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testRotateFast(300, 30, 3, 3, 1, 0);
//	TestScheme::testConjugate(300, 30, 2, 2);


//----------------------------------------------------------------------------------
//   POWER & PRODUCT TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testPowerOf2(300, 30, 2, 2, 4);
//	TestScheme::testPower(300, 30, 2, 2, 13);


//----------------------------------------------------------------------------------
//   FUNCTION TESTS
//----------------------------------------------------------------------------------


//	TestScheme::testInverse(300, 25, 2, 2, 5);
//	TestScheme::testLogarithm(300, 30, 2, 2, 7);
//	TestScheme::testExponent(300, 30, 2, 2, 7);
//	TestScheme::testSigmoid(300, 30, 2, 2, 7);
//	TestScheme::testSigmoidLazy(300, 30, 2, 2, 7);


//----------------------------------------------------------------------------------
//   MATRIX TESTS
//----------------------------------------------------------------------------------

//	TestScheme::testTranspose(65, 30, 6);
//	TestScheme::testSqrMatMult(300, 30, 6);
//	TestScheme::testSqrMatPow(300, 30, 6, 4);
//	TestScheme::testMatInv(300, 30, 6, 4);


//----------------------------------------------------------------------------------
//   OTHER TESTS
//----------------------------------------------------------------------------------


	TestScheme::testBootstrap(40, 35, 7, 8, 3, 4);
//	TestScheme::testCiphertextWriteAndRead(10, 65, 30, 2);
//	TestScheme::test();

	return 0;
}
