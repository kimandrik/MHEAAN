/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "TestScheme.h"

#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include "Ciphertext.h"
#include "EvaluatorUtils.h"
#include "Ring.h"
#include "Scheme.h"
#include "SchemeAlgo.h"
#include "SecretKey.h"
#include "StringUtils.h"
#include "TimeUtils.h"

using namespace std;
using namespace NTL;


//----------------------------------------------------------------------------------
//   STANDARD TESTS
//----------------------------------------------------------------------------------


void TestScheme::testEncrypt(long logq, long logp, long logn0, long logn1) {
	cout << "!!! START TEST ENCRYPT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);

	timeutils.start("Encode matrix");
	Plaintext* msg = scheme.encode(mmat, n0, n1, logp, logQ);
	timeutils.stop("Encode matrix");

	timeutils.start("Encrypt msg");
	Ciphertext* cipher = scheme.encryptMsg(msg);
	timeutils.stop("Encrypt msg");

	timeutils.start("Decrypt msg");
	Plaintext* dsg = scheme.decryptMsg(secretKey, cipher);
	timeutils.stop("Decrypt msg");

	timeutils.start("Decode matrix");
	complex<double>* dmat = scheme.decode(dsg);
	timeutils.stop("Decode matrix");

	StringUtils::compare(mmat, dmat, n, "val");

	cout << "!!! END TEST ENCRYPT !!!" << endl;
}

void TestScheme::testEncryptSingle(long logq, long logp) {
	cout << "!!! START TEST ENCRYPT SINGLE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	complex<double> mval = EvaluatorUtils::randomComplexSigned();

	timeutils.start("Encrypt single");
	Ciphertext* cipher = scheme.encryptSingle(mval, logp, logQ);
	timeutils.stop("Encrypt single");

	timeutils.start("Decrypt single");
	complex<double> dval = scheme.decryptSingle(secretKey, cipher);
	timeutils.stop("Decrypt single");

	StringUtils::compare(mval, dval, "val");

	cout << "!!! END TEST ENCRYPT SINGLE !!!" << endl;
}

void TestScheme::testStandard(long logq, long logp, long logn0, long logn1) {
	cout << "!!! START TEST STANDARD !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat1 = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* mmat2 = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* madd = new complex<double>[n];
	complex<double>* mmult = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mmult[i] = mmat1[i] * mmat2[i];
		madd[i] = mmat1[i] + mmat2[i];
	}
	Ciphertext* cipher1 = scheme.encrypt(mmat1, n0, n1, logp, logq);
	Ciphertext* cipher2 = scheme.encrypt(mmat2, n0, n1, logp, logq);

	Ciphertext* cadd = new Ciphertext(cipher1);
	timeutils.start("add matrix");
	scheme.addAndEqual(cadd, cipher2);
	timeutils.stop("add matrix");

	Ciphertext* cmult = new Ciphertext(cipher1);
	timeutils.start("mult matrix");
	scheme.multAndEqual(cmult, cipher2);
	scheme.reScaleByAndEqual(cmult, logp);
	timeutils.stop("mult matrix");

	complex<double>* dadd = scheme.decrypt(secretKey, cadd);
	complex<double>* dmult = scheme.decrypt(secretKey, cmult);

	StringUtils::compare(madd, dadd, n, "add");
	StringUtils::compare(mmult, dmult, n, "mult");

	cout << "!!! END TEST STANDARD !!!" << endl;
}

void TestScheme::testimult(long logq, long logp, long logn0, long logn1) {
	cout << "!!! START TEST i MULTIPLICATION !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* mmatimult = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mmatimult[i].real(-mmat[i].imag());
		mmatimult[i].imag(mmat[i].real());
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);

	timeutils.start("Multiplication by i");
	scheme.imultAndEqual(cipher);
	timeutils.stop("Multiplication by i");

	complex<double>* dmatimult = scheme.decrypt(secretKey, cipher);

	StringUtils::compare(mmatimult, dmatimult, n, "imult");

	cout << "!!! END TEST i MULTIPLICATION !!!" << endl;
}


//----------------------------------------------------------------------------------
//   ROTATION & CONJUGATION & TRANSPOSITION TESTS
//----------------------------------------------------------------------------------


void TestScheme::testRotateFast(long logq, long logp, long logn0, long logn1, long r0, long r1) {
	cout << "!!! START TEST ROTATE FAST !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addLeftRotKey(secretKey, r0, r1);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);
	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);

	timeutils.start("Left rotate fast");
	scheme.leftRotateAndEqual(cipher, r0, r1);
	timeutils.stop("Left rotate fast");

	complex<double>* dmat = scheme.decrypt(secretKey, cipher);
	EvaluatorUtils::leftRotateAndEqual(mmat, n0, n1, r0, r1);
	StringUtils::compare(mmat, dmat, n, "val");

	cout << "!!! END TEST ROTATE FAST !!!" << endl;
}

void TestScheme::testConjugate(long logq, long logp, long logn0, long logn1) {
	cout << "!!! START TEST CONJUGATE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addConjKey(secretKey);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* mmatconj = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mmatconj[i] = conj(mmat[i]);
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);

	timeutils.start("Conjugate");
	scheme.conjugateAndEqual(cipher);
	timeutils.stop("Conjugate");

	complex<double>* dmatconj = scheme.decrypt(secretKey, cipher);
	StringUtils::compare(mmatconj, dmatconj, n, "conj");

	cout << "!!! END TEST CONJUGATE !!!" << endl;
}


//----------------------------------------------------------------------------------
//   POWER & PRODUCT TESTS
//----------------------------------------------------------------------------------


void TestScheme::testPowerOf2(long logq, long logp, long logn0, long logn1, long logDegree) {
	cout << "!!! START TEST POWER OF 2 !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;
	long degree = 1 << logDegree;

	complex<double>* mmat = EvaluatorUtils::randomCircleArray(n);
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mpow[i] = pow(mmat[i], degree);
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* cpow = new Ciphertext();

	timeutils.start("Power of 2");
	algo.powerOf2(cpow, cipher, logp, logDegree);
	timeutils.stop("Power of 2");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	StringUtils::compare(mpow, dpow, n, "pow");

	cout << "!!! END TEST POWER OF 2 !!!" << endl;
}

void TestScheme::testPower(long logq, long logp, long logn0, long logn1, long degree) {
	cout << "!!! START TEST POWER !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomCircleArray(n);
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mpow[i] = pow(mmat[i], degree);
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* cpow = new Ciphertext();
	timeutils.start("Power");
	algo.power(cpow, cipher, logp, degree);
	timeutils.stop("Power");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	StringUtils::compare(mpow, dpow, n, "pow");

	cout << "!!! END TEST POWER !!!" << endl;
}


//----------------------------------------------------------------------------------
//   FUNCTION TESTS
//----------------------------------------------------------------------------------


void TestScheme::testInverse(long logq, long logp, long logn0, long logn1, long steps) {
	cout << "!!! START TEST INVERSE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomCircleArray(n, 0.1);
	complex<double>* minv = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		minv[i] = 1. / mmat[i];
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* cinv = new Ciphertext();
	timeutils.start("Inverse");
	algo.inverse(cinv, cipher, logp, steps);
	timeutils.stop("Inverse");

	complex<double>* dinv = scheme.decrypt(secretKey, cinv);
	StringUtils::compare(minv, dinv, n, "inv");

	cout << "!!! END TEST INVERSE !!!" << endl;
}

//-----------------------------------------

void TestScheme::testLogarithm(long logq, long logp, long logn0, long logn1, long degree) {
	cout << "!!! START TEST LOGARITHM !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n, 0.1);
	complex<double>* mlog = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mlog[i] = log(mmat[i] + 1.);
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* clog = new Ciphertext();
	timeutils.start(LOGARITHM);
	algo.function(clog, cipher, LOGARITHM, logp, degree);
	timeutils.stop(LOGARITHM);

	complex<double>* dlog = scheme.decrypt(secretKey, clog);
	StringUtils::compare(mlog, dlog, n, LOGARITHM);

	cout << "!!! END TEST LOGARITHM !!!" << endl;
}

void TestScheme::testExponent(long logq, long logp, long logn0, long logn1, long degree) {
	cout << "!!! START TEST EXPONENT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mmat[i]);
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* cexp = new Ciphertext();
	timeutils.start(EXPONENT);
	algo.function(cexp, cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT);

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	StringUtils::compare(mexp, dexp, n, EXPONENT);

	cout << "!!! END TEST EXPONENT !!!" << endl;
}

void TestScheme::testExponentLazy(long logq, long logp, long logn0, long logn1, long degree) {
	cout << "!!! START TEST EXPONENT LAZY !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mmat[i]);
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* cexp = new Ciphertext();
	timeutils.start(EXPONENT + " lazy");
	algo.functionLazy(cexp, cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT + " lazy");

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	StringUtils::compare(mexp, dexp, n, EXPONENT);

	cout << "!!! END TEST EXPONENT LAZY !!!" << endl;
}

void TestScheme::testSigmoid(long logq, long logp, long logn0, long logn1, long degree) {
	cout << "!!! START TEST SIGMOID !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mmat[i]) / (1. + exp(mmat[i]));
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* csig = new Ciphertext();
	timeutils.start(SIGMOID);
	algo.function(csig, cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID);

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	StringUtils::compare(msig, dsig, n, SIGMOID);

	cout << "!!! END TEST SIGMOID !!!" << endl;
}

void TestScheme::testSigmoidLazy(long logq, long logp, long logn0, long logn1, long degree) {
	cout << "!!! START TEST SIGMOID LAZY !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mmat[i]) / (1. + exp(mmat[i]));
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);
	Ciphertext* csig = new Ciphertext();
	timeutils.start(SIGMOID + " lazy");
	algo.functionLazy(csig, cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID + " lazy");

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	StringUtils::compare(msig, dsig, n, SIGMOID);

	cout << "!!! END TEST SIGMOID LAZY !!!" << endl;
}


//----------------------------------------------------------------------------------
//   MATRIX TESTS
//----------------------------------------------------------------------------------

void TestScheme::testTranspose(long logq, long logp, long logn) {
	cout << "!!! START TEST TRANSPOSE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	scheme.addTransposeKeys(secretKey, logn, logp);

	long n = (1 << logn);
	long n2 = n * n;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n2);
	complex<double>* mt = EvaluatorUtils::transpose(mmat, n);
	Ciphertext* cipher = scheme.encrypt(mmat, n, n, logp, logq);
	Ciphertext* ct = new Ciphertext();
	timeutils.start("Transpose");
	algo.transpose(ct, cipher, logp, n);
	timeutils.stop("Transpose");

	complex<double>* dt = scheme.decrypt(secretKey, ct);
	StringUtils::compare(mt, dt, n2, "matrix");

	cout << ct->logq << endl;

	cout << "!!! END TEST TRANSPOSE !!!" << endl;
}

void TestScheme::testSqrMatMult(long logq, long logp, long logn) {
	cout << "!!! START TEST SQUARE MATRIX !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	scheme.addSqrMatKeys(secretKey, logn, logp);

	long n = (1 << logn);
	long n2 = n * n;

	complex<double>* mmat1 = EvaluatorUtils::randomComplexArray(n2);
	complex<double>* mmat2 = EvaluatorUtils::randomComplexArray(n2);
	complex<double>* mmatmult = EvaluatorUtils::squareMatMult(mmat1, mmat2, n);

	Ciphertext* cipher1 = scheme.encrypt(mmat1, n, n, logp, logq);
	Ciphertext* cipher2 = scheme.encrypt(mmat2, n, n, logp, logq);
	Ciphertext* cmatmult = new Ciphertext();
	timeutils.start("Square Matrix Mult");
	algo.sqrMatMult(cmatmult, cipher1, cipher2, logp, n);
	timeutils.stop("Square Matrix Mult");

	complex<double>* dmatmult = scheme.decrypt(secretKey, cmatmult);
	StringUtils::compare(mmatmult, dmatmult, n2, "matrix");

	cout << cmatmult->logq << endl;
	cout << "!!! END TEST SQUARE MATRIX !!!" << endl;
}

void TestScheme::testSqrMatPow(long logq, long logp, long logn, long logDegree) {
	cout << "!!! START TEST SQUARE MATRIX POW!!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);
	scheme.addSqrMatKeys(secretKey, logn, logp);

	long n = (1 << logn);
	long n2 = n * n;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n2);

	Ciphertext* cipher = scheme.encrypt(mmat, n, n, logp, logq);
	Ciphertext* tmp = new Ciphertext();

	timeutils.start("Square Matrix Mult");
	for (long i = 0; i < logDegree; ++i) {
		algo.sqrMatSqr(tmp, cipher, logp, n);
		cipher->copy(tmp);
	}
	timeutils.stop("Square Matrix Mult");

	for (long i = 0; i < logDegree; ++i) {
		EvaluatorUtils::squareMatSquareAndEqual(mmat, n);
	}
	complex<double>* dmat = scheme.decrypt(secretKey, cipher);
	StringUtils::compare(mmat, dmat, n2, "matrix");
	cout << cipher->logq << endl;
	cout << "!!! END TEST SQUARE MATRIX POW !!!" << endl;
}

void TestScheme::testMatInv(long logq, long logp, long logn, long steps) {
	cout << "!!! START TEST MATRIX INV !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);
	scheme.addSqrMatKeys(secretKey, logn, logp);

	long n = (1 << logn);
	long n2 = n * n;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n2, 0.01/n);
	for (long i = 0; i < n; ++i) {
		mmat[i + i * n] += 1.0 - EvaluatorUtils::randomReal(0.3);
	}

	Ciphertext* cipher = scheme.encrypt(mmat, n, n, logp, logq);
	Ciphertext* cmatinv = new Ciphertext();

	timeutils.start("Matrix Inv");
	algo.matInv(cmatinv, cipher, logp, n, steps);
	timeutils.stop("Matrix Inv");

	complex<double>* dmatinv = scheme.decrypt(secretKey, cmatinv);
	complex<double>* imat = EvaluatorUtils::squareMatMult(mmat, dmatinv, n);
	StringUtils::showMat(imat, n, n);

	cout << cmatinv->logq << endl;
	cout << "!!! END TEST MATRIX INV !!!" << endl;
}


//----------------------------------------------------------------------------------
//   OTHER TESTS
//----------------------------------------------------------------------------------


void TestScheme::testBootstrap(long logq, long logp, long logn0, long logn1, long logT, long logI) {
	cout << "!!! START TEST BOOTSTRAP !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	timeutils.start("Scheme generating");
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	timeutils.stop("Scheme generated");

	timeutils.start("BootKey generating");
	scheme.addBootKey(secretKey, logn0, logn1, logq + logI);
	timeutils.stop("BootKey generated");

	long n0 = (1 << logn0);
	long n1 = (1 << logn1);
	long n = n0 * n1;

	complex<double>* mmat = EvaluatorUtils::randomComplexSignedArray(n);

	Ciphertext* cipher = scheme.encrypt(mmat, n0, n1, logp, logq);

	cout << "cipher logq before: " << cipher->logq << endl;
	scheme.normalizeAndEqual(cipher);

	cipher->logq = logQ;
	cipher->logp = logq + logI;

	timeutils.start("Coeff to Slot");
	scheme.coeffToSlotAndEqual(cipher);
	timeutils.stop("Coeff to Slot");

	timeutils.start("Remove I Part");
	scheme.removeIPartAndEqual(cipher, logT, logI);
	timeutils.stop("Remove I Part");

	timeutils.start("Slot to Coeff");
	scheme.slotToCoeffAndEqual(cipher);
	timeutils.stop("Slot to Coeff");

	cout << "cipher logp after: " << cipher->logp << endl;
	cout << "cipher logq after: " << cipher->logq << endl;

	cipher->logp = logp;

	complex<double>* dmat = scheme.decrypt(secretKey, cipher);
	StringUtils::compare(mmat, dmat, 10, "boot");

	cout << "!!! END TEST BOOTSRTAP !!!" << endl;
}

void TestScheme::testCiphertextWriteAndRead(long logq, long logp, long logn0, long logn1) {
	cout << "!!! START TEST WRITE AND READ !!!" << endl;
	cout << "!!! END TEST WRITE AND READ !!!" << endl;
}

void TestScheme::test() {
}
