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

#include "Common.h"
#include "Ciphertext.h"
#include "EvaluatorUtils.h"
#include "Ring2XY.h"
#include "Scheme.h"
#include "SchemeAlgo.h"
#include "SecretKey.h"
#include "StringUtils.h"
#include "TimeUtils.h"
#include "Numb.h"

using namespace std;
using namespace NTL;


//----------------------------------------------------------------------------------
//   STANDARD TESTS
//----------------------------------------------------------------------------------


void TestScheme::testEncode(long logNx, long logQ, long logp, long lognx) {
	long Ny = 256;
	Ring2XY ring(logNx, logQ);
	uint64_t* mx1 = new uint64_t[Ny]();
	uint64_t* mx2 = new uint64_t[Ny]();
	for (long i = 0; i < Ny; ++i) {
		mx1[i] = 1;
		mx2[i] = 1;
	}
	ring.multiplier.NTTY(mx1, 0);
	ring.multiplier.NTTY(mx2, 0);

	uint64_t p = ring.multiplier.pVec[0];
	for (long i = 0; i < Ny; ++i) {
		mulMod(mx1[i], mx1[i], mx2[i], p);
	}
	ring.multiplier.INTTY(mx1, 0);
	ring.multiplier.INTTY(mx2, 0);

	for (long i = 0; i < Ny; ++i) {
		cout << mx1[i] << " ";
	}
	cout << endl;

	for (long i = 0; i < Ny; ++i) {
		cout << mx2[i] << " ";
	}
	cout << endl;

}

void TestScheme::testEncrypt(long logNx, long logNy, long logQ, long logp, long lognx, long logny) {
	cout << "!!! START TEST ENCRYPT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ring.Ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);

	timeutils.start("Encode matrix");
	Plaintext msg = scheme.encode(mmat, nx, ny, logp, logQ);
	timeutils.stop("Encode matrix");

	timeutils.start("Encrypt msg");
	Ciphertext cipher = scheme.encryptMsg(msg);
	timeutils.stop("Encrypt msg");

	timeutils.start("Decrypt msg");
	Plaintext dsg = scheme.decryptMsg(secretKey, cipher);
	timeutils.stop("Decrypt msg");

	timeutils.start("Decode matrix");
	complex<double>* dmat = scheme.decode(dsg);
	timeutils.stop("Decode matrix");

	StringUtils::compare(mmat, dmat, n, "val");

	cout << "!!! END TEST ENCRYPT !!!" << endl;
}

void TestScheme::testEncryptSingle(long logNx, long logNy, long logQ, long logp) {
	cout << "!!! START TEST ENCRYPT SINGLE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy , logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	complex<double> mval = EvaluatorUtils::randomComplex();

	timeutils.start("Encrypt single");
	Ciphertext cipher = scheme.encryptSingle(mval, logp, logQ);
	timeutils.stop("Encrypt single");

	timeutils.start("Decrypt single");
	complex<double> dval = scheme.decryptSingle(secretKey, cipher);
	timeutils.stop("Decrypt single");

	StringUtils::compare(mval, dval, "val");

	cout << "!!! END TEST ENCRYPT SINGLE !!!" << endl;
}

void TestScheme::testStandard(long logNx, long logNy, long logQ, long logp, long lognx, long logny) {
	cout << "!!! START TEST STANDARD !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy , logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat1 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmat2 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* madd = new complex<double>[n];
	complex<double>* mmult = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mmult[i] = mmat1[i] * mmat2[i];
		madd[i] = mmat1[i] + mmat2[i];
	}
	Ciphertext cipher1 = scheme.encrypt(mmat1, nx, ny, logp, logQ);
	Ciphertext cipher2 = scheme.encrypt(mmat2, nx, ny, logp, logQ);

	timeutils.start("add matrix");
	Ciphertext cadd = scheme.add(cipher1, cipher2);
	timeutils.stop("add matrix");

	timeutils.start("mult matrix");
	Ciphertext cmult = scheme.mult(cipher1, cipher2);
	scheme.reScaleByAndEqual(cmult, logp);
	timeutils.stop("mult matrix");

	complex<double>* dadd = scheme.decrypt(secretKey, cadd);
	complex<double>* dmult = scheme.decrypt(secretKey, cmult);

	StringUtils::compare(madd, dadd, n, "add");
	StringUtils::compare(mmult, dmult, n, "mult");

	cout << "!!! END TEST STANDARD !!!" << endl;
}

void TestScheme::testimult(long logNx, long logNy, long logQ, long logp, long lognx, long logny) {
	cout << "!!! START TEST i MULTIPLICATION !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmatimult = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mmatimult[i].real(-mmat[i].imag());
		mmatimult[i].imag(mmat[i].real());
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Multiplication by i");
	Ciphertext cimult = scheme.imult(cipher);
	timeutils.stop("Multiplication by i");

	complex<double>* dmatimult = scheme.decrypt(secretKey, cimult);

	StringUtils::compare(mmatimult, dmatimult, n, "imult");

	cout << "!!! END TEST i MULTIPLICATION !!!" << endl;
}


//----------------------------------------------------------------------------------
//   ROTATION & CONJUGATION & TRANSPOSITION TESTS
//----------------------------------------------------------------------------------


void TestScheme::testRotateFast(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long rx, long ry) {
	cout << "!!! START TEST ROTATE FAST !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addLeftRotKey(secretKey, rx, ry);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Left rotate fast");
	scheme.leftRotateFastAndEqual(cipher, rx, ry);
	timeutils.stop("Left rotate fast");

	complex<double>* dmat = scheme.decrypt(secretKey, cipher);
	EvaluatorUtils::leftRotateAndEqual(mmat, nx, ny, rx, ry);
	StringUtils::compare(mmat, dmat, n, "val");

	cout << "!!! END TEST ROTATE FAST !!!" << endl;
}

void TestScheme::testRotate(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long rx, long ry) {
	cout << "!!! START TEST ROTATE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addLeftXRotKeys(secretKey);
	scheme.addLeftYRotKeys(secretKey);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Left rotate");
	scheme.leftRotateAndEqual(cipher, rx, ry);
	timeutils.stop("Left rotate");

	complex<double>* dmat = scheme.decrypt(secretKey, cipher);
	EvaluatorUtils::leftRotateAndEqual(mmat, nx, ny, rx, ry);
	StringUtils::compare(mmat, dmat, n, "val");

	cout << "!!! END TEST ROTATE !!!" << endl;
}

void TestScheme::testConjugate(long logNx, long logNy, long logQ, long logp, long lognx, long logny) {
	cout << "!!! START TEST CONJUGATE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addConjKey(secretKey);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmatconj = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mmatconj[i] = conj(mmat[i]);
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Conjugate");
	Ciphertext cconj = scheme.conjugate(cipher);
	timeutils.stop("Conjugate");

	complex<double>* dmatconj = scheme.decrypt(secretKey, cconj);
	StringUtils::compare(mmatconj, dmatconj, n, "conj");

	cout << "!!! END TEST CONJUGATE !!!" << endl;
}

void TestScheme::testTranspose(long logNx, long logNy, long logQ, long logp, long lognx, long logny) {
	cout << "!!! START TEST TRANSPOSE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addTranspKey(secretKey);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmattr = new complex<double>[n];
	for (long j = 0; j < nx; ++j) {
		for (long i = 0; i < ny; ++i) {
			mmattr[i + j * ny] = mmat[j + i * nx];
		}
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Transpose");
	Ciphertext ctr = scheme.transpose(cipher);
	timeutils.stop("Transpose");

	complex<double>* dvectr = scheme.decrypt(secretKey, ctr);
	StringUtils::compare(mmattr, dvectr, n, "tr");

	cout << "!!! END TEST TRANSPOSE !!!" << endl;
}


//----------------------------------------------------------------------------------
//   POWER & PRODUCT TESTS
//----------------------------------------------------------------------------------


void TestScheme::testPowerOf2(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long logDegree) {
	cout << "!!! START TEST POWER OF 2 !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = 1 << lognx;
	long ny = 1 << logny;
	long n = nx * ny;
	long degree = 1 << logDegree;

	complex<double>* mmat = EvaluatorUtils::randomCircleArray(n);
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mpow[i] = pow(mmat[i], degree);
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Power of 2");
	Ciphertext cpow = algo.powerOf2(cipher, logp, logDegree);
	timeutils.stop("Power of 2");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	StringUtils::compare(mpow, dpow, n, "pow");

	cout << "!!! END TEST POWER OF 2 !!!" << endl;
}

void TestScheme::testPower(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree) {
	cout << "!!! START TEST POWER !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomCircleArray(n);
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mpow[i] = pow(mmat[i], degree);
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Power");
	Ciphertext cpow = algo.power(cipher, logp, degree);
	timeutils.stop("Power");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	StringUtils::compare(mpow, dpow, n, "pow");

	cout << "!!! END TEST POWER !!!" << endl;
}

//-----------------------------------------

void TestScheme::testProdOfPo2(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long logDegree) {
	cout << "!!! START TEST PROD OF POWER OF 2 !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;
	long degree = 1 << logDegree;

	complex<double>** mmatvec = new complex<double>*[degree];
	for (long i = 0; i < degree; ++i) {
		mmatvec[i] = EvaluatorUtils::randomCircleArray(n);
	}

	complex<double>* pmat = new complex<double>[n];
	for (long j = 0; j < n; ++j) {
		pmat[j] = mmatvec[0][j];
		for (long i = 1; i < degree; ++i) {
			pmat[j] *= mmatvec[i][j];
		}
	}

	Ciphertext* cvec = new Ciphertext[degree];
	for (long i = 0; i < degree; ++i) {
		cvec[i] = scheme.encrypt(mmatvec[i], nx, ny, logp, logQ);
	}

	timeutils.start("Product of power of 2");
	Ciphertext cprod = algo.prodOfPo2(cvec, logp, logDegree);
	timeutils.stop("Product of power of 2");

	complex<double>* dmat = scheme.decrypt(secretKey, cprod);
	StringUtils::compare(pmat, dmat, n, "prod");

	cout << "!!! END TEST PROD OF POWER OF 2 !!!" << endl;
}

void TestScheme::testProd(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree) {
	cout << "!!! START TEST PROD !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>** mmatvec = new complex<double>*[degree];
	for (long i = 0; i < degree; ++i) {
		mmatvec[i] = EvaluatorUtils::randomCircleArray(n);
	}

	complex<double>* pmat = new complex<double>[n];
	for (long j = 0; j < n; ++j) {
		pmat[j] = mmatvec[0][j];
		for (long i = 1; i < degree; ++i) {
			pmat[j] *= mmatvec[i][j];
		}
	}

	Ciphertext* cvec = new Ciphertext[degree];
	for (long i = 0; i < degree; ++i) {
		cvec[i] = scheme.encrypt(mmatvec[i], nx, ny, logp, logQ);
	}

	timeutils.start("Product");
	Ciphertext cprod = algo.prod(cvec, logp, degree);
	timeutils.stop("Product");

	complex<double>* dmat = scheme.decrypt(secretKey, cprod);
	StringUtils::compare(pmat, dmat, n, "prod");

	cout << "!!! END TEST PROD !!!" << endl;
}


//----------------------------------------------------------------------------------
//   FUNCTION TESTS
//----------------------------------------------------------------------------------


void TestScheme::testInverse(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long steps) {
	cout << "!!! START TEST INVERSE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomCircleArray(n, 0.1);
	complex<double>* minv = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		minv[i] = 1. / mmat[i];
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start("Inverse");
	Ciphertext cinv = algo.inverse(cipher, logp, steps);
	timeutils.stop("Inverse");

	complex<double>* dinv = scheme.decrypt(secretKey, cinv);
	StringUtils::compare(minv, dinv, n, "inv");

	cout << "!!! END TEST INVERSE !!!" << endl;
}

//-----------------------------------------

void TestScheme::testLogarithm(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree) {
	cout << "!!! START TEST LOGARITHM !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n, 0.1);
	complex<double>* mlog = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mlog[i] = log(mmat[i] + 1.);
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start(LOGARITHM);
	Ciphertext clog = algo.function(cipher, LOGARITHM, logp, degree);
	timeutils.stop(LOGARITHM);

	complex<double>* dlog = scheme.decrypt(secretKey, clog);
	StringUtils::compare(mlog, dlog, n, LOGARITHM);

	cout << "!!! END TEST LOGARITHM !!!" << endl;
}

void TestScheme::testExponent(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree) {
	cout << "!!! START TEST EXPONENT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mmat[i]);
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start(EXPONENT);
	Ciphertext cexp = algo.function(cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT);

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	StringUtils::compare(mexp, dexp, n, EXPONENT);

	cout << "!!! END TEST EXPONENT !!!" << endl;
}

void TestScheme::testExponentLazy(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree) {
	cout << "!!! START TEST EXPONENT LAZY !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mmat[i]);
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start(EXPONENT + " lazy");
	Ciphertext cexp = algo.functionLazy(cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT + " lazy");

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	StringUtils::compare(mexp, dexp, n, EXPONENT);

	cout << "!!! END TEST EXPONENT LAZY !!!" << endl;
}

void TestScheme::testSigmoid(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree) {
	cout << "!!! START TEST SIGMOID !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mmat[i]) / (1. + exp(mmat[i]));
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start(SIGMOID);
	Ciphertext csig = algo.function(cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID);

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	StringUtils::compare(msig, dsig, n, SIGMOID);

	cout << "!!! END TEST SIGMOID !!!" << endl;
}

void TestScheme::testSigmoidLazy(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long degree) {
	cout << "!!! START TEST SIGMOID LAZY !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mmat[i]) / (1. + exp(mmat[i]));
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, ny, logp, logQ);

	timeutils.start(SIGMOID + " lazy");
	Ciphertext csig = algo.functionLazy(cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID + " lazy");

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	StringUtils::compare(msig, dsig, n, SIGMOID);

	cout << "!!! END TEST SIGMOID LAZY !!!" << endl;
}


//----------------------------------------------------------------------------------
//   MATRIX TESTS
//----------------------------------------------------------------------------------


void TestScheme::testSquareMatMult(long logNx, long logNy, long logQ, long logp, long lognx) {
	cout << "!!! START TEST SQUARE MATRIX !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	scheme.addSquareMatrixKeys(secretKey, lognx);

	long nx = (1 << lognx);
	long n = nx * nx;

	complex<double>* mmat1 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmat2 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmatmult = new complex<double>[n];
	EvaluatorUtils::squareMatMult(mmatmult, mmat1, mmat2, nx);

	Ciphertext cipher1 = scheme.encrypt(mmat1, nx, nx, logp, logQ);
	Ciphertext cipher2 = scheme.encrypt(mmat2, nx, nx, logp, logQ);

	timeutils.start("Square Matrix Mult");
	Ciphertext cmatmult = algo.squareMatMult(cipher1, cipher2, logp, nx);
	timeutils.stop("Square Matrix Mult");

	complex<double>* dmatmult = scheme.decrypt(secretKey, cmatmult);
	StringUtils::compare(mmatmult, dmatmult, n, "matrix");

	cout << "!!! END TEST SQUARE MATRIX !!!" << endl;
}

void TestScheme::testSquareMatPow(long logNx, long logNy, long logQ, long logp, long lognx, long logDegree) {
	cout << "!!! START TEST SQUARE MATRIX POW!!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);
	scheme.addSquareMatrixKeys(secretKey, lognx);

	long nx = (1 << lognx);
	long n = nx * nx;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n, 1.0/nx);

	Ciphertext cipher = scheme.encrypt(mmat, nx, nx, logp, logQ);

	timeutils.start("Square Matrix Mult");
	for (long i = 0; i < logDegree; ++i) {
		algo.squareMatMultAndEqual(cipher, logp, nx);
	}
	timeutils.stop("Square Matrix Mult");

	for (long i = 0; i < logDegree; ++i) {
		EvaluatorUtils::squareMatSquareAndEqual(mmat, nx);
	}
	complex<double>* dmat = scheme.decrypt(secretKey, cipher);
	StringUtils::compare(mmat, dmat, n, "matrix");

	cout << "!!! END TEST SQUARE MATRIX POW !!!" << endl;
}

void TestScheme::testSquareMatInv(long logNx, long logNy, long logQ, long logp, long lognx, long steps) {
	cout << "!!! START TEST MATRIX INV !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);
	scheme.addSquareMatrixKeys(secretKey, lognx);

	long nx = (1 << lognx);
	long n = nx * nx;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n, 0.01/nx);
	for (long i = 0; i < nx; ++i) {
		mmat[i + i * nx] += 1.0 - EvaluatorUtils::randomReal(0.3);
	}

	Ciphertext cipher = scheme.encrypt(mmat, nx, nx, logp, logQ);

	timeutils.start("Matrix Inv");
	Ciphertext cmatinv = algo.matInv(cipher, logp, nx, steps);
	timeutils.stop("Matrix Inv");

	complex<double>* dmatinv = scheme.decrypt(secretKey, cmatinv);
	complex<double>* imat = new complex<double>[n];
	EvaluatorUtils::squareMatMult(imat, mmat, dmatinv, nx);
	StringUtils::showMat(imat, nx, nx);

	cout << "!!! END TEST MATRIX INV !!!" << endl;
}

void TestScheme::testMatMult(long logNx, long logNy, long logQ, long logp, long lognx, long logny, long lognz) {

	//TODO implement the right method

	//	long logNx = 7;
	//	long logNy = 7;
	//	long logQ = 200;
	//	long logp = 30;
	//	long logSlotx = 2;
	//	long logSloty = 3;
	//	cout << "!!! START TEST NMN MATRIX !!!" << endl;
	//	//-----------------------------------------
	//	TimeUtils timeutils;
	//	Ring2XY ring(logNx, logNy , logQ);
	//	SecretKey secretKey(ring);
	//	Scheme scheme(secretKey, ring);
	//	SchemeAlgo algo(scheme);
	//	//-----------------------------------------
	//	srand(time(NULL));
	//	//-----------------------------------------
	//	long slotx = (1 << logSlotx);
	//	long sloty = (1 << logSloty);
	//	long slots = slotx * sloty;
	//	scheme.addSquareMatrixKeys(secretKey, slotx);
	//	scheme.addSquareMatrixKeys(secretKey, sloty);
	//
	//	complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
	//	complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots);
	//
	//	complex<double>* mvecM = new complex<double>[slots];
	//	for (long i = 0; i < slotx; ++i) {
	//		for (long j = 0; j < sloty; ++j) {
	//			for (long k = 0; k < slotx; ++k) {
	//				mvecM[i + j * slotx] += mvec1[k + j * slotx] * mvec2[i + k * sloty];
	//			}
	//		}
	//	}
	//
	//	Ciphertext cipher1 = scheme.encrypt(mvec1, slotx, sloty, logp, logQ);
	//	Ciphertext cipher2 = scheme.encrypt(mvec2, sloty, slotx, logp, logQ);
	//
	//	ZZ** ptmpm = new ZZ*[slotx];
	//
	//	for (long i = 0; i < slotx; ++i) {
	//		complex<double>* tmp = new complex<double>[slots];
	//		for (long j = 0; j < i; ++j) {
	//			tmp[j + (j - i + slotx) * slotx] = 1;
	//		}
	//		for (long j = i; j < sloty; ++j) {
	//			tmp[j + (j - i) * slotx] = 1;
	//		}
	//		ring.encode(ptmpm[i], tmp, slotx, slotx, logp);
	//		delete[] tmp;
	//	}
	//
	//	timeutils.start("NMN Matrix Mult");
	//	Ciphertext* cipherP = new Ciphertext[slotx];
	//	Ciphertext* cipherR = new Ciphertext[slotx];
	//
	//	cipherR[0] = cipher1;
	//	scheme.modDownByAndEqual(cipherR[0], logp);
	//
	//	cipherP[0] = scheme.multByPoly(cipher2, ptmpm[0], logp);
	//	scheme.reScaleByAndEqual(cipherP[0], logp);
	//	for (long j = 0; j < logSlotx; ++j) {
	//		Ciphertext rot = scheme.leftYRotateByPo2(cipherP[0], j);
	//		scheme.addAndEqual(cipherP[0], rot);
	//	}
	//	scheme.multAndEqual(cipherR[0], cipherP[0]);
	//
	//	for (long i = 1; i < slotx; ++i) {
	//		cipherR[i] = scheme.rightXRotateFast(cipher1, i);
	//		scheme.modDownByAndEqual(cipherR[i], logp);
	//
	//		cipherP[i] = scheme.multByPoly(cipher2, ptmpm[i], logp);
	//		scheme.reScaleByAndEqual(cipherP[i], logp);
	//		for (long j = 0; j < logSlotx; ++j) {
	//			Ciphertext rot = scheme.leftYRotateByPo2(cipherP[i], j);
	//			scheme.addAndEqual(cipherP[i], rot);
	//		}
	//		scheme.multAndEqual(cipherR[i], cipherP[i]);
	//		scheme.addAndEqual(cipherR[0], cipherR[i]);
	//	}
	//
	//	scheme.reScaleByAndEqual(cipherR[0], logp);
	//	timeutils.stop("NMN Matrix Mult");
	//
	//	complex<double>* dvecM = scheme.decrypt(secretKey, cipherR[0]);
	//
	//	StringUtils::compare(mvecM, dvecM, slots, "matrix");
	//
	//	cout << "!!! END TEST NMN MATRIX !!!" << endl;
}


//----------------------------------------------------------------------------------
//   OTHER TESTS
//----------------------------------------------------------------------------------


void TestScheme::testCiphertextWriteAndRead(long logNx, long logNy, long logQ, long logp, long logSlotx, long logSloty) {
	cout << "!!! START TEST WRITE AND READ !!!" << endl;
	cout << "!!! END TEST WRITE AND READ !!!" << endl;
}

void TestScheme::test(long logNx, long logNy, long logQ, long logp, long lognx, long logny) {
	cout << "!!! START TEST ENCRYPT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);

	TimeUtils timeutils;
	Ring2XY ring(logNx, logNy , logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long nx = (1 << lognx);
	long ny = (1 << logny);
	long n = nx * ny;

	complex<double>* mmat = EvaluatorUtils::randomComplexArray(n);

	timeutils.start("Encode matrix");
	Plaintext msg = scheme.encode(mmat, nx, ny, logp, logQ);
	timeutils.stop("Encode matrix");

	timeutils.start("Encrypt msg");
	Ciphertext cipher = scheme.encryptMsg(msg);
	timeutils.stop("Encrypt msg");

	scheme.modDownToAndEqual(cipher, 300);
	timeutils.start("Decrypt matrix");
	complex<double>* dmat = scheme.decrypt(secretKey, cipher);
	timeutils.stop("Decrypt matrix");

	StringUtils::compare(mmat, dmat, n, "val");

	cout << "!!! END TEST ENCRYPT !!!" << endl;

}

