/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef MPHEAAN_EVALUATORUTILS_H_
#define MPHEAAN_EVALUATORUTILS_H_

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <complex>

using namespace std;
using namespace NTL;

class EvaluatorUtils {
public:


	//----------------------------------------------------------------------------------
	//   RANDOM REAL AND COMPLEX NUMBERS
	//----------------------------------------------------------------------------------


	static double randomReal(double bound = 1.0);

	static double randomRealSigned(double bound = 1.0);

	static complex<double> randomComplex(double bound = 1.0);

	static complex<double> randomComplexSigned(double bound = 1.0);

	static complex<double> randomCircle(double anglebound = 1.0);

	static double* randomRealArray(long n, double bound = 1.0);

	static double* randomRealSignedArray(long n, double bound = 1.0);

	static complex<double>* randomComplexArray(long n, double bound = 1.0);

	static complex<double>* randomComplexSignedArray(long n, double bound = 1.0);

	static complex<double>* randomCircleArray(long n, double bound = 1.0);


	//----------------------------------------------------------------------------------
	//   DOUBLE & RR <-> ZZ
	//----------------------------------------------------------------------------------


	static double scaleDownToReal(const ZZ& x, const long logp);

	static ZZ scaleUpToZZ(const double x, const long logp);

	static ZZ scaleUpToZZ(const RR& x, const long logp);


	//----------------------------------------------------------------------------------
	//   ROTATIONS
	//----------------------------------------------------------------------------------


	static void leftRotateAndEqual(complex<double>* vals, const long nx, const long ny, const long rx, const long ry);

	static void rightRotateAndEqual(complex<double>* vals, const long nx, const long ny, const long rx, const long ry);


	//----------------------------------------------------------------------------------
	//   MATRIX
	//----------------------------------------------------------------------------------


	static void squareMatMult(complex<double>* res, complex<double>* vals1, complex<double>* vals2, const long nx);

	static void squareMatSquareAndEqual(complex<double>* vals, const long nx);
};

#endif
