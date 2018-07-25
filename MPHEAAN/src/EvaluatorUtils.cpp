/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "EvaluatorUtils.h"


//----------------------------------------------------------------------------------
//   RANDOM REAL & COMPLEX
//----------------------------------------------------------------------------------


double EvaluatorUtils::randomReal(double bound)  {
	return (double) rand()/(RAND_MAX) * bound;
}

complex<double> EvaluatorUtils::randomComplex(double bound) {
	complex<double> res;
	res.real(randomReal(bound));
	res.imag(randomReal(bound));
	return res;
}

complex<double> EvaluatorUtils::randomCircle(double anglebound) {
	double angle = randomReal(anglebound);
	complex<double> res;
	res.real(cos(angle * 2 * M_PI));
	res.imag(sin(angle * 2 * M_PI));
	return res;
}

double* EvaluatorUtils::randomRealArray(long n, double bound) {
	double* res = new double[n];
	for (long i = 0; i < n; ++i) {
		res[i] = randomReal(bound);
	}
	return res;
}

complex<double>* EvaluatorUtils::randomComplexArray(long n, double bound) {
	complex<double>* res = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		res[i] = randomComplex(bound);
	}
	return res;
}

complex<double>* EvaluatorUtils::randomCircleArray(long n, double bound) {
	complex<double>* res = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		res[i] = randomCircle(bound);
	}
	return res;
}


//----------------------------------------------------------------------------------
//   DOUBLE & RR <-> ZZ
//----------------------------------------------------------------------------------


double EvaluatorUtils::scaleDownToReal(const ZZ& x, const long logp) {
	RR xp = to_RR(x);
	xp.e -= logp;
	return to_double(xp);
}

ZZ EvaluatorUtils::scaleUpToZZ(const double x, const long logp) {
	return scaleUpToZZ(to_RR(x), logp);
}

ZZ EvaluatorUtils::scaleUpToZZ(const RR& x, const long logp) {
	RR xp = MakeRR(x.x, x.e + logp);
	return RoundToZZ(xp);
}


//----------------------------------------------------------------------------------
//   ROTATIONS
//----------------------------------------------------------------------------------


void EvaluatorUtils::leftRotateAndEqual(complex<double>* vals, const long nx, const long ny, const long rx, const long ry) {
	long rxRem = rx % nx;
	if(rxRem != 0) {
		long divisor = GCD(rxRem, nx);
		long steps = nx / divisor;
		for (long iy = 0; iy < ny; ++iy) {
			for (long ix = 0; ix < divisor; ++ix) {
				complex<double> tmp = vals[ix + iy * nx];
				long idx = ix;
				for (long k = 0; k < steps - 1; ++k) {
					vals[idx + iy * nx] = vals[((idx + rxRem) % nx) + iy * nx];
					idx = (idx + rxRem) % nx;
				}
				vals[idx + iy * nx] = tmp;
			}
		}
	}
	long ryRem = ry % ny;
	if(ryRem != 0) {
		long divisor = GCD(ryRem, ny);
		long steps = ny / divisor;
		for (long ix = 0; ix < nx; ++ix) {
			for (long iy = 0; iy < divisor; ++iy) {
				complex<double> tmp = vals[ix + iy * nx];
				long idy = iy;
				for (long k = 0; k < steps - 1; ++k) {
					vals[ix + idy * nx] = vals[ix + ((idy + ryRem) % ny) * nx];
					idy = (idy + ryRem) % ny;
				}
				vals[ix + idy * nx] = tmp;
			}
		}
	}
}

void EvaluatorUtils::rightRotateAndEqual(complex<double>* vals, const long nx, const long ny, const long rx, const long ry) {
	long rxRem = rx % nx;
	rxRem = (nx - rxRem) % nx;
	long ryRem = ry % ny;
	ryRem = (ny - ryRem) % ny;
	leftRotateAndEqual(vals, nx, ny, rxRem, ryRem);
}


//----------------------------------------------------------------------------------
//   MATRIX
//----------------------------------------------------------------------------------


void EvaluatorUtils::squareMatMult(complex<double>* res, complex<double>* vals1, complex<double>* vals2, const long nx) {
	for (long i = 0; i < nx; ++i) {
		for (long j = 0; j < nx; ++j) {
			for (long k = 0; k < nx; ++k) {
				res[i + j * nx] += vals1[k + j * nx] * vals2[i + k * nx];
			}
		}
	}
}

void EvaluatorUtils::squareMatSquareAndEqual(complex<double>* vals, const long nx) {
	long n = nx * nx;
	complex<double>* res = new complex<double>[n];
	for (long ix = 0; ix < nx; ++ix) {
		for (long iy = 0; iy < nx; ++iy) {
			for (long k = 0; k < nx; ++k) {
				res[ix + iy * nx] += vals[k + iy * nx] * vals[ix + k * nx];
			}
		}
	}
	for (long i = 0; i < n; ++i) {
		vals[i] = res[i];
	}
	delete[] res;
}
