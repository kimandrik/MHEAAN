/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MHEAAN_STRINGUTILS_H_
#define MHEAAN_STRINGUTILS_H_

#include "Common.h"

#include <complex>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class StringUtils {
public:

	static void showVec(long* vals, long n);
	static void showVec(double* vals, long n);
	static void showVec(complex<double>* vals, long n);
	static void showVec(ZZ* vals, long n);

	static void showMat(long* vals, long nx, long ny);
	static void showMat(double* vals, long nx, long ny);
	static void showMat(complex<double>* vals, long nx, long ny);
	static void showMat(ZZ* vals, long nx, long ny);

	static void compare(complex<double> val1, complex<double> val2, string prefix);
	static void compare(complex<double>* vals1, complex<double>* vals2, long n, string prefix);
	static void compare(complex<double>* vals1, complex<double> val2, long n, string prefix);
	static void compare(complex<double> val1, complex<double>* vals2, long n, string prefix);

};

#endif
