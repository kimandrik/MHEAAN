/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MPHEAAN_CIPHERTEXT_H_
#define MPHEAAN_CIPHERTEXT_H_

#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class Ciphertext {
public:

	ZZ* axy;
	ZZ* bxy;

	long logp;
	long logq;

	long Nx;
	long Ny;

	long nx;
	long ny;

	bool isComplex;

	//-----------------------------------------

	Ciphertext(ZZ* axy = NULL, ZZ* bxy = NULL, long logp = 0, long logq = 0, long Nx = 1, long Ny = 1, long nx = 1, long ny = 1, bool isComplex = true);

	Ciphertext(const Ciphertext& o);

	Ciphertext& operator=(const Ciphertext &o);

	~Ciphertext();
};

#endif
