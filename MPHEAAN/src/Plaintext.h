/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef MPHEAAN_PLAINTEXT_H_
#define MPHEAAN_PLAINTEXT_H_

#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class Plaintext {
public:

	ZZ* mxy;

	long logp;
	long logq;

	long Nx;
	long Ny;

	long nx;
	long ny;

	bool isComplex;

	//-----------------------------------------

	Plaintext(ZZ* mxy = NULL, long logp = 0, long logq = 0, long Nx = 1, long Ny = 1, long nx = 1, long ny = 1, bool isComplex = true);

	Plaintext(const Plaintext& o);

	Plaintext& operator=(const Plaintext &o);

	~Plaintext();

};

#endif
